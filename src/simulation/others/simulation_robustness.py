"""Supplementary robustness simulation driver for traceCB.

Purpose
-------
This script starts from the base ``simulation.py`` workflow and adds robustness
switches that are not part of the main simulation driver:
``--gmm_propt_mode``, ``--causal_overlap``, and
``--causal_partition_mode pop2_a_shared_b``.

Mode details
------------
``--gmm_propt_mode`` changes only the cell-type proportion supplied to
``GMMtissue``. The simulated tissue data are still generated from the true
``--propt`` value, so this mode tests sensitivity to misspecifying the cell-type
mixture proportion during inference. ``exact`` supplies the true ``propt``.
``underestimate`` supplies ``propt - --gmm_propt_mode_scale`` after truncating
to ``[0, 1]``. ``overestimate`` supplies ``propt + --gmm_propt_mode_scale``
after the same truncation. ``normal`` draws the supplied value from
``Normal(propt, sqrt(--gmm_propt_normal_var))`` and truncates it to ``[0, 1]``.

``--causal_overlap`` changes how causal SNP sets are generated across the two
populations. When omitted, the base simulation is used: pop1 and pop2 share the
same causal SNPs, and their effects are drawn from a bivariate normal with
genetic correlation ``--gc``. When set to a value in ``[0, 1]``, each population
has ``pcausal * nsnp`` causal SNPs, and the requested fraction of pop1 causal
SNPs is also causal in pop2. Shared causal SNPs use aligned standardized
effects; population-specific causal SNPs get independent effects. In this mode,
``--gc`` is not used to generate effect-size correlation.

``--causal_partition_mode pop2_a_shared_b`` creates a segmented null design.
The default ``none`` is the original base simulation where causal SNPs are sampled
independently across the genome. The ``pop2_a_shared_b`` mode creates two regions
with different causal architectures.
The first ``--null_region_prop`` fraction of SNPs is region A, where pop1 is
null but pop2 can have causal effects. The remaining SNPs are region B, where
pop1 and pop2 share causal SNPs with effect correlation controlled by ``--gc``.
This mode is used to evaluate type I error in a pop1-null region that still has
pop2 signal.

Typical workflow
----------------
Run from the repository root in conda env ``py312``:

```
conda activate py312
python3 src/simulation/others/simulation_robustness.py ...
python3 src/simulation/visual_simulation.py ...
```

The robustness command grid is collected in
``src/simulation/others/run_robustness.sh``.

Outputs
-------
Each parameter setting is written under ``<out_dir>/<runname>/<setting>/`` as
``simulation_<rep>.csv`` files. The visualizer collapses those replicate files
into ``result_df.csv`` and writes paper figures into ``--img_dir`` or ``$IMG_DIR``.

Important modeling switches
---------------------------
``--estimate_omega`` estimates omega with cross-LDSC; otherwise the true
simulated omega is used. ``--gmm_propt_mode`` perturbs the mean cell-type
proportion supplied to GMM tissue. ``--causal_overlap`` and
``--causal_partition_mode pop2_a_shared_b`` implement shared-causal-SNP and
segmented-null architectures used in supplementary runs.
"""

import argparse
import time
import os
import sys
from pathlib import Path

import numpy as np

SIMULATION_DIR = Path(__file__).resolve().parents[1]
ROOT_DIR = SIMULATION_DIR.parents[1]
SRC_DIR = ROOT_DIR / "src"
for path in (SRC_DIR, SIMULATION_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from traceCB.ldsc import Run_Cross_LDSC
from traceCB.utils import z2p, MIN_HERITABILITY
from simulation import (
    MIN_FLOAT,
    P_VAL_THRED,
    cal_ld,
    calculate_sumstats,
    get_genotype,
    run_gmm_meta_kernel,
    should_run_true_omega_gmm,
)
from simulation_utils import (
    calculate_pi2_omega_sum_const,
    flatten_float_seq,
    iter_gmm_propt_subsettings,
    make_sim_seed,
    perturb_gmm_propt,
    sanitize_ld_scores,
    seed_random_component,
    standardize_genotype,
    tail_panel_start,
    unknown_cell_effect_scale,
    validate_nonnegative,
    validate_unit_interval,
)

GMM_PROPT_MODES = (
    "exact",
    "underestimate",
    "overestimate",
    "normal",
)
CAUSAL_PARTITION_MODES = (
    "none",
    "pop2_a_shared_b",
)


def parse_args():
    parser = argparse.ArgumentParser(description="Generate data and simulation")
    parser.add_argument(
        "--pop1_geno",
        default="data/simulation/EAS_n5000_chr22_loci29.npy",
        type=str,
        help="Path to population 1 genotype file",
    )
    parser.add_argument(
        "--pop2_geno",
        default="data/simulation/EUR_n20000_chr22_loci29.npy",
        type=str,
        help="Path to population 2 genotype file",
    )
    parser.add_argument(
        "--runname",
        default="power",
        type=str,
        help="runname of this simulation",
    )
    parser.add_argument(
        "--h1sq",
        default=0.1,
        type=float,
        nargs="+",
        help="Heritability of population 1",
    )
    parser.add_argument(
        "--h2sq",
        default=0.1,
        type=float,
        nargs="+",
        help="Heritability of population 2",
    )
    parser.add_argument(
        "--gc",
        default=0.2,
        type=float,
        nargs="+",
        help="Genetic correlation between population 1 and 2",
    )
    parser.add_argument(
        "--n1",
        default=100,
        type=int,
        nargs="+",
        help="Sample size of population 1",
    )
    parser.add_argument(
        "--n2",
        default=200,
        type=int,
        nargs="+",
        help="Sample size of population 2",
    )
    parser.add_argument(
        "--nt",
        default=1000,
        type=int,
        nargs="+",
        help="Sample size of tissue of population 2",
    )
    parser.add_argument(
        "--nsnp",
        default=1000,
        type=int,
        help="Number of SNPs",
    )
    parser.add_argument(
        "--propt",
        default=0.1,
        type=float,
        nargs="+",
        help="Proportion of cell type in tissue",
    )
    parser.add_argument(
        "--pcausal",
        default=0.1,
        type=float,
        nargs="+",
        help="Proportion of causal SNPs",
    )
    parser.add_argument(
        "--causal_overlap",
        default=None,
        type=float,
        nargs="+",
        help=(
            "Optional proportion of pop1 causal SNPs also causal in pop2. "
            "When set, shared causal SNPs use aligned standardized effects and "
            "--gc is not used to generate effect-size correlation. Omit to keep "
            "the original --gc-controlled logic."
        ),
    )
    parser.add_argument(
        "--causal_max_abs_cor",
        default=None,
        type=float,
        help=(
            "Optional upper bound on pairwise absolute SNP correlation among "
            "selected causal SNPs. Correlation is max(abs(cor_pop1), "
            "abs(cor_pop2)). Omit to keep unconstrained causal SNP sampling."
        ),
    )
    parser.add_argument(
        "--causal_partition_mode",
        default="none",
        type=str,
        choices=CAUSAL_PARTITION_MODES,
        help=(
            "Optional segmented causal architecture. pop2_a_shared_b uses the "
            "first --null_region_prop SNPs as pop1-null region A with pop2-only "
            "causal SNPs, and the remaining SNPs as shared non-null region B."
        ),
    )
    parser.add_argument(
        "--null_region_prop",
        default=0.8,
        type=float,
        nargs="+",
        help="Proportion of SNPs assigned to pop1-null region A.",
    )
    parser.add_argument(
        "--estimate_omega",
        action="store_true",
        help="Estimate omega from data, otherwise use true omega",
    )
    parser.add_argument(
        "--gmm_propt_mode",
        default=["exact"],
        type=str,
        nargs="+",
        choices=GMM_PROPT_MODES,
        help="How to perturb the mean cell-type proportion supplied to GMM tissue.",
    )
    parser.add_argument(
        "--gmm_propt_mode_scale",
        default=0.1,
        type=float,
        nargs="+",
        help="Scale for --gmm_propt_mode underestimate and overestimate.",
    )
    parser.add_argument(
        "--gmm_propt_normal_var",
        default=0.01,
        type=float,
        nargs="+",
        help="Variance for --gmm_propt_mode normal; draws are truncated to [0, 1].",
    )
    parser.add_argument(
        "--out_dir",
        default="bench/result",
        type=str,
        help="Output directory",
    )
    parser.add_argument(
        "--nrep",
        type=int,
        default=100,
        help="number of repetition for simulation",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help=(
            "Base random seed. When omitted, the current start time is used. "
            "Only the global replicate id is used in seed derivation, so the "
            "same replicate reuses component streams across simulation settings."
        ),
    )
    return parser.parse_args()


def safe_omega_corr(omega):
    denom = np.sqrt(
        max(float(omega[0, 0]), MIN_HERITABILITY)
        * max(float(omega[1, 1]), MIN_HERITABILITY)
    )
    return float(omega[0, 1] / (denom + MIN_FLOAT))


def calculate_true_weighted_omega_o(omega_cb, propt):
    """Return the proportion-weighted other-cell contribution in tissue."""
    pi2_omega_sum_const = calculate_pi2_omega_sum_const(propt)
    return float(
        omega_cb[2, 2]
        - propt**2 * omega_cb[1, 1]
        - pi2_omega_sum_const * (omega_cb[1, 2] - propt * omega_cb[1, 1])
    )


def calculate_estimated_weighted_omega_o(aux_omega_matrix, omega, gmm_propt):
    """Return the cross-LDSC plug-in estimate of weighted omega_o."""
    pi2_omega_sum_const = calculate_pi2_omega_sum_const(gmm_propt)
    return float(
        aux_omega_matrix[1, 1]
        - gmm_propt**2 * omega[1, 1]
        - pi2_omega_sum_const
        * np.maximum(aux_omega_matrix[0, 1] - gmm_propt * omega[1, 1], 0)
    )


def clip_weighted_omega_o(weighted_omega_o):
    return float(np.maximum(weighted_omega_o, MIN_HERITABILITY))


def save_omega_summary(
    out_dir,
    simulation_path,
    id_sim,
    omega_true,
    omega_est,
    omega_est_se,
    omega_est_p,
    pi2_omega_sum,
    weighted_omega_o_true_raw,
    weighted_omega_o_true_clipped,
    weighted_omega_o_est_raw,
    weighted_omega_o_est_clipped,
    pi_mean,
    propt,
    gmm_propt,
    run_gmm,
    run_gmm_tissue,
):
    columns = [
        "id_sim",
        "omega11_true",
        "omega22_true",
        "omega12_true",
        "omega11_est",
        "omega22_est",
        "omega12_est",
        "omega11_se",
        "omega22_se",
        "omega12_se",
        "omega11_p",
        "omega22_p",
        "omega12_p",
        "omega11_diff",
        "omega22_diff",
        "omega12_diff",
        "corr_true",
        "corr_est",
        "pi2_omega_sum",
        "weighted_omega_o_true_raw",
        "weighted_omega_o_true_clipped",
        "weighted_omega_o_est_raw",
        "weighted_omega_o_est_clipped",
        "weighted_omega_o_used",
        "pi_mean",
        "propt",
        "gmm_propt",
        "run_gmm",
        "run_gmm_tissue",
    ]
    values = np.array(
        [
            id_sim,
            omega_true[0, 0],
            omega_true[1, 1],
            omega_true[0, 1],
            omega_est[0, 0],
            omega_est[1, 1],
            omega_est[0, 1],
            omega_est_se[0, 0],
            omega_est_se[1, 1],
            omega_est_se[0, 1],
            omega_est_p[0, 0],
            omega_est_p[1, 1],
            omega_est_p[0, 1],
            omega_est[0, 0] - omega_true[0, 0],
            omega_est[1, 1] - omega_true[1, 1],
            omega_est[0, 1] - omega_true[0, 1],
            safe_omega_corr(omega_true),
            safe_omega_corr(omega_est),
            pi2_omega_sum,
            weighted_omega_o_true_raw,
            weighted_omega_o_true_clipped,
            weighted_omega_o_est_raw,
            weighted_omega_o_est_clipped,
            pi2_omega_sum,
            pi_mean,
            propt,
            gmm_propt,
            float(run_gmm),
            float(run_gmm_tissue),
        ],
        dtype=float,
    ).reshape(1, -1)
    np.savetxt(
        os.path.join(out_dir, simulation_path, f"omega_summary_{id_sim}.csv"),
        values,
        delimiter=",",
        header=",".join(columns),
        comments="",
    )


def count_gmm_propt_settings(
    gmm_propt_modes, gmm_propt_mode_scales, gmm_propt_normal_vars
):
    total = 0
    for mode in gmm_propt_modes:
        if mode in ("underestimate", "overestimate"):
            total += len(gmm_propt_mode_scales)
        elif mode == "normal":
            total += len(gmm_propt_normal_vars)
        else:
            total += 1
    return total


def build_causal_corr_matrix(G1, G2):
    cor1 = np.corrcoef(G1, rowvar=False)
    cor2 = np.corrcoef(G2, rowvar=False)
    cor1, cor2 = sanitize_ld_scores(cor1, cor2)
    return np.maximum(np.abs(cor1), np.abs(cor2))


def sample_causal_ids(
    nsnp, num_causal, causal_corr=None, max_abs_cor=None, max_attempts=100
):
    return sample_causal_ids_from_pool(
        np.arange(nsnp), num_causal, causal_corr, max_abs_cor, None, max_attempts
    )


def sample_causal_ids_from_pool(
    pool_ids,
    num_causal,
    causal_corr=None,
    max_abs_cor=None,
    existing_ids=None,
    max_attempts=100,
):
    """Sample causal SNP ids from a candidate pool.

    Parameters
    ----------
    pool_ids : array-like of int
        Candidate SNP ids to sample from.
    num_causal : int
        Number of causal SNPs to sample.
    causal_corr : array-like, optional
        Pairwise SNP correlation matrix indexed by SNP id. Required when
        ``max_abs_cor`` is set.
    max_abs_cor : float, optional
        Maximum allowed absolute correlation between any newly sampled SNP and
        the already selected or existing causal SNPs. If omitted, sampling is
        unconstrained by correlation.
    existing_ids : array-like of int, optional
        Previously selected causal SNP ids that new samples must also be
        compatible with under ``max_abs_cor``.
    max_attempts : int, optional
        Number of randomized greedy sampling attempts before raising an error.
    """
    pool_ids = np.asarray(pool_ids, dtype=int)
    if num_causal <= 0:
        return np.array([], dtype=int)
    if max_abs_cor is None:
        return np.random.choice(pool_ids, num_causal, replace=False)
    if not 0 <= max_abs_cor <= 1:
        raise ValueError("--causal_max_abs_cor must be between 0 and 1.")
    if causal_corr is None:
        raise ValueError("causal_corr is required when --causal_max_abs_cor is set.")

    existing_ids = np.asarray(
        existing_ids if existing_ids is not None else [], dtype=int
    )
    for _ in range(max_attempts):
        candidates = np.random.permutation(pool_ids)
        selected = []
        for snp_id in candidates:
            compatible_with_selected = not selected or np.all(
                causal_corr[snp_id, selected] <= max_abs_cor
            )
            compatible_with_existing = existing_ids.size == 0 or np.all(
                causal_corr[snp_id, existing_ids] <= max_abs_cor
            )
            if compatible_with_selected and compatible_with_existing:
                selected.append(snp_id)
                if len(selected) == num_causal:
                    return np.array(selected, dtype=int)
    raise ValueError(
        "Could not sample enough causal SNPs under --causal_max_abs_cor="
        f"{max_abs_cor}. Increase the threshold or reduce pcausal/nsnp."
    )


def generate_causal_effects_by_overlap(
    nsnp, pcausal, h1sq, h2sq, causal_overlap, causal_corr=None, max_abs_cor=None
):
    """Generate two causal sets whose expected genetic correlation is causal_overlap."""
    num_causal = int(pcausal * nsnp)
    if not 0 <= causal_overlap <= 1:
        raise ValueError("--causal_overlap must be between 0 and 1.")
    if num_causal <= 0:
        return (
            np.zeros(nsnp),
            np.zeros(nsnp),
            np.array([], dtype=int),
            np.array([], dtype=int),
        )

    num_shared = int(round(causal_overlap * num_causal))
    min_shared = max(0, 2 * num_causal - nsnp)
    if num_shared < min_shared:
        raise ValueError(
            "causal_overlap is too small for the requested pcausal/nsnp: "
            f"need at least {min_shared / num_causal:g}."
        )

    num_unique = num_causal - num_shared
    num_union = num_shared + 2 * num_unique

    causal_union = sample_causal_ids(nsnp, num_union, causal_corr, max_abs_cor)
    shared_ids = causal_union[:num_shared]
    pop1_unique_ids = causal_union[num_shared : num_shared + num_unique]
    pop2_unique_ids = causal_union[num_shared + num_unique :]

    causal_ids1 = np.concatenate([shared_ids, pop1_unique_ids])
    causal_ids2 = np.concatenate([shared_ids, pop2_unique_ids])
    beta1 = np.zeros(nsnp)
    beta2 = np.zeros(nsnp)

    sd1 = np.sqrt(h1sq / num_causal)
    sd2 = np.sqrt(h2sq / num_causal)
    if num_shared > 0:
        shared_z = np.random.normal(size=num_shared)
        beta1[shared_ids] = sd1 * shared_z
        beta2[shared_ids] = sd2 * shared_z
    if num_unique > 0:
        beta1[pop1_unique_ids] = np.random.normal(scale=sd1, size=num_unique)
        beta2[pop2_unique_ids] = np.random.normal(scale=sd2, size=num_unique)

    return beta1, beta2, causal_ids1, causal_ids2


def generate_partitioned_causal_effects(
    nsnp,
    pcausal,
    h1sq,
    h2sq,
    gc,
    null_region_prop,
    causal_corr=None,
    max_abs_cor=None,
):
    """A: pop1-null/pop2-causal; B: shared pop1-pop2 causal."""
    if not 0 < null_region_prop < 1:
        raise ValueError("--null_region_prop must be between 0 and 1.")
    split_idx = int(round(null_region_prop * nsnp))
    region_a_ids = np.arange(split_idx)
    region_b_ids = np.arange(split_idx, nsnp)
    num_a_causal = int(pcausal * len(region_a_ids))
    num_b_causal = int(pcausal * len(region_b_ids))
    region_a = np.zeros(nsnp, dtype=bool)
    region_a[region_a_ids] = True
    if num_a_causal <= 0 or num_b_causal <= 0:
        if pcausal == 0 or (num_a_causal == 0 and num_b_causal == 0):
            return (
                np.zeros(nsnp),
                np.zeros(nsnp),
                np.array([], dtype=int),
                np.array([], dtype=int),
                region_a,
            )
        raise ValueError(
            "pcausal must select at least one causal SNP in both region A and B."
        )
    if len(region_a_ids) < num_a_causal or len(region_b_ids) < num_b_causal:
        raise ValueError(
            "Not enough SNPs in region A or B for the requested pcausal/nsnp."
        )

    pop2_a_ids = sample_causal_ids_from_pool(
        region_a_ids, num_a_causal, causal_corr, max_abs_cor
    )
    shared_b_ids = sample_causal_ids_from_pool(
        region_b_ids,
        num_b_causal,
        causal_corr,
        max_abs_cor,
        existing_ids=pop2_a_ids,
    )

    beta1 = np.zeros(nsnp)
    beta2 = np.zeros(nsnp)

    pop2_num_causal = num_a_causal + num_b_causal
    pop1_b_var = h1sq / num_b_causal
    pop2_var = h2sq / pop2_num_causal
    omega_b = np.array(
        [
            [pop1_b_var, np.sqrt(pop1_b_var * pop2_var) * gc],
            [np.sqrt(pop1_b_var * pop2_var) * gc, pop2_var],
        ]
    )
    beta_shared_b = np.random.multivariate_normal(
        mean=np.zeros(2), cov=omega_b, size=num_b_causal
    )
    beta1[shared_b_ids] = beta_shared_b[:, 0]
    beta2[shared_b_ids] = beta_shared_b[:, 1]

    pop2_sd = np.sqrt(pop2_var)
    beta2[pop2_a_ids] = np.random.normal(scale=pop2_sd, size=num_a_causal)

    return beta1, beta2, shared_b_ids, pop2_a_ids, region_a


def generate_data(
    G1,
    G2,
    h1sq,
    h2sq,
    gc,
    n1,
    n2,
    nt,
    nsnp,
    propt,
    pcausal,
    causal_overlap=None,
    causal_corr=None,
    causal_max_abs_cor=None,
    causal_partition_mode="none",
    null_region_prop=0.8,
    tissue_start=None,
    seed_base=None,
    seed_parts=(),
):
    """
    Generate data for simulation.

    Args:
        G1 (np.ndarray): Genotype matrix of population 1 (n1, nsnp).
        G2 (np.ndarray): Genotype matrix of population 2 (n2+nt, nsnp).
        h1sq (float): Heritability of population 1.
        h2sq (float): Heritability of population 2.
        gc (float): Genetic correlation between population 1 and 2.
        n1 (int): Sample size of population 1.
        n2 (int): Sample size of population 2.
        nt (int): Sample size of tissue samples from population 2.
        nsnp (int): Number of SNPs.
        propt (float): Proportion of cell type in tissue.
        pcausal (float): Proportion of causal SNPs.
        causal_overlap (float | None): If set, proportion of causal SNPs shared
            by the two populations. If None, use the original gc-controlled
            shared-causal-set effect generation.
        causal_corr (np.ndarray | None): Pairwise max absolute SNP correlation
            used by causal_max_abs_cor.
        causal_max_abs_cor (float | None): Optional pairwise causal SNP
            correlation threshold.
        causal_partition_mode (str): Optional segmented causal architecture.
        null_region_prop (float): Proportion assigned to pop1-null region A.
        tissue_start (int | None): Row offset for the population 2 tissue panel.
            Defaults to the last nt rows of G2.
        seed_base (int | None): Base seed for component-specific random streams.
            When None, use the caller's current NumPy RNG state.
        seed_parts (tuple): Stable identifiers for the replicate random streams.

    Returns:
        tuple: A tuple containing:
            - Omega_causal (np.ndarray): Causal SNPs covariance matrix.
            - b1_hat (np.ndarray): Effect size estimates for population 1.
            - se1_hat (np.ndarray): Standard errors for population 1.
            - b2_hat (np.ndarray): Effect size estimates for population 2.
            - se2_hat (np.ndarray): Standard errors for population 2.
            - bt_hat (np.ndarray): Effect size estimates for tissue.
            - se_t_hat (np.ndarray): Standard errors for tissue.
            - sig1 (np.ndarray): Boolean array of significant SNPs in population 1.
            - sig2 (np.ndarray): Boolean array of significant SNPs in population 2.
            - sigt (np.ndarray): Boolean array of significant SNPs in tissue.
            - causal_ids (np.ndarray): Indices of population 1 causal SNPs.
            - region_a (np.ndarray): Boolean indicator for pop1-null region A.
            - pi_mean (float): Mean individual cell type proportion in tissue.
    """
    Omega_causal = np.array(
        [[h1sq, np.sqrt(h1sq * h2sq) * gc], [np.sqrt(h1sq * h2sq) * gc, h2sq]]
    )
    # print("Omega:", Omega_causal / nsnp)
    if tissue_start is None:
        tissue_start = tail_panel_start(G2.shape[0], n2, nt, "population 2")
    if n1 > G1.shape[0]:
        raise ValueError(f"n1={n1} exceeds population 1 genotype rows={G1.shape[0]}")
    if n2 > G2.shape[0]:
        raise ValueError(f"n2={n2} exceeds population 2 genotype rows={G2.shape[0]}")
    if tissue_start < n2:
        raise ValueError(
            f"tissue_start={tissue_start} must be >= n2={n2} to keep panels disjoint"
        )
    if tissue_start + nt > G2.shape[0]:
        raise ValueError(
            f"tissue_start + nt = {tissue_start + nt} exceeds "
            f"population 2 genotype rows={G2.shape[0]}"
        )

    G1c = G1[:n1, :]
    X1 = standardize_genotype(G1c, MIN_FLOAT)
    G2c = G2[:n2, :]
    X2 = standardize_genotype(G2c, MIN_FLOAT)
    G2t = G2[tissue_start : tissue_start + nt, :]
    Xt = standardize_genotype(G2t, MIN_FLOAT)
    # cell type data
    num_causal = int(pcausal * nsnp)
    region_a = np.zeros(nsnp, dtype=bool)
    if causal_partition_mode == "pop2_a_shared_b":
        seed_random_component(seed_base, seed_parts, "causal_partition")
        beta1, beta2, causal_ids, _pop2_a_ids, region_a = (
            generate_partitioned_causal_effects(
                nsnp,
                pcausal,
                h1sq,
                h2sq,
                gc,
                null_region_prop,
                causal_corr,
                causal_max_abs_cor,
            )
        )
    elif causal_overlap is None:
        seed_random_component(seed_base, seed_parts, "causal_ids")
        causal_ids = sample_causal_ids(
            nsnp, num_causal, causal_corr, causal_max_abs_cor
        )
        beta1 = np.zeros(nsnp)
        beta2 = np.zeros(nsnp)
        if num_causal > 0:
            seed_random_component(seed_base, seed_parts, "causal_effects")
            beta_causal = np.random.multivariate_normal(
                mean=np.zeros(2), cov=Omega_causal / num_causal, size=num_causal
            )  # (M, <c11, c12>)
            beta1[causal_ids] = beta_causal[:, 0]
            beta2[causal_ids] = beta_causal[:, 1]
    else:
        seed_random_component(seed_base, seed_parts, "causal_overlap")
        beta1, beta2, causal_ids, _causal_ids2 = generate_causal_effects_by_overlap(
            nsnp,
            pcausal,
            h1sq,
            h2sq,
            causal_overlap,
            causal_corr,
            causal_max_abs_cor,
        )
    seed_random_component(seed_base, seed_parts, "pop1_noise")
    y1 = X1 @ beta1.T + np.sqrt(1 - h1sq) * np.random.randn(n1)
    seed_random_component(seed_base, seed_parts, "pop2_noise")
    y2 = X2 @ beta2.T + np.sqrt(1 - h2sq) * np.random.randn(n2)

    # tissue data
    delta = 5  # control the variance of pi
    seed_random_component(seed_base, seed_parts, "cell_proportion")
    if tissue_start + nt == G2.shape[0]:
        pi_all = np.random.beta(
            (propt + MIN_FLOAT) * delta,
            (1 - propt + MIN_FLOAT) * delta,
            G2.shape[0],
        )
        pi_ind = pi_all[tissue_start : tissue_start + nt]
    else:
        pi_ind = np.random.beta(
            (propt + MIN_FLOAT) * delta,
            (1 - propt + MIN_FLOAT) * delta,
            nt,
        )
    pi_mean = np.mean(pi_ind)
    ## define unknown cell type
    beta_unknown = np.zeros(nsnp)
    num_unknown_celltype = 1
    for celltype_id in range(num_unknown_celltype):
        num_unknown_causal = int(pcausal * nsnp)
        if num_unknown_causal <= 0:
            continue
        seed_random_component(seed_base, seed_parts, f"unknown_ids_{celltype_id}")
        causal_unknown_id = sample_causal_ids(
            nsnp, num_unknown_causal, causal_corr, causal_max_abs_cor
        )
        # causal_unknown_id = causal_ids  #! share causal SNPs
        seed_random_component(seed_base, seed_parts, f"unknown_effects_{celltype_id}")
        beta_causal_unknown = np.random.normal(
            loc=0,
            scale=unknown_cell_effect_scale(
                h2sq, num_unknown_causal, num_unknown_celltype
            ),
            size=num_unknown_causal,
        )
        # beta_causal_unknown = (
        #     beta_causal[:, 0] / 3 + beta_causal[:, 1] / 3 + beta_causal_unknown / 3
        # ) #! make unknown cell type correlated with known cell types
        beta_unknown[causal_unknown_id] += beta_causal_unknown
        # beta_unknown[causal_unknown_id] += beta_causal_unknown / 2
        # beta_unknown[causal_ids] += (
        #     beta_causal[:, 1] / 2
        # )  #! share effect with known cell type
    seed_random_component(seed_base, seed_parts, "tissue_noise")
    if tissue_start + nt == G2.shape[0]:
        tissue_noise = np.random.randn(G2.shape[0])[tissue_start : tissue_start + nt]
    else:
        tissue_noise = np.random.randn(nt)
    yt = (
        pi_ind * (Xt @ beta2.T)
        + (1 - pi_ind) * (Xt @ beta_unknown.T)
        + np.sqrt(
            np.maximum(1 - (pi_ind**2 + (1 - pi_ind) ** 2) * h2sq, MIN_FLOAT)
        )
        * tissue_noise
    )

    # sumstats
    b1_hat, se1_hat = calculate_sumstats(X1, y1, n1)
    b2_hat, se2_hat = calculate_sumstats(X2, y2, n2)
    bt_hat, se_t_hat = calculate_sumstats(Xt, yt, nt)
    z1 = b1_hat / (se1_hat + MIN_FLOAT)
    z2 = b2_hat / (se2_hat + MIN_FLOAT)
    zt = bt_hat / (se_t_hat + MIN_FLOAT)
    pval1 = z2p(z1)
    pval2 = z2p(z2)
    pvalt = z2p(zt)
    sig1 = pval1 < P_VAL_THRED
    sig2 = pval2 < P_VAL_THRED
    sigt = pvalt < P_VAL_THRED
    # breakpoint()
    # cal cov between target/auxiliary cell type and unknown cell type in tissue
    Omega = np.cov(
        np.stack([beta1, beta2, beta2 * pi_mean + beta_unknown * (1 - pi_mean)])
    )
    return (
        Omega,
        b1_hat,
        se1_hat,
        b2_hat,
        se2_hat,
        bt_hat,
        se_t_hat,
        sig1,
        sig2,
        sigt,
        causal_ids,
        region_a,
        pi_mean,
    )


def simulation(
    runname,
    G1,
    G2,
    ld1,
    ld2,
    ldx,
    h1sq,
    h2sq,
    gc,
    n1,
    n2,
    nt,
    nsnp,
    propt,
    gmm_propt_mode,
    gmm_propt_normal_var,
    gmm_propt_mode_scale,
    pcausal,
    causal_overlap,
    causal_max_abs_cor,
    causal_corr,
    causal_partition_mode,
    null_region_prop,
    out_dir,
    true_omega,
    id_sim,
    tissue_start=None,
    seed_base=None,
    seed_parts=(),
):
    """
    Run simulation for GMM, save results to out_dir
    :param
    runname: power or alpha
    G1, G2: genotype matrix of population 1 and 2 (n1, nsnp), (n2+nt, nsnp)
    ld1, ld2, ldx: LD between snp j and the rest of the snps of population 1 and 2
    h1sq, h2sq: heritability of population 1 and 2
    gc: genetic correlation between population 1 and 2
    n1, n2, nt: sample size of population 1, 2 and tissue
    nsnp: number of SNPs
    propt: proportion of cell type in tissue
    gmm_propt_mode: perturbation mode for the proportion supplied to GMM tissue
    gmm_propt_normal_var: variance for the normal perturbation mode
    gmm_propt_mode_scale: additive scale used by underestimate/overestimate modes
    pcausal: proportion of causal SNPs
    causal_overlap: optional proportion of shared causal SNPs
    causal_max_abs_cor: optional maximum absolute correlation among causal SNPs
    causal_corr: pairwise SNP correlation matrix used by causal_max_abs_cor
    causal_partition_mode: optional segmented causal architecture
    null_region_prop: proportion assigned to pop1-null region A
    out_dir: output directory
    true_omega: true covariance matrix
    id_sim: simulation id
    tissue_start: row offset for the population 2 tissue panel
    seed_base, seed_parts: common-random-number seed controls
    :return
    None
    """
    gmm_propt_path = ""
    if gmm_propt_mode != "exact":
        gmm_propt_path = f"_gmmproptmode_{gmm_propt_mode}"
        if gmm_propt_mode in ("underestimate", "overestimate"):
            gmm_propt_path += f"_gmmproptmodescale_{gmm_propt_mode_scale}"
        elif gmm_propt_mode == "normal":
            gmm_propt_path += f"_gmmproptnormalvar_{gmm_propt_normal_var:g}"
    causal_overlap_path = ""
    if causal_overlap is not None:
        causal_overlap_path = f"_causaloverlap_{causal_overlap}"
    causal_corr_path = ""
    if causal_max_abs_cor is not None:
        causal_corr_path = f"_causalmaxabscor_{causal_max_abs_cor:g}"
    causal_partition_path = ""
    if causal_partition_mode != "none":
        causal_partition_path = (
            f"_causalpartition_{causal_partition_mode}"
            f"_nullregionprop_{null_region_prop:g}"
        )
    simulation_path = f"{runname}/h1sq_{h1sq}_h2sq_{h2sq}_gc_{gc}_n1_{n1}_n2_{n2}_nt_{nt}_nsnp_{nsnp}_propt_{propt}{gmm_propt_path}_pcausal_{pcausal}{causal_overlap_path}{causal_corr_path}{causal_partition_path}_omega_{true_omega}"
    os.makedirs(os.path.join(out_dir, simulation_path), exist_ok=True)
    simulation_name = f"{simulation_path}/simulation_{id_sim}"
    (
        OmegaCB,
        b1_hat,
        se1_hat,
        b2_hat,
        se2_hat,
        bt_hat,
        se_t_hat,
        sig1,
        sig2,
        sigt,
        causal_ids,
        region_a,
        pi_mean,
    ) = generate_data(
        G1,
        G2,
        h1sq,
        h2sq,
        gc,
        n1,
        n2,
        nt,
        nsnp,
        propt,
        pcausal,
        causal_overlap,
        causal_corr,
        causal_max_abs_cor,
        causal_partition_mode,
        null_region_prop,
        tissue_start=tissue_start,
        seed_base=seed_base,
        seed_parts=seed_parts,
    )
    seed_random_component(seed_base, seed_parts, "gmm_propt")
    gmm_propt = perturb_gmm_propt(
        propt, gmm_propt_mode, gmm_propt_normal_var, gmm_propt_mode_scale
    )
    weighted_omega_o_true_raw = calculate_true_weighted_omega_o(OmegaCB, propt)
    weighted_omega_o_true_clipped = clip_weighted_omega_o(weighted_omega_o_true_raw)
    weighted_omega_o_est_raw = np.nan
    weighted_omega_o_est_clipped = np.nan
    Omega = np.zeros((2, 2))
    omega_est = np.full((2, 2), np.nan)
    omega_est_se = np.full((2, 2), np.nan)
    omega_est_p = np.full((2, 2), np.nan)
    if not true_omega:  # estimate omega
        run_gmm = False  # default no gmm
        run_gmm_tissue = False  # default no gmm tissue
        pi2_omega_sum = 0.0  # weighted Sigma_o contribution for non-target cell types
        Omega, Omega_se = Run_Cross_LDSC(
            b1_hat / (se1_hat + MIN_FLOAT),
            n1,
            ld1,
            b2_hat / (se2_hat + MIN_FLOAT),
            n2,
            ld2,
            ldx,
            np.array([1, 1, 0]),
        )
        Omega_p = z2p(Omega / (Omega_se + MIN_FLOAT))
        omega_est = Omega.copy()
        omega_est_se = Omega_se.copy()
        omega_est_p = Omega_p.copy()
        p_thred = 0.10
        if np.all(Omega_p < p_thred):
            #! if np.all(Omega_p < P_VAL_THRED):
            run_gmm = True
            aux_Omega_matrix, aux_Omega_matrix_se = Run_Cross_LDSC(
                b2_hat / (se2_hat + MIN_FLOAT),
                n2,
                ld2,
                bt_hat / (se_t_hat + MIN_FLOAT),
                nt,
                ldx,
                ldx,
                np.array([1.0, 1.0, 0.0]),
            )
            if np.all(
                z2p(aux_Omega_matrix / (aux_Omega_matrix_se + MIN_FLOAT))
                < p_thred
                #! z2p(aux_Omega_matrix / aux_Omega_matrix_se) < P_VAL_THRED
            ):
                run_gmm_tissue = True
                # \text{LDSC}(z_t, z_t) - \pi_c^2 \omega_2 - (2\pi_c + \frac{\sum_{i \neq j} \pi_i \pi_j}{1-\pi_c}) (\text{LDSC}(z_2, z_t)-\pi_c \omega_2)
                weighted_omega_o_est_raw = calculate_estimated_weighted_omega_o(
                    aux_Omega_matrix, Omega, gmm_propt
                )
                weighted_omega_o_est_clipped = clip_weighted_omega_o(
                    weighted_omega_o_est_raw
                )
                pi2_omega_sum = weighted_omega_o_est_clipped

    else:  # true Omega
        Omega = OmegaCB[:2, :2]
        omega_est = Omega.copy()
        run_gmm = should_run_true_omega_gmm(h1sq, h2sq, gc, Omega)
        run_gmm_tissue = run_gmm
        pi2_omega_sum = weighted_omega_o_true_raw
    save_omega_summary(
        out_dir,
        simulation_path,
        id_sim,
        OmegaCB[:2, :2],
        omega_est,
        omega_est_se,
        omega_est_p,
        pi2_omega_sum,
        weighted_omega_o_true_raw,
        weighted_omega_o_true_clipped,
        weighted_omega_o_est_raw,
        weighted_omega_o_est_clipped,
        pi_mean,
        propt,
        gmm_propt,
        run_gmm,
        run_gmm_tissue,
    )
    (
        pop1_beta,
        pop1_se,
        pop2_beta,
        pop2_se,
        meta_beta,
        meta_se,
        meta_tissue_beta,
        meta_tissue_se,
    ) = run_gmm_meta_kernel(
        nsnp,
        run_gmm,
        run_gmm_tissue,
        Omega,
        pi2_omega_sum,
        gmm_propt,
        b1_hat,
        se1_hat,
        ld1,
        b2_hat,
        se2_hat,
        ld2,
        ldx,
        bt_hat,
        se_t_hat,
    )

    pop1_z = pop1_beta / (pop1_se + MIN_FLOAT)
    pop2_z = pop2_beta / (pop2_se + MIN_FLOAT)
    # save results
    all_results_columns = [
        "causal",
        "sign1",
        "sign2",
        "sign_t",
        "z1_sumstat",
        "z1_cross",
        "z1_tissue",
        "z2_sumstat",
        "z2_cross",
        "z2_tissue",
        "zt_sumstat",
        "z_meta",
        "z_metatissue",
        "gmm_propt",
        "region_a",
    ]
    all_results = np.zeros((nsnp, len(all_results_columns)))
    all_results[causal_ids, 0] = 1
    all_results[:, 1] = sig1
    all_results[:, 2] = sig2
    all_results[:, 3] = sigt
    all_results[:, 4:7] = pop1_z
    all_results[:, 7:10] = pop2_z
    all_results[:, 10] = bt_hat / (se_t_hat + MIN_FLOAT)
    all_results[:, 11] = meta_beta / (meta_se + MIN_FLOAT)
    all_results[:, 12] = meta_tissue_beta / (meta_tissue_se + MIN_FLOAT)
    all_results[:, 13] = gmm_propt
    all_results[:, 14] = region_a
    np.savetxt(
        os.path.join(
            out_dir,
            simulation_name + ".csv",
        ),
        all_results,
        delimiter=",",
        header=",".join(all_results_columns),
        comments="",
    )


def main():
    args = parse_args()
    ## data
    G1 = get_genotype(args.pop1_geno, args.nsnp)
    G2 = get_genotype(args.pop2_geno, args.nsnp)
    ld1, ld2, ldx = cal_ld(G1, G2)
    ld1, ld2, ldx = sanitize_ld_scores(ld1, ld2, ldx)
    causal_corr = None
    if args.causal_max_abs_cor is not None:
        print(
            "Calculating causal SNP correlation matrix with max abs cor "
            f"{args.causal_max_abs_cor}"
        )
        causal_corr = build_causal_corr_matrix(G1, G2)
    ## parse
    h1sq = args.h1sq if isinstance(args.h1sq, list) else [args.h1sq]
    h2sq = args.h2sq if isinstance(args.h2sq, list) else [args.h2sq]
    gc = args.gc if isinstance(args.gc, list) else [args.gc]
    n1 = args.n1 if isinstance(args.n1, list) else [args.n1]
    n2 = args.n2 if isinstance(args.n2, list) else [args.n2]
    nt = args.nt if isinstance(args.nt, list) else [args.nt]
    propt = flatten_float_seq(args.propt)
    gmm_propt_mode = (
        args.gmm_propt_mode
        if isinstance(args.gmm_propt_mode, list)
        else [args.gmm_propt_mode]
    )
    gmm_propt_mode_scale = flatten_float_seq(args.gmm_propt_mode_scale)
    gmm_propt_normal_var = flatten_float_seq(args.gmm_propt_normal_var)
    validate_unit_interval("--h1sq", h1sq)
    validate_unit_interval("--h2sq", h2sq)
    validate_unit_interval("--propt", propt)
    validate_nonnegative("--gmm_propt_mode_scale", gmm_propt_mode_scale)
    validate_nonnegative("--gmm_propt_normal_var", gmm_propt_normal_var)
    pcausal = args.pcausal if isinstance(args.pcausal, list) else [args.pcausal]
    validate_unit_interval("--pcausal", pcausal)
    causal_overlap = (
        args.causal_overlap
        if isinstance(args.causal_overlap, list)
        else [args.causal_overlap]
    )
    null_region_prop = (
        args.null_region_prop
        if isinstance(args.null_region_prop, list)
        else [args.null_region_prop]
    )
    nsnp = args.nsnp
    trueOmega = not args.estimate_omega

    ## process indicators
    gmm_propt_setting_count = count_gmm_propt_settings(
        gmm_propt_mode, gmm_propt_mode_scale, gmm_propt_normal_var
    )
    total_combinations = (
        len(h1sq)
        * len(h2sq)
        * len(gc)
        * len(n1)
        * len(n2)
        * len(nt)
        * len(propt)
        * gmm_propt_setting_count
        * len(pcausal)
        * len(causal_overlap)
        * len(null_region_prop)
    )
    have_run = 0

    start_time = time.time()
    base_seed = int(args.seed) if args.seed is not None else int(start_time)
    max_tissue_end = max(n2) + max(nt)
    if max_tissue_end > G2.shape[0]:
        raise ValueError(
            f"required n2/nt sample rows = {max_tissue_end} exceeds population 2 "
            f"genotype rows={G2.shape[0]}"
        )
    print("Simulation start at ", time.ctime())
    print("Simulation base seed: ", base_seed)
    for h1sqi in h1sq:
        for h2sqj in h2sq:
            for gck in gc:
                for n1l in n1:
                    for n2m in n2:
                        for ntn in nt:
                            for proptp in propt:
                                for gmm_propt_modeq in gmm_propt_mode:
                                    for (
                                        _sub_idx,
                                        gmm_propt_scaleq,
                                        gmm_propt_nvq,
                                    ) in iter_gmm_propt_subsettings(
                                        gmm_propt_modeq,
                                        gmm_propt_mode_scale,
                                        gmm_propt_normal_var,
                                    ):
                                        for pcausalr in pcausal:
                                            for causal_overlapo in causal_overlap:
                                                for null_region_propa in null_region_prop:
                                                    for s in range(args.nrep):
                                                        seed_parts = (s,)
                                                        simulation(
                                                            args.runname,
                                                            G1,
                                                            G2,
                                                            ld1,
                                                            ld2,
                                                            ldx,
                                                            h1sqi,
                                                            h2sqj,
                                                            gck,
                                                            n1l,
                                                            n2m,
                                                            ntn,
                                                            nsnp,
                                                            proptp,
                                                            gmm_propt_modeq,
                                                            gmm_propt_nvq,
                                                            gmm_propt_scaleq,
                                                            pcausalr,
                                                            causal_overlapo,
                                                            args.causal_max_abs_cor,
                                                            causal_corr,
                                                            args.causal_partition_mode,
                                                            null_region_propa,
                                                            args.out_dir,
                                                            trueOmega,
                                                            s,
                                                            tissue_start=None,
                                                            seed_base=base_seed,
                                                            seed_parts=seed_parts,
                                                        )
                                                    have_run += 1
                                                    if have_run % 5 == 0:
                                                        print(
                                                            f"Simulation {have_run}/{total_combinations} done, {time.ctime()}"
                                                        )
    print("Simulation end at ", time.ctime())
    print(
        f"Total time: {time.time() - start_time:.2f} s for {total_combinations} settings"
    )


if __name__ == "__main__":
    main()
