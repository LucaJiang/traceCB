"""traceCB^2 simulation driver with tissue panels in both populations.

This script extends ``simulation.py`` by adding a population 1 bulk/tissue
summary-statistic panel. The main comparison is whether combining pop1 sc,
pop1 tissue, pop2 sc, and pop2 tissue summaries improves target-population
eQTL discovery compared with traceCB and conventional meta-analysis.

Run the curated grid with:

```
bash src/simulation/others/run_tracecb2.sh
```

The paired visualization script is ``visual_tracecb2.py``. It reads replicate
CSV files from ``bench/result/<runname>/``, writes ``result_df_tracecb2*.csv``,
and saves figures under ``bench/result/img``.
"""

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
from numba import njit, prange

SIMULATION_DIR = Path(__file__).resolve().parents[1]
ROOT_DIR = SIMULATION_DIR.parents[1]
SRC_DIR = ROOT_DIR / "src"
for path in (SRC_DIR, SIMULATION_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from traceCB.gmm import GMM, GMMtissue, GMMtissueBoth
from traceCB.ldsc import Run_Cross_LDSC
from traceCB.utils import MIN_HERITABILITY, z2p
from simulation import (
    MAX_CORR,
    MIN_FLOAT,
    P_VAL_THRED,
    cal_ld,
    calculate_sumstats,
    get_genotype,
    re2_meta,
    should_run_true_omega_gmm,
)
from simulation_utils import (
    sanitize_ld_scores,
    seed_random_component,
    standardize_genotype,
    tail_panel_start,
    unknown_cell_effect_scale,
    validate_unit_interval,
)


@njit(nogil=True, parallel=True, cache=True)
def run_tracecb2_inference_kernel(
    nsnp,
    run_gmm,
    run_tracecb,
    run_tracecb2,
    omega,
    pi2_omega_sum,
    pi2_omega_sum1,
    pi2_omega_sum2,
    pi12_omega_sum12,
    gmm_propt,
    b1_hat,
    se1_hat,
    bt1_hat,
    se_t1_hat,
    b2_hat,
    se2_hat,
    bt2_hat,
    se_t2_hat,
    ld1,
    ld2,
    ldx,
):
    pop1_beta = np.empty((nsnp, 4))
    pop2_beta = np.empty((nsnp, 4))
    pop1_se = np.empty((nsnp, 4))
    pop2_se = np.empty((nsnp, 4))
    meta_beta = np.empty(nsnp)
    meta_se = np.empty(nsnp)
    meta_tissue_beta = np.empty(nsnp)
    meta_tissue_se = np.empty(nsnp)
    eye2 = np.eye(2)
    eye3 = np.eye(3)
    eye4 = np.eye(4)

    for j in prange(nsnp):
        pop1_beta[j, 0] = b1_hat[j]
        pop1_se[j, 0] = se1_hat[j]
        pop2_beta[j, 0] = b2_hat[j]
        pop2_se[j, 0] = se2_hat[j]

        if run_gmm:
            (
                pop1_beta[j, 1],
                pop1_se[j, 1],
                pop2_beta[j, 1],
                pop2_se[j, 1],
            ) = GMM(
                omega,
                eye2,
                b1_hat[j],
                se1_hat[j],
                ld1[j],
                b2_hat[j],
                se2_hat[j],
                ld2[j],
                ldx[j],
            )
        else:
            pop1_beta[j, 1] = b1_hat[j]
            pop1_se[j, 1] = se1_hat[j]
            pop2_beta[j, 1] = b2_hat[j]
            pop2_se[j, 1] = se2_hat[j]

        if run_tracecb:
            (
                pop1_beta[j, 2],
                pop1_se[j, 2],
                pop2_beta[j, 2],
                pop2_se[j, 2],
            ) = GMMtissue(
                omega,
                eye3,
                b1_hat[j],
                se1_hat[j],
                ld1[j],
                b2_hat[j],
                se2_hat[j],
                ld2[j],
                ldx[j],
                bt2_hat[j],
                se_t2_hat[j],
                pi2_omega_sum,
                gmm_propt,
            )
        else:
            pop1_beta[j, 2] = pop1_beta[j, 1]
            pop1_se[j, 2] = pop1_se[j, 1]
            pop2_beta[j, 2] = pop2_beta[j, 1]
            pop2_se[j, 2] = pop2_se[j, 1]

        if run_tracecb2:
            (
                pop1_beta[j, 3],
                pop1_se[j, 3],
                pop2_beta[j, 3],
                pop2_se[j, 3],
            ) = GMMtissueBoth(
                omega,
                eye4,
                b1_hat[j],
                se1_hat[j],
                ld1[j],
                bt1_hat[j],
                se_t1_hat[j],
                b2_hat[j],
                se2_hat[j],
                ld2[j],
                ldx[j],
                bt2_hat[j],
                se_t2_hat[j],
                pi2_omega_sum1,
                pi2_omega_sum2,
                pi12_omega_sum12,
                gmm_propt,
                gmm_propt,
            )
        else:
            pop1_beta[j, 3] = pop1_beta[j, 2]
            pop1_se[j, 3] = pop1_se[j, 2]
            pop2_beta[j, 3] = pop2_beta[j, 2]
            pop2_se[j, 3] = pop2_se[j, 2]

        meta_beta[j], meta_se[j] = re2_meta(
            np.array([b1_hat[j], b2_hat[j]]),
            np.array([se1_hat[j], se2_hat[j]]),
        )
        meta_tissue_beta[j], meta_tissue_se[j] = re2_meta(
            np.array([b1_hat[j], bt1_hat[j], b2_hat[j], bt2_hat[j]]),
            np.array([se1_hat[j], se_t1_hat[j], se2_hat[j], se_t2_hat[j]]),
        )

    return (
        pop1_beta,
        pop1_se,
        pop2_beta,
        pop2_se,
        meta_beta,
        meta_se,
        meta_tissue_beta,
        meta_tissue_se,
    )


def parse_args():
    parser = argparse.ArgumentParser(description="Run traceCB^2 simulation.")
    parser.add_argument(
        "--pop1_geno", default="data/simulation/EAS_n5000_chr22_loci29.npy"
    )
    parser.add_argument(
        "--pop2_geno", default="data/simulation/EUR_n20000_chr22_loci29.npy"
    )
    parser.add_argument("--runname", default="nt1_nt2_propt")
    parser.add_argument("--h1sq", default=[0.1], type=float, nargs="+")
    parser.add_argument("--h2sq", default=[0.1], type=float, nargs="+")
    parser.add_argument("--gc", default=[0.2], type=float, nargs="+")
    parser.add_argument("--n1", default=[100], type=int, nargs="+")
    parser.add_argument("--n2", default=[200], type=int, nargs="+")
    parser.add_argument("--nt", default=None, type=int, nargs="+")
    parser.add_argument("--nt1", default=None, type=int, nargs="+")
    parser.add_argument("--nt2", default=None, type=int, nargs="+")
    parser.add_argument("--nsnp", default=1000, type=int)
    parser.add_argument("--propt", default=[0.1], type=float, nargs="+")
    parser.add_argument("--pcausal", default=[0.1], type=float, nargs="+")
    parser.add_argument(
        "--estimate_omega",
        action="store_true",
        help="Estimate omega from data, otherwise use true simulated omega.",
    )
    parser.add_argument("--out_dir", default="bench/result")
    parser.add_argument("--nrep", default=100, type=int)
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


def build_result_path(
    runname,
    h1sq,
    h2sq,
    gc,
    n1,
    n2,
    nt1,
    nt2,
    nsnp,
    propt,
    pcausal,
    true_omega,
):
    return (
        f"{runname}/h1sq_{h1sq}_h2sq_{h2sq}_gc_{gc}_n1_{n1}_n2_{n2}"
        f"_nt1_{nt1}_nt2_{nt2}_nsnp_{nsnp}_propt_{propt}_pcausal_{pcausal}"
        f"_omega_{true_omega}"
    )


def sample_causal_ids(nsnp, num_causal):
    if num_causal <= 0:
        return np.array([], dtype=int)
    return np.random.choice(np.arange(nsnp), num_causal, replace=False)


def sample_unknown_effects(nsnp, pcausal, hsq):
    beta_unknown = np.zeros(nsnp)
    num_causal = int(pcausal * nsnp)
    if num_causal <= 0:
        return beta_unknown
    causal_unknown_id = sample_causal_ids(nsnp, num_causal)
    beta_unknown[causal_unknown_id] = np.random.normal(
        loc=0,
        scale=unknown_cell_effect_scale(hsq, num_causal, 1),
        size=num_causal,
    )
    return beta_unknown


def generate_data_tracecb2(
    G1,
    G2,
    h1sq,
    h2sq,
    gc,
    n1,
    n2,
    nt1,
    nt2,
    nsnp,
    propt,
    pcausal,
    pop2_tissue_start=None,
    seed_base=None,
    seed_parts=(),
):
    omega_causal = np.array(
        [[h1sq, np.sqrt(h1sq * h2sq) * gc], [np.sqrt(h1sq * h2sq) * gc, h2sq]]
    )
    if pop2_tissue_start is None:
        pop2_tissue_start = tail_panel_start(
            G2.shape[0], n2, nt2, "population 2", allow_overlap=True
        )
    pop1_tissue_start = tail_panel_start(
        G1.shape[0], n1, nt1, "population 1", allow_overlap=True
    )
    if n1 > G1.shape[0]:
        raise ValueError(f"n1={n1} exceeds population 1 genotype rows={G1.shape[0]}")
    if n2 > G2.shape[0]:
        raise ValueError(f"n2={n2} exceeds population 2 genotype rows={G2.shape[0]}")
    if nt1 > G1.shape[0]:
        raise ValueError(f"nt1={nt1} exceeds population 1 genotype rows={G1.shape[0]}")
    if pop2_tissue_start < 0:
        raise ValueError("pop2_tissue_start must be nonnegative")
    if pop2_tissue_start + nt2 > G2.shape[0]:
        raise ValueError(
            f"pop2_tissue_start + nt2 = {pop2_tissue_start + nt2} exceeds "
            f"population 2 genotype rows={G2.shape[0]}"
        )

    X1 = standardize_genotype(G1[:n1, :], MIN_FLOAT)
    X2 = standardize_genotype(G2[:n2, :], MIN_FLOAT)
    Xt1 = standardize_genotype(
        G1[pop1_tissue_start : pop1_tissue_start + nt1, :], MIN_FLOAT
    )
    Xt2 = standardize_genotype(G2[pop2_tissue_start : pop2_tissue_start + nt2, :], MIN_FLOAT)

    num_causal = int(pcausal * nsnp)
    seed_random_component(seed_base, seed_parts, "causal_ids")
    causal_ids = sample_causal_ids(nsnp, num_causal)
    beta1 = np.zeros(nsnp)
    beta2 = np.zeros(nsnp)
    if num_causal > 0:
        seed_random_component(seed_base, seed_parts, "causal_effects")
        beta_causal = np.random.multivariate_normal(
            mean=np.zeros(2), cov=omega_causal / num_causal, size=num_causal
        )
        beta1[causal_ids] = beta_causal[:, 0]
        beta2[causal_ids] = beta_causal[:, 1]

    seed_random_component(seed_base, seed_parts, "pop1_noise")
    y1 = X1 @ beta1.T + np.sqrt(1 - h1sq) * np.random.randn(n1)
    seed_random_component(seed_base, seed_parts, "pop2_noise")
    y2 = X2 @ beta2.T + np.sqrt(1 - h2sq) * np.random.randn(n2)

    delta = 5
    seed_random_component(seed_base, seed_parts, "cell_proportion_pop1")
    pi_all1 = np.random.beta(
        (propt + MIN_FLOAT) * delta,
        (1 - propt + MIN_FLOAT) * delta,
        G1.shape[0],
    )
    pi_ind1 = pi_all1[pop1_tissue_start : pop1_tissue_start + nt1]
    seed_random_component(seed_base, seed_parts, "cell_proportion_pop2")
    pi_all2 = np.random.beta(
        (propt + MIN_FLOAT) * delta,
        (1 - propt + MIN_FLOAT) * delta,
        G2.shape[0],
    )
    pi_ind2 = pi_all2[pop2_tissue_start : pop2_tissue_start + nt2]
    pi_mean1 = float(np.mean(pi_ind1))
    pi_mean2 = float(np.mean(pi_ind2))
    seed_random_component(seed_base, seed_parts, "unknown_pop1")
    beta_unknown1 = sample_unknown_effects(nsnp, pcausal, h1sq)
    seed_random_component(seed_base, seed_parts, "unknown_pop2")
    beta_unknown2 = sample_unknown_effects(nsnp, pcausal, h2sq)
    seed_random_component(seed_base, seed_parts, "tissue_noise_pop1")
    tissue_noise1 = np.random.randn(G1.shape[0])[
        pop1_tissue_start : pop1_tissue_start + nt1
    ]
    yt1 = (
        pi_ind1 * (Xt1 @ beta1.T)
        + (1 - pi_ind1) * (Xt1 @ beta_unknown1.T)
        + np.sqrt(
            np.maximum(1 - (pi_ind1**2 + (1 - pi_ind1) ** 2) * h1sq, MIN_FLOAT)
        )
        * tissue_noise1
    )
    seed_random_component(seed_base, seed_parts, "tissue_noise_pop2")
    tissue_noise2 = np.random.randn(G2.shape[0])[
        pop2_tissue_start : pop2_tissue_start + nt2
    ]
    yt2 = (
        pi_ind2 * (Xt2 @ beta2.T)
        + (1 - pi_ind2) * (Xt2 @ beta_unknown2.T)
        + np.sqrt(
            np.maximum(1 - (pi_ind2**2 + (1 - pi_ind2) ** 2) * h2sq, MIN_FLOAT)
        )
        * tissue_noise2
    )

    b1_hat, se1_hat = calculate_sumstats(X1, y1, n1)
    b2_hat, se2_hat = calculate_sumstats(X2, y2, n2)
    bt1_hat, se_t1_hat = calculate_sumstats(Xt1, yt1, nt1)
    bt2_hat, se_t2_hat = calculate_sumstats(Xt2, yt2, nt2)

    omega_cb = np.cov(
        np.stack([beta1, beta2, beta2 * pi_mean2 + beta_unknown2 * (1 - pi_mean2)])
    )
    omega_cb2 = np.cov(
        np.stack(
            [
                beta1,
                beta1 * pi_mean1 + beta_unknown1 * (1 - pi_mean1),
                beta2,
                beta2 * pi_mean2 + beta_unknown2 * (1 - pi_mean2),
            ]
        )
    )
    return (
        omega_cb,
        omega_cb2,
        b1_hat,
        se1_hat,
        b2_hat,
        se2_hat,
        bt1_hat,
        se_t1_hat,
        bt2_hat,
        se_t2_hat,
        z2p(b1_hat / (se1_hat + MIN_FLOAT)) < P_VAL_THRED,
        z2p(b2_hat / (se2_hat + MIN_FLOAT)) < P_VAL_THRED,
        z2p(bt1_hat / (se_t1_hat + MIN_FLOAT)) < P_VAL_THRED,
        z2p(bt2_hat / (se_t2_hat + MIN_FLOAT)) < P_VAL_THRED,
        causal_ids,
        pi_mean1,
        pi_mean2,
    )


def estimate_tracecb2_moments(
    b1_hat,
    se1_hat,
    b2_hat,
    se2_hat,
    bt1_hat,
    se_t1_hat,
    bt2_hat,
    se_t2_hat,
    n1,
    n2,
    nt1,
    nt2,
    ld1,
    ld2,
    ldx,
    gmm_propt,
):
    z1 = b1_hat / (se1_hat + MIN_FLOAT)
    z2 = b2_hat / (se2_hat + MIN_FLOAT)
    zt1 = bt1_hat / (se_t1_hat + MIN_FLOAT)
    zt2 = bt2_hat / (se_t2_hat + MIN_FLOAT)
    intercept = np.array([1.0, 1.0, 0.0])
    omega, omega_se = Run_Cross_LDSC(z1, n1, ld1, z2, n2, ld2, ldx, intercept)
    omega_p = z2p(omega / (omega_se + MIN_FLOAT))

    run_gmm = False
    run_tracecb = False
    run_tracecb2 = False
    pi2_omega_sum = 0.0
    pi2_omega_sum1 = MIN_HERITABILITY
    pi2_omega_sum2 = MIN_HERITABILITY
    pi12_omega_sum12 = 0.0

    p_threshold = 0.10
    if np.all(omega_p < p_threshold):
        run_gmm = True
        aux1, aux1_se = Run_Cross_LDSC(
            z1, n1, ld1, zt1, nt1, ld1, ld1, intercept
        )
        aux2, aux2_se = Run_Cross_LDSC(
            z2, n2, ld2, zt2, nt2, ld2, ld2, intercept
        )
        aux1_p = z2p(aux1 / (aux1_se + MIN_FLOAT))
        aux2_p = z2p(aux2 / (aux2_se + MIN_FLOAT))

        if np.all(aux2_p < p_threshold):
            run_tracecb = True
            pi2_omega_sum = (
                aux2[1, 1]
                - gmm_propt**2 * omega[1, 1]
                - 2 * gmm_propt * (aux2[0, 1] - gmm_propt * omega[1, 1])
            )
            pi2_omega_sum = max(float(pi2_omega_sum), MIN_HERITABILITY)

        if np.all(aux1_p < p_threshold) and np.all(aux2_p < p_threshold):
            run_tracecb2 = True
            tissue_cross, _tissue_cross_se = Run_Cross_LDSC(
                zt1, nt1, ld1, zt2, nt2, ld2, ldx, intercept
            )
            pi2_omega_sum1 = max(
                float(aux1[1, 1] - gmm_propt**2 * omega[0, 0]),
                MIN_HERITABILITY,
            )
            pi2_omega_sum2 = max(
                float(aux2[1, 1] - gmm_propt**2 * omega[1, 1]),
                MIN_HERITABILITY,
            )
            pi12_omega_sum12 = float(tissue_cross[0, 1] - gmm_propt**2 * omega[0, 1])
            pi12_bound = np.sqrt(pi2_omega_sum1 * pi2_omega_sum2) * MAX_CORR
            pi12_omega_sum12 = float(
                np.clip(pi12_omega_sum12, -pi12_bound, pi12_bound)
            )

    return (
        omega,
        omega_se,
        omega_p,
        run_gmm,
        run_tracecb,
        run_tracecb2,
        pi2_omega_sum,
        pi2_omega_sum1,
        pi2_omega_sum2,
        pi12_omega_sum12,
    )


def run_one(
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
    nt1,
    nt2,
    nsnp,
    propt,
    pcausal,
    out_dir,
    true_omega,
    id_sim,
    pop2_tissue_start=None,
    seed_base=None,
    seed_parts=(),
):
    path = build_result_path(
        runname,
        h1sq,
        h2sq,
        gc,
        n1,
        n2,
        nt1,
        nt2,
        nsnp,
        propt,
        pcausal,
        true_omega,
    )
    os.makedirs(os.path.join(out_dir, path), exist_ok=True)
    (
        omega_cb,
        omega_cb2,
        b1_hat,
        se1_hat,
        b2_hat,
        se2_hat,
        bt1_hat,
        se_t1_hat,
        bt2_hat,
        se_t2_hat,
        sig1,
        sig2,
        sigt1,
        sigt2,
        causal_ids,
        pi_mean1,
        pi_mean2,
    ) = generate_data_tracecb2(
        G1,
        G2,
        h1sq,
        h2sq,
        gc,
        n1,
        n2,
        nt1,
        nt2,
        nsnp,
        propt,
        pcausal,
        pop2_tissue_start=pop2_tissue_start,
        seed_base=seed_base,
        seed_parts=seed_parts,
    )
    gmm_propt = propt
    if true_omega:
        omega = omega_cb[:2, :2]
        omega_se = np.full((2, 2), np.nan)
        omega_p = np.full((2, 2), np.nan)
        run_gmm = should_run_true_omega_gmm(h1sq, h2sq, gc, omega)
        run_tracecb = run_gmm
        run_tracecb2 = run_gmm
        pi2_omega_sum = (
            omega_cb[2, 2]
            - gmm_propt**2 * omega_cb[1, 1]
            - 2 * gmm_propt * (omega_cb[1, 2] - gmm_propt * omega_cb[1, 1])
        )
        pi2_omega_sum = max(float(pi2_omega_sum), MIN_HERITABILITY)
        pi2_omega_sum1 = max(
            float(omega_cb2[1, 1] - gmm_propt**2 * omega[0, 0]),
            MIN_HERITABILITY,
        )
        pi2_omega_sum2 = max(
            float(omega_cb2[3, 3] - gmm_propt**2 * omega[1, 1]),
            MIN_HERITABILITY,
        )
        pi12_omega_sum12 = float(omega_cb2[1, 3] - gmm_propt**2 * omega[0, 1])
        pi12_bound = np.sqrt(pi2_omega_sum1 * pi2_omega_sum2) * MAX_CORR
        pi12_omega_sum12 = float(np.clip(pi12_omega_sum12, -pi12_bound, pi12_bound))
    else:
        (
            omega,
            omega_se,
            omega_p,
            run_gmm,
            run_tracecb,
            run_tracecb2,
            pi2_omega_sum,
            pi2_omega_sum1,
            pi2_omega_sum2,
            pi12_omega_sum12,
        ) = estimate_tracecb2_moments(
            b1_hat,
            se1_hat,
            b2_hat,
            se2_hat,
            bt1_hat,
            se_t1_hat,
            bt2_hat,
            se_t2_hat,
            n1,
            n2,
            nt1,
            nt2,
            ld1,
            ld2,
            ldx,
            gmm_propt,
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
    ) = run_tracecb2_inference_kernel(
        nsnp,
        run_gmm,
        run_tracecb,
        run_tracecb2,
        omega,
        pi2_omega_sum,
        pi2_omega_sum1,
        pi2_omega_sum2,
        pi12_omega_sum12,
        gmm_propt,
        b1_hat,
        se1_hat,
        bt1_hat,
        se_t1_hat,
        b2_hat,
        se2_hat,
        bt2_hat,
        se_t2_hat,
        ld1,
        ld2,
        ldx,
    )

    pop1_z = pop1_beta / (pop1_se + MIN_FLOAT)
    pop2_z = pop2_beta / (pop2_se + MIN_FLOAT)
    columns = [
        "causal",
        "sign1",
        "sign2",
        "sign_t1",
        "sign_t2",
        "z1_sumstat",
        "z1_cross",
        "z1_tissue",
        "z1_tracecb2",
        "z2_sumstat",
        "z2_cross",
        "z2_tissue",
        "z2_tracecb2",
        "zt1_sumstat",
        "zt2_sumstat",
        "z_meta",
        "z_metatissue",
        "gmm_propt",
    ]
    out = np.zeros((nsnp, len(columns)))
    out[causal_ids, 0] = 1
    out[:, 1] = sig1
    out[:, 2] = sig2
    out[:, 3] = sigt1
    out[:, 4] = sigt2
    out[:, 5:9] = pop1_z
    out[:, 9:13] = pop2_z
    out[:, 13] = bt1_hat / (se_t1_hat + MIN_FLOAT)
    out[:, 14] = bt2_hat / (se_t2_hat + MIN_FLOAT)
    out[:, 15] = meta_beta / (meta_se + MIN_FLOAT)
    out[:, 16] = meta_tissue_beta / (meta_tissue_se + MIN_FLOAT)
    out[:, 17] = gmm_propt
    np.savetxt(
        os.path.join(out_dir, path, f"simulation_{id_sim}.csv"),
        out,
        delimiter=",",
        header=",".join(columns),
        comments="",
    )

    omega_columns = [
        "id_sim",
        "omega11_true",
        "omega22_true",
        "omega12_true",
        "omega11_est",
        "omega22_est",
        "omega12_est",
        "omega11_p",
        "omega22_p",
        "omega12_p",
        "pi_mean1",
        "pi_mean2",
        "gmm_propt",
        "run_gmm",
        "run_tracecb",
        "run_tracecb2",
        "pi2_omega_sum",
        "pi2_omega_sum1",
        "pi2_omega_sum2",
        "pi12_omega_sum12",
    ]
    omega_values = np.array(
        [
            id_sim,
            omega_cb[0, 0],
            omega_cb[1, 1],
            omega_cb[0, 1],
            omega[0, 0],
            omega[1, 1],
            omega[0, 1],
            omega_p[0, 0],
            omega_p[1, 1],
            omega_p[0, 1],
            pi_mean1,
            pi_mean2,
            gmm_propt,
            float(run_gmm),
            float(run_tracecb),
            float(run_tracecb2),
            pi2_omega_sum,
            pi2_omega_sum1,
            pi2_omega_sum2,
            pi12_omega_sum12,
        ]
    ).reshape(1, -1)
    np.savetxt(
        os.path.join(out_dir, path, f"omega_summary_{id_sim}.csv"),
        omega_values,
        delimiter=",",
        header=",".join(omega_columns),
        comments="",
    )


def main():
    args = parse_args()
    if args.nt is not None:
        args.nt1 = args.nt if args.nt1 is None else args.nt1
        args.nt2 = args.nt if args.nt2 is None else args.nt2
    if args.nt1 is None or args.nt2 is None:
        raise ValueError("Provide --nt1 and --nt2, or provide --nt for both.")

    G1 = get_genotype(args.pop1_geno, args.nsnp)
    G2 = get_genotype(args.pop2_geno, args.nsnp)
    ld1, ld2, ldx = cal_ld(G1, G2)
    ld1, ld2, ldx = sanitize_ld_scores(ld1, ld2, ldx)

    validate_unit_interval("--h1sq", args.h1sq)
    validate_unit_interval("--h2sq", args.h2sq)
    validate_unit_interval("--propt", args.propt)
    validate_unit_interval("--pcausal", args.pcausal)
    trueOmega = not args.estimate_omega

    total = (
        len(args.h1sq)
        * len(args.h2sq)
        * len(args.gc)
        * len(args.n1)
        * len(args.n2)
        * len(args.nt1)
        * len(args.nt2)
        * len(args.propt)
        * len(args.pcausal)
    )
    done = 0
    start_time = time.time()
    base_seed = int(args.seed) if args.seed is not None else int(start_time)
    if max(args.n1) > G1.shape[0]:
        raise ValueError(
            f"max(n1)={max(args.n1)} exceeds population 1 genotype rows={G1.shape[0]}"
        )
    if max(args.n2) > G2.shape[0]:
        raise ValueError(
            f"max(n2)={max(args.n2)} exceeds population 2 genotype rows={G2.shape[0]}"
        )
    if max(args.nt1) > G1.shape[0]:
        raise ValueError(
            f"max(nt1)={max(args.nt1)} exceeds population 1 genotype rows={G1.shape[0]}"
        )
    if max(args.nt2) > G2.shape[0]:
        raise ValueError(
            f"max(nt2)={max(args.nt2)} exceeds population 2 genotype rows={G2.shape[0]}"
        )
    print("traceCB^2 simulation start at ", time.ctime())
    print("traceCB^2 simulation base seed: ", base_seed)
    for i, h1sq in enumerate(args.h1sq):
        for j, h2sq in enumerate(args.h2sq):
            for k, gc in enumerate(args.gc):
                for l, n1 in enumerate(args.n1):
                    for m, n2 in enumerate(args.n2):
                        for n, nt1 in enumerate(args.nt1):
                            for o, nt2 in enumerate(args.nt2):
                                for p, propt in enumerate(args.propt):
                                    for q, pcausal in enumerate(args.pcausal):
                                        for rep in range(args.nrep):
                                            seed_parts = (rep,)
                                            run_one(
                                                args.runname,
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
                                                nt1,
                                                nt2,
                                                args.nsnp,
                                                propt,
                                                pcausal,
                                                args.out_dir,
                                                trueOmega,
                                                rep,
                                                pop2_tissue_start=None,
                                                seed_base=base_seed,
                                                seed_parts=seed_parts,
                                            )
                                        done += 1
                                        if done % 5 == 0:
                                            print(
                                                f"traceCB^2 simulation {done}/{total} done, {time.ctime()}"
                                            )
    print("traceCB^2 simulation end at ", time.ctime())
    print(f"Total time: {time.time() - start_time:.2f} s for {total} settings")


if __name__ == "__main__":
    main()
