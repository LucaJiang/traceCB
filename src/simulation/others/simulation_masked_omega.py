"""Simulate and compare original traceC/traceCB with masked-input variants.

This follows the original simulation workflow:

1. Generate pop1 sc, pop2 sc, and pop2 bulk summary statistics.
2. Use true omega from the simulated effects by default, or estimate omega if
   --estimate_omega is set.
3. Run original traceC and traceCB with pop1 as the target.
4. Run two traceCB-style variants by masking one sc input:
   - pop1 sc + pop2 bulk: mask pop2 sc.
   - pop2 sc + pop2 bulk: mask pop1 sc, but still estimate the pop1 target.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from numba import njit, prange

SIMULATION_DIR = Path(__file__).resolve().parents[1]
ROOT_DIR = SIMULATION_DIR.parents[1]
SRC_DIR = ROOT_DIR / "src"
for path in (SRC_DIR, SIMULATION_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from simulation import cal_ld, generate_data, get_genotype  # noqa: E402
from simulation_utils import (  # noqa: E402
    sanitize_ld_scores,
    validate_unit_interval,
)
from traceCB.gmm import GMM, GMMtissue  # noqa: E402
from traceCB.ldsc import Run_Cross_LDSC  # noqa: E402
from traceCB.utils import MIN_FLOAT, MIN_HERITABILITY, z2p  # noqa: E402

P_THRESHOLD_Z = 1.959963984540054
DEFAULT_OMEGA_P_THRESHOLD = 0.10
MAX_CORR = 0.99
METHODS = (
    "original_pop1",
    "original_pop2",
    "traceC_original",
    "traceCB_original",
    "pop1sc_pop2bulk",
    "pop2sc_pop2bulk",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare original traceC/traceCB with masked omega variants."
    )
    parser.add_argument(
        "--pop1_geno",
        default="data/simulation/EAS_n5000_chr22_loci29.npy",
        help="Population 1 genotype .npy file.",
    )
    parser.add_argument(
        "--pop2_geno",
        default="data/simulation/EUR_n20000_chr22_loci29.npy",
        help="Population 2 genotype .npy file.",
    )
    parser.add_argument("--runname", default="masked_omega_compare")
    parser.add_argument("--h1sq", type=float, nargs="+", default=[0.1])
    parser.add_argument("--h2sq", type=float, nargs="+", default=[0.1])
    parser.add_argument("--gc", type=float, nargs="+", default=[0.7])
    parser.add_argument("--n1", type=int, nargs="+", default=[100])
    parser.add_argument("--n2", type=int, nargs="+", default=[200])
    parser.add_argument("--nt", type=int, nargs="+", default=[1000])
    parser.add_argument("--nsnp", type=int, default=1000)
    parser.add_argument("--propt", type=float, nargs="+", default=[0.2, 0.5, 0.8])
    parser.add_argument("--pcausal", type=float, nargs="+", default=[0.01])
    parser.add_argument("--nrep", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260425)
    parser.add_argument(
        "--estimate_omega",
        action="store_true",
        help=(
            "Estimate omega from summary statistics. By default this script uses "
            "the true omega from simulated effects, matching simulation.py."
        ),
    )
    parser.add_argument(
        "--omega_p_threshold",
        type=float,
        default=DEFAULT_OMEGA_P_THRESHOLD,
        help="P-value threshold used with --use_original_gate.",
    )
    parser.add_argument(
        "--use_original_gate",
        dest="use_original_gate",
        action="store_true",
        default=True,
        help=(
            "Apply the original estimated-omega significance gate. By default, "
            "this is enabled to match src/simulation/simulation.py."
        ),
    )
    parser.add_argument(
        "--no_original_gate",
        dest="use_original_gate",
        action="store_false",
        help=(
            "Estimate omega but force all methods to run. This is mainly for "
            "diagnostics and can create tissue effects when rho is near zero."
        ),
    )
    parser.add_argument(
        "--missing_se",
        type=float,
        default=1e6,
        help="SE assigned to a masked single-cell summary statistic.",
    )
    parser.add_argument(
        "--out_dir",
        default="bench/result/img",
        help="Output directory.",
    )
    return parser.parse_args()


def clip_correlation(var1: float, var2: float, cov: float) -> float:
    var1 = max(float(var1), MIN_HERITABILITY)
    var2 = max(float(var2), MIN_HERITABILITY)
    denominator = np.sqrt(var1 * var2) + MIN_FLOAT
    cor = cov / denominator
    if abs(cor) > MAX_CORR:
        cov = np.sign(cor) * MAX_CORR * denominator
    return float(cov)


def sanitize_omega(omega: np.ndarray) -> np.ndarray:
    omega = np.array(omega, dtype=np.float64, copy=True)
    omega[0, 0] = max(omega[0, 0], MIN_HERITABILITY)
    omega[1, 1] = max(omega[1, 1], MIN_HERITABILITY)
    omega[0, 1] = clip_correlation(omega[0, 0], omega[1, 1], omega[0, 1])
    omega[1, 0] = omega[0, 1]
    return omega


def estimate_omega_terms(
    b1_hat: np.ndarray,
    se1_hat: np.ndarray,
    b2_hat: np.ndarray,
    se2_hat: np.ndarray,
    bt_hat: np.ndarray,
    se_t_hat: np.ndarray,
    ld1: np.ndarray,
    ld2: np.ndarray,
    ldx: np.ndarray,
    n1: int,
    n2: int,
    nt: int,
    propt: float,
    omega_p_threshold: float,
    use_original_gate: bool,
) -> tuple[np.ndarray, float, bool, bool, bool, bool, float, float]:
    """Estimate the omega terms used by the original simulation code."""
    omega_raw, omega_se = Run_Cross_LDSC(
        b1_hat / (se1_hat + MIN_FLOAT),
        n1,
        ld1,
        b2_hat / (se2_hat + MIN_FLOAT),
        n2,
        ld2,
        ldx,
        np.array([1.0, 1.0, 0.0]),
    )
    omega_p = z2p(omega_raw / (omega_se + MIN_FLOAT))
    omega = sanitize_omega(omega_raw)
    run_tracec = (
        bool(np.all(omega_p < omega_p_threshold)) if use_original_gate else True
    )
    run_pop1sc_bulk = bool(omega_p[0, 0] < omega_p_threshold)
    run_pop2sc_bulk = bool(omega_p[1, 1] < omega_p_threshold)

    aux_cov_p = np.nan
    pi2_omega_sum = MIN_HERITABILITY
    run_tracecb = False
    if use_original_gate and not run_tracec:
        return (
            omega,
            pi2_omega_sum,
            run_tracec,
            run_tracecb,
            run_pop1sc_bulk,
            run_pop2sc_bulk,
            float(omega_p[0, 1]),
            aux_cov_p,
        )

    aux_omega, aux_omega_se = Run_Cross_LDSC(
        b2_hat / (se2_hat + MIN_FLOAT),
        n2,
        ld2,
        bt_hat / (se_t_hat + MIN_FLOAT),
        nt,
        ldx,
        ldx,
        np.array([1.0, 1.0, 0.0]),
    )
    aux_omega_p = z2p(aux_omega / (aux_omega_se + MIN_FLOAT))
    aux_cov_p = float(aux_omega_p[0, 1])
    pi2_omega_raw = (
        aux_omega[1, 1]
        - propt**2 * omega[1, 1]
        - 2 * propt * np.maximum(aux_omega[0, 1] - propt * omega[1, 1], 0)
    )
    pi2_omega_sum = max(float(pi2_omega_raw), MIN_HERITABILITY)

    if use_original_gate:
        run_tracecb = bool(run_tracec and np.all(aux_omega_p < omega_p_threshold))
    else:
        run_tracecb = True

    return (
        omega,
        pi2_omega_sum,
        run_tracec,
        run_tracecb,
        run_pop1sc_bulk,
        run_pop2sc_bulk,
        float(omega_p[0, 1]),
        aux_cov_p,
    )


def safe_corr(omega: np.ndarray) -> float:
    denominator = np.sqrt(
        max(omega[0, 0], MIN_HERITABILITY) * max(omega[1, 1], MIN_HERITABILITY)
    )
    return float(omega[0, 1] / (denominator + MIN_FLOAT))


def true_omega_terms(omega_true: np.ndarray, propt: float) -> tuple[np.ndarray, float]:
    omega = sanitize_omega(omega_true[:2, :2])
    pi2_omega_sum = (
        omega_true[2, 2]
        - propt**2 * omega_true[1, 1]
        - 2 * propt * (omega_true[1, 2] - propt * omega_true[1, 1])
    )
    return omega, float(pi2_omega_sum)


@njit(nogil=True, parallel=True, cache=True)
def run_pop1_target_methods(
    b1_hat: np.ndarray,
    se1_hat: np.ndarray,
    b2_hat: np.ndarray,
    se2_hat: np.ndarray,
    bt_hat: np.ndarray,
    se_t_hat: np.ndarray,
    ld1: np.ndarray,
    ld2: np.ndarray,
    ldx: np.ndarray,
    omega: np.ndarray,
    pi2_omega_sum: float,
    propt: float,
    run_tracec: bool,
    run_tracecb: bool,
    run_pop1sc_bulk: bool,
    run_pop2sc_bulk: bool,
    missing_se: float,
) -> np.ndarray:
    z = np.zeros((6, b1_hat.shape[0]))
    eye2 = np.eye(2)
    eye3 = np.eye(3)
    dummy_beta = 0.0
    dummy_se = missing_se
    dummy_ld = MIN_FLOAT
    dummy_cross_ld = 0.0
    dummy_pi2_omega = MIN_HERITABILITY

    for j in prange(b1_hat.shape[0]):
        z_sumstat = b1_hat[j] / (se1_hat[j] + MIN_FLOAT)
        z[0, j] = z_sumstat
        z[1, j] = b2_hat[j] / (se2_hat[j] + MIN_FLOAT)

        if run_tracec:
            b1_c, se1_c, _b2_c, _se2_c = GMM(
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
            z[2, j] = b1_c / (se1_c + MIN_FLOAT)
        else:
            z[2, j] = z_sumstat

        if run_tracecb:
            b1_t, se1_t, _b2_t, _se2_t = GMMtissue(
                omega,
                eye3,
                b1_hat[j],
                se1_hat[j],
                ld1[j],
                b2_hat[j],
                se2_hat[j],
                ld2[j],
                ldx[j],
                bt_hat[j],
                se_t_hat[j],
                pi2_omega_sum,
                propt,
            )
            z[3, j] = b1_t / (se1_t + MIN_FLOAT)
        else:
            z[3, j] = z[2, j]

        if run_pop1sc_bulk:
            omega_pop1sc_bulk = np.array([[MIN_HERITABILITY, 0.0], [0.0, omega[0, 0]]])
            (
                _dummy_b,
                _dummy_se,
                b_pop1sc_bulk,
                se_pop1sc_bulk,
            ) = GMMtissue(
                omega_pop1sc_bulk,
                eye3,
                dummy_beta,
                dummy_se,
                dummy_ld,
                b1_hat[j],
                se1_hat[j],
                ld2[j],
                dummy_cross_ld,
                bt_hat[j],
                se_t_hat[j],
                dummy_pi2_omega,
                propt,
            )
            z[4, j] = b_pop1sc_bulk / (se_pop1sc_bulk + MIN_FLOAT)
        else:
            z[4, j] = z_sumstat

        if run_pop2sc_bulk:
            omega_pop2sc_bulk = np.array([[MIN_HERITABILITY, 0.0], [0.0, omega[1, 1]]])
            (
                _dummy_b,
                _dummy_se,
                b_pop2sc_bulk,
                se_pop2sc_bulk,
            ) = GMMtissue(
                omega_pop2sc_bulk,
                eye3,
                dummy_beta,
                dummy_se,
                dummy_ld,
                b2_hat[j],
                se2_hat[j],
                ld2[j],
                dummy_cross_ld,
                bt_hat[j],
                se_t_hat[j],
                dummy_pi2_omega,
                propt,
            )
            z[5, j] = b_pop2sc_bulk / (se_pop2sc_bulk + MIN_FLOAT)
        else:
            z[5, j] = 0.0

    return z


@njit(nogil=True, parallel=True, cache=True)
def calculate_metrics_all(causal: np.ndarray, z_scores: np.ndarray) -> np.ndarray:
    metrics = np.zeros((z_scores.shape[0], 2))
    causal_count = 0
    noncausal_count = 0
    for j in range(causal.shape[0]):
        if causal[j]:
            causal_count += 1
        else:
            noncausal_count += 1

    for method_id in prange(z_scores.shape[0]):
        causal_hits = 0
        noncausal_hits = 0
        for j in range(z_scores.shape[1]):
            z_abs = abs(z_scores[method_id, j])
            if z_abs > P_THRESHOLD_Z:
                if causal[j]:
                    causal_hits += 1
                else:
                    noncausal_hits += 1
        if causal_count > 0:
            metrics[method_id, 0] = causal_hits / causal_count
        else:
            metrics[method_id, 0] = np.nan
        if noncausal_count > 0:
            metrics[method_id, 1] = noncausal_hits / noncausal_count
        else:
            metrics[method_id, 1] = np.nan
    return metrics


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir) / args.runname
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading genotypes with nsnp={args.nsnp}")
    G1 = get_genotype(args.pop1_geno, args.nsnp)
    G2 = get_genotype(args.pop2_geno, args.nsnp)
    print("Calculating LD scores")
    ld1, ld2, ldx = cal_ld(G1, G2)
    ld1, ld2, ldx = sanitize_ld_scores(ld1, ld2, ldx)

    validate_unit_interval("--h1sq", args.h1sq)
    validate_unit_interval("--h2sq", args.h2sq)
    validate_unit_interval("--propt", args.propt)
    validate_unit_interval("--pcausal", args.pcausal)

    total_settings = (
        len(args.h1sq)
        * len(args.h2sq)
        * len(args.gc)
        * len(args.n1)
        * len(args.n2)
        * len(args.nt)
        * len(args.propt)
        * len(args.pcausal)
    )
    rows: list[dict] = []
    setting_id = 0
    start = time.time()
    pop2_tissue_start = None
    max_tissue_end = max(args.n2) + max(args.nt)
    if max_tissue_end > G2.shape[0]:
        raise ValueError(
            f"required n2/nt sample rows = {max_tissue_end} exceeds population 2 "
            f"genotype rows={G2.shape[0]}"
        )

    for h1sq in args.h1sq:
        for h2sq in args.h2sq:
            for gc in args.gc:
                for n1 in args.n1:
                    for n2 in args.n2:
                        for nt in args.nt:
                            for propt in args.propt:
                                for pcausal in args.pcausal:
                                    setting_id += 1
                                    print(
                                        f"Setting {setting_id}/{total_settings}: "
                                        f"h1sq={h1sq}, h2sq={h2sq}, gc={gc}, "
                                        f"n1={n1}, n2={n2}, nt={nt}, "
                                        f"propt={propt}, pcausal={pcausal}"
                                    )
                                    for rep in range(args.nrep):
                                        seed_parts = (
                                            h1sq,
                                            h2sq,
                                            gc,
                                            n1,
                                            n2,
                                            nt,
                                            args.nsnp,
                                            propt,
                                            pcausal,
                                            rep,
                                        )
                                        (
                                            _omega_true,
                                            b1_hat,
                                            se1_hat,
                                            b2_hat,
                                            se2_hat,
                                            bt_hat,
                                            se_t_hat,
                                            _sig1,
                                            _sig2,
                                            _sigt,
                                            causal_ids,
                                            _pi_mean,
                                        ) = generate_data(
                                            G1,
                                            G2,
                                            h1sq,
                                            h2sq,
                                            gc,
                                            n1,
                                            n2,
                                            nt,
                                            args.nsnp,
                                            propt,
                                            pcausal,
                                            tissue_start=pop2_tissue_start,
                                            seed_base=args.seed,
                                            seed_parts=seed_parts,
                                        )
                                        omega_true_corr = safe_corr(_omega_true)
                                        if args.estimate_omega:
                                            (
                                                omega,
                                                pi2_omega_sum,
                                                run_tracec,
                                                run_tracecb,
                                                run_pop1sc_bulk,
                                                run_pop2sc_bulk,
                                                omega_cov_p,
                                                aux_cov_p,
                                            ) = estimate_omega_terms(
                                                b1_hat,
                                                se1_hat,
                                                b2_hat,
                                                se2_hat,
                                                bt_hat,
                                                se_t_hat,
                                                ld1,
                                                ld2,
                                                ldx,
                                                n1,
                                                n2,
                                                nt,
                                                propt,
                                                args.omega_p_threshold,
                                                args.use_original_gate,
                                            )
                                            omega_est_corr = safe_corr(omega)
                                        else:
                                            omega, pi2_omega_sum = true_omega_terms(
                                                _omega_true, propt
                                            )
                                            run_tracec = True
                                            run_tracecb = True
                                            run_pop1sc_bulk = True
                                            run_pop2sc_bulk = True
                                            omega_cov_p = np.nan
                                            aux_cov_p = np.nan
                                            omega_est_corr = np.nan
                                        z_scores = run_pop1_target_methods(
                                            b1_hat,
                                            se1_hat,
                                            b2_hat,
                                            se2_hat,
                                            bt_hat,
                                            se_t_hat,
                                            ld1,
                                            ld2,
                                            ldx,
                                            omega,
                                            pi2_omega_sum,
                                            propt,
                                            run_tracec,
                                            run_tracecb,
                                            run_pop1sc_bulk,
                                            run_pop2sc_bulk,
                                            args.missing_se,
                                        )
                                        causal = np.zeros(args.nsnp, dtype=bool)
                                        causal[causal_ids] = True
                                        metrics = calculate_metrics_all(
                                            causal, z_scores
                                        )

                                        context = {
                                            "rep": rep,
                                            "h1sq": h1sq,
                                            "h2sq": h2sq,
                                            "gc": gc,
                                            "n1": n1,
                                            "n2": n2,
                                            "nt": nt,
                                            "nsnp": args.nsnp,
                                            "propt": propt,
                                            "pcausal": pcausal,
                                            "omega_true_corr": omega_true_corr,
                                            "omega_est_corr": omega_est_corr,
                                            "estimate_omega": args.estimate_omega,
                                            "omega_cov_p": omega_cov_p,
                                            "aux_cov_p": aux_cov_p,
                                            "run_tracec": run_tracec,
                                            "run_tracecb": run_tracecb,
                                            "run_pop1sc_bulk": run_pop1sc_bulk,
                                            "run_pop2sc_bulk": run_pop2sc_bulk,
                                            "use_original_gate": args.use_original_gate,
                                        }
                                        for method_id, method in enumerate(METHODS):
                                            rows.append(
                                                {
                                                    **context,
                                                    "method": method,
                                                    "power": metrics[method_id, 0],
                                                    "alpha": metrics[method_id, 1],
                                                }
                                            )

    result_df = pd.DataFrame(rows)
    result_csv = out_dir / "replicate_metrics.csv"
    result_df.to_csv(result_csv, index=False)
    summary_df = (
        result_df.groupby(
            [
                "method",
                "h1sq",
                "h2sq",
                "gc",
                "n1",
                "n2",
                "nt",
                "nsnp",
                "propt",
                "pcausal",
                "estimate_omega",
            ],
            dropna=False,
        )
        .agg(
            power_mean=("power", "mean"),
            power_sd=("power", "std"),
            alpha_mean=("alpha", "mean"),
            alpha_sd=("alpha", "std"),
            omega_true_corr_mean=("omega_true_corr", "mean"),
            omega_true_corr_sd=("omega_true_corr", "std"),
            omega_est_corr_mean=("omega_est_corr", "mean"),
            omega_est_corr_sd=("omega_est_corr", "std"),
            run_tracec_rate=("run_tracec", "mean"),
            run_tracecb_rate=("run_tracecb", "mean"),
            run_pop1sc_bulk_rate=("run_pop1sc_bulk", "mean"),
            run_pop2sc_bulk_rate=("run_pop2sc_bulk", "mean"),
            nrep=("rep", "nunique"),
        )
        .reset_index()
    )
    summary_csv = out_dir / "summary_metrics.csv"
    summary_df.to_csv(summary_csv, index=False)

    print(f"Saved replicate metrics to {result_csv}")
    print(f"Saved summary metrics to {summary_csv}")
    print(f"Finished in {time.time() - start:.2f} seconds")


if __name__ == "__main__":
    main()
