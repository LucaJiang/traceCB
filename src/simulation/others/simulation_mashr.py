"""Generate mashr benchmark simulations with true betas and estimated SEs."""

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

SIMULATION_DIR = Path(__file__).resolve().parents[1]
ROOT_DIR = SIMULATION_DIR.parents[1]
SRC_DIR = ROOT_DIR / "src"
for path in (SRC_DIR, SIMULATION_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from simulation import MIN_FLOAT, calculate_sumstats, get_genotype
from simulation_common import unknown_cell_effect_scale
from traceCB.utils import z2p


P_VAL_THRESHOLD = 0.05


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate simulation files with true betas and bhat/se for mashr."
    )
    parser.add_argument(
        "--pop1_geno",
        default="data/simulation/EAS_n5000_chr22_loci29.npy",
        type=str,
    )
    parser.add_argument(
        "--pop2_geno",
        default="data/simulation/EUR_n20000_chr22_loci29.npy",
        type=str,
    )
    parser.add_argument("--runname", default="nt_n2_propt_mashr", type=str)
    parser.add_argument("--out_dir", default="bench/result_mashr", type=str)
    parser.add_argument("--h1sq", default=0.1, type=float)
    parser.add_argument("--h2sq", default=0.1, type=float)
    parser.add_argument("--gc", default=0.7, type=float)
    parser.add_argument("--n1", default=100, type=int)
    parser.add_argument("--n2", default=400, type=int)
    parser.add_argument("--nt", default=5000, type=int)
    parser.add_argument("--nsnp", default=2000, type=int)
    parser.add_argument("--propt", default=[0.4, 0.8], type=float, nargs="+")
    parser.add_argument("--pcausal", default=0.005, type=float)
    parser.add_argument("--nrep", default=100, type=int)
    parser.add_argument("--seed", default=20260604, type=int)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def setting_path(runname, h1sq, h2sq, gc, n1, n2, nt, nsnp, propt, pcausal):
    return (
        f"{runname}/h1sq_{h1sq}_h2sq_{h2sq}_gc_{gc}_n1_{n1}_n2_{n2}_"
        f"nt_{nt}_nsnp_{nsnp}_propt_{propt}_pcausal_{pcausal}_omega_False"
    )


def generate_mashr_data(G1, G2, h1sq, h2sq, gc, n1, n2, nt, nsnp, propt, pcausal):
    omega_causal = np.array(
        [[h1sq, np.sqrt(h1sq * h2sq) * gc], [np.sqrt(h1sq * h2sq) * gc, h2sq]]
    )
    x1 = (G1[:n1, :] - np.mean(G1[:n1, :], axis=0)) / (
        np.std(G1[:n1, :], axis=0) + MIN_FLOAT
    )
    x2 = (G2[:n2, :] - np.mean(G2[:n2, :], axis=0)) / (
        np.std(G2[:n2, :], axis=0) + MIN_FLOAT
    )
    xt_raw = G2[n2 : n2 + nt, :]
    xt = (xt_raw - np.mean(xt_raw, axis=0)) / (np.std(xt_raw, axis=0) + MIN_FLOAT)

    num_causal = int(pcausal * nsnp)
    causal_ids = np.array([], dtype=int)
    beta1_true = np.zeros(nsnp)
    beta2_true = np.zeros(nsnp)
    if num_causal > 0:
        causal_ids = np.random.choice(np.arange(nsnp), num_causal, replace=False)
        beta_causal = np.random.multivariate_normal(
            mean=np.zeros(2),
            cov=omega_causal / (pcausal * nsnp),
            size=num_causal,
        )
        beta1_true[causal_ids] = beta_causal[:, 0]
        beta2_true[causal_ids] = beta_causal[:, 1]

    y1 = x1 @ beta1_true.T + np.sqrt(1 - h1sq) * np.random.randn(n1)
    y2 = x2 @ beta2_true.T + np.sqrt(1 - h2sq) * np.random.randn(n2)

    pi_ind = np.random.beta(
        (propt + MIN_FLOAT) * 5, (1 - propt + MIN_FLOAT) * 5, nt
    )
    pi_mean = np.mean(pi_ind)
    beta_unknown_true = np.zeros(nsnp)
    num_unknown_causal = int(pcausal * nsnp)
    if num_unknown_causal > 0:
        causal_unknown_id = np.random.choice(
            np.arange(nsnp), num_unknown_causal, replace=False
        )
        beta_unknown_true[causal_unknown_id] += np.random.normal(
            loc=0,
            scale=unknown_cell_effect_scale(h2sq, num_unknown_causal, 1),
            size=num_unknown_causal,
        )
    yt = (
        pi_ind * (xt @ beta2_true.T)
        + (1 - pi_ind) * (xt @ beta_unknown_true.T)
        + np.sqrt(
            np.maximum(1 - (pi_ind**2 + (1 - pi_ind) ** 2) * h2sq, MIN_FLOAT)
        )
        * np.random.randn(nt)
    )
    beta_tissue_mean_true = pi_mean * beta2_true + (1 - pi_mean) * beta_unknown_true

    b1_hat, se1_hat = calculate_sumstats(x1, y1, n1)
    b2_hat, se2_hat = calculate_sumstats(x2, y2, n2)
    bt_hat, se_t_hat = calculate_sumstats(xt, yt, nt)
    z1 = b1_hat / (se1_hat + MIN_FLOAT)
    z2 = b2_hat / (se2_hat + MIN_FLOAT)
    zt = bt_hat / (se_t_hat + MIN_FLOAT)

    return pd.DataFrame(
        {
            "causal": np.isin(np.arange(nsnp), causal_ids).astype(int),
            "beta1_true": beta1_true,
            "beta2_true": beta2_true,
            "beta_unknown_true": beta_unknown_true,
            "beta_tissue_mean_true": beta_tissue_mean_true,
            "pi_mean": pi_mean,
            "b1_hat": b1_hat,
            "se1_hat": se1_hat,
            "b2_hat": b2_hat,
            "se2_hat": se2_hat,
            "bt_hat": bt_hat,
            "se_t_hat": se_t_hat,
            "sign1": (z2p(z1) < P_VAL_THRESHOLD).astype(int),
            "sign2": (z2p(z2) < P_VAL_THRESHOLD).astype(int),
            "sign_t": (z2p(zt) < P_VAL_THRESHOLD).astype(int),
            "z1_sumstat": z1,
            "z2_sumstat": z2,
            "zt_sumstat": zt,
        }
    )


def main():
    args = parse_args()
    start_time = time.time()
    print("Simulation start at", time.ctime())
    print("Simulation base seed:", args.seed)
    g1 = get_genotype(args.pop1_geno, args.nsnp)
    g2 = get_genotype(args.pop2_geno, args.nsnp)

    written = 0
    for p_idx, propt in enumerate(args.propt):
        sim_path = setting_path(
            args.runname,
            args.h1sq,
            args.h2sq,
            args.gc,
            args.n1,
            args.n2,
            args.nt,
            args.nsnp,
            propt,
            args.pcausal,
        )
        output_dir = os.path.join(args.out_dir, sim_path)
        os.makedirs(output_dir, exist_ok=True)
        for rep in range(args.nrep):
            output_file = os.path.join(output_dir, f"simulation_{rep}.csv")
            if os.path.exists(output_file) and not args.force:
                continue
            np.random.seed(args.seed + p_idx * args.nrep + rep)
            df = generate_mashr_data(
                g1,
                g2,
                args.h1sq,
                args.h2sq,
                args.gc,
                args.n1,
                args.n2,
                args.nt,
                args.nsnp,
                propt,
                args.pcausal,
            )
            df.to_csv(output_file, index=False)
            written += 1
        print(f"propt={propt} done at {time.ctime()}")
    print("Simulation end at", time.ctime())
    print(f"Written files: {written}; total time: {time.time() - start_time:.2f}s")


if __name__ == "__main__":
    main()
