# python -m src.simulation.simulation_double_gmm --pop1_geno data/simulation/EAS_n5000_chr22_loci29.npy --pop2_geno data/simulation/EUR_n20000_chr22_loci29.npy --runname nt_n2_propt --h1sq 0.1 --h2sq 0.1 --gc 0.7 --n1 100 --n2 100 200 400 --nt 500 10000 --nsnp 2000 --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 --out_dir bench/result_double --nrep 50 --estimate_omega
# python -m src.simulation.simulation_double_gmm --pop1_geno data/simulation/EAS_n5000_chr22_loci29.npy --pop2_geno data/simulation/EUR_n20000_chr22_loci29.npy --runname h2sq_gc_propt --h1sq 0.1 --h2sq 0.1 0.2 --gc 0.01 0.5 0.9 --n1 100 --n2 200 --nt 500 --nsnp 2000 --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 --out_dir bench/result_double --nrep 50 --estimate_omega
# python -m src.simulation.simulation_double_gmm --pop1_geno data/simulation/EAS_n5000_chr22_loci29.npy --pop2_geno data/simulation/EUR_n20000_chr22_loci29.npy --runname n1_pcausal_propt --h1sq 0.1 --h2sq 0.1 --gc 0.7 --n1 50 100 --n2 200 --nt 500 --nsnp 2000 --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 0.01 0.02 --out_dir bench/result_double --nrep 50 --estimate_omega
# python -m src.simulation.simulation_double_gmm --pop1_geno data/simulation/EAS_n5000_chr22_loci29.npy --pop2_geno data/simulation/EUR_n20000_chr22_loci29.npy --runname alpha_h2sq_pcausal_propt --h1sq 0.000000000001 --h2sq 0.1 0.2 --gc 0 --n1 100 --n2 400 --nt 500 --nsnp 2000 --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 0.01 0.02 --out_dir bench/result_double --nrep 50 --estimate_omega


import argparse
import time
import os
import numpy as np
import pandas as pd
from numba import njit, prange

from ..traceCB.ldsc import Run_cross_LDSC
from ..traceCB.gmm import GMM, GMMtissue
from ..traceCB.utils import (
    make_pd_shrink,
    z2p,
    MIN_HERITABILITY,
    MIN_FLOAT,
    P_VAL_THRED,
    MAX_CORR,
    eSNP_THRESHOLD,
)
from .simulation import (
    generate_data,
    get_genotype,
    cal_ld,
    calculate_sumstats,
    re2_meta,
)

SNP_START = 1000  # skip first part of chr


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate data and simulation for double GMM"
    )
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
        default="double_gmm",
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
        "--estimate_omega",
        action="store_true",
        help="Estimate omega from data, otherwise use true omega",
    )
    parser.add_argument(
        "--out_dir",
        default="bench/result_double",
        type=str,
        help="Output directory",
    )
    parser.add_argument(
        "--nrep",
        type=int,
        default=100,
        help="number of repetition for simulation",
    )
    return parser.parse_args()


def calculate_Neff(z: np.ndarray, omega: float, ld: np.ndarray) -> float:
    """
    Calculate effective sample size from z-scores, heritability omega, and LD vector.
    """
    return (np.mean(z**2) - 1) / (omega + MIN_FLOAT) / (np.mean(ld).item() + MIN_FLOAT)


def count_eSNPs(pval_list: np.ndarray) -> int:
    """
    Count the number of eSNPs based on p-values.
    """
    return np.sum(pval_list < eSNP_THRESHOLD).item()


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
    pcausal,
    out_dir,
    true_omega,
    id_sim,
):
    """
    Run simulation for Double GMM, save results to out_dir
    """
    simulation_path = f"{runname}/h1sq_{h1sq}_h2sq_{h2sq}_gc_{gc}_n1_{n1}_n2_{n2}_nt_{nt}_nsnp_{nsnp}_propt_{propt}_pcausal_{pcausal}_omega_{true_omega}"
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
        pi_mean,
    ) = generate_data(G1, G2, h1sq, h2sq, gc, n1, n2, nt, nsnp, propt, pcausal)

    # Z-scores
    z1 = b1_hat / se1_hat
    z2 = b2_hat / se2_hat
    zt = bt_hat / se_t_hat

    # Sample sizes (constant for all SNPs in simulation)
    n1_vec = np.full(nsnp, n1)
    n2_vec = np.full(nsnp, n2)
    nt_vec = np.full(nsnp, nt)

    # 1. Run TAR + AUX to get traceC (OmegaC)
    run_trace_c = False
    OmegaC, OmegaC_se = Run_cross_LDSC(
        z1,
        n1_vec,
        ld1,
        z2,
        n2_vec,
        ld2,
        ldx,
        np.array([1.0, 1.0, 0.0]),  # Intercepts fixed
    )
    OmegaC_p = z2p(OmegaC / OmegaC_se)
    if np.all(OmegaC_p < P_VAL_THRED):
        run_trace_c = True

    # Run GMM (TraceC)
    gmm_b_tar = np.full(nsnp, np.nan)
    gmm_se_tar = np.full(nsnp, np.nan)
    gmm_b_aux = np.full(nsnp, np.nan)
    gmm_se_aux = np.full(nsnp, np.nan)

    if run_trace_c:
        for i in range(nsnp):
            (gmm_b_tar[i], gmm_se_tar[i], gmm_b_aux[i], gmm_se_aux[i]) = GMM(
                OmegaC,
                np.eye(2),
                b1_hat[i],
                se1_hat[i],
                ld1[i],
                b2_hat[i],
                se2_hat[i],
                ld2[i],
                ldx[i],
            )
    else:
        gmm_b_tar = b1_hat
        gmm_se_tar = se1_hat
        gmm_b_aux = b2_hat
        gmm_se_aux = se2_hat

    gmm_z_tar = gmm_b_tar / gmm_se_tar
    gmm_z_aux = gmm_b_aux / gmm_se_aux

    # 2. Run AUX + Tissue to get updated AUX
    run_aux_tissue = False
    OmegaCB_aux_tissue = None

    Omega_aux_tissue, Omega_aux_tissue_se = Run_cross_LDSC(
        z2,
        n2_vec,
        ld2,
        zt,
        nt_vec,
        ld2,
        ld2,
        np.array([1.0, 1.0, 0.0]),
    )
    Omega_aux_tissue_p = z2p(Omega_aux_tissue / Omega_aux_tissue_se)

    if np.all(Omega_aux_tissue_p < P_VAL_THRED):
        run_aux_tissue = True
        # Construct OmegaCB for GMMtissue
        # This is for updating Aux using Tissue. Target is dummy.
        OmegaCB_aux_tissue = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, Omega_aux_tissue[0, 0], Omega_aux_tissue[0, 1]],
                [0.0, Omega_aux_tissue[0, 1], Omega_aux_tissue[1, 1]],
            ]
        )

    # Update AUX summary statistics with GMMtissue
    gmm_b_aux_t = np.full(nsnp, np.nan)
    gmm_se_aux_t = np.full(nsnp, np.nan)

    # Initialize TraceCB results
    gmm_b_tar_cb = np.full(nsnp, np.nan)
    gmm_se_tar_cb = np.full(nsnp, np.nan)
    gmm_b_aux_cb = np.full(nsnp, np.nan)
    gmm_se_aux_cb = np.full(nsnp, np.nan)

    run_trace_cb = False
    if run_aux_tissue:
        for i in range(nsnp):
            (
                _,
                _,
                gmm_b_aux_t[i],
                gmm_se_aux_t[i],
            ) = GMMtissue(
                OmegaCB_aux_tissue,
                np.eye(3),
                b1_hat[i],  # Dummy
                se1_hat[i],  # Dummy
                ld1[i],  # Dummy
                b2_hat[i],
                se2_hat[i],
                ld2[i],
                ldx[i],
                bt_hat[i],
                se_t_hat[i],
                propt,
            )

        # 3. Run TAR + updated AUX to get final TAR (TraceCB)
        # Use gmm_b_aux_t / gmm_se_aux_t as z-scores for updated Aux
        z2_updated = gmm_b_aux_t / gmm_se_aux_t
        Omega_final, Omega_final_se = Run_cross_LDSC(
            z1,
            n1_vec,
            ld1,
            z2_updated,
            n2_vec,
            ld2,
            ldx,
            np.array([1.0, 1.0, 0.0]),
        )
        Omega_final_p = z2p(Omega_final / Omega_final_se)
        if np.all(Omega_final_p < P_VAL_THRED):
            run_trace_cb = True
            for i in range(nsnp):
                (
                    gmm_b_tar_cb[i],
                    gmm_se_tar_cb[i],
                    gmm_b_aux_cb[i],
                    gmm_se_aux_cb[i],
                ) = GMM(
                    Omega_final,
                    np.eye(2),
                    b1_hat[i],
                    se1_hat[i],
                    ld1[i],
                    gmm_b_aux_t[i],
                    gmm_se_aux_t[i],
                    ld2[i],
                    ldx[i],
                )

    if not run_trace_cb:
        gmm_b_tar_cb = gmm_b_tar
        gmm_se_tar_cb = gmm_se_tar
        gmm_b_aux_cb = gmm_b_aux
        gmm_se_aux_cb = gmm_se_aux
    gmm_z_tar_cb = gmm_b_tar_cb / gmm_se_tar_cb
    gmm_z_aux_cb = gmm_b_aux_cb / gmm_se_aux_cb

    ## meta-analysis
    meta_beta = np.zeros((nsnp,))
    meta_se = np.zeros((nsnp,))
    meta_tissue_beta = np.zeros((nsnp,))
    meta_tissue_se = np.zeros((nsnp,))
    for j in prange(nsnp):
        ## meta-analysis
        meta_beta[j], meta_se[j] = re2_meta(
            np.array([b1_hat[j], b2_hat[j]]),
            np.array([se1_hat[j], se2_hat[j]]),
        )
        meta_tissue_beta[j], meta_tissue_se[j] = re2_meta(
            np.array([b1_hat[j], b2_hat[j], bt_hat[j]]),
            np.array([se1_hat[j], se2_hat[j], se_t_hat[j]]),
        )

    # save results
    # Columns: causal,sign1,sign2,sign_t,z1_sumstat,z1_cross,z1_tissue,z2_sumstat,z2_cross,z2_tissue,zt_sumstat,z_meta,z_metatissue
    # Mapping:
    # z1_sumstat: z1
    # z1_cross: gmm_z_tar (TraceC)
    # z1_tissue: gmm_z_tar_cb (TraceCB)
    # z2_sumstat: z2
    # z2_cross: gmm_z_aux (TraceC)
    # z2_tissue: gmm_z_aux_cb (TraceCB) or z2_updated?
    # In simulation.py, z2_tissue is from GMMtissue.
    # Here we use TraceCB which is GMM(Tar, UpdatedAux).
    # So gmm_z_aux_cb is the corresponding Aux estimate in the final step.

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
    ]
    all_results = np.zeros((nsnp, len(all_results_columns)))
    all_results[causal_ids, 0] = 1
    all_results[:, 1] = sig1
    all_results[:, 2] = sig2
    all_results[:, 3] = sigt
    all_results[:, 4] = z1
    all_results[:, 5] = gmm_z_tar
    all_results[:, 6] = gmm_z_tar_cb
    all_results[:, 7] = z2
    all_results[:, 8] = gmm_z_aux
    all_results[:, 9] = gmm_z_aux_cb
    all_results[:, 10] = zt
    all_results[:, 11] = meta_beta / meta_se
    all_results[:, 12] = meta_tissue_beta / meta_tissue_se

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
    ## parse
    h1sq = args.h1sq if isinstance(args.h1sq, list) else [args.h1sq]
    h2sq = args.h2sq if isinstance(args.h2sq, list) else [args.h2sq]
    gc = args.gc if isinstance(args.gc, list) else [args.gc]
    n1 = args.n1 if isinstance(args.n1, list) else [args.n1]
    n2 = args.n2 if isinstance(args.n2, list) else [args.n2]
    nt = args.nt if isinstance(args.nt, list) else [args.nt]
    propt = args.propt if isinstance(args.propt, list) else [args.propt]
    pcausal = args.pcausal if isinstance(args.pcausal, list) else [args.pcausal]
    nsnp = args.nsnp
    trueOmega = not args.estimate_omega

    ## process indicators
    total_combinations = (
        len(h1sq)
        * len(h2sq)
        * len(gc)
        * len(n1)
        * len(n2)
        * len(nt)
        * len(propt)
        * len(pcausal)
    )
    have_run = 0

    start_time = time.time()
    print("Simulation start at ", time.ctime())
    for i, h1sqi in enumerate(h1sq):
        for j, h2sqj in enumerate(h2sq):
            for k, gck in enumerate(gc):
                for l, n1l in enumerate(n1):
                    for m, n2m in enumerate(n2):
                        for n, ntn in enumerate(nt):
                            for p, proptp in enumerate(propt):
                                for q, pcausalq in enumerate(pcausal):
                                    for r in range(args.nrep):
                                        np.random.seed(
                                            int(start_time)
                                            + i
                                            + j
                                            + k
                                            + l
                                            + m
                                            + n
                                            + p
                                            + q
                                            + r
                                        )
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
                                            pcausalq,
                                            args.out_dir,
                                            trueOmega,
                                            r,
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
