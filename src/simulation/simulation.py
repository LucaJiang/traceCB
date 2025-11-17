# generate test data and run simulation for gmm
# python -m src.simulation.simulation --pop1_geno data/simulation/EAS_n5000_chr22_loci29.npy --pop2_geno data/simulation/EUR_n20000_chr22_loci29.npy --runname nt_n2_propt --h1sq 0.1 --h2sq 0.1 --gc 0.7 --n1 100 --n2 100 200 400 --nt 500 10000 --nsnp 2000 --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 --out_dir bench/result --nrep 100 --estimate_omega
# python -m src.simulation.simulation --pop1_geno data/simulation/EAS_n5000_chr22_loci29.npy --pop2_geno data/simulation/EUR_n20000_chr22_loci29.npy --runname h2sq_gc_propt --h1sq 0.1 --h2sq 0.1 0.2 --gc 0.01 0.5 0.9 --n1 100 --n2 200 --nt 500 --nsnp 2000 --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 --out_dir bench/result --nrep 100 --estimate_omega
# python -m src.simulation.simulation --pop1_geno data/simulation/EAS_n5000_chr22_loci29.npy --pop2_geno data/simulation/EUR_n20000_chr22_loci29.npy --runname n1_pcausal_propt --h1sq 0.1 --h2sq 0.1 --gc 0.7 --n1 50 100 --n2 200 --nt 500 --nsnp 2000 --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 0.01 0.02 --out_dir bench/result --nrep 100 --estimate_omega
# python -m src.simulation.simulation --pop1_geno data/simulation/EAS_n5000_chr22_loci29.npy --pop2_geno data/simulation/EUR_n20000_chr22_loci29.npy --runname alpha_h2sq_pcausal_propt --h1sq 0.000000000001 --h2sq 0.1 0.2 --gc 0 --n1 100 --n2 400 --nt 500 --nsnp 2000 --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 0.01 0.02 --out_dir bench/result --nrep 100 --estimate_omega

# python -m src.simulation.simulation --pop1_geno data/simulation/EAS_n5000_chr22_loci29.npy --pop2_geno data/simulation/EUR_n20000_chr22_loci29.npy --runname nt_n2_propt --h1sq 0.3 --h2sq 0.3 --gc 0.7 --n1 100 --n2 400 --nt 1000 --nsnp 2000 --propt 0.01 0.1 0.2 0.5 0.8 --pcausal 0.005 --out_dir bench/result --nrep 20 --estimate_omega

import argparse
import time
import os
import numpy as np
from numba import njit, prange

from ..traceCB.ldsc import Run_cross_LDSC
from ..traceCB.gmm import GMM, GMMtissue
from ..traceCB.utils import make_pd_shrink, z2p

MIN_FLOAT = 1e-32
P_VAL_THRED = 0.05  # for h2 and cov in omega
SNP_START = 1000  # skip first part of chr
MAX_CORR = 0.99


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
        "--estimate_omega",
        action="store_true",
        help="Estimate omega from data, otherwise use true omega",
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
    return parser.parse_args()


get_genotype = lambda geno_file, Nsnp: np.load(geno_file, mmap_mode="r")[
    :, SNP_START : Nsnp + SNP_START
]


@njit(nogil=True, parallel=True)
def cal_ld(G1, G2):
    cor1 = np.corrcoef(G1, rowvar=False)
    cor2 = np.corrcoef(G2, rowvar=False)
    ld1 = np.sum(cor1**2, axis=0)
    ld2 = np.sum(cor2**2, axis=0)
    ldx = np.sum(cor1 * cor2, axis=0)
    return ld1, ld2, ldx


@njit(nogil=True, parallel=True)
def calculate_sumstats(X, y, n):
    nsnp = X.shape[1]
    b_hat = np.zeros(nsnp)
    se_hat = np.zeros(nsnp)

    for j in range(nsnp):
        x_inner = np.dot(X[:, j], X[:, j]) + MIN_FLOAT
        b_hat[j] = np.dot(X[:, j], y) / x_inner
        e = y - X[:, j] * b_hat[j]
        se_hat[j] = np.sqrt(np.dot(e, e) / (n * x_inner))

    return b_hat, se_hat


def generate_data(G1, G2, h1sq, h2sq, gc, n1, n2, nt, nsnp, propt, pcausal):
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
            - causal_ids (np.ndarray): Indices of causal SNPs.
            - pi_mean (float): Mean individual cell type proportion in tissue.
    """
    Omega_causal = np.array(
        [[h1sq, np.sqrt(h1sq * h2sq) * gc], [np.sqrt(h1sq * h2sq) * gc, h2sq]]
    )
    # print("Omega:", Omega_causal / nsnp)
    G1c = G1[:n1, :]
    X1 = (G1c - np.mean(G1c, axis=0)) / (np.std(G1c, axis=0) + MIN_FLOAT)
    G2c = G2[:n2, :]
    X2 = (G2c - np.mean(G2c, axis=0)) / (np.std(G2c, axis=0) + MIN_FLOAT)
    G2t = G2[n2 : n2 + nt, :]
    Xt = (G2t - np.mean(G2t, axis=0)) / (np.std(G2t, axis=0) + MIN_FLOAT)
    # cell type data
    num_causal = int(pcausal * nsnp)
    causal_ids = np.random.choice(np.arange(nsnp), num_causal, replace=False)
    beta_causal = np.random.multivariate_normal(
        mean=np.zeros(2), cov=Omega_causal / (pcausal * nsnp), size=num_causal
    )  # (M, <c11, c12>)
    beta1 = np.zeros(nsnp)
    beta1[causal_ids] = beta_causal[:, 0]
    beta2 = np.zeros(nsnp)
    beta2[causal_ids] = beta_causal[:, 1]
    y1 = X1 @ beta1.T + np.sqrt(1 - h1sq) * np.random.randn(n1)
    y2 = X2 @ beta2.T + np.sqrt(1 - h2sq) * np.random.randn(n2)

    # tissue data
    delta = 5  # control the variance of pi
    pi_ind = np.random.beta(propt * delta, (1 - propt + MIN_FLOAT) * delta, nt)
    pi_mean = np.mean(pi_ind)
    ## define unknown cell type
    beta_unknown = np.zeros(nsnp)
    num_unknown_celltype = 10
    for _ in range(num_unknown_celltype):
        causal_unknown_id = np.random.choice(
            np.arange(nsnp), int(pcausal * nsnp), replace=False
        )
        # causal_unknown_id = causal_ids  #! share causal SNPs
        beta_causal_unknown = np.random.normal(
            loc=0,
            scale=h2sq / (pcausal * nsnp) / num_unknown_celltype,
            size=int(pcausal * nsnp),
        )
        # beta_causal_unknown = (
        #     beta_causal[:, 0] / 3 + beta_causal[:, 1] / 3 + beta_causal_unknown / 3
        # ) #! make unknown cell type correlated with known cell types
        # beta_unknown[causal_unknown_id] += beta_causal_unknown
        beta_unknown[causal_unknown_id] += beta_causal_unknown / 2
        beta_unknown[causal_ids] += (
            beta_causal[:, 1] / 2
        )  #! share effect with known cell type
    yt = (
        pi_ind * (Xt @ beta2.T)
        + (1 - pi_ind) * (Xt @ beta_unknown.T)
        + np.sqrt(1 - (pi_ind**2 + (1 - pi_ind) ** 2) * h2sq) * np.random.randn(nt)
    )

    # sumstats
    b1_hat, se1_hat = calculate_sumstats(X1, y1, n1)
    b2_hat, se2_hat = calculate_sumstats(X2, y2, n2)
    bt_hat, se_t_hat = calculate_sumstats(Xt, yt, nt)
    z1 = b1_hat / se1_hat
    z2 = b2_hat / se2_hat
    zt = bt_hat / se_t_hat
    pval1 = z2p(z1)
    pval2 = z2p(z2)
    pvalt = z2p(zt)
    sig1 = pval1 < P_VAL_THRED
    sig2 = pval2 < P_VAL_THRED
    sigt = pvalt < P_VAL_THRED
    # breakpoint()
    return (
        Omega_causal,
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
    )


@njit(nogil=True)
def re2_meta(betas, ses, converge_tol=1e-6, max_iter=100):
    """
    Calculate RE2 meta-analysis for multiple populations using random-effects model with iterative tau**2 estimation

    :param betas: array of effect sizes (beta coefficients) for each study (n_studies,)
    :param ses: array of standard errors for each study (n_studies,)
    :param converge_tol: convergence tolerance for tau**2 and beta estimates
    :param max_iter: maximum number of iterations
    :return: tuple (beta: combined effect size, se: standard error of combined effect)
    """
    # \widehat{\mu}_{(n+1)} & =\frac{\sum \frac{X_i}{V_i+\hat{\tau}_{(n)}^2}}{\sum \frac{1}{V_i+\hat{\tau}_{(n)}^2}} \\
    # \hat{\tau}_{(n+1)}^2 & =\frac{\sum \frac{\left(X_i-\hat{\mu}_{(n+1)}\right)^2-V_i}{\left(V_i+\hat{\tau}_{(n)}^2\right)^2}}{\sum \frac{1}{\left(V_i+\hat{\tau}_{(n)}^2\right)^2}}
    variances = ses**2

    # Step 1: Initialize
    tau2 = 0.0
    weights = 1.0 / variances
    sum_weights = np.sum(weights)
    combined_beta = np.sum(weights * betas) / sum_weights

    # Initialize variables for iteration
    converged = False

    # Step 2: Iterate to estimate tau
    for i in range(max_iter):
        # Calculate weights for the current tau
        denom = variances + tau2
        weights = 1.0 / denom
        sum_weights = np.sum(weights)

        # Update combined beta using Han & Eskin formula
        new_combined_beta = np.sum(betas * weights) / sum_weights

        # Calculate new tau using Han & Eskin formula
        num_numerator = (betas - new_combined_beta) ** 2 - variances
        num_denom = denom**2
        numerator = np.sum(num_numerator / num_denom)
        denominator = np.sum(1.0 / num_denom)
        new_tau2 = numerator / denominator

        if new_tau2 < 0:  # If tau2 is negative, switch to fixed effects model
            tau2 = 0.0
            weights = 1.0 / variances
            sum_weights = np.sum(weights)
            combined_beta = np.sum(betas * weights) / sum_weights
            converged = True
            break

        beta_diff = abs(new_combined_beta - combined_beta)
        tau_diff = abs(new_tau2 - tau2)

        combined_beta = new_combined_beta
        tau2 = new_tau2

        if max(beta_diff, tau_diff) < converge_tol:
            converged = True
            break

    # Step 3: Calculate final combined beta and standard error
    if tau2 > 0:  # RE2 model
        weights = 1.0 / (variances + tau2)
    else:  # FE model
        weights = 1.0 / variances

    sum_weights = np.sum(weights)
    combined_beta = np.sum(betas * weights) / sum_weights
    combined_se = np.sqrt(1.0 / sum_weights)

    return combined_beta, combined_se


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
    pcausal: proportion of causal SNPs
    out_dir: output directory
    true_omega: true covariance matrix
    id_sim: simulation id
    :return
    None
    """
    simulation_path = f"{runname}/h1sq_{h1sq}_h2sq_{h2sq}_gc_{gc}_n1_{n1}_n2_{n2}_nt_{nt}_nsnp_{nsnp}_propt_{propt}_pcausal_{pcausal}_omega_{true_omega}"
    os.makedirs(os.path.join(out_dir, simulation_path), exist_ok=True)
    simulation_name = f"{simulation_path}/simulation_{id_sim}"
    (
        Omega_causal,
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
    if not true_omega:
        Omega, Omega_se = Run_cross_LDSC(
            b1_hat / se1_hat,
            n1,
            ld1,
            b2_hat / se2_hat,
            n2,
            ld2,
            ldx,
            np.array([1, 1, 0]),
        )
        Omega_p = z2p(Omega / Omega_se)
        if Omega_p[0, 0] >= P_VAL_THRED:
            run_gmm = False  # skip genes with unsignificant tar heritability
        else:
            run_gmm = True
        if Omega_p[0, 1] >= P_VAL_THRED:
            Omega[0, 1] = 0.0
            Omega[1, 0] = 0.0
        if run_gmm:
            ## ! get omega_co, omega_oo
            aux_Omega_matrix, aux_Omega_matrix_se = Run_cross_LDSC(
                b2_hat / se2_hat,
                n2,
                ld2,
                bt_hat / se_t_hat,
                nt,
                ldx,
                ldx,
                np.array([1.0, 1.0, 0.0]),
            )
            ## omega_co: var between target celltype in aux pop and tissue
            # aux_Omega_matrix[0, 1] = propt * omega_2 + (1-propt) * omega_co
            # if z2p(aux_Omega_matrix[0, 1] / aux_Omega_matrix_se[0, 1]) < P_VAL_THRED:
            # else:
            #     omega_co = MIN_FLOAT
            omega_co = (aux_Omega_matrix[0, 1] - propt * Omega[1, 1]) / (1 - propt)
            omega_co_se = (
                aux_Omega_matrix_se[0, 1] ** 2 + (propt**2) * (Omega_se[1, 1] ** 2)
            ) ** 0.5 / (1 - propt)
            omega_co_p = z2p(omega_co / np.maximum(omega_co_se, MIN_FLOAT))
            if omega_co_p >= P_VAL_THRED:
                omega_co = 0.0
            ## omega_oo: var of non-target celltype in tissue
            ## aux_Omega_matrix[1, 1] = propt^2 * omega_2 + 2*propt*(1-propt)*omega_co + (1-propt)^2 * omega_oo
            omega_oo = (
                aux_Omega_matrix[1, 1]
                - propt**2 * Omega[1, 1]
                - 2 * propt * (1 - propt) * omega_co
            ) / ((1 - propt) ** 2)
            omega_oo = np.maximum(omega_oo, 1e-12)

            omega_co_cor = omega_co / (Omega[1, 1] ** 0.5 * omega_oo**0.5 + 1e-12)
            if omega_co_cor > MAX_CORR:
                omega_co = MAX_CORR * (Omega[1, 1] ** 0.5 * omega_oo**0.5 + 1e-12)
            elif omega_co_cor < -MAX_CORR:
                omega_co = -MAX_CORR * (Omega[1, 1] ** 0.5 * omega_oo**0.5 + 1e-12)
                # if z2p(aux_Omega_matrix[1, 1] / aux_Omega_matrix_se[1, 1]) < P_VAL_THRED:
                # else:
                #     omega_oo = 1e-12
            # omega_oo_se = (aux_Omega_matrix_se[1, 1]**2 + (propt**4) * (Omega_se[1, 1]**2) + (2*propt*(1-propt))**2 * (omega_co_se**2))**0.5 / ((1-propt)**2)
            # omega_oo_p = z2p(omega_oo / omega_oo_se)
            # if omega_oo_p > P_VAL_THRED:
            #     omega_oo = MIN_FLOAT

            ## ! get omega_xo
            tar_tissue_Omega_matrix, tar_tissue_Omega_matrix_se = Run_cross_LDSC(
                b1_hat / se1_hat,
                n1,
                ld1,
                b2_hat / se2_hat,
                n2,
                ld2,
                ldx,
                np.array([1.0, 1.0, 0.0]),
            )
            ## omega_xo: var between target celltype in target pop and non-target celltype in aux tissue
            ## tar_tissue_Omega_matrix[0, 1] = propt * omega_x + (1-propt) * omega_xo
            # if z2p(tar_tissue_Omega_matrix[0, 1] / tar_tissue_Omega_matrix_se[0, 1]) < P_VAL_THRED:
            # else:
            #     omega_xo = MIN_FLOAT
            omega_xo = (tar_tissue_Omega_matrix[0, 1] - propt * Omega[0, 1]) / (
                1 - propt
            )
            omega_xo_se = (
                tar_tissue_Omega_matrix_se[0, 1] ** 2
                + (propt**2) * (Omega_se[0, 1] ** 2)
            ) ** 0.5 / (1 - propt)
            omega_xo_p = z2p(omega_xo / np.maximum(omega_xo_se, MIN_FLOAT))
            if omega_xo_p > P_VAL_THRED:
                omega_xo = MIN_FLOAT
            omega_xo_cor = omega_xo / (Omega[0, 0] ** 0.5 * omega_oo**0.5 + 1e-12)
            if omega_xo_cor > MAX_CORR:
                omega_xo = MAX_CORR * (Omega[0, 0] ** 0.5 * omega_oo**0.5 + 1e-12)
            elif omega_xo_cor < -MAX_CORR:
                omega_xo = -MAX_CORR * (Omega[0, 0] ** 0.5 * omega_oo**0.5 + 1e-12)

            ## construct OmegaCB for GMMtissue
            OmegaCB = np.array(
                [
                    [Omega[0, 0], Omega[0, 1], omega_xo],
                    [Omega[0, 1], Omega[1, 1], omega_co],
                    [omega_xo, omega_co, omega_oo],
                ]
            )
            # if np.all(Omega > MIN_FLOAT):
            #     print("True Omega:\n", Omega_causal / nsnp)
            #     print("Estimated Omega:\n", Omega)
            #     print("Estimated OmegaCB:\n", OmegaCB)
            #     print(f"xo_cor: {omega_xo_cor:.4f}, co_cor: {omega_co_cor:.4f}")
            #     print("=" * 10, "\n")

    else:  # true omega
        raise NotImplementedError("True omega not implemented yet")
        # Omega = Omega_causal / nsnp
    # GMM
    pop1_beta = np.zeros((nsnp, 3))  # sumstat, cross, tissue
    pop2_beta = np.zeros((nsnp, 3))
    pop1_se = np.zeros((nsnp, 3))
    pop2_se = np.zeros((nsnp, 3))
    pop1_beta[:, 0] = b1_hat
    pop1_se[:, 0] = se1_hat
    pop2_beta[:, 0] = b2_hat
    pop2_se[:, 0] = se2_hat

    if run_gmm:
        for j in prange(nsnp):
            pop1_beta[j, 1], pop1_se[j, 1], pop2_beta[j, 1], pop2_se[j, 1] = GMM(
                Omega,
                np.eye(2),
                b1_hat[j],
                se1_hat[j],
                ld1[j],
                b2_hat[j],
                se2_hat[j],
                ld2[j],
                ldx[j],
            )
            pop1_beta[j, 2], pop1_se[j, 2], pop2_beta[j, 2], pop2_se[j, 2] = GMMtissue(
                OmegaCB,
                np.eye(3),
                b1_hat[j],
                se1_hat[j],
                ld1[j],
                b2_hat[j],
                se2_hat[j],
                ld2[j],
                ldx[j],
                bt_hat[j],
                se_t_hat[j],
                propt,
            )
        pop1_z = pop1_beta / pop1_se
        pop2_z = pop2_beta / pop2_se
    else:  # no run gmm, keep sumstat results
        pop1_z = np.zeros((nsnp, 3))
        pop2_z = np.zeros((nsnp, 3))
        pop1_z[:, 0] = b1_hat / se1_hat
        pop1_z[:, 1] = b1_hat / se1_hat
        pop1_z[:, 2] = b1_hat / se1_hat
        pop2_z[:, 0] = b2_hat / se2_hat
        pop2_z[:, 1] = b2_hat / se2_hat
        pop2_z[:, 2] = b2_hat / se2_hat

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
    all_results[:, 4:7] = pop1_z
    all_results[:, 7:10] = pop2_z
    all_results[:, 10] = bt_hat / se_t_hat
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
