# Apply GMM (Generalized Method of Moments)
# GMM, GMMtissue, GMMtissueTarget, GMMtissueBoth examples are provided in the main function
from numba import njit
import numpy as np


@njit(nogil=True)
def GMM(
    Omega: np.ndarray,
    C: np.ndarray,
    beta1: float,
    se1: float,
    ld1: float,
    beta2: float,
    se2: float,
    ld2: float,
    ldx: float,
):
    """Apply cross population GMM (without tissue) to j th SNP

    Parameters
    ----------
    Omega : np.ndarray
        (2, 2) per-snp covariance matrix
    C : np.ndarray
        (2, 2) genetic drift matrix estimated from the LDSC
    beta1 : float
        beta for snp j in population 1
    se1 : float
        standard error for snp j in population 1
    ld1 : float
        LD between snp j and the rest of the snps of target gene in population 1
    beta2 : float
        beta for snp j in population 2
    se2 : float
        standard error for snp j in population 2
    ld2 : float
        LD between snp j and the rest of the snps of target gene in population 2
    ldx : float
        LD of snp j between population 1 and population 2

    Returns
    -------
    beta1_blue, se1_blue, beta2_blue, se2_blue : float
        GMM estimates for population 1 and population 2
    """
    A = np.eye(2)
    Betas = np.array([beta1, beta2]).reshape(-1, 1)
    # weight omega by ld
    Omegaj = np.array(
        [[Omega[0, 0] * ld1, Omega[0, 1] * ldx], [Omega[0, 1] * ldx, Omega[1, 1] * ld2]]
    )
    ## check if Omega is positive definite
    # Omegaj = make_pd_shrink(Omegaj)
    # weight C by se
    Sj = np.array([[se1, 0], [0, se2]])
    SCSj = Sj @ C @ Sj
    # lambda
    lambda1 = np.array([1, Omegaj[0, 1] / Omegaj[0, 0]]).reshape(-1, 1)
    Lambda1_inv = A @ Omegaj @ A.T + SCSj - Omegaj[0, 0] * lambda1 @ lambda1.T
    Lambda1 = np.linalg.inv(Lambda1_inv)
    lambda2 = np.array([Omegaj[0, 1] / Omegaj[1, 1], 1]).reshape(-1, 1)
    Lambda2_inv = A @ Omegaj @ A.T + SCSj - Omegaj[1, 1] * lambda2 @ lambda2.T
    Lambda2 = np.linalg.inv(Lambda2_inv)
    # blue estimates
    var1_blue = 1 / (lambda1.T @ Lambda1 @ lambda1).item()
    beta1_blue = (var1_blue * lambda1.T @ Lambda1 @ Betas).item()
    var2_blue = 1 / (lambda2.T @ Lambda2 @ lambda2).item()
    beta2_blue = (var2_blue * lambda2.T @ Lambda2 @ Betas).item()
    if var1_blue <= 0:
        beta1_blue = beta1
        se1_blue = se1
    else:
        se1_blue = np.sqrt(var1_blue)
    if var2_blue <= 0:
        beta2_blue = beta2
        se2_blue = se2
    else:
        se2_blue = np.sqrt(var2_blue)

    return beta1_blue, se1_blue, beta2_blue, se2_blue


@njit(nogil=True)
def GMMtissue(
    Omega: np.ndarray,
    C: np.ndarray,
    beta1: float,
    se1: float,
    ld1: float,
    beta2: float,
    se2: float,
    ld2: float,
    ldx: float,
    beta_t: float,
    se_t: float,
    pi2_omega_o: float,
    propt: float,
):
    """Apply cross population GMM (with tissue) to j th SNP
    Parameters
    ----------
    Omega : np.ndarray
        (2, 2) per-snp covariance matrix
    C : np.ndarray
        (3, 3) genetic drift matrix estimated from the LDSC
    beta1 : float
        beta for snp j in population 1
    se1 : float
        standard error for snp j in population 1
    ld1 : float
        LD between snp j and the rest of the snps of target gene in population 1
    beta2 : float
        beta for snp j in population 2
    se2 : float
        standard error for snp j in population 2
    ld2 : float
        LD between snp j and the rest of the snps of target gene in population 2
    ldx : float
        LD of snp j between population 1 and population 2
    beta_t : float
        beta for snp j in tissue
    se_t : float
        standard error for snp j in tissue
    pi2_omega_o : float
        mixture-proportion-weighted variance contribution from other cell types
        in tissue, e.g. (1 - pi)^2 * Var(beta_other) in a two-cell mixture
    propt : float
        proportion of cell type in tissue

    Returns
    -------
    beta1_blue, se1_blue, beta2_blue, se2_blue: float
        GMM estimates for population 1 and population 2
    """
    Betas = np.array([beta1, beta2, beta_t]).reshape(-1, 1)
    A = np.array([[1.0, 0.0], [0.0, 1.0], [0.0, propt]])
    # weight omega by ld
    Omegaj = np.array(
        [[Omega[0, 0] * ld1, Omega[0, 1] * ldx], [Omega[0, 1] * ldx, Omega[1, 1] * ld2]]
    )
    # weight C by se
    Sj = np.array([[se1, 0, 0], [0, se2, 0], [0, 0, se_t]])
    SCSj = Sj @ C @ Sj
    # lambda
    lambda1 = (
        np.array([Omegaj[0, 0], Omegaj[0, 1], propt * Omegaj[0, 1]]).reshape(-1, 1)
        / Omegaj[0, 0]
    )
    Lambda1_inv = A @ Omegaj @ A.T + SCSj - Omegaj[0, 0] * lambda1 @ lambda1.T
    Lambda1_inv[2, 2] += pi2_omega_o * ld2  # add weighted other-cell contribution
    Lambda1 = np.linalg.inv(Lambda1_inv)
    lambda2 = np.array([Omegaj[1, 0] / Omegaj[1, 1], 1, propt]).reshape(-1, 1)
    Lambda2_inv = A @ Omegaj @ A.T + SCSj - Omegaj[1, 1] * lambda2 @ lambda2.T
    Lambda2_inv[2, 2] += pi2_omega_o * ld2  # add weighted other-cell contribution
    Lambda2 = np.linalg.inv(Lambda2_inv)
    # blue estimates
    var1_blue = 1 / (lambda1.T @ Lambda1 @ lambda1).item()
    beta1_blue = (var1_blue * lambda1.T @ Lambda1 @ Betas).item()
    var2_blue = 1 / (lambda2.T @ Lambda2 @ lambda2).item()
    beta2_blue = (var2_blue * lambda2.T @ Lambda2 @ Betas).item()
    if var1_blue <= 0:
        beta1_blue = beta1
        se1_blue = se1
    else:
        se1_blue = np.sqrt(var1_blue)
    if var2_blue <= 0:
        beta2_blue = beta2
        se2_blue = se2
    else:
        se2_blue = np.sqrt(var2_blue)
    return beta1_blue, se1_blue, beta2_blue, se2_blue


@njit(nogil=True)
def GMMtissueBoth(
    Omega: np.ndarray,
    C: np.ndarray,
    beta1: float,
    se1: float,
    ld1: float,
    beta_t1: float,
    se_t1: float,
    beta2: float,
    se2: float,
    ld2: float,
    ldx: float,
    beta_t2: float,
    se_t2: float,
    pi2_omega_o1: float,
    pi2_omega_o2: float,
    pi12_omega_o12: float,
    propt1: float,
    propt2: float,
):
    """Apply cross population GMM with tissue in both populations to j th SNP

    Parameters
    ----------
    Omega : np.ndarray
        (2, 2) per-snp covariance matrix
    C : np.ndarray
        (4, 4) genetic drift matrix estimated from the LDSC
    beta1 : float
        beta for snp j in population 1
    se1 : float
        standard error for snp j in population 1
    ld1 : float
        LD between snp j and the rest of the snps of target gene in population 1
    beta_t1 : float
        beta for snp j in tissue from population 1
    se_t1 : float
        standard error for snp j in tissue from population 1
    beta2 : float
        beta for snp j in population 2
    se2 : float
        standard error for snp j in population 2
    ld2 : float
        LD between snp j and the rest of the snps of target gene in population 2
    ldx : float
        LD of snp j between population 1 and population 2
    beta_t2 : float
        beta for snp j in tissue from population 2
    se_t2 : float
        standard error for snp j in tissue from population 2
    pi2_omega_o1 : float
        mixture-proportion-weighted variance contribution from other cell types
        in tissue from population 1
    pi2_omega_o2 : float
        mixture-proportion-weighted variance contribution from other cell types
        in tissue from population 2
    pi12_omega_o12 : float
        mixture-proportion-weighted covariance contribution from other cell types
        between the two tissue studies
    propt1 : float
        proportion of cell type in tissue from population 1
    propt2 : float
        proportion of cell type in tissue from population 2

    Returns
    -------
    beta1_blue, se1_blue, beta2_blue, se2_blue: float
        GMM estimates for population 1 and population 2
    """
    Betas = np.array([beta1, beta_t1, beta2, beta_t2]).reshape(-1, 1)
    A = np.array([[1.0, 0.0], [propt1, 0.0], [0.0, 1.0], [0.0, propt2]])
    # weight omega by ld
    Omegaj = np.array(
        [[Omega[0, 0] * ld1, Omega[0, 1] * ldx], [Omega[0, 1] * ldx, Omega[1, 1] * ld2]]
    )
    # weight C by se
    Sj = np.array([[se1, 0, 0, 0], [0, se_t1, 0, 0], [0, 0, se2, 0], [0, 0, 0, se_t2]])
    SCSj = Sj @ C @ Sj
    Sigma_oj = np.zeros((4, 4))
    Sigma_oj[1, 1] = pi2_omega_o1 * ld1
    Sigma_oj[3, 3] = pi2_omega_o2 * ld2
    Sigma_oj[1, 3] = pi12_omega_o12 * ldx
    Sigma_oj[3, 1] = pi12_omega_o12 * ldx
    # lambda
    lambda1 = (
        np.array(
            [Omegaj[0, 0], propt1 * Omegaj[0, 0], Omegaj[0, 1], propt2 * Omegaj[0, 1]]
        ).reshape(-1, 1)
        / Omegaj[0, 0]
    )
    Lambda1_inv = (
        A @ Omegaj @ A.T + Sigma_oj + SCSj - Omegaj[0, 0] * lambda1 @ lambda1.T
    )
    Lambda1 = np.linalg.inv(Lambda1_inv)
    lambda2 = (
        np.array(
            [Omegaj[1, 0], propt1 * Omegaj[1, 0], Omegaj[1, 1], propt2 * Omegaj[1, 1]]
        ).reshape(-1, 1)
        / Omegaj[1, 1]
    )
    Lambda2_inv = (
        A @ Omegaj @ A.T + Sigma_oj + SCSj - Omegaj[1, 1] * lambda2 @ lambda2.T
    )
    Lambda2 = np.linalg.inv(Lambda2_inv)
    # blue estimates
    var1_blue = 1 / (lambda1.T @ Lambda1 @ lambda1).item()
    beta1_blue = (var1_blue * lambda1.T @ Lambda1 @ Betas).item()
    var2_blue = 1 / (lambda2.T @ Lambda2 @ lambda2).item()
    beta2_blue = (var2_blue * lambda2.T @ Lambda2 @ Betas).item()
    if var1_blue <= 0:
        beta1_blue = beta1
        se1_blue = se1
    else:
        se1_blue = np.sqrt(var1_blue)
    if var2_blue <= 0:
        beta2_blue = beta2
        se2_blue = se2
    else:
        se2_blue = np.sqrt(var2_blue)
    return beta1_blue, se1_blue, beta2_blue, se2_blue


if __name__ == "__main__":
    # test data for illustration
    h1sq = 7e-3
    h2sq = 5e-3
    cor = 0.8
    Omega = np.array([[h1sq, cor], [cor, h2sq]])  # covariance matrix
    # pi2_omega_o = 0.0001  # weighted other-cell contribution in tissue
    C2 = np.eye(2)  # genetic drift matrix
    C3 = np.eye(3)  # genetic drift matrix for GMMtissue
    C4 = np.eye(4)  # genetic drift matrix for GMMtissueBoth
    propt = 0.4  # proportion of cell type in tissue
    # summary statistics
    beta1 = -0.51
    se1 = 0.05
    beta2 = -0.31
    se2 = 0.10
    beta_t = -0.77
    se_t = 0.04
    beta_t1 = -0.58
    se_t1 = 0.06
    beta_t2 = -0.77
    se_t2 = 0.04
    ld1 = 1.168  # LD between snp j and the rest of the snps of target gene in population 1
    ld2 = 2.356  # LD between snp j and the rest of the snps of target gene in population 2
    ldx = 1.334  # LD of snp j between population 1 and population 2
    pi2_omega_o = 0.001  # constant for weighted other-cell contribution in tissue
    pi2_omega_o1 = (
        0.001  # constant for weighted other-cell contribution in population 1 tissue
    )
    pi2_omega_o2 = (
        0.001  # constant for weighted other-cell contribution in population 2 tissue
    )
    pi12_omega_o12 = (
        0.0001  # constant for weighted other-cell covariance across tissues
    )
    propt1 = 0.3  # proportion of cell type in population 1 tissue
    propt2 = 0.4  # proportion of cell type in population 2 tissue
    # example usage
    beta1_blue, se1_blue, beta2_blue, se2_blue = GMM(
        Omega=Omega,
        C=C2,
        beta1=beta1,
        se1=se1,
        ld1=ld1,
        beta2=beta2,
        se2=se2,
        ld2=ld2,
        ldx=ldx,
    )
    print("GMM estimates:")
    print(f"beta1_blue = {beta1_blue:.6f}, se1_blue = {se1_blue:.6f}")
    print(f"beta2_blue = {beta2_blue:.6f}, se2_blue = {se2_blue:.6f}")

    beta1_blue, se1_blue, beta2_blue, se2_blue = GMMtissue(
        Omega=Omega,
        C=C3,
        beta1=beta1,
        se1=se1,
        ld1=ld1,
        beta2=beta2,
        se2=se2,
        ld2=ld2,
        ldx=ldx,
        beta_t=beta_t,
        se_t=se_t,
        pi2_omega_o=pi2_omega_o,
        propt=propt,
    )
    print("\nGMMtissue estimates:")
    print(f"beta1_blue = {beta1_blue:.6f}, se1_blue = {se1_blue:.6f}")
    print(f"beta2_blue = {beta2_blue:.6f}, se2_blue = {se2_blue:.6f}")

    beta1_blue, se1_blue, beta2_blue, se2_blue = GMMtissueBoth(
        Omega=Omega,
        C=C4,
        beta1=beta1,
        se1=se1,
        ld1=ld1,
        beta_t1=beta_t1,
        se_t1=se_t1,
        beta2=beta2,
        se2=se2,
        ld2=ld2,
        ldx=ldx,
        beta_t2=beta_t2,
        se_t2=se_t2,
        pi2_omega_o1=pi2_omega_o1,
        pi2_omega_o2=pi2_omega_o2,
        pi12_omega_o12=pi12_omega_o12,
        propt1=propt1,
        propt2=propt2,
    )
    print("\nGMMtissueBoth estimates:")
    print(f"beta1_blue = {beta1_blue:.6f}, se1_blue = {se1_blue:.6f}")
    print(f"beta2_blue = {beta2_blue:.6f}, se2_blue = {se2_blue:.6f}")

# GMM estimates:
# beta1_blue = -0.733634, se1_blue = 0.060011
# beta2_blue = -2.021197, se2_blue = 0.257334

# GMMtissue estimates:
# beta1_blue = -0.732058, se1_blue = 0.060011
# beta2_blue = -1.951153, se2_blue = 0.134176

# GMMtissueBoth estimates:
# beta1_blue = -0.805281, se1_blue = 0.058068
# beta2_blue = -1.929212, se2_blue = 0.134127
