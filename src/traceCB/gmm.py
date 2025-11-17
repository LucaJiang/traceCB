# Apply GMM (Generalized Method of Moments)
# GMM, GMMtissue examples are provided in the main function
from numba import njit
import numpy as np

from .utils import make_pd_shrink_numba as make_pd_shrink


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
    se1_blue = np.sqrt(var1_blue)
    beta1_blue = (var1_blue * lambda1.T @ Lambda1 @ Betas).item()
    var2_blue = 1 / (lambda2.T @ Lambda2 @ lambda2).item()
    se2_blue = np.sqrt(var2_blue)
    beta2_blue = (var2_blue * lambda2.T @ Lambda2 @ Betas).item()

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
    sigma_o : float
        variance of other cell types in tissue
    propt : float
        proportion of cell type in tissue

    Returns
    -------
    beta1_blue, se1_blue, beta2_blue, se2_blue: float
        GMM estimates for population 1 and population 2
    """
    # breakpoint()
    Betas = np.array([beta1, beta2, beta_t]).reshape(-1, 1)
    A = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, propt, 1.0 - propt]])
    # weight omega by ld
    Omegaj = np.array(
        [
            [Omega[0, 0] * ld1, Omega[0, 1] * ldx, Omega[0, 2] * ldx],
            [Omega[0, 1] * ldx, Omega[1, 1] * ld2, Omega[2, 1] * ld2],
            [Omega[0, 2] * ldx, Omega[2, 1] * ld2, Omega[2, 2] * ld2],
        ]
    )

    ## ! check if Omega is positive definite
    # Omegaj = make_pd_shrink(Omegaj)
    # if not np.all(np.linalg.eigvals(Omegaj) > 0):
    #     print("Warning: Omega matrix is not positive definite in GMMtissue.")
    #     raise ValueError("Omega matrix is not positive definite in GMMtissue.")
    # weight C by se
    Sj = np.array([[se1, 0, 0], [0, se2, 0], [0, 0, se_t]])
    SCSj = Sj @ C @ Sj
    # lambda
    lambda1 = (
        np.array(
            [
                Omegaj[0, 0],
                Omegaj[0, 1],
                propt * Omegaj[0, 1] + (1 - propt) * Omegaj[0, 2],
            ]
        ).reshape(-1, 1)
        / Omegaj[0, 0]
    )
    Lambda1_inv = A @ Omegaj @ A.T + SCSj - Omegaj[0, 0] * lambda1 @ lambda1.T
    Lambda1 = np.linalg.inv(Lambda1_inv)
    lambda2 = np.array(
        [
            Omegaj[1, 0] / Omegaj[1, 1],
            1.0,
            propt + (1 - propt) * Omegaj[1, 2] / Omegaj[1, 1],
        ]
    ).reshape(-1, 1)
    Lambda2_inv = A @ Omegaj @ A.T + SCSj - Omegaj[1, 1] * lambda2 @ lambda2.T
    Lambda2 = np.linalg.inv(Lambda2_inv)
    # blue estimates
    var1_blue = 1 / (lambda1.T @ Lambda1 @ lambda1).item()
    se1_blue = np.sqrt(var1_blue)
    beta1_blue = (var1_blue * lambda1.T @ Lambda1 @ Betas).item()
    var2_blue = 1 / (lambda2.T @ Lambda2 @ lambda2).item()
    se2_blue = np.sqrt(var2_blue)
    beta2_blue = (var2_blue * lambda2.T @ Lambda2 @ Betas).item()
    return beta1_blue, se1_blue, beta2_blue, se2_blue


if __name__ == "__main__":
    # test data for illustration
    h1sq = 0.00005005371452165598
    h2sq = 0.00004256252586151709
    cor = 0.000022333797435697717
    Omega = np.array([[h1sq, cor], [cor, h2sq]])  # covariance matrix
    # sigma_o = 0.0001  # variance of other cell types in tissue
    C2 = np.eye(2)  # genetic drift matrix
    C3 = np.eye(3)  # genetic drift matrix for GMMtissue
    propt = 0.4  # proportion of cell type in tissue
    # summary statistics
    beta1 = 0.153536154422526
    se1 = 0.1645668866928277
    beta2 = -0.0279196
    se2 = 0.0982573
    beta_t = -0.477724
    se_t = 0.0352819
    ld1 = 0.0160532668232917  # LD between snp j and the rest of the snps of target gene in population 1
    ld2 = 2.35622501373291  # LD between snp j and the rest of the snps of target gene in population 2
    ldx = 1.3348826169967651  # LD of snp j between population 1 and population 2
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
        # sigma_o=sigma_o,
        propt=propt,
    )
    print(f"beta1_blue = {beta1_blue:.6f}, se1_blue = {se1_blue:.6f}")
    print(f"beta2_blue = {beta2_blue:.6f}, se2_blue = {se2_blue:.6f}")

# beta1_blue = -0.002158, se1_blue = 0.009367
# beta2_blue = -0.023266, se2_blue = 0.098134
# beta1_blue = -0.061843, se1_blue = 0.006401
# beta2_blue = -0.648565, se2_blue = 0.066993
