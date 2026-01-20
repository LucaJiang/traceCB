# Functions for running single and cross population LD score regression (LDSC)
import numpy as np
from traceCB.utils import MIN_HERITABILITY, MIN_FLOAT

# Set Threshold
MIN_WEIGHT = 1e-12  # Minimum weight for IRWLS


############################################
#      Functions for IRWLS in LDSC         #
############################################
def get_coef_raw(
    ldscore: np.ndarray,
    sqrt_nprod: np.ndarray,
    zprod: np.ndarray,
    max_intercept: float = 1.0,
) -> tuple[float, float]:
    """Obtain coefficient for heritability / genetic covariance

    Parameters
    ----------
    ldscore : np.array
        ldscore, (N, 1)
    sqrt_nprod : np.array
        sqrt product of n1 and n2, (N, 1)
    zprod : np.array
        product of z1 and z2, (N, 1)
    max_intercept : float, optional
        maximum intercept of regression, default 1

    Returns
    -------
    tuple[float, float]
        A tuple containing (intercept, coef)
    """
    # get averages
    mean_ldscore = np.mean(ldscore)
    mean_nprod = np.mean(sqrt_nprod)
    mean_zprod = np.mean(zprod)
    # estimate intercept
    zprod_sort = np.sort(zprod)
    idx = int(len(zprod) * 0.95)
    mean_val = np.mean(zprod_sort[0:idx])
    intercept = mean_val if mean_val < max_intercept else max_intercept
    # get raw coefficients
    coef = (mean_zprod - intercept) / (mean_ldscore * mean_nprod + MIN_FLOAT)
    return float(intercept), float(coef)


def get_pred(
    coef: float, ldscore: np.ndarray, sqrt_nprod: np.ndarray, intercept: float
) -> np.ndarray:
    """Get prediction from LDSC regression.

    Parameters
    ----------
    coef : float
        coefficient of regression
    ldscore : np.ndarray
        LD scores, shape (N,)
    sqrt_nprod : np.ndarray
        sqrt root of n1*n2, shape (N,)
    intercept : float
        intercept of regression

    Returns
    -------
    np.ndarray
        prediction values, shape (N,)
    """
    return coef * ldscore * sqrt_nprod + intercept


def regression_jk(
    x: np.ndarray,
    y: np.ndarray,
    intercept: float = np.nan,  # fixed intercept if not NaN (R: constrain_intercept with `subtract`)
    estimate_se: bool = True,
    jackknife: bool = True,
    nblocks: int = 200,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Least squares with optional fixed intercept and block jackknife SE
    replicating ldsc.R::regression().

    x: (N, p) design matrix (first column should be the intercept column if present)
    y: (N,) response
    intercept : fixed intercept value; if NaN, intercept is estimated
    estimate_se : whether to compute SEs
    jackknife : whether to use block jackknife
    nblocks : number of blocks for jackknife
    """
    y = y.reshape(-1)
    n, p = x.shape

    # precompute cross-products
    xtx = x.T @ x
    xty = x.T @ y

    # coefficients
    if not np.isnan(intercept):  # constrained intercept
        coefs = np.zeros(p)
        # solve for slopes using the submatrix (drop intercept column/row)
        if p == 1:
            # degenerate case: only intercept column present
            coefs[0] = intercept
        else:
            coefs[1:] = np.linalg.solve(xtx[1:, 1:], xty[1:])
            coefs[0] = intercept
    else:
        coefs = np.linalg.solve(xtx, xty)

    if not estimate_se:
        return coefs, np.full(p, np.nan)

    if not jackknife:
        # fallback: classic OLS SEs
        k = p if np.isnan(intercept) else (p - 1)  # number of estimated params
        resid = y - x @ coefs
        sigsq = (resid @ resid) / max(n - k, 1)
        cov_beta = sigsq * np.linalg.inv(xtx if np.isnan(intercept) else xtx)
        return coefs, np.sqrt(np.maximum(np.diag(cov_beta), 0.0))

    # --- block jackknife ---
    edges = np.floor(np.linspace(0, n, num=nblocks + 1)).astype(int)
    # drop empty blocks (can happen when N < nblocks)
    blocks = [
        (edges[i], edges[i + 1])
        for i in range(len(edges) - 1)
        if edges[i + 1] > edges[i]
    ]
    B = len(blocks)
    if B < 2:
        # not enough blocks to jackknife — return NaNs (or OLS as fallback if prefer)
        return coefs, np.full(p, np.nan)

    coefs_jk = np.zeros((B, p))
    for b, (lo, hi) in enumerate(blocks):
        # leave-one-block-out cross-products
        xtx_blk = x[lo:hi, :].T @ x[lo:hi, :]
        xty_blk = x[lo:hi, :].T @ y[lo:hi]
        xtx_loo = xtx - xtx_blk
        xty_loo = xty - xty_blk

        # solve on the reduced system
        if not np.isnan(intercept):  # constrained intercept: only (p-1) unknowns
            coef_b = np.zeros(p)
            if p > 1:
                coef_b[1:] = np.linalg.solve(xtx_loo[1:, 1:], xty_loo[1:])
            coef_b[0] = intercept
        else:
            coef_b = np.linalg.solve(xtx_loo, xty_loo)
        coefs_jk[b, :] = coef_b

    # np.cov with rowvar=False, ddof=1 == sample covariance (divide by B-1)
    jk_cov = np.cov(coefs_jk, rowvar=False, ddof=1) * (B - 1)
    jk_se = np.sqrt(np.maximum(np.diag(jk_cov), 0.0))

    return coefs, jk_se


def regression(
    x: np.ndarray,
    y: np.ndarray,
    intercept: float = np.nan,
    estimate_se: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Perform least square regression.

    Parameters
    ----------
    x : np.ndarray
        Design matrix, shape (N, 2)
    y : np.ndarray
        Response variable, shape (N,)
    intercept : float, optional
        Fixed value for intercept. If NaN, intercept is estimated.
    estimate_se : bool, optional
        Whether to estimate standard errors, default True

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Coefficients and standard errors. If estimate_se is False, the second element
        is NaN array.
    """
    # perform regression
    xtx = x.T @ x
    xty = x.T @ y
    n, p = x.shape  # p should be 2 in this study
    coefs = np.zeros((p, 1))
    if not np.isnan(intercept):
        coefs[1] = xty[1] / xtx[1, 1]
        coefs[0] = intercept
        k = p - 1
    else:
        coefs = np.linalg.solve(xtx, xty)
        k = p
    if not estimate_se:
        return coefs, np.full(p, np.nan)

    ## calculate standard errors
    err = y - x @ coefs
    sigsq = err.T @ err / (n - k)
    cov_beta = sigsq * np.linalg.inv(xtx)
    coefs_se = np.sqrt(np.diag(cov_beta))
    return coefs, coefs_se


def get_coef(
    ldscore: np.ndarray,
    sqrt_nprod: np.ndarray,
    zprod: np.ndarray,
    weight: np.ndarray,
    intercept: float = np.nan,
    estimate_se: bool = True,
) -> np.ndarray:
    """Obtain coefficient for heritability / genetic covariance

    Parameters
    ----------
    ldscore : np.ndarray
        LD scores, shape (N,)
    sqrt_nprod : np.ndarray
        Sample size values (n1, n2, or sqrt(n1 * n2)), shape (N,)
    zprod : np.ndarray
        Product of Z-scores, shape (N,)
    weight : np.ndarray
        Regression weights, shape (N,)
    intercept : float, optional
        Fixed value for intercept. If NaN, intercept is estimated.
    estimate_se : bool, optional
        Whether to return standard error, default True

    Returns
    -------
    np.ndarray
        If estimate_se=True: [intercept, coef, intercept_se, coef_se]
        If estimate_se=False: [intercept, coef, NaN, NaN]
    """
    zprod_adj = zprod.copy()
    if not np.isnan(intercept):
        zprod_adj -= intercept
    zprod_adj = zprod_adj.reshape(-1, 1)

    # Rescale to improve matrix conditioning
    nbar = np.mean(sqrt_nprod)
    ld_rescale = ldscore * sqrt_nprod / nbar
    ld_rescale_add_column = np.ones((len(ld_rescale), 2))
    ld_rescale_add_column[:, 1] = ld_rescale
    weight = np.maximum(weight, MIN_WEIGHT)
    weight_sqrt = np.sqrt(weight).reshape(-1, 1)

    # Weighted regression
    coefs, coefs_se = regression(
        ld_rescale_add_column * weight_sqrt,
        zprod_adj * weight_sqrt,
        intercept=intercept,
        estimate_se=estimate_se,
    )

    # Rescale coefficient back
    coefs[1] /= nbar
    coefs_se[1] /= nbar  # ok for nan
    return np.append(coefs.flatten(), coefs_se.flatten())


def update_weights(weights: np.ndarray, preds: np.ndarray) -> np.ndarray:
    """Update weights by variance of prediction.

    Parameters
    ----------
    weights : np.ndarray
        Weights array, shape (N, 3)
    preds : np.ndarray
        Predictions array, shape (N, 3)

    Returns
    -------
    np.ndarray
        Updated weights, shape (N, 3)
    """
    var_all = np.zeros((len(preds), 3))
    var_all[:, :2] = 2.0 * np.square(preds[:, :2])
    var_all[:, 2] = preds[:, 0] * preds[:, 1] + np.square(
        preds[:, 2]
    )  #! test，its wrong
    return weights / var_all


##################################################
#     LDSC Regression for Cross Population       #
##################################################
def estimate_gc(
    zscore: np.ndarray,
    n: np.ndarray,
    ldscore: np.ndarray,
    crossld: np.ndarray,
    filter_idx: np.ndarray,
    reg_w1: np.ndarray = np.array([1.0]),
    reg_w2: np.ndarray = np.array([1.0]),
    intercept: np.ndarray = np.array([np.nan, np.nan, np.nan]),
    estimate_se: bool = True,
) -> np.ndarray:
    """Estimate genetic covariance from LDSC regression.

    Parameters
    ----------
    zscore : np.ndarray
        Z-score matrix, shape (num_SNP, 2)
    n : np.ndarray
        Sample size matrix, shape (num_SNP, 2)
    ldscore : np.ndarray
        LD Score matrix, shape (num_SNP, 2)
    crossld : np.ndarray
        Cross population LD Score, shape (num_SNP,)
    filter_idx : np.ndarray
        Boolean array indicating SNPs with Z-score**2 < 30, if use filtered data.
        Else, all one array.
    reg_w1 : np.ndarray, optional
        Regularization weight for population 1, default 1.0
    reg_w2 : np.ndarray, optional
        Regularization weight for population 2, default 1.0
    intercept : np.ndarray, optional
        Intercepts [h11, h22, h12], default NaN array
    estimate_se : bool, optional
        Whether to estimate standard errors, default True

    Returns
    -------
    np.ndarray
        Coefficients array, always shape (3, 4) with [pop1,
        pop2, popx] x [intercept, coef, intercept_se, coef_se].
        If estimate_se=False, standard errors are set to NaN.
    """
    num_snp = len(zscore)
    # Prepare data arrays
    z1 = zscore[:, 0]
    z2 = zscore[:, 1]
    n1 = n[:, 0]
    n2 = n[:, 1]

    sqrt_nprods = np.zeros((num_snp, 3))
    sqrt_nprods[:, 0] = n1
    sqrt_nprods[:, 1] = n2
    sqrt_nprods[:, 2] = np.sqrt(n1 * n2)

    zprods = np.zeros((num_snp, 3))
    zprods[:, 0] = z1**2
    zprods[:, 1] = z2**2
    zprods[:, 2] = z1 * z2

    ldscores = np.zeros((num_snp, 3))
    ldscores[:, 0] = ldscore[:, 0]
    ldscores[:, 1] = ldscore[:, 1]
    ldscores[:, 2] = crossld

    # Select data based on filter_idx
    sqrt_nprods_use = sqrt_nprods[filter_idx]
    zprods_use = zprods[filter_idx]
    ldscores_use = ldscores[filter_idx]
    # Get initial predictions for weights
    preds = np.zeros_like(ldscores_use)
    max_ints = np.array([1.0, 1.0, 0.0])
    for i in range(3):
        int_i, tau_i = get_coef_raw(
            ldscores_use[:, i], sqrt_nprods_use[:, i], zprods_use[:, i], max_ints[i]
        )
        preds[:, i] = get_pred(tau_i, ldscores_use[:, i], sqrt_nprods_use[:, i], int_i)
    # Update weights
    weights = np.zeros((len(preds), 3))
    weights[:, 0] = np.maximum(reg_w1[filter_idx], MIN_WEIGHT)
    weights[:, 1] = np.maximum(reg_w2[filter_idx], MIN_WEIGHT)
    weights[:, 2] = np.sqrt(weights[:, 0] * weights[:, 1])
    new_weights = update_weights(weights, preds)
    # Get final estimates
    taus = np.zeros((3, 4))
    for i in range(3):
        result = get_coef(
            ldscores_use[:, i],
            sqrt_nprods_use[:, i],
            zprods_use[:, i],
            new_weights[:, i],
            intercept[i],
            estimate_se,
        )
        taus[i, :] = result

    return taus


def ldscore_regression_gc(
    zscore: np.ndarray,
    n: np.ndarray,
    ldscore: np.ndarray,
    crossld: np.ndarray,
    filter_idx: np.ndarray,
    intercept: np.ndarray = np.array([np.nan, np.nan, np.nan]),
    estimate_se: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """LDSC regression for genetic covariance between two populations.

    Parameters
    ----------
    zscore : np.ndarray
        Z-score matrix of 2 pop, shape (num_SNP, 2)
    n : np.ndarray
        Sample size matrix of 2 pop, shape (num_SNP, 2)
    ldscore : np.ndarray
        LD Score matrix of 2 pop, shape (num_SNP, 2)
    crossld : np.ndarray
        Cross population LD Score, shape (num_SNP,)
    filter_idx : np.ndarray
        Boolean array indicating SNPs with Z-score**2 < 30
    intercept : np.ndarray, optional
        Intercept values [I1, I2, Ix], default NaN array
    estimate_se : bool, optional
        Whether to estimate standard errors, default True

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        Arrays for population 1, population 2, and cross-population results,
        each with shape (4,) containing [intercept, coef, intercept_se, coef_se]
    """
    reg_w1 = 1.0 / np.maximum(ldscore[:, 0], 1.0)
    reg_w2 = 1.0 / np.maximum(ldscore[:, 1], 1.0)

    # If any intercept is NaN, estimate intercept
    if np.any(np.isnan(intercept)):
        # coefs shape (3, 4) for [pop1, pop2, popx]
        # x [intercept, coef, intercept_se, coef_se]
        coefs = estimate_gc(
            zscore,
            n,
            ldscore,
            crossld,
            filter_idx,
            reg_w1,
            reg_w2,
            intercept,
            False,
        )
        for i in range(3):
            if np.isnan(intercept[i]):
                intercept[i] = coefs[i, 0]

    # Estimate genetic covariance with full data
    coefs = estimate_gc(
        zscore,
        n,
        ldscore,
        crossld,
        np.full(len(zscore), True, dtype=np.bool_),
        reg_w1,
        reg_w2,
        intercept,
        estimate_se,
    )
    return coefs[0], coefs[1], coefs[2]  # shape (4, ) for pop1, pop2, cross-population


###################################################
#         LDSC Regression Main Function           #
###################################################
def Run_Single_LDSC(
    zscore: np.ndarray,
    n: np.ndarray,
    ldscore: np.ndarray,
    intercept: float = np.nan,
) -> tuple[float, float]:
    """Run LDSC regression for heritability in a single population.

    Parameters
    ----------
    zscore : np.ndarray
        Z-score of SNPs, shape (num_SNP,)
    n : np.ndarray
        Sample size of Z-score, shape (num_SNP,)
    ldscore : np.ndarray
        LD Score of SNPs, shape (num_SNP,)
    intercept : float, optional
        Intercept value, default NaN

    Returns
    -------
    tuple[float, float]
        Heritability and its standard error.
    """
    # num_snp = len(zscore)

    # Filter SNPs with Z-score squared > 30
    filter_idx = zscore**2 < 30

    # Prepare data
    zprod = zscore**2
    sqrt_nprod = n

    # Get intercept if not provided
    if np.isnan(intercept):
        # Get initial prediction for weights
        max_int = 1.0
        int_init, tau_init = get_coef_raw(
            ldscore[filter_idx], sqrt_nprod[filter_idx], zprod[filter_idx], max_int
        )
        pred = get_pred(tau_init, ldscore[filter_idx], sqrt_nprod[filter_idx], int_init)

        # Calculate initial weights
        reg_w = np.maximum(1.0 / ldscore[filter_idx], MIN_WEIGHT)
        var = 2.0 * np.square(pred)
        weights = reg_w / var
        # Initial estimation of intercept
        result_init = get_coef(
            ldscore[filter_idx],
            sqrt_nprod[filter_idx],
            zprod[filter_idx],
            weights,
            intercept=np.nan,
            estimate_se=False,
        )
        intercept = result_init[0]
    # Final estimation with all data
    int_all, tau_all = get_coef_raw(ldscore, sqrt_nprod, zprod, 1.0)
    pred_all = get_pred(tau_all, ldscore, sqrt_nprod, int_all)
    reg_w_all = np.maximum(1.0 / ldscore, MIN_WEIGHT)
    var_all = 2.0 * np.square(pred_all)
    weights = reg_w_all / var_all
    result = get_coef(
        ldscore,
        sqrt_nprod,
        zprod,
        weights,
        intercept=intercept,
        estimate_se=True,
    )

    h2 = max(result[1], MIN_HERITABILITY)
    h2_se = result[3]
    return h2, h2_se


def Run_Cross_LDSC(
    zscore1: np.ndarray,
    n1: np.ndarray,
    ldscore1: np.ndarray,
    zscore2: np.ndarray,
    n2: np.ndarray,
    ldscore2: np.ndarray,
    crossld: np.ndarray,
    intercept: np.ndarray = np.array([np.nan, np.nan, np.nan]),
) -> tuple[np.ndarray, np.ndarray]:
    """Run LDSC regression for genetic covariance between two populations.
    Note: heritability have been clipped to >0, while genetic correlation maintains original estimates.

    Parameters
    ----------
    zscore1 : np.ndarray
        Z-score of SNPs in population 1, shape (num_SNP,)
    n1 : np.ndarray
        Sample size of Z-score in population 1, shape (num_SNP,)
    ldscore1 : np.ndarray
        LD Score of SNPs in population 1, shape (num_SNP,)
    zscore2 : np.ndarray
        Z-score of SNPs in population 2, shape (num_SNP,)
    n2 : np.ndarray
        Sample size of Z-score in population 2, shape (num_SNP,)
    ldscore2 : np.ndarray
        LD Score of SNPs in population 2, shape (num_SNP,)
    crossld : np.ndarray
        Cross population LD Score, shape (num_SNP,)
    intercept : np.ndarray, optional
        Intercept values [h11, h22, h12], default NaN array

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Genetic covariance matrix Omega and its standard error Omega_se,
        both shape (2, 2) indicating [[pop1, popx], [popx, pop2]].
    """
    num_snp = len(zscore1)
    zscore = np.zeros((num_snp, 2))
    zscore[:, 0] = zscore1
    zscore[:, 1] = zscore2

    n = np.zeros((num_snp, 2))
    n[:, 0] = n1
    n[:, 1] = n2

    ldscore = np.zeros((num_snp, 2))
    ldscore[:, 0] = ldscore1
    ldscore[:, 1] = ldscore2

    # Filter SNPs with Z-score squared > 30
    z1_filter = zscore[:, 0] ** 2 < 30
    z2_filter = zscore[:, 1] ** 2 < 30
    filter_idx = z1_filter & z2_filter

    pop1, pop2, cross = ldscore_regression_gc(
        zscore, n, ldscore, crossld, filter_idx, intercept, estimate_se=True
    )  # each shape (4,) [intercept, coef, intercept_se, coef_se]

    # Format results
    Omega = np.zeros((2, 2))
    Omega_se = np.zeros((2, 2))
    Omega[0, 0] = pop1[1]
    Omega_se[0, 0] = pop1[3]
    Omega[1, 1] = pop2[1]
    Omega_se[1, 1] = pop2[3]
    Omega[0, 1] = cross[1]
    Omega_se[0, 1] = cross[3]
    Omega[1, 0] = cross[1]
    Omega_se[1, 0] = cross[3]

    # Clip to suitable range, optional
    # Omega = np.maximum(Omega, MIN_HERITABILITY)
    Omega[0, 0] = max(Omega[0, 0], MIN_HERITABILITY)
    Omega[1, 1] = max(Omega[1, 1], MIN_HERITABILITY)
    # thred = np.sqrt(Omega[0, 0] * Omega[1, 1] + MIN_FLOAT) * MAX_CORR
    # Omega[0, 1] = np.minimum(np.maximum(Omega[0, 1], -thred), thred)
    # Omega[1, 0] = Omega[0, 1]
    return Omega, Omega_se


if __name__ == "__main__":
    # Minimum example usage
    ## Attention: inputs MUST be numpy array
    zscore1 = np.array([-1.2, 5, -0.3])
    n1 = np.array([1000, 1500, 1200])
    ldscore1 = np.array([0.2, 0.3, 0.4])
    zscore2 = np.array([-1.5, 4, -0.6])
    n2 = np.array([1100, 1400, 1300])
    ldscore2 = np.array([0.25, 0.35, 0.45])
    crossld = np.array([0.15, 0.25, 0.35])
    intercept = np.array([np.nan, np.nan, np.nan])

    h2, h2_se = Run_Single_LDSC(zscore1, n1, ldscore1, np.nan)
    print("Single-pop h2, h2_se:", h2, h2_se)

    print("\nCross-pop LDSC:")
    Omega, Omega_se = Run_Cross_LDSC(
        zscore1, n1, ldscore1, zscore2, n2, ldscore2, crossld, intercept
    )
    print("Omega:", Omega)
    print("Omega_se:", Omega_se)

    print("\nCross-pop LDSC with fixed intercept:")
    intercept_fixed = np.array([1.0, 1.0, 0.0])
    Omega_fixed, Omega_se_fixed = Run_Cross_LDSC(
        zscore1, n1, ldscore1, zscore2, n2, ldscore2, crossld, intercept_fixed
    )
    print("Omega with fixed intercept:", Omega_fixed)
    print("Omega_se with fixed intercept:", Omega_se_fixed)
