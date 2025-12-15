# Useful functions
import numpy as np
from numba import njit
from scipy.stats import norm


## convert p-value to z-score and vice versa
p2z = lambda p: np.abs(norm.ppf(p / 2))
z2p = lambda z: 2 * norm.sf(abs(z))

## Constants
MIN_FLOAT = np.finfo(float).eps
MIN_HERITABILITY = 1e-12
MAX_CORR = 0.99
P_VAL_THRED = 0.05
P_VAL_THRED_Z = p2z(P_VAL_THRED)  # convert p-value threshold to z-score threshold
eSNP_THRESHOLD = 1e-5

## Function to check if a matrix is positive definite
is_pd = lambda x: np.all(np.linalg.eigvals(x) > -MIN_FLOAT)


def is_pd(x):
    """Check if a matrix is positive definite"""
    eigns = np.linalg.eigvals(x)
    return np.all(eigns > -MIN_FLOAT)


@njit
def is_pd_numba(x):
    """Check if a matrix is positive definite using Numba"""
    eigns = np.linalg.eigvals(x)
    for val in eigns:
        if val <= -MIN_FLOAT:
            return False
    return True


## Force a matrix to be positive definite
def make_pd_shrink(Omega, shrink=0.9):
    """shrink off-diagonal elements towards 0
    return Omega
    !NOTE: if the diagonal elements of Omega are negative, return a matrix with MIN_FLOAT on the diagonal.
    """
    if np.any(np.diag(Omega) < 0):
        return np.eye(Omega.shape[0]) * MIN_HERITABILITY
    if is_pd(Omega):
        return Omega
    diag_Omega = np.eye(Omega.shape[0]) * Omega
    off_Omega = Omega - diag_Omega
    j = 0
    while not is_pd(Omega):
        # shrink off-diagonal elements towards 0
        off_Omega *= shrink
        Omega = diag_Omega + off_Omega
        j += 1
        if j == 100:
            return np.eye(Omega.shape[0]) * MIN_HERITABILITY
    return Omega


@njit
def make_pd_shrink_numba(Omega, shrink=0.9):
    """Check if a matrix is positive definite using Numba
    !NOTE: if the diagonal elements of Omega are negative, return a matrix with MIN_FLOAT on the diagonal.
    """
    if np.any(np.diag(Omega) < 0):
        return np.eye(Omega.shape[0]) * MIN_HERITABILITY
    if is_pd_numba(Omega):
        return Omega
    diag_Omega = np.eye(Omega.shape[0]) * Omega
    off_Omega = Omega - diag_Omega
    j = 0
    while not is_pd_numba(Omega):
        # shrink off-diagonal elements towards 0
        off_Omega *= shrink
        Omega = diag_Omega + off_Omega
        j += 1
        if j == 100:
            return np.eye(Omega.shape[0]) * MIN_HERITABILITY
    return Omega
