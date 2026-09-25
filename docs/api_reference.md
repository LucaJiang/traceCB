# API Reference

The functions below operate on already aligned summary statistics. Use the
same SNP order, allele orientation, and genome build across populations and LD
scores. Pass finite NumPy arrays; sample sizes and standard errors must be
positive. These functions do not perform allele harmonization or file loading.

```python
from traceCB import GMM, GMMtissue, Run_Cross_LDSC
from traceCB.ldsc import Run_Single_LDSC
```

The GMM routines use Numba and compile on their first call. The command-line
runner is available as `python -m traceCB.run_gmm --help`; full-data preparation
is described in the [pipeline guide](pipeline.md).

- [API Reference](#api-reference)
  - [traceCB.gmm](#tracecbgmm)
    - [`GMM`](#gmm)
    - [`GMMtissue`](#gmmtissue)
  - [traceCB.ldsc](#tracecbldsc)
    - [`Run_Single_LDSC`](#run_single_ldsc)
    - [`Run_Cross_LDSC`](#run_cross_ldsc)


## traceCB.gmm

Core functions for the Generalized Method of Moments (GMM) estimation.

### `GMM`

Apply cross-population GMM (without tissue-specific information) to estimate effect sizes.

```python
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
) -> tuple[float, float, float, float]
```

**Parameters**

- **Omega** (`np.ndarray`): A (2, 2) per-SNP covariance matrix.
- **C** (`np.ndarray`): A (2, 2) sampling-error scaling matrix, for example from LDSC intercepts.
- **beta1** (`float`): Effect size (beta) for the SNP in population 1.
- **se1** (`float`): Standard error for the SNP in population 1.
- **ld1** (`float`): LD score between the SNP and the rest of the SNPs in the target gene in population 1.
- **beta2** (`float`): Effect size (beta) for the SNP in population 2.
- **se2** (`float`): Standard error for the SNP in population 2.
- **ld2** (`float`): LD score between the SNP and the rest of the SNPs in the target gene in population 2.
- **ldx** (`float`): Cross-population LD score for the SNP between population 1 and population 2.

**Returns**

- **beta1_blue** (`float`): GMM estimate (BLUE) for population 1.
- **se1_blue** (`float`): Standard error of the GMM estimate for population 1.
- **beta2_blue** (`float`): GMM estimate (BLUE) for population 2.
- **se2_blue** (`float`): Standard error of the GMM estimate for population 2.

---

### `GMMtissue`

Apply cross-population GMM using an additional bulk-tissue eQTL estimate from population 2.

```python
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
) -> tuple[float, float, float, float]
```

**Parameters**

- **Omega** (`np.ndarray`): A (2, 2) per-SNP covariance matrix.
- **C** (`np.ndarray`): A (3, 3) sampling-error scaling matrix, for example from LDSC intercepts.
- **beta1** (`float`): Effect size (beta) for the SNP in population 1.
- **se1** (`float`): Standard error for the SNP in population 1.
- **ld1** (`float`): LD score for population 1.
- **beta2** (`float`): Effect size (beta) for the SNP in population 2.
- **se2** (`float`): Standard error for the SNP in population 2.
- **ld2** (`float`): LD score for population 2.
- **ldx** (`float`): Cross-population LD score.
- **beta_t** (`float`): Bulk-tissue effect size for the SNP in population 2.
- **se_t** (`float`): Standard error of the bulk-tissue effect size.
- **pi2_omega_o** (`float`): Mixture-proportion-weighted per-SNP variance contribution of other cell types. For a two-cell mixture this is `(1 − propt)² × Var(beta_other)`. The function multiplies it by `ld2`.
- **propt** (`float`): Fraction of the focal cell type in population 2 tissue, between 0 and 1. Convert percentages to fractions first.

**Returns**

- **beta1_blue** (`float`): GMM estimate for population 1.
- **se1_blue** (`float`): Standard error for population 1.
- **beta2_blue** (`float`): GMM estimate for population 2.
- **se2_blue** (`float`): Standard error for population 2.

---

## traceCB.ldsc

Functions for running Single and Cross-Population LD Score Regression (LDSC).

### `Run_Single_LDSC`

Estimate the per-SNP genetic variance coefficient from single-population LDSC.

```python
def Run_Single_LDSC(
    zscore: np.ndarray,
    n: np.ndarray,
    ldscore: np.ndarray,
    intercept: float = np.nan,
) -> tuple[float, float]
```

**Parameters**

- **zscore** (`np.ndarray`): Array of Z-scores for SNPs, shape `(num_SNP,)`.
- **n** (`np.ndarray`): Array of sample sizes for each SNP, shape `(num_SNP,)`.
- **ldscore** (`np.ndarray`): Array of LD scores, shape `(num_SNP,)`.
- **intercept** (`float`, optional): Fixed intercept value. If `np.nan` (default), the intercept is estimated from the data.

**Returns**

- **h2** (`float`): Estimated per-SNP variance coefficient, bounded below by `MIN_HERITABILITY = 1e-12`. This is the slope in `E[z²] = intercept + n × LD × h2`; it is not summed over the SNPs in the locus.
- **h2_se** (`float`): Standard error of the per-SNP variance estimate.

---

### `Run_Cross_LDSC`

Estimate the per-SNP genetic covariance matrix Ω between two populations.

```python
def Run_Cross_LDSC(
    zscore1: np.ndarray,
    n1: np.ndarray,
    ldscore1: np.ndarray,
    zscore2: np.ndarray,
    n2: np.ndarray,
    ldscore2: np.ndarray,
    crossld: np.ndarray,
    intercept: np.ndarray = np.array([np.nan, np.nan, np.nan]),
) -> tuple[np.ndarray, np.ndarray]
```

**Parameters**

- **zscore1** (`np.ndarray`): Z-scores for population 1.
- **n1** (`np.ndarray`): Sample sizes for population 1.
- **ldscore1** (`np.ndarray`): LD scores for population 1.
- **zscore2** (`np.ndarray`): Z-scores for population 2.
- **n2** (`np.ndarray`): Sample sizes for population 2.
- **ldscore2** (`np.ndarray`): LD scores for population 2.
- **crossld** (`np.ndarray`): Cross-population LD scores.
- **intercept** (`np.ndarray`, optional): Array of intercept values `[I1, I2, Ix]`. Default is `[nan, nan, nan]`, which estimates all intercepts independently for each call. `[1.0, 1.0, 0.0]` fixes the two within-population intercepts to 1 and the cross-population intercept to 0. The input array is not modified.

**Returns**

- **Omega** (`np.ndarray`): Estimated per-SNP genetic covariance matrix of shape `(2, 2)`.
  - `Omega[0, 0]`: Per-SNP genetic variance in population 1
  - `Omega[1, 1]`: Per-SNP genetic variance in population 2
  - `Omega[0, 1]` / `Omega[1, 0]`: Genetic covariance
- **Omega_se** (`np.ndarray`): Standard error matrix for `Omega`, shape `(2, 2)`.


The diagonal elements of Ω are bounded below by `1e-12`; the off-diagonal
covariance is not clipped by `Run_Cross_LDSC`. Before calling GMM, check that
Ω and the LD-weighted covariance matrices are suitable for the intended model.
The tutorial and full-data runner apply covariance significance checks and
correlation clipping. Singular inputs can raise a numerical linear algebra
error.

For `GMM`, rows and columns of `C` correspond to populations 1 and 2. For
`GMMtissue`, they correspond to population 1, population 2, and bulk tissue.
The sampling covariance used by the estimator is `diag(SE) @ C @ diag(SE)`.
An identity matrix represents independent sampling errors with unit LDSC
intercepts. This matrix describes sampling error, separately from Ω.
