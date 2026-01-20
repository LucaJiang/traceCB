# API Reference

This page provides detailed documentation for the core functions in the `traceCB` package.

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
- **C** (`np.ndarray`): A (2, 2) genetic drift matrix estimated from LDSC.
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

Apply cross-population GMM incorporating tissue-specific information (e.g., from scRNA-seq derived eQTLs).

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
- **C** (`np.ndarray`): A (3, 3) genetic drift matrix estimated from LDSC.
- **beta1** (`float`): Effect size (beta) for the SNP in population 1.
- **se1** (`float`): Standard error for the SNP in population 1.
- **ld1** (`float`): LD score for population 1.
- **beta2** (`float`): Effect size (beta) for the SNP in population 2.
- **se2** (`float`): Standard error for the SNP in population 2.
- **ld2** (`float`): LD score for population 2.
- **ldx** (`float`): Cross-population LD score.
- **beta_t** (`float`): Effect size (beta) for the SNP in the specific tissue/cell-type.
- **se_t** (`float`): Standard error for the SNP in the specific tissue/cell-type.
- **pi2_omega_o** (`float`): Variance component of other cell types within the tissue.
- **propt** (`float`): Proportion of the specific cell type in the tissue.

**Returns**

- **beta1_blue** (`float`): GMM estimate for population 1.
- **se1_blue** (`float`): Standard error for population 1.
- **beta2_blue** (`float`): GMM estimate for population 2.
- **se2_blue** (`float`): Standard error for population 2.

---

## traceCB.ldsc

Functions for running Single and Cross-Population LD Score Regression (LDSC).

### `Run_Single_LDSC`

Run LDSC regression to estimate heritability ($h^2$) in a single population.

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

- **h2** (`float`): Estimated heritability. (Clipped to be $\ge$ MIN_HERITABILITY).
- **h2_se** (`float`): Standard error of the heritability estimate.

---

### `Run_Cross_LDSC`

Run LDSC regression to estimate genetic covariance ($\Omega$) between two populations.

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
- **intercept** (`np.ndarray`, optional): Array of intercept values `[h11, h22, h12]`. Default is `[nan, nan, nan]`, which estimates all intercepts.

**Returns**

- **Omega** (`np.ndarray`): Estimated genetic covariance matrix of shape `(2, 2)`.
  - `Omega[0, 0]`: Heritability pop 1
  - `Omega[1, 1]`: Heritability pop 2
  - `Omega[0, 1]` / `Omega[1, 0]`: Genetic covariance
- **Omega_se** (`np.ndarray`): Standard error matrix for `Omega`, shape `(2, 2)`.
