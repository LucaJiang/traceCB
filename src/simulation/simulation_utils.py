"""Shared helpers for simulation drivers."""

import hashlib

import numpy as np

MISSING_GENOTYPE_VALUE = 3


def impute_missing_genotype(
    genotype: np.ndarray, missing_value: int | float = MISSING_GENOTYPE_VALUE
) -> np.ndarray:
    """Return float genotype dosage with missing values imputed by SNP mean."""
    geno = np.asarray(genotype, dtype=np.float32).copy()
    valid = np.isfinite(geno) & (geno != missing_value)
    if np.all(valid):
        return geno

    counts = np.sum(valid, axis=0)
    sums = np.sum(np.where(valid, geno, 0.0), axis=0, dtype=np.float64)
    col_mean = np.divide(
        sums,
        counts,
        out=np.zeros(geno.shape[1], dtype=np.float64),
        where=counts > 0,
    ).astype(np.float32)
    missing = ~valid
    geno[missing] = np.take(col_mean, np.where(missing)[1])
    return geno


def load_genotype_window(
    geno_file: str,
    nsnp: int,
    snp_start: int,
    missing_value: int | float = MISSING_GENOTYPE_VALUE,
) -> np.ndarray:
    """Load a SNP window from a simulation .npy genotype file."""
    geno = np.load(geno_file, mmap_mode="r")[:, snp_start : nsnp + snp_start]
    return impute_missing_genotype(geno, missing_value=missing_value)


def standardize_genotype(genotype: np.ndarray, min_float: float) -> np.ndarray:
    x = np.asarray(genotype, dtype=float)
    return (x - np.mean(x, axis=0)) / (np.std(x, axis=0) + min_float)


def flatten_float_seq(x):
    """Scalar, list, tuple, or array -> flat list of floats."""
    return [float(v) for v in np.asarray(x, dtype=float).ravel()]


def list_arg(value):
    return value if isinstance(value, list) else [value]


def validate_unit_interval(name, values):
    for value in values:
        if not np.isfinite(value) or value < 0 or value > 1:
            raise ValueError(f"{name} values must be in [0, 1], got {value}")


def validate_nonnegative(name, values):
    for value in values:
        if not np.isfinite(value) or value < 0:
            raise ValueError(f"{name} values must be nonnegative, got {value}")


def tail_panel_start(
    n_rows, n_singlecell, n_tissue, panel_name="tissue", allow_overlap=False
):
    """Return the start row for a tail tissue panel, optionally requiring no overlap."""
    n_rows = int(n_rows)
    n_singlecell = int(n_singlecell)
    n_tissue = int(n_tissue)
    if n_singlecell < 0 or n_tissue < 0:
        raise ValueError(f"{panel_name} sample sizes must be nonnegative")
    if n_singlecell > n_rows:
        raise ValueError(
            f"single-cell n={n_singlecell} exceeds {panel_name} rows={n_rows}"
        )
    if n_tissue > n_rows:
        raise ValueError(f"tissue n={n_tissue} exceeds {panel_name} rows={n_rows}")
    tissue_start = n_rows - n_tissue
    if not allow_overlap and tissue_start < n_singlecell:
        raise ValueError(
            f"{panel_name} tail tissue panel overlaps prefix single-cell panel: "
            f"n_singlecell={n_singlecell}, n_tissue={n_tissue}, rows={n_rows}"
        )
    return tissue_start


def make_sim_seed(base_seed, *parts):
    """Create a stable component seed from a base seed and global replicate parts."""
    tokens = []
    for part in parts:
        if part is None:
            tokens.append("None")
        elif isinstance(part, (float, np.floating)):
            tokens.append(f"{part:.17g}")
        else:
            tokens.append(str(part))
    seed_key = "|".join([str(int(base_seed)), *tokens])
    return int.from_bytes(
        hashlib.blake2s(seed_key.encode("utf-8"), digest_size=4).digest(),
        "little",
    )


def seed_random_component(seed_base, seed_parts, component):
    """Reset NumPy's legacy RNG for one simulation random component."""
    if seed_base is None:
        return
    np.random.seed(make_sim_seed(seed_base, component, *seed_parts))


def unknown_cell_effect_scale(h2sq, num_causal, num_unknown_celltype):
    if num_causal <= 0:
        return 0.0
    return np.sqrt(h2sq / num_causal / num_unknown_celltype)


def calculate_pi2_omega_sum_const(proportion):
    return 2 * proportion  # only for two cell types scenario


def perturb_gmm_propt(propt, mode="exact", normal_var=0.01, scale=0.1):
    """Return the cell-type proportion supplied to GMM tissue."""
    if mode == "exact":
        gmm_propt = propt
    elif mode == "underestimate":
        gmm_propt = propt - scale
    elif mode == "overestimate":
        gmm_propt = propt + scale
    elif mode == "normal":
        gmm_propt = np.random.normal(propt, np.sqrt(float(normal_var)))
    else:
        raise ValueError(f"Unknown GMM propt mode: {mode}")
    return float(np.clip(gmm_propt, 0, 1))


def iter_gmm_propt_subsettings(mode, mode_scales, normal_vars):
    """Yield (sub_idx, mode_scale, normal_var) for one GMM propt mode."""
    if mode in ("underestimate", "overestimate"):
        for sub_idx, sc in enumerate(mode_scales):
            yield sub_idx, float(sc), float(normal_vars[0])
    elif mode == "normal":
        for sub_idx, nv in enumerate(normal_vars):
            yield sub_idx, 0.0, float(nv)
    else:
        yield 0, 0.0, float(normal_vars[0])


def sanitize_ld_scores(*ld_scores):
    """Replace non-finite LD-score entries caused by constant genotype columns."""
    return tuple(
        np.nan_to_num(
            np.asarray(ld_score, dtype=float),
            nan=0.0,
            posinf=0.0,
            neginf=0.0,
        )
        for ld_score in ld_scores
    )
