"""Shared helpers for simulation drivers."""

import hashlib

import numpy as np


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


def make_sim_seed(base_seed, *parts):
    """Create a stable per-replicate seed from data-generating settings."""
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


def unknown_cell_effect_scale(h2sq, num_causal, num_unknown_celltype):
    if num_causal <= 0:
        return 0.0
    return np.sqrt(h2sq / num_causal / num_unknown_celltype)


def calculate_pi2_omega_sum_const(proportion):
    return 2 * proportion


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
