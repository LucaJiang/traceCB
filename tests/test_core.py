from __future__ import annotations

import numpy as np
import pytest

import traceCB
from traceCB.gmm import GMM, GMMtissue
from traceCB.run_gmm import clip_correlation
from traceCB.utils import is_pd, make_pd_shrink


def test_public_version_matches_release() -> None:
    assert traceCB.__version__ == "1.0"


def test_gmm_estimators_return_finite_standard_errors() -> None:
    omega = np.array([[0.20, 0.05], [0.05, 0.25]])
    gmm_result = GMM(
        omega,
        np.eye(2),
        0.20,
        0.10,
        1.20,
        0.15,
        0.12,
        1.10,
        0.80,
    )
    tissue_result = GMMtissue(
        omega,
        np.eye(3),
        0.20,
        0.10,
        1.20,
        0.15,
        0.12,
        1.10,
        0.80,
        0.10,
        0.08,
        0.02,
        0.30,
    )

    assert np.isfinite(gmm_result).all()
    assert np.isfinite(tissue_result).all()
    assert gmm_result[1] > 0 and gmm_result[3] > 0
    assert tissue_result[1] > 0 and tissue_result[3] > 0


def test_matrix_and_correlation_safety_bounds() -> None:
    repaired = make_pd_shrink(np.array([[1.0, 2.0], [2.0, 1.0]]))
    assert is_pd(repaired)

    covariance, correlation = clip_correlation(1.0, 1.0, 2.0)
    assert correlation == pytest.approx(0.99)
    assert covariance == pytest.approx(0.99)
