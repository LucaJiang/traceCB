import math
import sys
from pathlib import Path

import numpy as np

SIMULATION_DIR = Path(__file__).resolve().parents[1] / "src" / "simulation"
EXPERIMENTS_DIR = SIMULATION_DIR / "experiments"
for path in (SIMULATION_DIR, EXPERIMENTS_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import simulation
import simulate_masked_omega as simulation_masked_omega
import simulation_utils
import simulate_robustness as simulation_robustness
import simulate_tracecb2 as simulation_tracecb2


def test_unknown_cell_effect_scale_uses_standard_deviation():
    scale = simulation.unknown_cell_effect_scale(
        h2sq=0.1, num_causal=10, num_unknown_celltype=1
    )

    assert math.isclose(scale, math.sqrt(0.1 / 10), rel_tol=0, abs_tol=1e-12)


def test_true_weighted_omega_o_uses_data_generating_proportion():
    omega_cb = np.array(
        [
            [0.2, 0.01, 0.03],
            [0.01, 0.4, 0.12],
            [0.03, 0.12, 0.2],
        ]
    )

    weighted_omega_o = simulation_robustness.calculate_true_weighted_omega_o(
        omega_cb, propt=0.25
    )

    expected = 0.2 - 0.25**2 * 0.4 - 2 * 0.25 * (0.12 - 0.25 * 0.4)
    assert math.isclose(weighted_omega_o, expected, rel_tol=0, abs_tol=1e-12)


def test_pi2_omega_sum_const_matches_simulation_formula():
    assert math.isclose(
        simulation_robustness.calculate_pi2_omega_sum_const(0.25),
        2 * 0.25 + 0 / (1 - 0.25),
        rel_tol=0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        simulation_robustness.calculate_pi2_omega_sum_const(1.0),
        2.0,
        rel_tol=0,
        abs_tol=1e-12,
    )


def test_gmm_propt_perturbation_is_clipped_to_unit_interval():
    assert simulation_robustness.perturb_gmm_propt(
        0.8, mode="overestimate", scale=0.5
    ) == 1.0
    assert simulation_robustness.perturb_gmm_propt(
        0.1, mode="underestimate", scale=0.5
    ) == 0.0


def test_unit_interval_validation_rejects_invalid_propt():
    try:
        simulation_robustness.validate_unit_interval("--propt", [-0.1])
    except ValueError as exc:
        assert "--propt" in str(exc)
    else:
        raise AssertionError("Expected invalid propt to raise ValueError")


def test_sim_seed_changes_with_replicate_id():
    base_parts = (
        20260525,
        1e-12,
        0.1,
        0.0,
        100,
        400,
        5000,
        2000,
        0.01,
        0.005,
        None,
        None,
        "none",
        0.8,
    )

    seeds = [
        simulation_robustness.make_sim_seed(*base_parts, replicate_id)
        for replicate_id in range(100)
    ]

    assert len(set(seeds)) == len(seeds)


def test_estimated_weighted_omega_o_matches_ldsc_formula():
    omega = np.array([[0.1, 0.02], [0.02, 0.4]])
    aux_omega = np.array([[0.4, 0.13], [0.13, 0.22]])

    weighted_omega_o = simulation_robustness.calculate_estimated_weighted_omega_o(
        aux_omega, omega, gmm_propt=0.25
    )

    expected = 0.22 - 0.25**2 * 0.4 - 2 * 0.25 * max(0.13 - 0.25 * 0.4, 0)
    assert math.isclose(weighted_omega_o, expected, rel_tol=0, abs_tol=1e-12)


def test_base_generate_data_allows_zero_pcausal():
    rng = np.random.default_rng(123)
    nsnp = 12
    n1 = 8
    n2 = 8
    nt = 8
    G1 = rng.normal(size=(n1, nsnp))
    G2 = rng.normal(size=(n2 + nt, nsnp))

    np.random.seed(456)
    result = simulation.generate_data(
        G1,
        G2,
        h1sq=0.1,
        h2sq=0.1,
        gc=0.0,
        n1=n1,
        n2=n2,
        nt=nt,
        nsnp=nsnp,
        propt=0.2,
        pcausal=0.0,
    )

    omega, b1_hat, se1_hat, b2_hat, se2_hat, bt_hat, se_t_hat = result[:7]
    causal_ids = result[10]
    assert causal_ids.size == 0
    for arr in (omega, b1_hat, se1_hat, b2_hat, se2_hat, bt_hat, se_t_hat):
        assert np.all(np.isfinite(arr))


def test_zero_pcausal_overlap_and_partition_helpers_return_zero_effects():
    beta1, beta2, causal_ids1, causal_ids2 = (
        simulation_robustness.generate_causal_effects_by_overlap(
            nsnp=20,
            pcausal=0.0,
            h1sq=0.1,
            h2sq=0.1,
            causal_overlap=0.5,
        )
    )
    assert causal_ids1.size == 0
    assert causal_ids2.size == 0
    assert np.all(beta1 == 0)
    assert np.all(beta2 == 0)

    beta1, beta2, causal_ids, pop2_a_ids, region_a = (
        simulation_robustness.generate_partitioned_causal_effects(
            nsnp=20,
            pcausal=0.0,
            h1sq=0.1,
            h2sq=0.1,
            gc=0.0,
            null_region_prop=0.5,
        )
    )
    assert causal_ids.size == 0
    assert pop2_a_ids.size == 0
    assert np.all(beta1 == 0)
    assert np.all(beta2 == 0)
    assert region_a.dtype == bool


def test_tracecb2_unknown_effects_allow_zero_pcausal():
    beta_unknown = simulation_tracecb2.sample_unknown_effects(
        nsnp=20,
        pcausal=0.0,
        hsq=0.1,
    )

    assert np.all(beta_unknown == 0)


def test_sanitize_ld_scores_removes_nonfinite_values():
    ld1, ld2, ldx = simulation_utils.sanitize_ld_scores(
        np.array([np.nan, np.inf, -np.inf, 1.0]),
        np.array([2.0, np.nan]),
        np.array([np.inf, 3.0]),
    )

    for arr in (ld1, ld2, ldx):
        assert np.all(np.isfinite(arr))


def test_masked_methods_fallback_to_their_source_sumstats_when_not_run():
    b1_hat = np.array([0.1, -0.2, 0.05])
    se1_hat = np.array([0.2, 0.4, 0.1])
    b2_hat = np.array([1.0, -1.5, 0.8])
    se2_hat = np.array([0.5, 0.5, 0.4])
    bt_hat = np.array([0.7, -0.9, 0.3])
    se_t_hat = np.array([0.3, 0.3, 0.2])
    ld1 = np.ones(3)
    ld2 = np.ones(3)
    ldx = np.ones(3)
    omega = np.array([[0.1, 0.02], [0.02, 0.2]])

    z_scores = simulation_masked_omega.run_pop1_target_methods(
        b1_hat,
        se1_hat,
        b2_hat,
        se2_hat,
        bt_hat,
        se_t_hat,
        ld1,
        ld2,
        ldx,
        omega,
        pi2_omega_sum=0.05,
        propt=0.2,
        run_tracec=False,
        run_tracecb=False,
        run_pop1sc_bulk=False,
        run_pop2sc_bulk=False,
        missing_se=1e6,
    )

    expected_pop1 = b1_hat / se1_hat
    assert np.allclose(z_scores[0], expected_pop1)
    assert np.allclose(z_scores[2], expected_pop1)
    assert np.allclose(z_scores[3], expected_pop1)
    assert np.allclose(z_scores[4], expected_pop1)
    assert np.allclose(z_scores[5], b2_hat / se2_hat)
