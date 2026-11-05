"""Behavioral tests for fixed-support, relative-only phase contrasts."""

import json

import numpy as np
import pytest

from peatguard.timeseries.local_contrast import extract_local_contrast, solve_contrast_network


DATES = ["20260622", "20260704", "20260716"]
PAIRS = [(DATES[0], DATES[1]), (DATES[1], DATES[2]), (DATES[0], DATES[2])]


def _inputs(pair_count=3):
    phase = np.zeros((pair_count, 4, 8))
    coherence = np.ones_like(phase)
    components = np.ones_like(phase)
    target = np.zeros(phase.shape[1:], dtype=bool)
    target[:, :3] = True
    background = np.zeros_like(target)
    background[:, 5:] = True
    return phase, coherence, components, target, background


def test_arbitrary_shared_pair_offsets_and_integer_cycles_cancel():
    phase, coh, cc, target, background = _inputs()
    expected = np.array([0.2, -0.5, -0.3])
    phase[:, target] = expected[:, None]
    before = extract_local_contrast(phase, coh, cc, target, background)
    phase += np.array([12 * np.pi, -101 * np.pi, 123.456])[:, None, None]
    after = extract_local_contrast(phase, coh, cc, target, background)
    assert before["available"] and after["available"]
    np.testing.assert_allclose(after["contrast_rad"], expected, atol=1e-12)
    assert after["target_count"] == after["background_count"] == 12
    assert after["component_labels"] == [1, 1, 1]
    assert after["support_shape"] == [4, 8]
    assert after["target_support_flat_indices"] == np.flatnonzero(target).tolist()
    assert after["background_support_flat_indices"] == np.flatnonzero(background).tolist()
    assert not after["release_ready"] and not after["annual_velocity_ready"]
    json.dumps(after, allow_nan=False)


def test_background_motion_is_relative_signal_not_declared_target_motion():
    phase, coh, cc, target, background = _inputs()
    phase[:, background] = np.array([1, 2, 3])[:, None]
    result = extract_local_contrast(phase, coh, cc, target, background)
    assert result["contrast_rad"] == [-1, -2, -3]
    assert "may move" in result["limitations"][0]


def test_finite_zero_phase_is_valid_and_scatter_is_not_uncertainty():
    result = extract_local_contrast(*_inputs())
    assert result["available"]
    assert result["contrast_rad"] == [0, 0, 0]
    assert result["per_pair"][0]["target_spatial_std_rad"] == 0
    assert not result["spatial_scatter_is_uncertainty"]
    assert "observation_covariance" not in result


def test_support_is_frozen_across_every_pair_and_fraction_uses_original_mask():
    phase, coh, cc, target, background = _inputs()
    phase[:, 0, 0] = 100
    coh[1, 0, 0] = 0.1
    phase[:, 1, 0] = 50
    phase[2, 1, 0] = np.nan
    result = extract_local_contrast(phase, coh, cc, target, background)
    assert result["available"]
    assert result["contrast_rad"] == [0, 0, 0]
    assert result["target_count"] == 10
    assert result["target_fraction"] == pytest.approx(10 / 12)
    assert all(pair["target_count"] == 10 for pair in result["per_pair"])
    strict = extract_local_contrast(
        phase, coh, cc, target, background, min_target_fraction=0.9
    )
    assert not strict["available"]
    assert strict["reasons"] == ["insufficient_target_fraction"]
    assert strict["contrast_rad"] is None
    assert strict["target_support_flat_indices"] == result["target_support_flat_indices"]
    assert strict["background_support_flat_indices"] == result["background_support_flat_indices"]


def test_component_labels_must_match_in_every_pair():
    phase, coh, cc, target, background = _inputs()
    cc[2, background] = 2
    result = extract_local_contrast(phase, coh, cc, target, background)
    assert not result["available"]
    assert "no_matching_background_support" in result["reasons"]
    assert result["component_labels"] == [1, 1, 1]


def test_same_component_label_across_invalid_spatial_gap_is_not_a_shared_island():
    phase, coh, cc, target, background = _inputs()
    phase[1, :, 4] = np.nan
    result = extract_local_contrast(phase, coh, cc, target, background)
    assert not result["available"]
    assert "no_matching_background_support" in result["reasons"]
    assert result["complete_background_candidate_count"] == 12


def test_diagonal_only_contact_does_not_connect_spatial_islands():
    phase = np.array([[[0, np.nan], [np.nan, 0]]])
    target = np.array([[1, 0], [0, 0]], dtype=bool)
    background = np.array([[0, 0], [0, 1]], dtype=bool)
    result = extract_local_contrast(
        phase, np.ones_like(phase), np.ones_like(phase), target, background, min_pixels=1
    )
    assert not result["available"]
    assert "no_matching_background_support" in result["reasons"]


def test_labels_can_change_between_pairs_while_original_ids_are_preserved():
    phase, coh, cc, target, background = _inputs()
    cc[1] = 17
    cc[2] = 2
    result = extract_local_contrast(phase, coh, cc, target, background)
    assert result["available"]
    assert result["component_labels"] == [1, 17, 2]
    assert result["spatial_island_ids"] == [1, 1, 1]
    assert result["per_pair"][1]["component_label"] == 17


def test_signature_selection_is_support_based_and_does_not_switch_for_background():
    phase, coh, cc, target, background = _inputs()
    cc[:, 3, :3] = 2  # Minority target signature has all available background.
    cc[:, background] = 2
    phase[:, 3, :3] = 1234
    result = extract_local_contrast(
        phase, coh, cc, target, background, min_pixels=1, min_target_fraction=0.1
    )
    assert result["component_labels"] == [1, 1, 1]
    assert result["target_count"] == 9
    assert not result["available"]
    assert "no_matching_background_support" in result["reasons"]


def test_signature_tie_uses_geometry_not_phase_values_or_arbitrary_component_numbers():
    phase, coh, cc, target, background = _inputs()
    cc[:, 2:, :] = 2
    phase[:, :2] = 1000
    result = extract_local_contrast(
        phase, coh, cc, target, background, min_pixels=1
    )
    assert result["available"]
    assert result["component_labels"] == [1, 1, 1]
    assert result["target_count"] == result["background_count"] == 6
    assert result["contrast_rad"] == [0, 0, 0]
    cc[:, :2] = 100
    cc[:, 2:] = 7
    renamed = extract_local_contrast(phase, coh, cc, target, background, min_pixels=1)
    assert renamed["component_labels"] == [100, 100, 100]
    assert renamed["target_support_flat_indices"] == result["target_support_flat_indices"]
    assert renamed["background_support_flat_indices"] == result["background_support_flat_indices"]
    assert renamed["contrast_rad"] == result["contrast_rad"]


@pytest.mark.parametrize("label_value", [0, -1, 1.5, np.nan, np.inf])
def test_invalid_component_values_cannot_supply_support(label_value):
    phase, coh, cc, target, background = _inputs()
    cc[0, target] = label_value
    result = extract_local_contrast(phase, coh, cc, target, background)
    assert result["reasons"] == ["no_complete_target_support"]
    assert result["target_support_flat_indices"] == []
    assert result["background_support_flat_indices"] == []
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("coherence", [-0.1, 0.2, 1.01, np.nan, np.inf])
def test_bad_coherence_cannot_supply_support(coherence):
    phase, coh, cc, target, background = _inputs()
    coh[0, target] = coherence
    result = extract_local_contrast(phase, coh, cc, target, background)
    assert not result["available"]
    assert result["reasons"] == ["no_complete_target_support"]


def test_coherence_at_threshold_and_one_are_valid():
    phase, coh, cc, target, background = _inputs()
    coh[:, target] = 0.3
    assert extract_local_contrast(phase, coh, cc, target, background)["available"]


@pytest.mark.parametrize(
    "kwargs", [
        {"min_pixels": 0}, {"min_pixels": True}, {"min_pixels": 1.5},
        {"min_coherence": np.nan}, {"min_coherence": -0.1},
        {"min_target_fraction": 1.1}, {"min_background_fraction": -0.1},
    ]
)
def test_malformed_thresholds_are_rejected(kwargs):
    with pytest.raises(ValueError):
        extract_local_contrast(*_inputs(), **kwargs)


def test_mismatching_shapes_overlap_empty_and_nonbinary_masks_are_rejected():
    phase, coh, cc, target, background = _inputs()
    with pytest.raises(ValueError, match="identical shapes"):
        extract_local_contrast(phase, coh[:1], cc, target, background)
    with pytest.raises(ValueError, match="disjoint"):
        extract_local_contrast(phase, coh, cc, target, target)
    with pytest.raises(ValueError, match="nonempty"):
        extract_local_contrast(phase, coh, cc, target & False, background)
    with pytest.raises(ValueError, match="binary"):
        extract_local_contrast(phase, coh, cc, target * np.nan, background)


def test_known_redundant_network_and_cycle_residuals():
    result = solve_contrast_network(PAIRS, [2, 3, 5])
    np.testing.assert_allclose(result["epoch_contrast_rad"], [0, 2, 5], atol=1e-12)
    np.testing.assert_allclose(result["residual_rad"], 0, atol=1e-12)
    assert result["rank"] == 2
    assert result["cycle_count"] == result["redundancy"] == 1
    assert result["cycle_matrix"] == [[-1, -1, 1]]
    assert result["cycle_closure_rad"] == [0]
    assert result["epoch_covariance_rad2"] is None
    assert not result["release_ready"] and not result["annual_velocity_ready"]
    noisy = solve_contrast_network(PAIRS, [2, 3, 6])
    assert noisy["cycle_closure_rad"] == [1]
    np.testing.assert_allclose(noisy["residual_rad"], [-1 / 3, -1 / 3, 1 / 3])
    json.dumps(noisy, allow_nan=False)


def test_gls_uses_full_correlated_covariance_and_actual_linear_operator():
    covariance = np.array([[4, 1, 0.5], [1, 9, 2], [0.5, 2, 16]])
    observations = np.array([2, 3, 6])
    result = solve_contrast_network(PAIRS, observations, observation_covariance=covariance)
    design = np.array([[1, 0], [-1, 1], [0, 1]])
    precision = np.linalg.inv(covariance)
    analytic_covariance = np.linalg.inv(design.T @ precision @ design)
    expected_operator = analytic_covariance @ design.T @ precision
    operator = np.asarray(result["linear_operator"])
    np.testing.assert_allclose(operator[1:], expected_operator, atol=1e-12)
    np.testing.assert_allclose(result["epoch_contrast_rad"], operator @ observations)
    propagated = np.asarray(result["epoch_covariance_rad2"])
    np.testing.assert_allclose(propagated, operator @ covariance @ operator.T, atol=1e-12)
    np.testing.assert_allclose(propagated[1:, 1:], analytic_covariance, atol=1e-12)
    np.testing.assert_array_equal(propagated[0], 0)
    assert result["solver"] == "generalized_least_squares"
    json.dumps(result, allow_nan=False)


def test_tree_has_no_redundancy_or_inferred_precision_despite_exact_fit():
    result = solve_contrast_network(PAIRS[:2], [4, -3])
    np.testing.assert_allclose(result["epoch_contrast_rad"], [0, 4, 1], atol=1e-12)
    np.testing.assert_allclose(result["residual_rad"], 0, atol=1e-12)
    assert result["cycle_count"] == result["redundancy"] == 0
    assert result["cycle_closure_rad"] == []
    assert not result["has_redundant_closure_test"]
    assert result["epoch_covariance_rad2"] is None


def test_tree_covariance_is_cumulative_not_mean_interferogram_sigma():
    result = solve_contrast_network(PAIRS[:2], [4, -3], observation_covariance=np.diag([4, 9]))
    np.testing.assert_allclose(result["epoch_covariance_rad2"], [[0, 0, 0], [0, 4, 4], [0, 4, 13]])


def test_gls_pair_order_does_not_change_epochs_or_propagated_covariance():
    covariance = np.array([[4, 1, 0.5], [1, 9, 2], [0.5, 2, 16]])
    observations = np.array([2, 3, 6])
    order = np.array([2, 0, 1])
    original = solve_contrast_network(PAIRS, observations, observation_covariance=covariance)
    permuted = solve_contrast_network(
        [PAIRS[index] for index in order],
        observations[order],
        observation_covariance=covariance[np.ix_(order, order)],
    )
    assert permuted["dates"] == original["dates"]
    np.testing.assert_allclose(permuted["epoch_contrast_rad"], original["epoch_contrast_rad"])
    np.testing.assert_allclose(
        permuted["epoch_covariance_rad2"], original["epoch_covariance_rad2"]
    )


def test_every_fundamental_cycle_annihilates_the_date_incidence_matrix():
    dates = DATES + ["20260728", "20260809"]
    pairs = [
        (dates[1], dates[3]), (dates[0], dates[4]), (dates[2], dates[3]),
        (dates[0], dates[1]), (dates[0], dates[2]), (dates[3], dates[4]),
        (dates[1], dates[2]),
    ]
    incidence = np.zeros((len(pairs), len(dates)))
    for row, (first, second) in enumerate(pairs):
        incidence[row, dates.index(first)] = -1
        incidence[row, dates.index(second)] = 1
    observations = incidence @ np.array([0, 5, 1, -2, 8])
    result = solve_contrast_network(pairs, observations)
    assert result["cycle_count"] == 3
    cycles = np.asarray(result["cycle_matrix"])
    assert np.linalg.matrix_rank(cycles) == 3
    np.testing.assert_array_equal(cycles @ incidence, 0)
    np.testing.assert_array_equal(result["cycle_closure_rad"], 0)


def test_one_edge_network_is_allowed_but_remains_relative_only():
    result = solve_contrast_network([PAIRS[0]], [0])
    assert result["epoch_contrast_rad"] == [0, 0]
    assert result["rank"] == 1 and result["redundancy"] == 0
    assert not result["annual_velocity_ready"]


@pytest.mark.parametrize(
    "pairs, observations, message", [
        (None, [], "sequence"),
        ([], [], "at least one"),
        ([PAIRS[0], PAIRS[0]], [1, 1], "duplicate"),
        ([PAIRS[0], PAIRS[0][::-1]], [1, -1], "earlier-to-later"),
        ([(DATES[0], DATES[0])], [0], "earlier-to-later"),
        ([("20260230", DATES[1])], [0], "invalid calendar"),
        ([("2026062", DATES[1])], [0], "YYYYMMDD"),
        ([(b"20260622", DATES[1])], [0], "YYYYMMDD"),
        ([PAIRS[0], ("20260801", "20260813")], [1, 2], "disconnected"),
        ([PAIRS[0]], [np.nan], "finite vector"),
        ([PAIRS[0]], [1, 2], "finite vector"),
        (["2026062220260704"], [1], "exactly two"),
    ]
)
def test_invalid_networks_are_rejected(pairs, observations, message):
    with pytest.raises(ValueError, match=message):
        solve_contrast_network(pairs, observations)


@pytest.mark.parametrize(
    "covariance", [
        np.eye(2), np.full((3, 3), np.nan), np.zeros((3, 3)),
        np.diag([1, -1, 1]), np.array([[1, 0.1, 0], [0, 1, 0], [0, 0, 1]]),
    ]
)
def test_invalid_covariances_are_rejected(covariance):
    with pytest.raises(ValueError, match="observation_covariance"):
        solve_contrast_network(PAIRS, [2, 3, 5], observation_covariance=covariance)
