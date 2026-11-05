"""Tests for historical and short-span NISAR evidence fusion."""

from __future__ import annotations

import numpy as np
import pytest

from peatguard.analysis.historical_nisar import (
    DAYS_PER_YEAR,
    EVIDENCE_CONCORDANT_DOWNWARD,
    EVIDENCE_HISTORICAL_OTHER_NISAR_INCONCLUSIVE,
    EVIDENCE_HISTORICAL_SUBSIDENCE_ONLY,
    EVIDENCE_NISAR_ONLY_DOWNWARD,
    EVIDENCE_NISAR_ONLY_UPWARD,
    EVIDENCE_POSSIBLE_CHANGED_DOWNWARD,
    annualize_short_span_displacement,
    bayesian_rate_sensitivity,
    classify_historical_nisar_evidence,
)


def test_annualize_short_span_displacement_is_unit_conversion() -> None:
    displacement = np.array([-10.0, 5.0])
    uncertainty = np.array([20.0, 10.0])
    rate, sigma = annualize_short_span_displacement(
        displacement,
        uncertainty,
        DAYS_PER_YEAR,
    )
    assert np.array_equal(rate, displacement)
    assert np.array_equal(sigma, uncertainty)


def test_bayesian_fusion_exposes_small_nisar_weight() -> None:
    result = bayesian_rate_sensitivity(
        np.array([-20.0]),
        np.array([-200.0]),
        np.array([800.0]),
        historical_measurement_sigma_mm_per_year=3.6,
        historical_change_sigma_mm_per_year=50.0,
    )
    expected_prior_variance = 3.6**2 + 50.0**2
    expected_weight = expected_prior_variance / (
        expected_prior_variance + 800.0**2
    )
    assert result.nisar_weight_upper_bound[0] == pytest.approx(expected_weight)
    assert result.posterior_rate_mm_per_year[0] == pytest.approx(
        -20.0 * (1.0 - expected_weight) - 200.0 * expected_weight
    )
    assert result.nisar_weight_upper_bound[0] < 0.01


def test_evidence_classification_keeps_epochs_explicit() -> None:
    historical = np.array(
        [
            [-20.0, -2.0, -20.0],
            [np.nan, np.nan, np.nan],
        ]
    )
    displacement = np.array(
        [
            [-30.0, -30.0, np.nan],
            [-30.0, 30.0, np.nan],
        ]
    )
    uncertainty = np.array(
        [
            [10.0, 10.0, np.nan],
            [10.0, 10.0, np.nan],
        ]
    )
    evidence, z_score = classify_historical_nisar_evidence(
        historical,
        displacement,
        uncertainty,
    )
    assert evidence[0, 0] == EVIDENCE_CONCORDANT_DOWNWARD
    assert evidence[0, 1] == EVIDENCE_POSSIBLE_CHANGED_DOWNWARD
    assert evidence[0, 2] == EVIDENCE_HISTORICAL_SUBSIDENCE_ONLY
    assert evidence[1, 0] == EVIDENCE_NISAR_ONLY_DOWNWARD
    assert evidence[1, 1] == EVIDENCE_NISAR_ONLY_UPWARD
    assert evidence[1, 2] == 0
    assert z_score[0, 0] == pytest.approx(-3.0)


def test_inconclusive_nisar_does_not_become_current_subsidence() -> None:
    evidence, _ = classify_historical_nisar_evidence(
        np.array([[-2.0]]),
        np.array([[-5.0]]),
        np.array([[50.0]]),
    )
    assert evidence[0, 0] == EVIDENCE_HISTORICAL_OTHER_NISAR_INCONCLUSIVE
