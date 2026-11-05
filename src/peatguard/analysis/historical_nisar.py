"""Uncertainty-aware fusion of historical peat motion and NISAR evidence.

This module deliberately separates a historical long-term rate from a short
NISAR commissioning displacement. The quantitative Bayesian result is a
sensitivity experiment, not a current velocity product, because annualizing a
24-day displacement assumes linear motion and ignores seasonal peat motion.

The historical source used by PeatGuard is the ALOS-2 long-term-rate product
from Tay et al. (2025), doi:10.1038/s43247-025-02926-6. Tay remains an
independent benchmark. No bias fit, scale fit, or datum fit to Tay is allowed.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

DAYS_PER_YEAR = 365.2425

# Historical and current-evidence classification codes.
EVIDENCE_NONE = 0
EVIDENCE_HISTORICAL_SUBSIDENCE_ONLY = 1
EVIDENCE_HISTORICAL_OTHER_ONLY = 2
EVIDENCE_HISTORICAL_SUBSIDENCE_NISAR_INCONCLUSIVE = 3
EVIDENCE_HISTORICAL_OTHER_NISAR_INCONCLUSIVE = 4
EVIDENCE_CONCORDANT_DOWNWARD = 5
EVIDENCE_POSSIBLE_CHANGED_DOWNWARD = 6
EVIDENCE_SIGNIFICANT_UPWARD_WITH_HISTORY = 7
EVIDENCE_NISAR_ONLY_INCONCLUSIVE = 8
EVIDENCE_NISAR_ONLY_DOWNWARD = 9
EVIDENCE_NISAR_ONLY_UPWARD = 10

EVIDENCE_LABELS = {
    EVIDENCE_NONE: "No historical or current evidence",
    EVIDENCE_HISTORICAL_SUBSIDENCE_ONLY: (
        "Historical subsidence; no valid NISAR observation"
    ),
    EVIDENCE_HISTORICAL_OTHER_ONLY: (
        "Historical non-subsiding or rebound; no valid NISAR observation"
    ),
    EVIDENCE_HISTORICAL_SUBSIDENCE_NISAR_INCONCLUSIVE: (
        "Historical subsidence; current NISAR inconclusive"
    ),
    EVIDENCE_HISTORICAL_OTHER_NISAR_INCONCLUSIVE: (
        "Historical non-subsiding or rebound; current NISAR inconclusive"
    ),
    EVIDENCE_CONCORDANT_DOWNWARD: (
        "Historical subsidence and significant current downward motion"
    ),
    EVIDENCE_POSSIBLE_CHANGED_DOWNWARD: (
        "No historical subsidence; significant current downward motion"
    ),
    EVIDENCE_SIGNIFICANT_UPWARD_WITH_HISTORY: (
        "Significant current upward motion or rebound with historical coverage"
    ),
    EVIDENCE_NISAR_ONLY_INCONCLUSIVE: (
        "Current NISAR inconclusive; no historical coverage"
    ),
    EVIDENCE_NISAR_ONLY_DOWNWARD: (
        "Significant current downward motion; no historical coverage"
    ),
    EVIDENCE_NISAR_ONLY_UPWARD: (
        "Significant current upward motion; no historical coverage"
    ),
}


@dataclass(frozen=True)
class BayesianRateFusion:
    """Arrays and assumptions from a historical-prior sensitivity test."""

    posterior_rate_mm_per_year: np.ndarray
    posterior_sigma_mm_per_year: np.ndarray
    nisar_weight_upper_bound: np.ndarray
    historical_total_sigma_mm_per_year: float
    historical_measurement_sigma_mm_per_year: float
    historical_change_sigma_mm_per_year: float

    def assumptions(self) -> dict[str, float]:
        """Return serializable scalar assumptions."""
        return {
            "historical_total_sigma_mm_per_year": (
                self.historical_total_sigma_mm_per_year
            ),
            "historical_measurement_sigma_mm_per_year": (
                self.historical_measurement_sigma_mm_per_year
            ),
            "historical_change_sigma_mm_per_year": (
                self.historical_change_sigma_mm_per_year
            ),
        }


def annualize_short_span_displacement(
    displacement_mm: np.ndarray,
    uncertainty_mm: np.ndarray,
    span_days: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert a short-span displacement and uncertainty to equivalent rates.

    This is a diagnostic unit conversion. It does not make the short
    displacement an annual velocity estimate.
    """
    if not np.isfinite(span_days) or span_days <= 0:
        raise ValueError("span_days must be positive")
    displacement = np.asarray(displacement_mm, dtype=np.float64)
    uncertainty = np.asarray(uncertainty_mm, dtype=np.float64)
    if displacement.shape != uncertainty.shape:
        raise ValueError("displacement and uncertainty shapes differ")
    if np.any(uncertainty[np.isfinite(uncertainty)] < 0):
        raise ValueError("uncertainty cannot be negative")
    scale = DAYS_PER_YEAR / span_days
    return displacement * scale, uncertainty * scale


def bayesian_rate_sensitivity(
    historical_rate_mm_per_year: np.ndarray,
    nisar_equivalent_rate_mm_per_year: np.ndarray,
    nisar_sigma_mm_per_year: np.ndarray,
    *,
    historical_measurement_sigma_mm_per_year: float = 3.6,
    historical_change_sigma_mm_per_year: float = 50.0,
) -> BayesianRateFusion:
    """Run an inverse-variance historical-prior sensitivity experiment.

    The process-change term acknowledges that a 2018-2024 rate need not equal
    a 2026 rate. NISAR uncertainty is currently an ionosphere-plus-reference
    lower bound, so the returned NISAR weight is an upper bound.
    """
    historical = np.asarray(historical_rate_mm_per_year, dtype=np.float64)
    nisar = np.asarray(nisar_equivalent_rate_mm_per_year, dtype=np.float64)
    nisar_sigma = np.asarray(nisar_sigma_mm_per_year, dtype=np.float64)
    if historical.shape != nisar.shape or historical.shape != nisar_sigma.shape:
        raise ValueError("historical, NISAR rate, and NISAR sigma shapes differ")
    if (
        not np.isfinite(historical_measurement_sigma_mm_per_year)
        or historical_measurement_sigma_mm_per_year <= 0
    ):
        raise ValueError("historical measurement sigma must be positive")
    if (
        not np.isfinite(historical_change_sigma_mm_per_year)
        or historical_change_sigma_mm_per_year < 0
    ):
        raise ValueError("historical change sigma cannot be negative")

    prior_variance = (
        historical_measurement_sigma_mm_per_year**2
        + historical_change_sigma_mm_per_year**2
    )
    historical_total_sigma = float(np.sqrt(prior_variance))
    posterior = np.full(historical.shape, np.nan, dtype=np.float64)
    posterior_sigma = np.full(historical.shape, np.nan, dtype=np.float64)
    nisar_weight = np.full(historical.shape, np.nan, dtype=np.float64)
    valid = (
        np.isfinite(historical)
        & np.isfinite(nisar)
        & np.isfinite(nisar_sigma)
        & (nisar_sigma > 0)
    )
    observation_variance = nisar_sigma[valid] ** 2
    weight = prior_variance / (prior_variance + observation_variance)
    posterior[valid] = historical[valid] * (1.0 - weight) + nisar[valid] * weight
    posterior_sigma[valid] = np.sqrt(
        1.0 / (1.0 / prior_variance + 1.0 / observation_variance)
    )
    nisar_weight[valid] = weight
    return BayesianRateFusion(
        posterior_rate_mm_per_year=posterior,
        posterior_sigma_mm_per_year=posterior_sigma,
        nisar_weight_upper_bound=nisar_weight,
        historical_total_sigma_mm_per_year=historical_total_sigma,
        historical_measurement_sigma_mm_per_year=(
            historical_measurement_sigma_mm_per_year
        ),
        historical_change_sigma_mm_per_year=(
            historical_change_sigma_mm_per_year
        ),
    )


def classify_historical_nisar_evidence(
    historical_rate_mm_per_year: np.ndarray,
    nisar_displacement_mm: np.ndarray,
    nisar_uncertainty_mm: np.ndarray,
    *,
    historical_subsidence_threshold_mm_per_year: float = -5.0,
    significance_sigma: float = 2.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Classify historical state and short-span NISAR directional evidence.

    A significant NISAR displacement is still described as downward or upward
    motion, not long-term subsidence or uplift.

    The z statistic divides by the supplied PARTIAL uncertainty (ionosphere
    only: troposphere, unwrapping, orbit, DEM, and datum error are excluded),
    so significance here is an upper bound on real significance and must be
    labeled as such wherever it is reported.
    """
    historical = np.asarray(historical_rate_mm_per_year, dtype=np.float64)
    displacement = np.asarray(nisar_displacement_mm, dtype=np.float64)
    uncertainty = np.asarray(nisar_uncertainty_mm, dtype=np.float64)
    if historical.shape != displacement.shape or historical.shape != uncertainty.shape:
        raise ValueError("historical, displacement, and uncertainty shapes differ")
    if not np.isfinite(significance_sigma) or significance_sigma <= 0:
        raise ValueError("significance_sigma must be positive")

    evidence = np.full(historical.shape, EVIDENCE_NONE, dtype=np.uint8)
    z_score = np.full(historical.shape, np.nan, dtype=np.float64)
    historical_valid = np.isfinite(historical)
    historical_subsidence = (
        historical_valid
        & (historical < historical_subsidence_threshold_mm_per_year)
    )
    historical_other = historical_valid & ~historical_subsidence
    nisar_valid = (
        np.isfinite(displacement)
        & np.isfinite(uncertainty)
        & (uncertainty > 0)
    )
    z_score[nisar_valid] = displacement[nisar_valid] / uncertainty[nisar_valid]
    downward = nisar_valid & (z_score <= -significance_sigma)
    upward = nisar_valid & (z_score >= significance_sigma)
    inconclusive = nisar_valid & ~downward & ~upward

    evidence[historical_subsidence & ~nisar_valid] = (
        EVIDENCE_HISTORICAL_SUBSIDENCE_ONLY
    )
    evidence[historical_other & ~nisar_valid] = EVIDENCE_HISTORICAL_OTHER_ONLY
    evidence[historical_subsidence & inconclusive] = (
        EVIDENCE_HISTORICAL_SUBSIDENCE_NISAR_INCONCLUSIVE
    )
    evidence[historical_other & inconclusive] = (
        EVIDENCE_HISTORICAL_OTHER_NISAR_INCONCLUSIVE
    )
    evidence[historical_subsidence & downward] = EVIDENCE_CONCORDANT_DOWNWARD
    evidence[historical_other & downward] = EVIDENCE_POSSIBLE_CHANGED_DOWNWARD
    evidence[historical_valid & upward] = EVIDENCE_SIGNIFICANT_UPWARD_WITH_HISTORY
    evidence[~historical_valid & inconclusive] = EVIDENCE_NISAR_ONLY_INCONCLUSIVE
    evidence[~historical_valid & downward] = EVIDENCE_NISAR_ONLY_DOWNWARD
    evidence[~historical_valid & upward] = EVIDENCE_NISAR_ONLY_UPWARD
    return evidence, z_score
