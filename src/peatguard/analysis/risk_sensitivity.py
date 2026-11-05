"""Sensitivity analysis for PeatGuard's project-priority risk heuristic.

The canal influence radius and proximity/subsidence weights are transparent
project assumptions, not calibrated physical parameters. This module measures
how strongly those choices affect the resulting priority surface.
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
from scipy.ndimage import distance_transform_edt
from scipy.stats import spearmanr


def _validate_inputs(
    canal_mask: np.ndarray,
    subsidence_mm_yr: np.ndarray,
    pixel_size_m: tuple[float, float],
) -> None:
    if canal_mask.ndim != 2 or subsidence_mm_yr.ndim != 2:
        raise ValueError("canal_mask and subsidence_mm_yr must be 2D")
    if canal_mask.shape != subsidence_mm_yr.shape:
        raise ValueError("canal_mask and subsidence_mm_yr must share a grid")
    if len(pixel_size_m) != 2 or any(value <= 0 for value in pixel_size_m):
        raise ValueError("pixel_size_m must contain positive (y, x) spacing")


def _top_fraction_mask(values: np.ndarray, fraction: float) -> np.ndarray:
    if not 0 < fraction < 1:
        raise ValueError("top fraction must be between 0 and 1")
    threshold = float(np.quantile(values, 1.0 - fraction))
    return values >= threshold


def _jaccard(first: np.ndarray, second: np.ndarray) -> float:
    union = int(np.count_nonzero(first | second))
    return (
        float(np.count_nonzero(first & second) / union)
        if union
        else float("nan")
    )


def audit_priority_sensitivity(
    canal_mask: np.ndarray,
    subsidence_mm_yr: np.ndarray,
    *,
    pixel_size_m: tuple[float, float],
    influence_radii_m: Iterable[float] = (300.0, 600.0, 900.0, 1200.0, 1500.0),
    proximity_weights: Iterable[float] = (0.30, 0.45, 0.60),
    reference_radius_m: float = 1200.0,
    reference_proximity_weight: float = 0.45,
    minimum_severe_scale_mm_yr: float = 5.0,
) -> dict:
    """Quantify ranking sensitivity to the project risk parameters.

    The production Tay builder normalizes negative rates by the magnitude of
    their 5th percentile, floored at ``minimum_severe_scale_mm_yr``. This audit
    mirrors that behavior and changes only the influence radius and component
    weights.
    """
    canal_mask = np.asarray(canal_mask, dtype=bool)
    subsidence_mm_yr = np.asarray(subsidence_mm_yr, dtype=np.float64)
    _validate_inputs(canal_mask, subsidence_mm_yr, pixel_size_m)

    radii = sorted({float(value) for value in influence_radii_m})
    weights = sorted({float(value) for value in proximity_weights})
    if any(value <= 0 for value in radii):
        raise ValueError("influence radii must be positive")
    if any(value < 0 or value > 1 for value in weights):
        raise ValueError("proximity weights must be between 0 and 1")
    if reference_radius_m not in radii:
        radii.append(float(reference_radius_m))
        radii.sort()
    if reference_proximity_weight not in weights:
        weights.append(float(reference_proximity_weight))
        weights.sort()

    valid = np.isfinite(subsidence_mm_yr)
    subsiding = valid & (subsidence_mm_yr < 0)
    if np.count_nonzero(valid) < 100:
        raise ValueError("at least 100 finite subsidence pixels are required")
    if np.count_nonzero(subsiding) < 100:
        raise ValueError("at least 100 subsiding pixels are required")

    severe_scale = max(
        abs(float(np.percentile(subsidence_mm_yr[subsiding], 5))),
        float(minimum_severe_scale_mm_yr),
    )
    subsidence_risk = np.zeros(subsidence_mm_yr.shape, dtype=np.float32)
    subsidence_risk[subsiding] = np.clip(
        -subsidence_mm_yr[subsiding] / severe_scale,
        0.0,
        1.0,
    )
    distance_m = distance_transform_edt(
        ~canal_mask,
        sampling=pixel_size_m,
    ).astype(np.float32)

    scenarios: dict[tuple[float, float], np.ndarray] = {}
    for radius in radii:
        proximity = np.clip(1.0 - distance_m / radius, 0.0, 1.0)
        for proximity_weight in weights:
            scenarios[(radius, proximity_weight)] = (
                proximity_weight * proximity[valid]
                + (1.0 - proximity_weight) * subsidence_risk[valid]
            ).astype(np.float32)

    reference_key = (
        float(reference_radius_m),
        float(reference_proximity_weight),
    )
    reference = scenarios[reference_key]
    reference_top_10 = _top_fraction_mask(reference, 0.10)
    reference_top_05 = _top_fraction_mask(reference, 0.05)

    scenario_reports = []
    for (radius, proximity_weight), values in scenarios.items():
        correlation = spearmanr(reference, values).statistic
        scenario_reports.append(
            {
                "influence_radius_m": radius,
                "proximity_weight": proximity_weight,
                "subsidence_weight": 1.0 - proximity_weight,
                "spearman_vs_reference": float(correlation),
                "mean_absolute_score_change": float(
                    np.mean(np.abs(values - reference))
                ),
                "p95_absolute_score_change": float(
                    np.percentile(np.abs(values - reference), 95)
                ),
                "top_10_percent_jaccard": _jaccard(
                    reference_top_10,
                    _top_fraction_mask(values, 0.10),
                ),
                "top_5_percent_jaccard": _jaccard(
                    reference_top_05,
                    _top_fraction_mask(values, 0.05),
                ),
            }
        )

    non_reference = [
        item
        for item in scenario_reports
        if not (
            item["influence_radius_m"] == reference_key[0]
            and item["proximity_weight"] == reference_key[1]
        )
    ]
    return {
        "status": "project_parameter_sensitivity_not_field_calibration",
        "interpretation": (
            "Ranking robustness across declared project priors only. This does "
            "not validate the hydrological response shape or radius."
        ),
        "coverage": {
            "grid_pixels": int(canal_mask.size),
            "canal_pixels": int(np.count_nonzero(canal_mask)),
            "finite_subsidence_pixels": int(np.count_nonzero(valid)),
            "subsiding_pixels": int(np.count_nonzero(subsiding)),
        },
        "normalization": {
            "method": "magnitude of 5th percentile among negative rates",
            "severe_scale_mm_yr": severe_scale,
            "minimum_scale_mm_yr": float(minimum_severe_scale_mm_yr),
        },
        "reference_scenario": {
            "influence_radius_m": reference_key[0],
            "proximity_weight": reference_key[1],
            "subsidence_weight": 1.0 - reference_key[1],
        },
        "scenario_count": len(scenario_reports),
        "scenarios": scenario_reports,
        "robustness_envelope_excluding_reference": {
            "minimum_spearman": float(
                min(item["spearman_vs_reference"] for item in non_reference)
            ),
            "minimum_top_10_percent_jaccard": float(
                min(item["top_10_percent_jaccard"] for item in non_reference)
            ),
            "minimum_top_5_percent_jaccard": float(
                min(item["top_5_percent_jaccard"] for item in non_reference)
            ),
            "maximum_mean_absolute_score_change": float(
                max(item["mean_absolute_score_change"] for item in non_reference)
            ),
            "maximum_p95_absolute_score_change": float(
                max(item["p95_absolute_score_change"] for item in non_reference)
            ),
        },
    }
