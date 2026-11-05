"""Tests for project-priority parameter sensitivity reporting."""

from __future__ import annotations

import numpy as np
import pytest

from peatguard.analysis.risk_sensitivity import audit_priority_sensitivity


def test_reference_scenario_is_exact_and_variants_are_reported() -> None:
    canal = np.zeros((30, 30), dtype=bool)
    canal[:, 2] = True
    subsidence = np.linspace(-5.0, -50.0, canal.size).reshape(canal.shape)

    report = audit_priority_sensitivity(
        canal,
        subsidence,
        pixel_size_m=(10.0, 10.0),
        influence_radii_m=(100.0, 200.0),
        proximity_weights=(0.25, 0.45),
        reference_radius_m=200.0,
        reference_proximity_weight=0.45,
    )

    assert report["scenario_count"] == 4
    reference = next(
        item
        for item in report["scenarios"]
        if item["influence_radius_m"] == 200.0
        and item["proximity_weight"] == 0.45
    )
    assert reference["spearman_vs_reference"] == pytest.approx(1.0)
    assert reference["mean_absolute_score_change"] == pytest.approx(0.0)
    assert reference["top_10_percent_jaccard"] == pytest.approx(1.0)


def test_rejects_insufficient_measured_subsidence() -> None:
    canal = np.zeros((10, 10), dtype=bool)
    subsidence = np.full((10, 10), np.nan)

    with pytest.raises(ValueError, match="finite subsidence"):
        audit_priority_sensitivity(
            canal,
            subsidence,
            pixel_size_m=(10.0, 10.0),
        )
