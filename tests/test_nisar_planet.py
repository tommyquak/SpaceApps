"""Tests for NISAR and PlanetScope complementarity diagnostics."""

from __future__ import annotations

import numpy as np
import pytest

from peatguard.analysis.nisar_planet import (
    block_bootstrap_near_far_difference,
    canal_distance_m,
    distance_band_statistics,
)


def test_distance_bands_report_coherence_threshold_coverage() -> None:
    canals = np.zeros((100, 100), dtype=bool)
    canals[:, 0] = True
    distance = canal_distance_m(canals, (20.0, 20.0))
    coherence = np.where(distance < 100.0, 0.3, 0.8)
    statistics = distance_band_statistics(coherence, distance)
    near = statistics[0]
    far = next(
        item for item in statistics if item.lower_distance_m == 300.0
    )
    assert near.median_coherence == pytest.approx(0.3)
    assert near.fraction_ge_0_50 == 0.0
    assert far.median_coherence == pytest.approx(0.8)
    assert far.fraction_ge_0_70 == 1.0


def test_block_bootstrap_is_deterministic() -> None:
    canals = np.zeros((400, 400), dtype=bool)
    canals[:, ::60] = True
    distance = canal_distance_m(canals, (20.0, 20.0))
    coherence = np.where(distance < 100.0, 0.4, 0.6)
    first = block_bootstrap_near_far_difference(
        coherence,
        distance,
        block_size_pixels=100,
        minimum_pixels_per_group=10,
        bootstrap_iterations=100,
    )
    second = block_bootstrap_near_far_difference(
        coherence,
        distance,
        block_size_pixels=100,
        minimum_pixels_per_group=10,
        bootstrap_iterations=100,
    )
    assert first == second
    assert first["block_difference_mean"] < 0
