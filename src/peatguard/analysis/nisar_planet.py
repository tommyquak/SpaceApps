"""Quantify how NISAR coherence complements PlanetScope canal geometry."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Optional, Sequence

import numpy as np
from scipy.ndimage import distance_transform_edt


@dataclass(frozen=True)
class DistanceBandStatistic:
    """Coherence distribution within one distance band from mapped canals."""

    lower_distance_m: float
    upper_distance_m: Optional[float]
    pixel_count: int
    mean_coherence: float
    median_coherence: float
    percentile_25: float
    percentile_75: float
    fraction_ge_0_35: float
    fraction_ge_0_50: float
    fraction_ge_0_70: float

    def to_dict(self) -> dict:
        return asdict(self)


def canal_distance_m(
    canal_mask: np.ndarray,
    pixel_spacing_m: tuple[float, float],
) -> np.ndarray:
    """Return Euclidean distance to the nearest mapped canal pixel."""
    canals = np.asarray(canal_mask, dtype=bool)
    if canals.ndim != 2:
        raise ValueError("Canal mask must be two-dimensional")
    if not canals.any():
        raise ValueError("Canal mask contains no canal pixels")
    spacing_y, spacing_x = pixel_spacing_m
    if spacing_y <= 0 or spacing_x <= 0:
        raise ValueError("Pixel spacing must be positive")
    return distance_transform_edt(
        ~canals, sampling=(spacing_y, spacing_x)
    )


def distance_band_statistics(
    coherence: np.ndarray,
    distance_m: np.ndarray,
    *,
    bands_m: Sequence[tuple[float, Optional[float]]] = (
        (0.0, 25.0),
        (25.0, 100.0),
        (100.0, 300.0),
        (300.0, 1_000.0),
        (1_000.0, None),
    ),
    mask: Optional[np.ndarray] = None,
) -> tuple[DistanceBandStatistic, ...]:
    """Summarize coherence in predeclared canal-distance bands."""
    values = np.asarray(coherence, dtype=np.float64)
    distance = np.asarray(distance_m, dtype=np.float64)
    if values.shape != distance.shape:
        raise ValueError("Coherence and canal-distance grids differ")
    valid = np.isfinite(values) & np.isfinite(distance)
    if mask is not None:
        if np.asarray(mask).shape != values.shape:
            raise ValueError("Analysis mask and coherence grids differ")
        valid &= np.asarray(mask, dtype=bool)

    output: list[DistanceBandStatistic] = []
    prior_upper = 0.0
    for lower, upper in bands_m:
        if lower < prior_upper or lower < 0:
            raise ValueError("Distance bands must be ordered and non-negative")
        selected = valid & (distance >= lower)
        if upper is not None:
            if upper <= lower:
                raise ValueError("Distance-band upper bound must exceed lower")
            selected &= distance < upper
            prior_upper = upper
        else:
            prior_upper = lower
        data = values[selected]
        if data.size == 0:
            continue
        p25, median, p75 = np.percentile(data, [25.0, 50.0, 75.0])
        output.append(
            DistanceBandStatistic(
                lower_distance_m=float(lower),
                upper_distance_m=(
                    float(upper) if upper is not None else None
                ),
                pixel_count=int(data.size),
                mean_coherence=float(np.mean(data)),
                median_coherence=float(median),
                percentile_25=float(p25),
                percentile_75=float(p75),
                fraction_ge_0_35=float(np.mean(data >= 0.35)),
                fraction_ge_0_50=float(np.mean(data >= 0.50)),
                fraction_ge_0_70=float(np.mean(data >= 0.70)),
            )
        )
    return tuple(output)


def block_bootstrap_near_far_difference(
    coherence: np.ndarray,
    distance_m: np.ndarray,
    *,
    mask: Optional[np.ndarray] = None,
    near_range_m: tuple[float, float] = (0.0, 100.0),
    far_range_m: tuple[float, float] = (300.0, 1_000.0),
    block_size_pixels: int = 200,
    minimum_pixels_per_group: int = 25,
    bootstrap_iterations: int = 5_000,
    random_seed: int = 20260723,
) -> dict[str, object]:
    """Block-bootstrap the near-minus-far median coherence difference.

    This is a spatial sensitivity diagnostic, not a causal canal-impact test.
    Blocks reduce, but do not eliminate, spatial autocorrelation.
    """
    values = np.asarray(coherence, dtype=np.float64)
    distance = np.asarray(distance_m, dtype=np.float64)
    if values.shape != distance.shape:
        raise ValueError("Coherence and canal-distance grids differ")
    if block_size_pixels < 1 or bootstrap_iterations < 1:
        raise ValueError("Block size and bootstrap iterations must be positive")
    valid = np.isfinite(values) & np.isfinite(distance)
    if mask is not None:
        if np.asarray(mask).shape != values.shape:
            raise ValueError("Analysis mask and coherence grids differ")
        valid &= np.asarray(mask, dtype=bool)

    near = (
        valid
        & (distance >= near_range_m[0])
        & (distance < near_range_m[1])
    )
    far = (
        valid
        & (distance >= far_range_m[0])
        & (distance < far_range_m[1])
    )
    block_differences: list[float] = []
    height, width = values.shape
    for row in range(0, height, block_size_pixels):
        row_stop = min(height, row + block_size_pixels)
        for column in range(0, width, block_size_pixels):
            column_stop = min(width, column + block_size_pixels)
            window = np.s_[row:row_stop, column:column_stop]
            near_values = values[window][near[window]]
            far_values = values[window][far[window]]
            if (
                near_values.size < minimum_pixels_per_group
                or far_values.size < minimum_pixels_per_group
            ):
                continue
            block_differences.append(
                float(np.median(near_values) - np.median(far_values))
            )
    if len(block_differences) < 3:
        raise ValueError("Too few mixed near/far blocks for spatial bootstrap")

    differences = np.asarray(block_differences, dtype=np.float64)
    generator = np.random.default_rng(random_seed)
    indices = generator.integers(
        0,
        differences.size,
        size=(bootstrap_iterations, differences.size),
    )
    bootstrap_means = np.mean(differences[indices], axis=1)
    lower, upper = np.percentile(bootstrap_means, [2.5, 97.5])
    return {
        "status": "spatial_sensitivity_not_causal_inference",
        "near_range_m": list(near_range_m),
        "far_range_m": list(far_range_m),
        "block_size_pixels": block_size_pixels,
        "mixed_block_count": int(differences.size),
        "block_difference_mean": float(np.mean(differences)),
        "block_difference_median": float(np.median(differences)),
        "bootstrap_mean_95_percent_interval": [
            float(lower),
            float(upper),
        ],
        "bootstrap_iterations": bootstrap_iterations,
        "random_seed": random_seed,
    }
