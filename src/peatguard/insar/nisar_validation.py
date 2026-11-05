"""Validation for PeatGuard NISAR interferograms.

The fine GSLC workflow and mission GUNW workflow do not share a pixel grid.
Wrapped phase must not be averaged arithmetically as radians. The legacy
diagnostic averages unit phasors. An opt-in mode averages amplitude-preserving
complex cross-products instead. Both operate on the geocoded GSLC grid and are
not an exact reproduction of radar-coordinate multilooking then geocoding.
This module tests both phase signs after a single circular reference offset
and reports coherence agreement.

For a three-date GSLC network, this module also audits the native-grid wrapped
phase closure A-B + B-C - (A-C). The closure calculation is block-wise and
requires the acquisition-date tags to prove that the three supplied rasters
form the intended triangle.
"""

from __future__ import annotations

import json
import logging
import math
from contextlib import ExitStack
from datetime import date
from pathlib import Path
from typing import Optional

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject, transform_bounds
from rasterio.windows import Window, bounds as window_bounds, from_bounds
from scipy.stats import spearmanr

logger = logging.getLogger(__name__)

_WORK_NODATA = -9999.0


def wrap_phase(phase: np.ndarray) -> np.ndarray:
    """Wrap phase angles to the half-open interval [-pi, pi)."""
    return (phase + np.pi) % (2.0 * np.pi) - np.pi


def _same_grid(first: rasterio.DatasetReader, second: rasterio.DatasetReader) -> bool:
    return (
        first.crs == second.crs
        and first.width == second.width
        and first.height == second.height
        and first.transform.almost_equals(second.transform)
    )


def _valid_raster_values(
    data: np.ndarray,
    nodata: Optional[float],
) -> np.ndarray:
    valid = np.isfinite(data)
    if nodata is not None and np.isfinite(nodata):
        valid &= data != nodata
    return valid


def _source_window_for_destination(
    source: rasterio.DatasetReader,
    destination: rasterio.DatasetReader,
    destination_window: Window,
) -> Window:
    """Return a clipped source window covering one destination block."""
    left, bottom, right, top = window_bounds(
        destination_window, destination.transform
    )
    if source.crs != destination.crs:
        left, bottom, right, top = transform_bounds(
            destination.crs,
            source.crs,
            left,
            bottom,
            right,
            top,
            densify_pts=21,
        )
    raw = from_bounds(left, bottom, right, top, source.transform)
    # One source-pixel margin avoids edge loss from floating-point grid bounds.
    col_start = max(0, math.floor(raw.col_off) - 1)
    row_start = max(0, math.floor(raw.row_off) - 1)
    col_stop = min(source.width, math.ceil(raw.col_off + raw.width) + 1)
    row_stop = min(source.height, math.ceil(raw.row_off + raw.height) + 1)
    return Window(
        col_start,
        row_start,
        max(0, col_stop - col_start),
        max(0, row_stop - row_start),
    )


def _reproject_fine_block(
    fine_phase: rasterio.DatasetReader,
    fine_coherence: rasterio.DatasetReader,
    mission_phase: rasterio.DatasetReader,
    mission_window: Window,
    fine_cross_product: Optional[rasterio.DatasetReader] = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Average unit phasors or complex cross-products onto a mission block.

    Strength is abs(mean(z)) / mean(abs(z)); z is a unit phasor in legacy
    mode, or the original cross-product in complex-weighted mode. It is a
    normalized resultant, not interferometric coherence.
    """
    source_window = _source_window_for_destination(
        fine_phase, mission_phase, mission_window
    )
    output_shape = (int(mission_window.height), int(mission_window.width))
    if source_window.width <= 0 or source_window.height <= 0:
        empty = np.full(output_shape, np.nan, dtype=np.float32)
        return empty, empty.copy(), empty.copy()

    phase = fine_phase.read(1, window=source_window).astype(np.float32)
    coherence = fine_coherence.read(1, window=source_window).astype(np.float32)
    valid = _valid_raster_values(phase, fine_phase.nodata)
    valid &= _valid_raster_values(coherence, fine_coherence.nodata)
    valid &= (coherence >= 0) & (coherence <= 1)
    if fine_cross_product is None:
        real = np.cos(phase)
        imaginary = np.sin(phase)
        amplitude = np.ones(phase.shape, dtype=np.float32)
    else:
        real = fine_cross_product.read(1, window=source_window).astype(np.float32)
        imaginary = fine_cross_product.read(2, window=source_window).astype(np.float32)
        valid &= _valid_raster_values(real, fine_cross_product.nodata)
        valid &= _valid_raster_values(imaginary, fine_cross_product.nodata)
        amplitude = np.hypot(real, imaginary)
        valid &= np.isfinite(amplitude) & (amplitude > 0)
        if np.any(valid):
            disagreement = wrap_phase(np.arctan2(imaginary[valid], real[valid]) - phase[valid])
            if np.any(np.abs(disagreement) > 1e-4):
                raise ValueError("Fine cross-product phase does not match the phase raster")

    # NaN cannot collide with a legitimate signed cross-product component.
    cosine = np.full(phase.shape, np.nan, dtype=np.float32)
    sine = np.full(phase.shape, np.nan, dtype=np.float32)
    source_amplitude = np.full(phase.shape, np.nan, dtype=np.float32)
    source_coherence = np.full(phase.shape, np.nan, dtype=np.float32)
    cosine[valid] = real[valid]
    sine[valid] = imaginary[valid]
    source_amplitude[valid] = amplitude[valid]
    source_coherence[valid] = coherence[valid]

    destination_transform = rasterio.windows.transform(
        mission_window, mission_phase.transform
    )
    source_transform = rasterio.windows.transform(source_window, fine_phase.transform)
    averaged = []
    for source_data in (cosine, sine, source_coherence, source_amplitude):
        destination_data = np.full(output_shape, np.nan, dtype=np.float32)
        reproject(
            source=source_data,
            destination=destination_data,
            src_transform=source_transform,
            src_crs=fine_phase.crs,
            src_nodata=np.nan,
            dst_transform=destination_transform,
            dst_crs=mission_phase.crs,
            dst_nodata=np.nan,
            resampling=Resampling.average,
        )
        averaged.append(destination_data)

    average_cosine, average_sine, average_coherence, average_amplitude = averaged
    magnitude = np.hypot(average_cosine, average_sine)
    phasor_strength = np.full(output_shape, np.nan, dtype=np.float32)
    np.divide(magnitude, average_amplitude, out=phasor_strength, where=average_amplitude > 0)
    phasor_strength = np.clip(phasor_strength, 0.0, 1.0)
    aggregated_phase = np.arctan2(average_sine, average_cosine).astype(np.float32)
    # Reprojection may leave roundoff after exact cancellation. Such a
    # resultant has no numerically meaningful phase even with a zero gate.
    invalid = (
        ~np.isfinite(average_cosine)
        | ~np.isfinite(average_sine)
        | (phasor_strength <= np.finfo(np.float32).eps)
    )
    aggregated_phase[invalid] = np.nan
    phasor_strength[invalid] = np.nan
    return aggregated_phase, average_coherence, phasor_strength


def _validate_cross_product_provenance(
    cross: rasterio.DatasetReader,
    phase: rasterio.DatasetReader,
    coherence: rasterio.DatasetReader,
) -> dict[str, str]:
    """Fail closed on missing or conflicting complex-product provenance."""
    if not _same_grid(cross, phase):
        raise ValueError("Fine cross-product and phase grids do not match")
    if cross.count != 2 or cross.descriptions != (
        "cross_product_real", "cross_product_imaginary"
    ):
        raise ValueError("Fine cross-product requires identified real and imaginary bands")
    tags = cross.tags()
    required_values = {
        "complex_representation": "real_imaginary_bands",
        "amplitude_source": "unscaled_reference_times_conjugate_secondary",
        "phase_convention": "reference_times_conjugate_secondary",
        "spatial_aggregation": "none_native_geocoded_gslc_grid",
    }
    for key, expected in required_values.items():
        if tags.get(key) != expected:
            raise ValueError(f"Fine cross-product has missing or invalid {key} provenance")
    for key in (
        "reference_date", "secondary_date", "reference_granule", "secondary_granule",
        "frequency", "polarization", "phase_convention",
    ):
        values = (tags.get(key), phase.tags().get(key), coherence.tags().get(key))
        if not all(values) or len(set(values)) != 1:
            raise ValueError(f"Fine cross-product/phase/coherence {key} tags do not match")
    return tags


def _validate_mission_pair_dates(
    fine_phase: rasterio.DatasetReader,
    fine_coherence: rasterio.DatasetReader,
    mission_phase: rasterio.DatasetReader,
    mission_coherence: rasterio.DatasetReader,
) -> tuple[str, str]:
    """Require the same ordered acquisition days without fitting the sign.

    Product footprint start times can differ slightly, so acquisition-day
    agreement is the check here. This does not prove matching processing
    versions or identical spatial support.
    """
    mission_dates = _require_pair_dates(mission_phase, mission_coherence, "Mission")
    fine_dates = _require_pair_dates(fine_phase, fine_coherence, "Fine")
    try:
        fine_days = tuple(date.fromisoformat(value[:10]) for value in fine_dates)
        mission_days = tuple(date.fromisoformat(value[:10]) for value in mission_dates)
    except ValueError as exc:
        raise ValueError("Fine and mission acquisition dates must use ISO dates") from exc
    if fine_days != mission_days or fine_days[0] >= fine_days[1]:
        raise ValueError("Fine and mission ordered acquisition dates do not match")
    return mission_dates


def _pearson_from_sums(
    count: int,
    sum_x: float,
    sum_y: float,
    sum_x2: float,
    sum_y2: float,
    sum_xy: float,
) -> float:
    if count < 2:
        return float("nan")
    numerator = count * sum_xy - sum_x * sum_y
    denominator = math.sqrt(
        max(0.0, count * sum_x2 - sum_x * sum_x)
        * max(0.0, count * sum_y2 - sum_y * sum_y)
    )
    return numerator / denominator if denominator > 0 else float("nan")


def _circular_offset(complex_sum: complex) -> float:
    return float(np.angle(complex_sum)) if complex_sum != 0 else float("nan")


def compare_fine_to_mission(
    fine_phase_path: Path,
    fine_coherence_path: Path,
    mission_phase_path: Path,
    mission_coherence_path: Path,
    output_json: Optional[Path] = None,
    min_coherence: float = 0.25,
    min_phasor_strength: float = 0.25,
    samples_per_block: int = 5000,
    aggregation: str = "unit_phasor",
    fine_cross_product_path: Optional[Path] = None,
) -> dict:
    """Compare a fine GSLC pair against the matching mission GUNW pair.

    Only a constant circular phase reference is removed. No scale, plane,
    regression, or Tay-based fit is permitted.

    ``complex_weighted`` requires the original unscaled cross-product export.
    Coherence must not be used to reconstruct its amplitudes. The default
    ``unit_phasor`` preserves the historical diagnostic for existing callers.
    """
    if not 0 <= min_coherence <= 1:
        raise ValueError("min_coherence must be between 0 and 1")
    if not 0 <= min_phasor_strength <= 1:
        raise ValueError("min_phasor_strength must be between 0 and 1")
    if samples_per_block < 1:
        raise ValueError("samples_per_block must be positive")
    if aggregation not in {"unit_phasor", "complex_weighted"}:
        raise ValueError("aggregation must be unit_phasor or complex_weighted")
    if (aggregation == "complex_weighted") != (fine_cross_product_path is not None):
        raise ValueError(
            "fine_cross_product_path is required only for complex_weighted aggregation"
        )

    fine_phase_path = Path(fine_phase_path)
    fine_coherence_path = Path(fine_coherence_path)
    mission_phase_path = Path(mission_phase_path)
    mission_coherence_path = Path(mission_coherence_path)

    with (
        ExitStack() as optional_inputs,
        rasterio.open(fine_phase_path) as fine_phase,
        rasterio.open(fine_coherence_path) as fine_coherence,
        rasterio.open(mission_phase_path) as mission_phase,
        rasterio.open(mission_coherence_path) as mission_coherence,
    ):
        if not _same_grid(fine_phase, fine_coherence):
            raise ValueError("Fine phase and coherence grids do not match")
        if not _same_grid(mission_phase, mission_coherence):
            raise ValueError("Mission phase and coherence grids do not match")
        fine_cross_product = (
            optional_inputs.enter_context(rasterio.open(fine_cross_product_path))
            if fine_cross_product_path is not None else None
        )
        cross_provenance = (
            _validate_cross_product_provenance(fine_cross_product, fine_phase, fine_coherence)
            if fine_cross_product is not None else None
        )
        mission_pair_dates = (
            _validate_mission_pair_dates(fine_phase, fine_coherence, mission_phase, mission_coherence)
        )

        phase_complex_sums = {1: 0.0j, -1: 0.0j}
        phase_count = 0
        coherence_count = 0
        coherence_sums = {
            "fine": 0.0,
            "mission": 0.0,
            "fine2": 0.0,
            "mission2": 0.0,
            "cross": 0.0,
            "difference": 0.0,
            "absolute_difference": 0.0,
            "squared_difference": 0.0,
            "phasor_strength": 0.0,
        }
        thresholds = (0.25, 0.35, 0.5, 0.7)
        threshold_counts = {
            threshold: {"fine": 0, "mission": 0, "both": 0}
            for threshold in thresholds
        }
        coherence_samples_fine: list[np.ndarray] = []
        coherence_samples_mission: list[np.ndarray] = []

        windows = [window for _, window in mission_phase.block_windows(1)]
        for mission_window in windows:
            aggregated_phase, aggregated_coherence, phasor_strength = (
                _reproject_fine_block(
                    fine_phase,
                    fine_coherence,
                    mission_phase,
                    mission_window,
                    fine_cross_product,
                )
            )
            mission_phase_block = mission_phase.read(1, window=mission_window)
            mission_coherence_block = mission_coherence.read(
                1, window=mission_window
            )
            mission_phase_valid = _valid_raster_values(
                mission_phase_block, mission_phase.nodata
            )
            mission_coherence_valid = _valid_raster_values(
                mission_coherence_block, mission_coherence.nodata
            )
            mission_coherence_valid &= (mission_coherence_block >= 0) & (mission_coherence_block <= 1)
            common = (
                np.isfinite(aggregated_phase)
                & np.isfinite(aggregated_coherence)
                & np.isfinite(phasor_strength)
                & mission_phase_valid
                & mission_coherence_valid
            )
            if not np.any(common):
                continue

            fine_coh = aggregated_coherence[common].astype(np.float64)
            mission_coh = mission_coherence_block[common].astype(np.float64)
            strength = phasor_strength[common].astype(np.float64)
            difference = fine_coh - mission_coh
            count = fine_coh.size
            coherence_count += count
            coherence_sums["fine"] += float(np.sum(fine_coh))
            coherence_sums["mission"] += float(np.sum(mission_coh))
            coherence_sums["fine2"] += float(np.sum(fine_coh * fine_coh))
            coherence_sums["mission2"] += float(np.sum(mission_coh * mission_coh))
            coherence_sums["cross"] += float(np.sum(fine_coh * mission_coh))
            coherence_sums["difference"] += float(np.sum(difference))
            coherence_sums["absolute_difference"] += float(
                np.sum(np.abs(difference))
            )
            coherence_sums["squared_difference"] += float(
                np.sum(difference * difference)
            )
            coherence_sums["phasor_strength"] += float(np.sum(strength))
            for threshold in thresholds:
                fine_above = fine_coh >= threshold
                mission_above = mission_coh >= threshold
                threshold_counts[threshold]["fine"] += int(np.sum(fine_above))
                threshold_counts[threshold]["mission"] += int(
                    np.sum(mission_above)
                )
                threshold_counts[threshold]["both"] += int(
                    np.sum(fine_above & mission_above)
                )

            sample_stride = max(1, count // samples_per_block)
            coherence_samples_fine.append(fine_coh[::sample_stride][:samples_per_block])
            coherence_samples_mission.append(
                mission_coh[::sample_stride][:samples_per_block]
            )

            high_quality = common.copy()
            high_quality[common] &= (
                (fine_coh >= min_coherence)
                & (mission_coh >= min_coherence)
                & (strength >= min_phasor_strength)
            )
            if np.any(high_quality):
                fine_values = aggregated_phase[high_quality].astype(np.float64)
                mission_values = mission_phase_block[high_quality].astype(np.float64)
                for sign in (1, -1):
                    differences = wrap_phase(sign * fine_values - mission_values)
                    phase_complex_sums[sign] += complex(
                        np.sum(np.exp(1j * differences))
                    )
                phase_count += fine_values.size

        if coherence_count == 0:
            raise RuntimeError("No overlapping valid coherence pixels")
        if phase_count == 0:
            raise RuntimeError("No phase pixels passed the comparison thresholds")

        offsets = {
            sign: _circular_offset(phase_complex_sums[sign]) for sign in (1, -1)
        }
        phase_metrics = {
            sign: {
                "absolute": 0.0,
                "squared": 0.0,
                "complex_sum": 0.0j,
                "count": 0,
                "samples": [],
            }
            for sign in (1, -1)
        }

        # Second pass computes residuals after each candidate circular offset.
        for mission_window in windows:
            aggregated_phase, aggregated_coherence, phasor_strength = (
                _reproject_fine_block(
                    fine_phase,
                    fine_coherence,
                    mission_phase,
                    mission_window,
                    fine_cross_product,
                )
            )
            mission_phase_block = mission_phase.read(1, window=mission_window)
            mission_coherence_block = mission_coherence.read(
                1, window=mission_window
            )
            common = (
                np.isfinite(aggregated_phase)
                & np.isfinite(aggregated_coherence)
                & np.isfinite(phasor_strength)
                & _valid_raster_values(mission_phase_block, mission_phase.nodata)
                & _valid_raster_values(
                    mission_coherence_block, mission_coherence.nodata
                )
                & (mission_coherence_block >= 0)
                & (mission_coherence_block <= 1)
                & (aggregated_coherence >= min_coherence)
                & (mission_coherence_block >= min_coherence)
                & (phasor_strength >= min_phasor_strength)
            )
            if not np.any(common):
                continue
            fine_values = aggregated_phase[common].astype(np.float64)
            mission_values = mission_phase_block[common].astype(np.float64)
            for sign in (1, -1):
                residual = wrap_phase(
                    sign * fine_values - mission_values - offsets[sign]
                )
                metric = phase_metrics[sign]
                metric["absolute"] += float(np.sum(np.abs(residual)))
                metric["squared"] += float(np.sum(residual * residual))
                metric["complex_sum"] += complex(np.sum(np.exp(1j * residual)))
                metric["count"] += residual.size
                stride = max(1, residual.size // samples_per_block)
                metric["samples"].append(residual[::stride][:samples_per_block])

    coherence_sample_x = np.concatenate(coherence_samples_fine)
    coherence_sample_y = np.concatenate(coherence_samples_mission)
    if np.ptp(coherence_sample_x) == 0 or np.ptp(coherence_sample_y) == 0:
        spearman = float("nan")
    else:
        spearman = spearmanr(coherence_sample_x, coherence_sample_y).statistic
    coherence_report = {
        "valid_pixels": coherence_count,
        "fine_mean": coherence_sums["fine"] / coherence_count,
        "mission_mean": coherence_sums["mission"] / coherence_count,
        "bias_fine_minus_mission": coherence_sums["difference"] / coherence_count,
        "mae": coherence_sums["absolute_difference"] / coherence_count,
        "rmse": math.sqrt(
            coherence_sums["squared_difference"] / coherence_count
        ),
        "pearson": _pearson_from_sums(
            coherence_count,
            coherence_sums["fine"],
            coherence_sums["mission"],
            coherence_sums["fine2"],
            coherence_sums["mission2"],
            coherence_sums["cross"],
        ),
        "spearman_sample": float(spearman),
        "spearman_sample_pixels": int(coherence_sample_x.size),
        "mean_circular_aggregation_strength": (
            coherence_sums["phasor_strength"] / coherence_count
        ),
        "threshold_coverage": {
            str(threshold): {
                "fine_fraction": counts["fine"] / coherence_count,
                "mission_fraction": counts["mission"] / coherence_count,
                "both_fraction": counts["both"] / coherence_count,
            }
            for threshold, counts in threshold_counts.items()
        },
    }

    candidate_reports = {}
    for sign in (1, -1):
        metric = phase_metrics[sign]
        count = metric["count"]
        samples = np.concatenate(metric["samples"])
        candidate_reports["same_sign" if sign == 1 else "opposite_sign"] = {
            "fine_phase_multiplier": sign,
            "constant_circular_offset_rad": offsets[sign],
            "valid_pixels": count,
            "circular_mae_rad": metric["absolute"] / count,
            "circular_rmse_rad": math.sqrt(metric["squared"] / count),
            "circular_p95_abs_rad_sample": float(
                np.percentile(np.abs(samples), 95)
            ),
            "resultant_length": abs(metric["complex_sum"]) / count,
            "sample_pixels": int(samples.size),
        }
    selected_name = min(
        candidate_reports,
        key=lambda name: candidate_reports[name]["circular_mae_rad"],
    )
    selected_phase = candidate_reports[selected_name]
    coherence_gate_checks = {
        "at_least_1000000_valid_pixels": coherence_count >= 1_000_000,
        "absolute_bias_at_most_0.03": (
            abs(coherence_report["bias_fine_minus_mission"]) <= 0.03
        ),
        "mae_at_most_0.10": coherence_report["mae"] <= 0.10,
        "pearson_at_least_0.85": bool(
            np.isfinite(coherence_report["pearson"])
            and coherence_report["pearson"] >= 0.85
        ),
    }
    phase_gate_checks = {
        "expected_same_sign_selected": selected_name == "same_sign",
        "at_least_1000000_valid_pixels": (
            selected_phase["valid_pixels"] >= 1_000_000
        ),
        "circular_mae_at_most_0.10_rad": (
            selected_phase["circular_mae_rad"] <= 0.10
        ),
        "circular_rmse_at_most_0.15_rad": (
            selected_phase["circular_rmse_rad"] <= 0.15
        ),
        "resultant_length_at_least_0.98": (
            selected_phase["resultant_length"] >= 0.98
        ),
    }

    report = {
        "method": (
            f"fine phase aggregated using {aggregation}; both signs tested after "
            "one constant circular reference offset; no scale or planar fit"
        ),
        "aggregation": {
            "mode": aggregation,
            "operator": (
                "arg(area_mean(reference_times_conjugate_secondary))"
                if aggregation == "complex_weighted"
                else "arg(area_mean(exp(i * fine_phase)))"
            ),
            "strength_definition": "abs(area_mean(z)) / area_mean(abs(z))",
            "coherence_operator": (
                "area_mean(local_fine_coherence), not mission coherence reconstruction"
            ),
            "mission_operator_equivalence": False,
            "caveat": (
                "Aggregation occurs after GSLC geocoding, not mission radar-coordinate "
                "multilooking followed by geocoding. Improved agreement would not establish "
                "operator equivalence or independently validated displacement."
            ),
            "cross_product_provenance": cross_provenance,
            "verified_ordered_mission_acquisition_dates": mission_pair_dates,
        },
        "inputs": {
            "fine_phase": str(fine_phase_path),
            "fine_coherence": str(fine_coherence_path),
            "mission_phase": str(mission_phase_path),
            "mission_coherence": str(mission_coherence_path),
            "fine_cross_product": (
                str(fine_cross_product_path) if fine_cross_product_path is not None else None
            ),
        },
        "thresholds": {
            "minimum_coherence": min_coherence,
            "minimum_circular_aggregation_strength": min_phasor_strength,
        },
        "phase": {
            "selected_convention": selected_name,
            "selected_fine_phase_multiplier": candidate_reports[selected_name][
                "fine_phase_multiplier"
            ],
            "candidates": candidate_reports,
        },
        "coherence": coherence_report,
        "release_gates": {
            "supporting_coherence": {
                "passed": all(coherence_gate_checks.values()),
                "checks": coherence_gate_checks,
                "interpretation": (
                    "Passing permits a quality/observability layer only."
                ),
            },
            "displacement_phase": {
                "passed": all(phase_gate_checks.values()),
                "checks": phase_gate_checks,
                "interpretation": (
                    "Passing is necessary but not sufficient for fine "
                    "displacement; unwrapping, corrections, datum control, "
                    "temporal redundancy, and independent validation remain "
                    "required."
                ),
            },
        },
    }
    if output_json is not None:
        output_json = Path(output_json)
        output_json.parent.mkdir(parents=True, exist_ok=True)
        output_json.write_text(json.dumps(report, indent=2) + "\n")
        logger.info("Wrote NISAR validation report: %s", output_json)
    return report


def _require_pair_dates(
    phase: rasterio.DatasetReader,
    coherence: rasterio.DatasetReader,
    label: str,
) -> tuple[str, str]:
    """Return matching acquisition-date tags for a phase/coherence pair."""
    phase_tags = phase.tags()
    coherence_tags = coherence.tags()
    phase_dates = (
        phase_tags.get("reference_date"),
        phase_tags.get("secondary_date"),
    )
    coherence_dates = (
        coherence_tags.get("reference_date"),
        coherence_tags.get("secondary_date"),
    )
    if not all(phase_dates):
        raise ValueError(f"{label} phase raster lacks acquisition-date tags")
    if not all(coherence_dates):
        raise ValueError(f"{label} coherence raster lacks acquisition-date tags")
    if phase_dates != coherence_dates:
        raise ValueError(
            f"{label} phase/coherence acquisition-date tags do not match"
        )
    return str(phase_dates[0]), str(phase_dates[1])


def _closure_metrics(
    *,
    count: int,
    absolute_sum: float,
    squared_sum: float,
    complex_sum: complex,
    maximum_absolute: float,
    threshold_counts: dict[float, int],
    samples: list[np.ndarray],
) -> dict:
    """Finalize one set of streaming phase-closure statistics."""
    if count == 0:
        return {
            "valid_pixels": 0,
            "circular_mean_rad": None,
            "circular_mae_rad": None,
            "circular_rmse_rad": None,
            "circular_max_abs_rad": None,
            "resultant_length": None,
            "circular_p95_abs_rad_sample": None,
            "circular_p99_abs_rad_sample": None,
            "sample_pixels": 0,
            "fraction_within_rad": {
                str(threshold): None for threshold in threshold_counts
            },
        }

    sample = np.concatenate(samples) if samples else np.empty(0, dtype=np.float32)
    return {
        "valid_pixels": count,
        "circular_mean_rad": float(np.angle(complex_sum)),
        "circular_mae_rad": absolute_sum / count,
        "circular_rmse_rad": math.sqrt(squared_sum / count),
        "circular_max_abs_rad": maximum_absolute,
        "resultant_length": abs(complex_sum) / count,
        "circular_p95_abs_rad_sample": (
            float(np.percentile(np.abs(sample), 95)) if sample.size else None
        ),
        "circular_p99_abs_rad_sample": (
            float(np.percentile(np.abs(sample), 99)) if sample.size else None
        ),
        "sample_pixels": int(sample.size),
        "fraction_within_rad": {
            str(threshold): threshold_count / count
            for threshold, threshold_count in threshold_counts.items()
        },
    }


def audit_three_date_phase_closure(
    phase_ab_path: Path,
    coherence_ab_path: Path,
    phase_bc_path: Path,
    coherence_bc_path: Path,
    phase_ac_path: Path,
    coherence_ac_path: Path,
    output_json: Optional[Path] = None,
    min_coherence: float = 0.25,
    samples_per_block: int = 2000,
) -> dict:
    """Audit wrapped phase closure for a native-grid three-date GSLC triangle.

    With the PeatGuard convention ``earlier * conjugate(later)``, valid wrapped
    phases should satisfy ``wrap(phi_ab + phi_bc - phi_ac) == 0``. This is an
    internal algebra and grid-integrity test. It does not establish geophysical
    accuracy, remove atmospheric effects, or make a short commissioning stack
    suitable for velocity estimation.
    """
    if not 0 <= min_coherence <= 1:
        raise ValueError("min_coherence must be between 0 and 1")
    if samples_per_block < 1:
        raise ValueError("samples_per_block must be positive")

    input_paths = {
        "phase_ab": Path(phase_ab_path),
        "coherence_ab": Path(coherence_ab_path),
        "phase_bc": Path(phase_bc_path),
        "coherence_bc": Path(coherence_bc_path),
        "phase_ac": Path(phase_ac_path),
        "coherence_ac": Path(coherence_ac_path),
    }
    with (
        rasterio.open(input_paths["phase_ab"]) as phase_ab,
        rasterio.open(input_paths["coherence_ab"]) as coherence_ab,
        rasterio.open(input_paths["phase_bc"]) as phase_bc,
        rasterio.open(input_paths["coherence_bc"]) as coherence_bc,
        rasterio.open(input_paths["phase_ac"]) as phase_ac,
        rasterio.open(input_paths["coherence_ac"]) as coherence_ac,
    ):
        datasets = {
            "phase_ab": phase_ab,
            "coherence_ab": coherence_ab,
            "phase_bc": phase_bc,
            "coherence_bc": coherence_bc,
            "phase_ac": phase_ac,
            "coherence_ac": coherence_ac,
        }
        for label, dataset in datasets.items():
            if not _same_grid(phase_ab, dataset):
                raise ValueError(f"Closure raster grid mismatch: phase_ab != {label}")

        dates_ab = _require_pair_dates(phase_ab, coherence_ab, "A-B")
        dates_bc = _require_pair_dates(phase_bc, coherence_bc, "B-C")
        dates_ac = _require_pair_dates(phase_ac, coherence_ac, "A-C")
        if not (
            dates_ab[0] == dates_ac[0]
            and dates_ab[1] == dates_bc[0]
            and dates_bc[1] == dates_ac[1]
        ):
            raise ValueError(
                "Input acquisition-date tags do not form A-B, B-C, A-C: "
                f"{dates_ab}, {dates_bc}, {dates_ac}"
            )

        thresholds = (0.01, 0.05, 0.1, 0.25, 0.5)
        accumulators = {
            "all_overlap": {
                "count": 0,
                "absolute_sum": 0.0,
                "squared_sum": 0.0,
                "complex_sum": 0.0j,
                "maximum_absolute": 0.0,
                "threshold_counts": {threshold: 0 for threshold in thresholds},
                "samples": [],
            },
            "coherence_masked": {
                "count": 0,
                "absolute_sum": 0.0,
                "squared_sum": 0.0,
                "complex_sum": 0.0j,
                "maximum_absolute": 0.0,
                "threshold_counts": {threshold: 0 for threshold in thresholds},
                "samples": [],
            },
        }
        total_pixels = phase_ab.width * phase_ab.height

        for _, window in phase_ab.block_windows(1):
            phase_blocks = [
                dataset.read(1, window=window)
                for dataset in (phase_ab, phase_bc, phase_ac)
            ]
            coherence_blocks = [
                dataset.read(1, window=window)
                for dataset in (coherence_ab, coherence_bc, coherence_ac)
            ]
            valid = np.ones(phase_blocks[0].shape, dtype=bool)
            for data, dataset in zip(
                phase_blocks,
                (phase_ab, phase_bc, phase_ac),
            ):
                valid &= _valid_raster_values(data, dataset.nodata)
            for data, dataset in zip(
                coherence_blocks,
                (coherence_ab, coherence_bc, coherence_ac),
            ):
                valid &= _valid_raster_values(data, dataset.nodata)
                valid &= (data >= 0.0) & (data <= 1.0)
            if not np.any(valid):
                continue

            residual = wrap_phase(
                phase_blocks[0].astype(np.float64)
                + phase_blocks[1].astype(np.float64)
                - phase_blocks[2].astype(np.float64)
            )
            quality = valid.copy()
            for coherence in coherence_blocks:
                quality &= coherence >= min_coherence

            for name, mask in (
                ("all_overlap", valid),
                ("coherence_masked", quality),
            ):
                if not np.any(mask):
                    continue
                values = residual[mask]
                absolute = np.abs(values)
                accumulator = accumulators[name]
                count = int(values.size)
                accumulator["count"] += count
                accumulator["absolute_sum"] += float(np.sum(absolute))
                accumulator["squared_sum"] += float(np.sum(values * values))
                accumulator["complex_sum"] += complex(
                    np.sum(np.exp(1j * values))
                )
                accumulator["maximum_absolute"] = max(
                    accumulator["maximum_absolute"],
                    float(np.max(absolute)),
                )
                for threshold in thresholds:
                    accumulator["threshold_counts"][threshold] += int(
                        np.sum(absolute <= threshold)
                    )
                stride = max(1, count // samples_per_block)
                accumulator["samples"].append(
                    values[::stride][:samples_per_block].astype(np.float32)
                )

    metrics = {
        name: _closure_metrics(**accumulator)
        for name, accumulator in accumulators.items()
    }
    quality_metrics = metrics["coherence_masked"]
    gate_checks = {
        "at_least_10000_coherence_masked_pixels": (
            quality_metrics["valid_pixels"] >= 10_000
        ),
        "circular_mae_at_most_0.05_rad": (
            quality_metrics["circular_mae_rad"] is not None
            and quality_metrics["circular_mae_rad"] <= 0.05
        ),
        "circular_rmse_at_most_0.10_rad": (
            quality_metrics["circular_rmse_rad"] is not None
            and quality_metrics["circular_rmse_rad"] <= 0.10
        ),
        "resultant_length_at_least_0.99": (
            quality_metrics["resultant_length"] is not None
            and quality_metrics["resultant_length"] >= 0.99
        ),
        "at_least_99_percent_within_0.10_rad": (
            quality_metrics["fraction_within_rad"]["0.1"] is not None
            and quality_metrics["fraction_within_rad"]["0.1"] >= 0.99
        ),
    }
    report = {
        "method": (
            "native-grid wrapped closure wrap(phi_ab + phi_bc - phi_ac); "
            "no fitted offset, scale, plane, or external calibration"
        ),
        "interpretation_limit": (
            "internal interferogram/grid integrity only; not a geophysical "
            "accuracy, atmospheric-correction, unwrapping, or velocity test"
        ),
        "inputs": {name: str(path) for name, path in input_paths.items()},
        "dates": {
            "a": dates_ab[0],
            "b": dates_ab[1],
            "c": dates_bc[1],
        },
        "grid": {
            "crs": str(phase_ab.crs),
            "width": phase_ab.width,
            "height": phase_ab.height,
            "transform": tuple(phase_ab.transform),
            "total_pixels": total_pixels,
        },
        "thresholds": {"minimum_coherence_each_pair": min_coherence},
        "metrics": metrics,
        "gate": {
            "passed": all(gate_checks.values()),
            "checks": gate_checks,
        },
    }
    if output_json is not None:
        output_json = Path(output_json)
        output_json.parent.mkdir(parents=True, exist_ok=True)
        output_json.write_text(json.dumps(report, indent=2) + "\n")
        logger.info("Wrote NISAR phase-closure report: %s", output_json)
    return report
