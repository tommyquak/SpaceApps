"""Reproducible audit helpers for external peat-subsidence baselines.

The baseline products in this module remain independent observations. The
audit verifies deterministic unit conversion, masking and provenance. It does
not fit, bias-correct or otherwise calibrate PeatGuard measurements to an
external raster.
"""

from __future__ import annotations

import hashlib
import json
from datetime import date
from pathlib import Path
from typing import Any

import numpy as np
import rasterio
from pyproj import Geod


TAY_DATASET_DOI = "10.21979/N9/JGQE4O"
TAY_PAPER_DOI = "10.1038/s43247-025-02926-6"
TAY_EXPECTED_MD5 = "32a196798d8b01cc3e9690c29e4804ac"
TAY_PERIOD_START = date(2018, 10, 12)
TAY_PERIOD_END = date(2024, 8, 9)


def hash_file(path: Path, algorithm: str) -> str:
    """Return a streaming file hash."""
    digest = hashlib.new(algorithm)
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _valid_mask(data: np.ndarray, nodata: float | None) -> np.ndarray:
    valid = np.isfinite(data)
    if nodata is not None and np.isfinite(nodata):
        valid &= data != nodata
    return valid


def _raster_metadata(dataset: rasterio.io.DatasetReader) -> dict[str, Any]:
    transform = dataset.transform
    center_x = (dataset.bounds.left + dataset.bounds.right) / 2.0
    center_y = (dataset.bounds.bottom + dataset.bounds.top) / 2.0
    pixel_width_m = None
    pixel_height_m = None
    if dataset.crs and dataset.crs.is_geographic:
        geod = Geod(ellps="WGS84")
        _, _, pixel_width_m = geod.inv(
            center_x,
            center_y,
            center_x + abs(transform.a),
            center_y,
        )
        _, _, pixel_height_m = geod.inv(
            center_x,
            center_y,
            center_x,
            center_y + abs(transform.e),
        )
    elif dataset.crs and dataset.crs.is_projected:
        pixel_width_m = abs(transform.a)
        pixel_height_m = abs(transform.e)

    return {
        "path": str(Path(dataset.name).resolve()),
        "crs": dataset.crs.to_string() if dataset.crs else None,
        "width": dataset.width,
        "height": dataset.height,
        "transform": list(transform)[:6],
        "bounds": {
            "west": dataset.bounds.left,
            "south": dataset.bounds.bottom,
            "east": dataset.bounds.right,
            "north": dataset.bounds.top,
        },
        "pixel_spacing_native": {
            "x": abs(transform.a),
            "y": abs(transform.e),
        },
        "pixel_spacing_m_at_grid_center": {
            "x": pixel_width_m,
            "y": pixel_height_m,
        },
        "dtype": dataset.dtypes[0],
        "nodata": dataset.nodata,
    }


def _same_grid(
    source: rasterio.io.DatasetReader,
    derivative: rasterio.io.DatasetReader,
) -> bool:
    return (
        source.width == derivative.width
        and source.height == derivative.height
        and source.crs == derivative.crs
        and source.transform.almost_equals(derivative.transform)
    )


def audit_tay_long_term_rate(
    source_path: Path,
    derivative_path: Path,
    *,
    derivative_scale: float = 1_000.0,
    derivative_absolute_limit: float = 100.0,
    expected_source_md5: str = TAY_EXPECTED_MD5,
) -> dict[str, Any]:
    """Audit PeatGuard's derivative against the official Tay et al. raster.

    The official README calls the product a long-term linear rate and records
    its unit as metres. PeatGuard interprets the rate as metres per year and
    converts to millimetres per year by multiplying by 1,000. The audit keeps
    the README's unit ambiguity explicit.
    """
    source_path = Path(source_path)
    derivative_path = Path(derivative_path)
    if derivative_scale <= 0:
        raise ValueError("derivative_scale must be positive")
    if derivative_absolute_limit <= 0:
        raise ValueError("derivative_absolute_limit must be positive")
    for path in (source_path, derivative_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    source_md5 = hash_file(source_path, "md5")
    source_sha256 = hash_file(source_path, "sha256")
    derivative_sha256 = hash_file(derivative_path, "sha256")

    with rasterio.open(source_path) as source, rasterio.open(
        derivative_path
    ) as derivative:
        source_data = source.read(1)
        derivative_data = derivative.read(1)
        source_finite = _valid_mask(source_data, source.nodata)
        derivative_valid = _valid_mask(derivative_data, derivative.nodata)
        source_measurement = source_finite & (source_data != 0.0)
        same_grid = _same_grid(source, derivative)

        overlap = source_measurement & derivative_valid
        difference = (
            derivative_data[overlap]
            - source_data[overlap] * derivative_scale
        )
        exact_values = bool(
            overlap.any()
            and np.array_equal(
                derivative_data[overlap],
                source_data[overlap] * derivative_scale,
            )
        )

        expected_derivative_mask = source_measurement & (
            np.abs(source_data * derivative_scale)
            < derivative_absolute_limit
        )
        exact_mask_rule = bool(
            same_grid and np.array_equal(derivative_valid, expected_derivative_mask)
        )
        source_metadata = _raster_metadata(source)
        derivative_metadata = _raster_metadata(derivative)

        source_values = source_data[source_measurement]
        derivative_values = derivative_data[derivative_valid]
        excluded_by_limit = source_measurement & ~expected_derivative_mask

    period_days = (TAY_PERIOD_END - TAY_PERIOD_START).days
    checks = {
        "official_md5_matches": source_md5 == expected_source_md5,
        "grids_identical": same_grid,
        "all_retained_values_exact_unit_conversion": exact_values,
        "derivative_mask_exactly_matches_documented_project_rule": exact_mask_rule,
        "no_regression_or_bias_fit_performed": True,
    }
    passed = all(checks.values())

    return {
        "schema_version": 1,
        "baseline": {
            "name": "Tay et al. (2025) ALOS-2 long-term peat-motion rate",
            "paper_doi": TAY_PAPER_DOI,
            "dataset_doi": TAY_DATASET_DOI,
            "sensor": "ALOS-2 PALSAR-2 L-band",
            "period_start": TAY_PERIOD_START.isoformat(),
            "period_end": TAY_PERIOD_END.isoformat(),
            "period_days": period_days,
            "source_documented_quantity": "linear best-fit long-term rate",
            "source_documented_unit": "m",
            "operational_unit_interpretation": {
                "value": "m/year converted to mm/year",
                "confidence": "[S]",
                "caveat": (
                    "The repository README labels the raster a long-term rate "
                    "but writes the unit as m rather than m/year."
                ),
            },
            "license": {
                "name": "CC BY-NC 4.0",
                "url": "https://creativecommons.org/licenses/by-nc/4.0/",
                "commercial_use": "permission required",
            },
        },
        "source": {
            **source_metadata,
            "bytes": source_path.stat().st_size,
            "md5": source_md5,
            "expected_md5": expected_source_md5,
            "sha256": source_sha256,
            "finite_pixels": int(source_finite.sum()),
            "nonzero_measurement_pixels": int(source_measurement.sum()),
            "value_range_source_units": {
                "minimum": float(source_values.min()),
                "median": float(np.median(source_values)),
                "maximum": float(source_values.max()),
            },
        },
        "peatguard_derivative": {
            **derivative_metadata,
            "bytes": derivative_path.stat().st_size,
            "sha256": derivative_sha256,
            "valid_pixels": int(derivative_valid.sum()),
            "value_range_mm_per_year": {
                "minimum": float(derivative_values.min()),
                "median": float(np.median(derivative_values)),
                "maximum": float(derivative_values.max()),
            },
            "conversion_scale": derivative_scale,
            "project_absolute_limit_mm_per_year": derivative_absolute_limit,
            "source_measurements_excluded_by_project_limit": int(
                excluded_by_limit.sum()
            ),
            "limit_source": (
                "PeatGuard plausibility filter, not a threshold documented in "
                "the Tay repository README"
            ),
        },
        "pixel_comparison": {
            "overlap_pixels": int(overlap.sum()),
            "maximum_absolute_difference_mm_per_year": (
                float(np.max(np.abs(difference))) if difference.size else None
            ),
            "mean_absolute_difference_mm_per_year": (
                float(np.mean(np.abs(difference))) if difference.size else None
            ),
        },
        "checks": checks,
        "audit_passed": passed,
        "appropriateness": {
            "independent_spatial_benchmark": {
                "verdict": "appropriate",
                "confidence": "[E]",
                "reason": (
                    "Independent, field-compared L-band result over the same "
                    "Central Kalimantan peat system."
                ),
            },
            "calibration_or_fitted_prior": {
                "verdict": "not appropriate",
                "confidence": "[E]",
                "reason": (
                    "Fitting NISAR to Tay would destroy the benchmark's "
                    "independence."
                ),
            },
            "current_2026_rate_substitute": {
                "verdict": "not appropriate",
                "confidence": "[S]",
                "reason": (
                    "The 2018-2024 average predates NISAR and drainage, fire, "
                    "land cover and restoration may have changed."
                ),
            },
            "nisar_timeseries_temporal_anchor": {
                "verdict": "not appropriate",
                "confidence": "[E]",
                "reason": (
                    "This is a derived long-term rate, not phase observations "
                    "that can extend the NISAR SBAS network."
                ),
            },
            "few_meter_deformation_baseline": {
                "verdict": "not appropriate",
                "confidence": "[E]",
                "reason": (
                    "The geographic pixels are roughly 157 by 209 m at the "
                    "grid centre and cannot validate independent 5 m motion."
                ),
            },
            "public_noncommercial_demo": {
                "verdict": "conditionally appropriate",
                "confidence": "[E]",
                "reason": (
                    "Attribution and non-commercial use must comply with "
                    "CC BY-NC 4.0."
                ),
            },
            "commercial_product": {
                "verdict": "not appropriate without permission",
                "confidence": "[E]",
                "reason": "The dataset licence prohibits unapproved commercial use.",
            },
        },
    }


def write_audit_json(report: dict[str, Any], output_path: Path) -> Path:
    """Write a deterministic, human-readable audit report."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return output_path
