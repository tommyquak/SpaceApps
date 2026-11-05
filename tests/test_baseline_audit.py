"""Tests for independent alternate-baseline provenance."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin

from peatguard.analysis.baseline_audit import audit_tay_long_term_rate, hash_file


def _write_raster(path: Path, data: np.ndarray, nodata: float | None) -> None:
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=data.shape[1],
        height=data.shape[0],
        count=1,
        dtype="float32",
        crs="EPSG:4326",
        transform=from_origin(114.0, -2.0, 0.001, 0.001),
        nodata=nodata,
    ) as output:
        output.write(data.astype(np.float32), 1)


def test_tay_audit_verifies_exact_scale_and_mask(tmp_path: Path) -> None:
    source_path = tmp_path / "source.tif"
    derivative_path = tmp_path / "derivative.tif"
    source = np.array(
        [
            [np.nan, 0.0, -0.010],
            [0.015, 0.100, -0.101],
        ],
        dtype=np.float32,
    )
    derivative = np.array(
        [
            [-9999.0, -9999.0, -10.0],
            [15.0, -9999.0, -9999.0],
        ],
        dtype=np.float32,
    )
    _write_raster(source_path, source, nodata=None)
    _write_raster(derivative_path, derivative, nodata=-9999.0)

    report = audit_tay_long_term_rate(
        source_path,
        derivative_path,
        expected_source_md5=hash_file(source_path, "md5"),
    )

    assert report["audit_passed"]
    assert report["pixel_comparison"]["maximum_absolute_difference_mm_per_year"] == 0.0
    assert (
        report["peatguard_derivative"][
            "source_measurements_excluded_by_project_limit"
        ]
        == 2
    )
    assert (
        report["appropriateness"]["nisar_timeseries_temporal_anchor"]["verdict"]
        == "not appropriate"
    )


def test_tay_audit_rejects_modified_derivative(tmp_path: Path) -> None:
    source_path = tmp_path / "source.tif"
    derivative_path = tmp_path / "derivative.tif"
    _write_raster(source_path, np.array([[-0.01]], dtype=np.float32), nodata=None)
    _write_raster(derivative_path, np.array([[-9.0]], dtype=np.float32), nodata=-9999.0)

    report = audit_tay_long_term_rate(
        source_path,
        derivative_path,
        expected_source_md5=hash_file(source_path, "md5"),
    )

    assert not report["audit_passed"]
    assert not report["checks"]["all_retained_values_exact_unit_conversion"]
