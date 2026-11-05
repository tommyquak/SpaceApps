"""Tests for restartable NISAR GSLC AOI staging."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np

from peatguard.ingest.nisar_gslc import (
    COMPLEX_CACHE_BLOCK_BYTES,
    MASK_CACHE_BLOCK_BYTES,
    _dataset_resume_row,
    _is_expired_signed_url_error,
    _remote_cache_settings,
    subset_open_gslc,
)


def _write_source(path: Path) -> None:
    rows, columns = 100, 120
    with h5py.File(path, "w") as source:
        identification = source.require_group("science/LSAR/identification")
        identification.create_dataset(
            "zeroDopplerStartTime", data=np.bytes_("2026-06-22T22:07:11")
        )
        identification.create_dataset("trackNumber", data=np.uint32(125))
        identification.create_dataset("frameNumber", data=np.uint32(176))
        identification.create_dataset("granuleId", data=np.bytes_("TEST_GSLC"))

        grid = source.require_group("science/LSAR/GSLC/grids/frequencyA")
        grid.create_dataset(
            "xCoordinates", data=114.0 + (np.arange(columns) + 0.5) * 0.001
        )
        grid.create_dataset(
            "yCoordinates", data=-2.0 - (np.arange(rows) + 0.5) * 0.001
        )
        grid.create_dataset("projection", data=np.uint32(4326))
        grid.create_dataset("xCoordinateSpacing", data=0.001)
        grid.create_dataset("yCoordinateSpacing", data=-0.001)
        phase = np.arange(rows * columns, dtype=np.float32).reshape(rows, columns) / 100
        grid.create_dataset(
            "HH",
            data=np.exp(1j * phase).astype(np.complex64),
            chunks=(20, 20),
            compression="gzip",
        )
        grid.create_dataset(
            "mask",
            data=np.ones((rows, columns), dtype=np.uint8),
            chunks=(20, 20),
            compression="gzip",
        )


def test_subset_open_gslc_preserves_paths_and_values(tmp_path):
    source_path = tmp_path / "source.h5"
    output_path = tmp_path / "subset.h5"
    _write_source(source_path)

    with h5py.File(source_path, "r") as source:
        result = subset_open_gslc(
            source,
            output_path,
            source_name="source.h5",
            bbox_wgs84=(114.02, -2.08, 114.08, -2.02),
            block_rows=13,
        )

    assert result == output_path
    assert not output_path.with_suffix(".h5.part").exists()
    with h5py.File(output_path, "r") as subset:
        assert bool(subset.attrs["peatguard_complete"])
        assert int(subset.attrs["peatguard_rows_completed"]) == 60
        assert int(subset.attrs["peatguard_rows_completed__HH"]) == 60
        assert int(subset.attrs["peatguard_rows_completed__mask"]) == 60
        grid = subset["science/LSAR/GSLC/grids/frequencyA"]
        assert grid["HH"].shape == (60, 60)
        assert grid["mask"].shape == (60, 60)
        assert np.all(grid["mask"][:] == 1)
        with h5py.File(source_path, "r") as source:
            expected = source["science/LSAR/GSLC/grids/frequencyA/HH"][
                20:80, 20:80
            ]
        np.testing.assert_allclose(grid["HH"][:], expected)


def test_subset_open_gslc_reuses_valid_completed_output(tmp_path):
    source_path = tmp_path / "source.h5"
    output_path = tmp_path / "subset.h5"
    _write_source(source_path)
    bbox = (114.02, -2.08, 114.08, -2.02)

    with h5py.File(source_path, "r") as source:
        first = subset_open_gslc(
            source, output_path, "source.h5", bbox_wgs84=bbox
        )
    original_mtime = first.stat().st_mtime_ns
    with h5py.File(source_path, "r") as source:
        second = subset_open_gslc(
            source, output_path, "source.h5", bbox_wgs84=bbox
        )

    assert second == first
    assert second.stat().st_mtime_ns == original_mtime


def test_expired_signed_url_detection():
    class StatusError(Exception):
        status = 403

    assert _is_expired_signed_url_error(StatusError())
    assert _is_expired_signed_url_error(
        RuntimeError("403, message='Forbidden', url='https://example.invalid'")
    )
    assert not _is_expired_signed_url_error(RuntimeError("connection reset"))


def test_remote_cache_switches_after_complex_raster_checkpoint(tmp_path):
    output_path = tmp_path / "subset.h5"
    assert _remote_cache_settings(output_path, "HH") == (
        COMPLEX_CACHE_BLOCK_BYTES,
        "readahead",
        {},
    )

    part_path = output_path.with_suffix(".h5.part")
    with h5py.File(part_path, "w") as partial:
        partial.attrs["peatguard_source_window"] = np.asarray(
            (100, 200, 300, 400), dtype=np.int64
        )
        partial.attrs["peatguard_rows_completed__HH"] = 100

    assert _remote_cache_settings(output_path, "HH") == (
        MASK_CACHE_BLOCK_BYTES,
        "blockcache",
        {"maxblocks": 128},
    )


def test_legacy_checkpoint_applies_only_to_primary_complex_raster():
    attrs = {"peatguard_rows_completed": 2560}
    assert _dataset_resume_row(attrs, "HH", "HH") == 2560
    assert _dataset_resume_row(attrs, "mask", "HH") == 0
    assert _dataset_resume_row(
        attrs, "inputDataExceptionMask", "HH"
    ) == 0

    attrs["peatguard_rows_completed__mask"] = 4096
    assert _dataset_resume_row(attrs, "mask", "HH") == 4096
