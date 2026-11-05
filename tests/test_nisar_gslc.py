"""Tests for block-wise NISAR GSLC interferogram formation."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest
import rasterio

from peatguard.insar.nisar_gslc import form_gslc_interferogram, inspect_gslc


def _write_gslc(
    path: Path,
    acquisition: str,
    phase: np.ndarray,
    *,
    x_offset: float = 0.0,
    track: int = 125,
    frame: int = 176,
    mask: np.ndarray | None = None,
    exception_mask: np.ndarray | None = None,
) -> None:
    height, width = phase.shape
    x = 114.0 + x_offset + (np.arange(width) + 0.5) * 0.001
    y = -2.0 - (np.arange(height) + 0.5) * 0.001
    complex_data = np.exp(1j * phase).astype(np.complex64)

    with h5py.File(path, "w") as h5:
        grid = h5.require_group("science/LSAR/GSLC/grids/frequencyA")
        grid.create_dataset("xCoordinates", data=x)
        grid.create_dataset("yCoordinates", data=y)
        grid.create_dataset("projection", data=np.uint32(4326))
        grid.create_dataset("centerFrequency", data=np.float64(1_239_000_000.0))
        grid.create_dataset("numberOfSubSwaths", data=np.uint8(1))
        grid.create_dataset("HH", data=complex_data)
        grid.create_dataset(
            "mask",
            data=(
                np.ones((height, width), dtype=np.uint8)
                if mask is None
                else mask.astype(np.uint8)
            ),
        )
        grid.create_dataset(
            "inputDataExceptionMask",
            data=(
                np.zeros((height, width), dtype=np.uint8)
                if exception_mask is None
                else exception_mask.astype(np.uint8)
            ),
        )

        identification = h5.require_group("science/LSAR/identification")
        identification.create_dataset("zeroDopplerStartTime", data=np.bytes_(acquisition))
        identification.create_dataset("trackNumber", data=np.uint32(track))
        identification.create_dataset("frameNumber", data=np.uint32(frame))
        identification.create_dataset(
            "granuleId", data=np.bytes_(f"TEST_GSLC_{acquisition[:10]}")
        )


def test_inspect_gslc_reads_grid_without_loading_raster(tmp_path):
    path = tmp_path / "reference.h5"
    _write_gslc(path, "2026-06-22T22:07:11", np.zeros((20, 30), dtype=np.float32))

    info = inspect_gslc(path)

    assert info.shape == (20, 30)
    assert info.epsg == 4326
    assert info.x_spacing == pytest.approx(0.001)
    assert info.y_spacing == pytest.approx(-0.001)
    assert info.track_number == 125
    assert info.frame_number == 176
    assert info.center_frequency_hz == pytest.approx(1_239_000_000.0)
    assert info.number_of_subswaths == 1
    assert info.wavelength_m == pytest.approx(0.2419632429378531)


def test_form_gslc_interferogram_phase_sign_and_coherence(tmp_path):
    reference = tmp_path / "reference.h5"
    secondary = tmp_path / "secondary.h5"
    phase_difference = 0.4
    shape = (48, 64)
    _write_gslc(
        reference,
        "2026-06-22T22:07:11",
        np.zeros(shape, dtype=np.float32),
    )
    _write_gslc(
        secondary,
        "2026-07-04T22:07:10",
        np.full(shape, phase_difference, dtype=np.float32),
    )

    products = form_gslc_interferogram(
        reference,
        secondary,
        tmp_path / "outputs",
        coherence_window=5,
        block_size=16,
    )
    assert set(products) == {"wrapped_phase", "coherence"}

    with rasterio.open(products["wrapped_phase"]) as phase_raster:
        phase = phase_raster.read(1)
        assert phase_raster.tags()["phase_convention"] == (
            "reference_times_conjugate_secondary"
        )
        assert float(phase_raster.tags()["wavelength_m"]) == pytest.approx(
            0.2419632429378531
        )
        assert phase_raster.tags()["effective_support"] == (
            "native_mode_dependent"
        )
        np.testing.assert_allclose(phase, -phase_difference, atol=1e-6)

    with rasterio.open(products["coherence"]) as coherence_raster:
        coherence = coherence_raster.read(1)
        assert coherence_raster.tags()["effective_support"] == (
            "5_by_5_source_pixels"
        )
        valid = coherence != coherence_raster.nodata
        # A full coherence window is unavailable along the outer raster edge.
        assert valid[2:-2, 2:-2].all()
        np.testing.assert_allclose(coherence[valid], 1.0, atol=1e-6)


def test_optional_cross_product_preserves_amplitudes_masks_and_provenance(tmp_path):
    reference = tmp_path / "reference.h5"
    secondary = tmp_path / "secondary.h5"
    shape = (16, 16)
    phase = np.zeros(shape, dtype=np.float32)
    mask = np.ones(shape, dtype=np.uint8)
    mask[6, 7] = 0
    _write_gslc(reference, "2026-06-22T22:07:11", phase, mask=mask)
    _write_gslc(secondary, "2026-07-04T22:07:10", phase)
    ref_data = np.full(shape, 3 + 0j, dtype=np.complex64)
    sec_data = np.full(shape, 3 + 0j, dtype=np.complex64)
    ref_data[:, 1::2] = 1 + 0j
    sec_data[:, 1::2] = -1j
    ref_data[4, 4] = 0j
    sec_data[5, 5] = complex(np.nan, 0)
    for path, data in ((reference, ref_data), (secondary, sec_data)):
        with h5py.File(path, "r+") as h5:
            h5["science/LSAR/GSLC/grids/frequencyA/HH"][:] = data

    products = form_gslc_interferogram(
        reference, secondary, tmp_path / "outputs", block_size=7,
        export_cross_product=True,
    )
    with rasterio.open(products["complex_cross_product"]) as cross, rasterio.open(
        products["wrapped_phase"]
    ) as phase_raster:
        data = cross.read()
        assert cross.count == 2
        assert cross.transform == phase_raster.transform
        assert cross.crs == phase_raster.crs
        assert np.isnan(cross.nodata)
        assert cross.descriptions == ("cross_product_real", "cross_product_imaginary")
        assert data[:, 2, 2].tolist() == [9.0, 0.0]
        assert data[:, 2, 3].tolist() == [0.0, 1.0]
        for row, col in ((4, 4), (5, 5), (6, 7)):
            assert np.isnan(data[:, row, col]).all()
            assert phase_raster.read(1)[row, col] == phase_raster.nodata
        valid = np.isfinite(data).all(axis=0)
        np.testing.assert_allclose(
            np.angle(data[0][valid] + 1j * data[1][valid]),
            phase_raster.read(1)[valid], atol=1e-6,
        )
        tags = cross.tags()
        assert tags["amplitude_source"] == "unscaled_reference_times_conjugate_secondary"
        assert tags["spatial_aggregation"] == "none_native_geocoded_gslc_grid"
        for key in (
            "reference_granule", "secondary_granule", "reference_date",
            "secondary_date", "phase_convention",
        ):
            assert tags[key] == phase_raster.tags()[key]


def test_form_gslc_interferogram_crops_to_bbox(tmp_path):
    reference = tmp_path / "reference.h5"
    secondary = tmp_path / "secondary.h5"
    shape = (100, 100)
    _write_gslc(reference, "2026-06-22T22:07:11", np.zeros(shape, dtype=np.float32))
    _write_gslc(secondary, "2026-07-04T22:07:10", np.zeros(shape, dtype=np.float32))

    products = form_gslc_interferogram(
        reference,
        secondary,
        tmp_path / "outputs",
        bbox_wgs84=(114.02, -2.08, 114.08, -2.02),
        block_size=32,
    )

    with rasterio.open(products["wrapped_phase"]) as output:
        assert 55 <= output.width <= 62
        assert 55 <= output.height <= 62
        assert output.bounds.left >= 114.019
        assert output.bounds.right <= 114.081


def test_form_gslc_interferogram_rejects_grid_mismatch(tmp_path):
    reference = tmp_path / "reference.h5"
    secondary = tmp_path / "secondary.h5"
    phase = np.zeros((20, 20), dtype=np.float32)
    _write_gslc(reference, "2026-06-22T22:07:11", phase)
    _write_gslc(secondary, "2026-07-04T22:07:10", phase, x_offset=0.01)

    with pytest.raises(ValueError, match="x-coordinate grids"):
        form_gslc_interferogram(reference, secondary, tmp_path / "outputs")


@pytest.mark.parametrize("shift_axis", [None, "x", "y"])
def test_form_gslc_requires_exact_utm_grid_alignment(tmp_path, shift_axis):
    """A 5 m shift is a full pixel even when coordinates are millions of metres."""
    reference = tmp_path / "reference.h5"
    secondary = tmp_path / "secondary.h5"
    phase = np.zeros((20, 20), dtype=np.float32)
    _write_gslc(reference, "2026-06-22T22:07:11", phase)
    _write_gslc(secondary, "2026-07-04T22:07:10", phase)
    for path in (reference, secondary):
        with h5py.File(path, "r+") as h5:
            grid = h5["science/LSAR/GSLC/grids/frequencyA"]
            grid["projection"][()] = 32750
            grid["xCoordinates"][:] = 500_000.0 + np.arange(20) * 5.0
            grid["yCoordinates"][:] = 9_800_000.0 - np.arange(20) * 5.0
            if path == secondary and shift_axis is not None:
                coordinates = grid[f"{shift_axis}Coordinates"]
                coordinates[:] = coordinates[:] + 5.0

    if shift_axis is not None:
        with pytest.raises(ValueError, match=f"{shift_axis}-coordinate grids"):
            form_gslc_interferogram(reference, secondary, tmp_path / "outputs")
        assert not (tmp_path / "outputs").exists()
    else:
        products = form_gslc_interferogram(reference, secondary, tmp_path / "outputs")
        with rasterio.open(products["wrapped_phase"]) as output:
            assert output.crs.to_epsg() == 32750
            assert output.transform.a == 5.0
            assert output.transform.e == -5.0
            np.testing.assert_allclose(output.read(1), 0.0, atol=1e-6)


def test_form_gslc_interferogram_applies_mission_valid_data_mask(
    tmp_path,
):
    reference = tmp_path / "reference.h5"
    secondary = tmp_path / "secondary.h5"
    shape = (32, 32)
    reference_mask = np.ones(shape, dtype=np.uint8)
    secondary_mask = np.ones(shape, dtype=np.uint8)
    secondary_exception_mask = np.zeros(shape, dtype=np.uint8)
    reference_mask[8:16, 8:16] = 0
    secondary_mask[20:24, 20:24] = 255
    secondary_mask[24:28, 24:28] = 2
    secondary_exception_mask[4:8, 4:8] = 2
    phase = np.zeros(shape, dtype=np.float32)
    _write_gslc(
        reference,
        "2026-06-22T22:07:11",
        phase,
        mask=reference_mask,
    )
    _write_gslc(
        secondary,
        "2026-07-04T22:07:10",
        phase,
        mask=secondary_mask,
        exception_mask=secondary_exception_mask,
    )

    products = form_gslc_interferogram(
        reference,
        secondary,
        tmp_path / "outputs",
        block_size=16,
    )
    with rasterio.open(products["wrapped_phase"]) as output:
        data = output.read(1)
        assert np.all(data[8:16, 8:16] == output.nodata)
        assert np.all(data[20:24, 20:24] == output.nodata)
        assert np.all(data[24:28, 24:28] == output.nodata)
        assert np.all(data[4:8, 4:8] == output.nodata)
        assert output.tags()["mission_valid_data_mask"] == (
            "subswath_ids_1_to_product_number_applied"
        )
        assert output.tags()["mission_input_data_exception_mask"] == (
            "nonzero_values_excluded"
        )
        assert output.tags()["reference_number_of_subswaths"] == "1"
        assert output.tags()["secondary_number_of_subswaths"] == "1"
    with rasterio.open(products["coherence"]) as output:
        data = output.read(1)
        assert np.all(data[8:16, 8:16] == output.nodata)
        assert np.all(data[20:24, 20:24] == output.nodata)
        assert np.all(data[24:28, 24:28] == output.nodata)
        assert np.all(data[4:8, 4:8] == output.nodata)


def test_form_gslc_interferogram_requires_mission_mask_by_default(
    tmp_path,
):
    reference = tmp_path / "reference.h5"
    secondary = tmp_path / "secondary.h5"
    phase = np.zeros((16, 16), dtype=np.float32)
    _write_gslc(reference, "2026-06-22T22:07:11", phase)
    _write_gslc(secondary, "2026-07-04T22:07:10", phase)
    with h5py.File(secondary, "r+") as h5:
        del h5["science/LSAR/GSLC/grids/frequencyA/mask"]

    with pytest.raises(KeyError, match="valid-data mask"):
        form_gslc_interferogram(
            reference,
            secondary,
            tmp_path / "required",
        )
    assert not list((tmp_path / "required").glob("*.tif"))

    products = form_gslc_interferogram(
        reference,
        secondary,
        tmp_path / "diagnostic",
        require_quality_mask=False,
    )
    assert products["wrapped_phase"].exists()
