from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest
import rasterio

from peatguard.insar.nisar_gunw import (
    export_gunw_layers,
    inspect_gunw,
    inspect_gunw_grid,
)


def _scalar(group: h5py.Group, name: str, value: object) -> None:
    group.create_dataset(name, data=value)


def _make_gunw(path: Path) -> Path:
    with h5py.File(path, "w") as h5:
        identification = h5.require_group("science/LSAR/identification")
        values = {
            "granuleId": b"NISAR_L2_PR_GUNW_TEST_P05023",
            "compositeReleaseId": b"P05023",
            "referenceZeroDopplerStartTime": b"2026-06-22T22:07:10.000000000",
            "secondaryZeroDopplerStartTime": b"2026-07-04T22:07:10.000000000",
            "trackNumber": 125,
            "frameNumber": 176,
            "orbitPassDirection": b"Ascending",
            "productVersion": b"1.0.11",
            "productSpecificationVersion": b"1.5.0",
        }
        for name, value in values.items():
            _scalar(identification, name, value)

        frequency = h5.require_group("science/LSAR/GUNW/grids/frequencyA")
        frequency.create_dataset(
            "centerFrequency", data=np.float64(1_239_000_000.0)
        )
        wrapped = frequency.require_group("wrappedInterferogram")
        wrapped_hh = wrapped.require_group("HH")
        wrapped_hh.create_dataset(
            "xCoordinates", data=np.array([114.00, 114.01, 114.02, 114.03])
        )
        wrapped_hh.create_dataset(
            "yCoordinates", data=np.array([-2.00, -2.01, -2.02])
        )
        wrapped_hh.create_dataset("projection", data=4326)
        complex_data = np.exp(
            1j
            * np.array(
                [
                    [0.0, 0.2, 0.4, 0.6],
                    [0.8, 1.0, 1.2, 1.4],
                    [1.6, 1.8, 2.0, 2.2],
                ],
                dtype=np.float32,
            )
        ).astype(np.complex64)
        wrapped_hh.create_dataset("wrappedInterferogram", data=complex_data)
        wrapped_hh.create_dataset(
            "coherenceMagnitude",
            data=np.full((3, 4), 0.75, dtype=np.float32),
        )

        unwrapped = frequency.require_group("unwrappedInterferogram")
        unwrapped_hh = unwrapped.require_group("HH")
        unwrapped_hh.create_dataset(
            "xCoordinates", data=np.array([114.00, 114.02])
        )
        unwrapped_hh.create_dataset(
            "yCoordinates", data=np.array([-2.00, -2.02])
        )
        unwrapped_hh.create_dataset("projection", data=4326)
        unwrapped_hh.create_dataset(
            "unwrappedPhase",
            data=np.array([[1.0, 2.0], [3.0, np.nan]], dtype=np.float32),
        )
        unwrapped_hh.create_dataset(
            "coherenceMagnitude",
            data=np.full((2, 2), 0.6, dtype=np.float32),
        )
        connected = unwrapped_hh.create_dataset(
            "connectedComponents",
            data=np.array([[1, 1], [2, 65535]], dtype=np.uint16),
        )
        connected.attrs["_FillValue"] = np.uint16(65535)
        unwrapped_hh.create_dataset(
            "ionospherePhaseScreen",
            data=np.full((2, 2), 0.3, dtype=np.float32),
        )
        unwrapped_hh.create_dataset(
            "ionospherePhaseScreenUncertainty",
            data=np.full((2, 2), 0.1, dtype=np.float32),
        )
    return path


def test_inspect_gunw(tmp_path: Path) -> None:
    path = _make_gunw(tmp_path / "test.h5")
    info = inspect_gunw(path)
    assert info.reference_date == "20260622"
    assert info.secondary_date == "20260704"
    assert info.track_number == 125
    assert info.maturity == "PROVISIONAL"
    assert info.center_frequency_hz == pytest.approx(1_239_000_000.0)
    assert info.wavelength_m == pytest.approx(0.2419632429378531)

    grid = inspect_gunw_grid(path, "wrappedInterferogram")
    assert grid.epsg == 4326
    assert grid.shape == (3, 4)
    assert grid.x_spacing == pytest.approx(0.01)
    assert grid.y_spacing == pytest.approx(-0.01)


def test_export_gunw_layers(tmp_path: Path) -> None:
    path = _make_gunw(tmp_path / "test.h5")
    outputs = export_gunw_layers(
        path,
        tmp_path / "outputs",
        bbox_wgs84=None,
        layers=[
            "wrapped_phase_20m",
            "coherence_20m",
            "unwrapped_phase_80m",
            "connected_components_80m",
        ],
        block_size=2,
    )
    assert set(outputs) == {
        "wrapped_phase_20m",
        "coherence_20m",
        "unwrapped_phase_80m",
        "connected_components_80m",
    }

    with rasterio.open(outputs["wrapped_phase_20m"]) as dataset:
        phase = dataset.read(1)
        assert phase.shape == (3, 4)
        assert phase[1, 1] == pytest.approx(1.0)
        assert dataset.tags()["product_maturity"] == "PROVISIONAL"
        assert dataset.tags()["corrections_applied"] == "none"
        assert float(dataset.tags()["wavelength_m"]) == pytest.approx(
            0.2419632429378531
        )

    with rasterio.open(outputs["unwrapped_phase_80m"]) as dataset:
        data = dataset.read(1)
        assert data[0, 1] == pytest.approx(2.0)
        assert data[1, 1] == dataset.nodata

    with rasterio.open(outputs["connected_components_80m"]) as dataset:
        assert dataset.dtypes[0] == "uint16"
        assert dataset.nodata == 65535


def test_export_rejects_unknown_layer(tmp_path: Path) -> None:
    path = _make_gunw(tmp_path / "test.h5")
    with pytest.raises(ValueError, match="Unknown GUNW layers"):
        export_gunw_layers(path, tmp_path / "outputs", layers=["velocity"])
