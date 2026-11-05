"""Tests for reliable MintPy velocity export."""

from pathlib import Path

import h5py
import numpy as np
import pytest
import rasterio

from peatguard.timeseries.velocity import export_velocity


def _write_velocity_inputs(
    directory: Path,
    velocity: np.ndarray,
    velocity_std: np.ndarray,
    incidence: np.ndarray | None,
) -> None:
    directory.mkdir(parents=True)
    with h5py.File(directory / "velocity.h5", "w") as output:
        output.create_dataset("velocity", data=velocity.astype(np.float32))
        output.create_dataset("velocityStd", data=velocity_std.astype(np.float32))
        output.attrs.update(
            {
                "X_FIRST": 114.0,
                "Y_FIRST": -2.0,
                "X_STEP": 0.001,
                "Y_STEP": -0.001,
                "EPSG": 4326,
            }
        )
    with h5py.File(directory / "temporalCoherence.h5", "w") as output:
        output.create_dataset(
            "temporalCoherence",
            data=np.full(velocity.shape, 0.9, dtype=np.float32),
        )
    if incidence is not None:
        inputs = directory / "inputs"
        inputs.mkdir()
        with h5py.File(inputs / "geometryGeo.h5", "w") as output:
            output.create_dataset(
                "incidenceAngle", data=incidence.astype(np.float32)
            )


def test_export_velocity_uses_pixel_incidence_and_masks_bad_values(
    tmp_path: Path,
) -> None:
    mintpy = tmp_path / "mintpy"
    velocity = np.array([[0.001, 0.001], [0.200, -0.001]], dtype=np.float32)
    velocity_std = np.full((2, 2), 0.0001, dtype=np.float32)
    incidence = np.array([[60.0, 0.0], [60.0, 30.0]], dtype=np.float32)
    _write_velocity_inputs(mintpy, velocity, velocity_std, incidence)

    products = export_velocity(mintpy, tmp_path / "products")

    with rasterio.open(products["subsidence_velocity"]) as dataset:
        exported = dataset.read(1)
        nodata = dataset.nodata
    assert exported[0, 0] == pytest.approx(2.0, abs=1e-5)
    assert exported[0, 1] == nodata
    assert exported[1, 0] == nodata
    assert exported[1, 1] == pytest.approx(
        -1.0 / np.cos(np.deg2rad(30.0)), abs=1e-5
    )

    with rasterio.open(products["velocity_uncertainty"]) as dataset:
        uncertainty = dataset.read(1)
        nodata = dataset.nodata
    assert uncertainty[0, 0] == pytest.approx(0.2, abs=1e-5)
    assert uncertainty[0, 1] == nodata
    assert uncertainty[1, 0] == nodata
    assert uncertainty[1, 1] == pytest.approx(
        0.1 / np.cos(np.deg2rad(30.0)), abs=1e-5
    )


def test_export_velocity_uses_constant_fallback_without_geometry(
    tmp_path: Path,
) -> None:
    mintpy = tmp_path / "mintpy"
    velocity = np.full((2, 2), 0.001, dtype=np.float32)
    velocity_std = np.full((2, 2), 0.0001, dtype=np.float32)
    _write_velocity_inputs(mintpy, velocity, velocity_std, incidence=None)

    products = export_velocity(mintpy, tmp_path / "products")

    with rasterio.open(products["subsidence_velocity"]) as dataset:
        exported = dataset.read(1)
    expected = 1.0 / np.cos(np.deg2rad(37.0))
    np.testing.assert_allclose(exported, expected, atol=1e-5)
