"""Behavioral regression guards for PeatGuard's settled scientific release policy.

The analysis stage must preserve the observation-derived InSAR datum
(a published peat subsidence rate is not a geodetic reference, so no
shift toward a literature mean is permitted) and must not emit a
universal carbon-emission raster (Hooijer et al. (2012) requires
process partitioning, peat properties, and field evidence before any
subsidence-to-carbon conversion).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pytest
import rasterio

from peatguard.config import load_config
from peatguard.export.cog import read_raster, write_cog
from peatguard.pipeline.orchestrator import run_analysis_stage

CONFIG_PATH = Path(__file__).parent.parent / "config" / "default.yaml"
NODATA = -9999.0


@pytest.fixture
def analysis_config(tmp_path):
    """Default config trimmed to the offline analysis-stage core."""
    config = load_config(CONFIG_PATH)
    config.storage.output_dir = tmp_path / "output"
    config.storage.output_dir.mkdir()
    # Stages that download external data are out of scope: the policy
    # under test concerns the velocity datum and derived rasters only.
    config.peat_mask.enabled = False
    config.forest_change.enabled = False
    config.validation.enabled = False
    config.fusion.enabled = False
    return config


@pytest.fixture
def velocity_file(tmp_path):
    """Small synthetic velocity GeoTIFF with a known nonzero mean."""
    rng = np.random.default_rng(7)
    velocity = rng.normal(-12.0, 6.0, size=(24, 24)).astype(np.float32)
    velocity[:2, :2] = NODATA
    path = tmp_path / "subsidence_velocity.tif"
    write_cog(
        data=velocity,
        output_path=path,
        crs=32750,
        transform=rasterio.transform.from_origin(300000.0, 9720000.0, 10.0, 10.0),
        nodata=NODATA,
        band_names=["velocity_mm_yr"],
        dtype="float32",
    )
    return path


def test_analysis_stage_preserves_native_velocity_datum(
    analysis_config, velocity_file, caplog
):
    """The velocity field keeps its native InSAR reference, unshifted."""
    before, _ = read_raster(velocity_file)
    input_velocity = before[0]
    input_mean = input_velocity[input_velocity != NODATA].mean()

    with caplog.at_level(logging.INFO, logger="peatguard.pipeline.orchestrator"):
        products = run_analysis_stage(analysis_config, velocity_file)

    after, _ = read_raster(velocity_file)
    output_velocity = after[0]
    output_mean = output_velocity[output_velocity != NODATA].mean()

    assert np.array_equal(input_velocity, output_velocity), (
        "analysis stage must not rewrite or datum-shift the velocity raster"
    )
    assert output_mean == pytest.approx(input_mean), (
        "valid-pixel mean must be unchanged: a literature target mean is "
        "not a geodetic reference"
    )
    assert "Retaining native InSAR reference datum" in caplog.text
    assert not (analysis_config.storage.output_dir / "calibrated_velocity.tif").exists()
    assert products["subsidence_class"].exists()


def test_analysis_stage_emits_no_universal_carbon_conversion(
    analysis_config, velocity_file
):
    """Carbon accounting requires process partitioning and field evidence."""
    products = run_analysis_stage(analysis_config, velocity_file)

    output_dir = analysis_config.storage.output_dir
    assert not (output_dir / "carbon_loss.tif").exists()
    assert not any("carbon" in name for name in products)
