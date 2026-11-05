"""Pins for scientific constants and thresholds so silent edits fail loudly.

Authoritative values live in config/default.yaml; the pydantic model
defaults must agree so code paths that construct config models without
the YAML behave identically. Each assertion message names the source of
the pinned value; project-defined parameters without a literature
calibration are labeled as such rather than attributed to a paper.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from peatguard.config import (
    ClassificationConfig,
    FusionConfig,
    MintPyConfig,
    NISARConfig,
    PeatMaskConfig,
    ProcessingConfig,
    RiskScoreConfig,
    Sentinel1Config,
    WaterMaskConfig,
    load_config,
)
from peatguard.timeseries.mintpy_prep import _get_sar_params

CONFIG_PATH = Path(__file__).parent.parent / "config" / "default.yaml"

SPEED_OF_LIGHT_M_S = 299792458.0
# Sentinel-1 C-band center frequency, 5.405 GHz (ESA Sentinel-1 SAR
# User Guide / Sentinel-1 Product Definition).
SENTINEL1_CENTER_FREQUENCY_HZ = 5.405e9


@pytest.fixture
def config():
    return load_config(CONFIG_PATH)


def test_risk_score_weights_pinned(config):
    rs = config.risk_score
    assert rs.proximity_weight == 0.45, (
        "Project proximity weight (config/default.yaml risk_score); "
        "Hooijer et al. (2012) supports the drainage-subsidence-carbon "
        "process link but does not prescribe this weight"
    )
    assert rs.subsidence_weight == 0.55, (
        "Project subsidence weight (config/default.yaml risk_score); "
        "Hooijer et al. (2012) does not prescribe this weight"
    )
    assert rs.proximity_weight + rs.subsidence_weight == pytest.approx(1.0), (
        "Risk components must form a convex combination"
    )


def test_canal_influence_radius_pinned(config):
    assert config.risk_score.max_influence_m == 600.0, (
        "Project-assumed canal influence radius (config/default.yaml "
        "risk_score.max_influence_m); a sensitivity parameter pending "
        "field hydrology calibration, not a Hooijer et al. (2012) value"
    )


def test_severe_velocity_threshold_pinned(config):
    assert config.risk_score.severe_velocity_mm_yr == -40.0, (
        "Project severity saturation threshold (config/default.yaml "
        "risk_score.severe_velocity_mm_yr), within the drained-peatland "
        "subsidence range reported by Hooijer et al. (2012) but not a "
        "literature calibration"
    )


def test_classification_thresholds_pinned(config):
    cls = config.classification
    assert cls.severe_threshold == -50.0, (
        "[H] project classification threshold (config/default.yaml "
        "classification.severe_threshold); no literature calibration"
    )
    assert cls.active_drying_threshold == -20.0, (
        "[H] project classification threshold (config/default.yaml "
        "classification.active_drying_threshold); no literature calibration"
    )
    assert cls.moderate_drying_threshold == -5.0, (
        "[H] project threshold (config/default.yaml "
        "classification.moderate_drying_threshold): below -5 mm/yr "
        "exceeds the ~1 mm/yr natural peat accumulation envelope"
    )
    assert cls.stable_threshold == 5.0, (
        "[H] project threshold (config/default.yaml "
        "classification.stable_threshold): upper bound of the natural "
        "variability envelope"
    )
    assert (
        cls.severe_threshold
        < cls.active_drying_threshold
        < cls.moderate_drying_threshold
        < cls.stable_threshold
    ), "Classification thresholds must be strictly ordered"


def test_sentinel1_wavelength_derives_from_cband_center_frequency(config):
    params = _get_sar_params(config)
    wavelength = float(params["WAVELENGTH"])
    derived = SPEED_OF_LIGHT_M_S / SENTINEL1_CENTER_FREQUENCY_HZ
    assert wavelength == 0.05546576, (
        "Sentinel-1 wavelength pin used for MintPy phase-to-displacement "
        "scaling (mintpy_prep._get_sar_params)"
    )
    assert wavelength == pytest.approx(derived, abs=1e-8), (
        "Wavelength must equal c / 5.405 GHz, the Sentinel-1 C-band "
        "center frequency (ESA Sentinel-1 SAR User Guide), rounded to "
        "8 decimal places"
    )


def test_model_defaults_match_default_yaml(config):
    """Constructing a section model without the YAML must yield the same
    values as config/default.yaml, so a missing section cannot silently
    change the pipeline."""
    assert Sentinel1Config() == config.sentinel1
    assert NISARConfig() == config.nisar
    assert ProcessingConfig() == config.processing
    assert ClassificationConfig() == config.classification
    assert WaterMaskConfig() == config.water_mask
    assert PeatMaskConfig() == config.peat_mask
    assert RiskScoreConfig() == config.risk_score
    assert FusionConfig() == config.fusion
    assert MintPyConfig() == config.mintpy
