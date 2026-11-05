"""Shared test fixtures and utilities."""

import os
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture(autouse=True)
def _clear_peatguard_env(monkeypatch):
    """Remove PEATGUARD_ environment variables before every test.

    Config values under test must come from YAML files and model
    defaults, not from the developer's shell. Tests that exercise an
    override set it explicitly with monkeypatch.setenv.
    """
    for env_key in list(os.environ):
        if env_key.startswith("PEATGUARD_"):
            monkeypatch.delenv(env_key, raising=False)


@pytest.fixture
def tmp_output(tmp_path):
    """Temporary output directory."""
    output = tmp_path / "output"
    output.mkdir()
    return output


@pytest.fixture
def sample_raster_100x100():
    """Generate a 100x100 float32 test raster with known values."""
    np.random.seed(42)
    return np.random.randn(100, 100).astype(np.float32)


@pytest.fixture
def sample_complex_ifg():
    """Generate a synthetic wrapped interferogram."""
    np.random.seed(42)
    rows, cols = 100, 100
    # Create a linear phase ramp (simulating deformation)
    row_coords = np.arange(rows, dtype=np.float64)
    col_coords = np.arange(cols, dtype=np.float64)
    rr, cc = np.meshgrid(row_coords, col_coords, indexing="ij")
    phase = 0.05 * rr + 0.03 * cc + np.random.randn(rows, cols) * 0.3
    return np.exp(1j * phase).astype(np.complex64)


@pytest.fixture
def sample_coherence():
    """Generate a synthetic coherence map (0-1)."""
    np.random.seed(42)
    return np.clip(np.random.beta(5, 2, size=(100, 100)), 0, 1).astype(np.float32)


@pytest.fixture
def sample_velocity():
    """Generate a synthetic subsidence velocity map (mm/yr)."""
    np.random.seed(42)
    rows, cols = 100, 100
    rr, cc = np.meshgrid(np.arange(rows), np.arange(cols), indexing="ij")
    # Gradient from -60 to +10 mm/yr
    velocity = -60 + 70 * cc / cols + np.random.randn(rows, cols) * 5
    return velocity.astype(np.float32)


@pytest.fixture
def sample_vv_backscatter():
    """Generate synthetic VV backscatter with canal-like features."""
    np.random.seed(42)
    data = np.random.uniform(0.01, 0.5, size=(100, 100)).astype(np.float32)
    # Add dark linear canal features
    data[45:55, :] *= 0.05  # Horizontal canal
    data[:, 20:25] *= 0.05  # Vertical canal
    return data


@pytest.fixture
def default_config():
    """Load the default pipeline configuration."""
    from peatguard.config import load_config

    config_path = Path(__file__).parent.parent / "config" / "default.yaml"
    return load_config(config_path)
