"""Tests for the C-band velocity fusion quarantine gate.

Sentinel-1 C-band peat subsidence was falsified over this AOI (r ~ -0.03
against the independent ALOS-2 L-band benchmark of Tay et al. (2025)).
fuse_velocity_maps must refuse to blend C-band velocity into a product
unless PEATGUARD_ALLOW_CBAND_FUSION is set to exactly "1" for explicit
archival reproduction.
"""

from __future__ import annotations

import numpy as np
import pytest

from peatguard.analysis.fusion import fuse_velocity_maps

NODATA = -9999.0


@pytest.fixture
def tiny_inputs():
    """2x2 arrays covering both-valid, C-only, L-only, and low-coherence cases."""
    c_vel = np.array([[-10.0, -20.0], [NODATA, -40.0]], dtype=np.float32)
    c_coh = np.array([[0.8, 0.6], [0.9, 0.1]], dtype=np.float32)
    l_vel = np.array([[-30.0, NODATA], [-50.0, -60.0]], dtype=np.float32)
    l_coh = np.array([[0.4, 0.2], [0.7, 0.5]], dtype=np.float32)
    return c_vel, c_coh, l_vel, l_coh


def test_fusion_raises_without_override(tiny_inputs, monkeypatch):
    monkeypatch.delenv("PEATGUARD_ALLOW_CBAND_FUSION", raising=False)
    with pytest.raises(RuntimeError, match="quarantined"):
        fuse_velocity_maps(*tiny_inputs)


@pytest.mark.parametrize("value", ["0", "true", "yes", ""])
def test_fusion_raises_for_non_exact_override_values(tiny_inputs, monkeypatch, value):
    monkeypatch.setenv("PEATGUARD_ALLOW_CBAND_FUSION", value)
    with pytest.raises(RuntimeError, match="quarantined"):
        fuse_velocity_maps(*tiny_inputs)


def test_fusion_computes_coherence_weighted_result_when_unlocked(
    tiny_inputs, monkeypatch
):
    monkeypatch.setenv("PEATGUARD_ALLOW_CBAND_FUSION", "1")
    c_vel, c_coh, l_vel, l_coh = tiny_inputs

    fused = fuse_velocity_maps(
        c_vel, c_coh, l_vel, l_coh, nodata=NODATA, min_coherence=0.3
    )

    # Both sensors valid: coherence-weighted average
    expected_both = (0.8 * -10.0 + 0.4 * -30.0) / (0.8 + 0.4)
    assert fused[0, 0] == pytest.approx(expected_both, rel=1e-5)
    # L-band nodata: C-band value passes through
    assert fused[0, 1] == pytest.approx(-20.0)
    # C-band nodata: L-band value passes through
    assert fused[1, 0] == pytest.approx(-50.0)
    # C-band coherence below min_coherence: L-band value passes through
    assert fused[1, 1] == pytest.approx(-60.0)
