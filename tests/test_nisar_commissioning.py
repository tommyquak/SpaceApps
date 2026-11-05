"""Tests for NISAR commissioning algebra and uncertainty gates."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest

from peatguard.timeseries.nisar_commissioning import (
    build_commissioning_quality_mask,
    phase_to_los_displacement,
    propagate_ionosphere_uncertainty,
    sha256_file,
    validate_tree_inversion,
    verify_file_manifest,
)
from peatguard.timeseries.nisar_gunw import ReferenceSelection


def _write_tree_stack(path: Path) -> tuple[Path, np.ndarray, float]:
    phase = np.stack(
        [
            np.full((5, 6), 1.5, dtype=np.float32),
            np.full((5, 6), -0.4, dtype=np.float32),
        ]
    )
    wavelength = 0.2419632429378531
    with h5py.File(path, "w") as h5:
        h5.create_dataset("unwrapPhase", data=phase)
        h5.create_dataset(
            "connectComponent", data=np.ones(phase.shape, dtype=np.uint16)
        )
        h5.create_dataset(
            "date",
            data=np.array(
                [
                    [b"20260622", b"20260704"],
                    [b"20260704", b"20260716"],
                ]
            ),
        )
        h5.attrs["WAVELENGTH"] = wavelength
    return path, phase, wavelength


def _write_timeseries(
    path: Path, phase: np.ndarray, wavelength: float
) -> Path:
    pair_displacement = phase_to_los_displacement(phase, wavelength)
    timeseries = np.stack(
        [
            np.zeros(phase.shape[1:], dtype=np.float32),
            pair_displacement[0],
            pair_displacement[0] + pair_displacement[1],
        ]
    )
    with h5py.File(path, "w") as h5:
        h5.create_dataset("timeseries", data=timeseries)
        h5.create_dataset(
            "date",
            data=np.array([b"20260622", b"20260704", b"20260716"]),
        )
    return path


def _reference() -> ReferenceSelection:
    return ReferenceSelection(
        x=3,
        y=2,
        longitude=114.5,
        latitude=-2.7,
        patch_radius_pixels=1,
        pixel_spacing_m=(80.0, 80.0),
        candidate_count=1,
        minimum_pair_coherence=0.9,
        pair_coherence=(0.9, 0.9),
        maximum_local_phase_std_rad=0.1,
        distance_to_mapped_canal_m=2_000.0,
        distance_to_mapped_peat_m=2_000.0,
        dominant_connected_components=(1, 1),
        reference_status="relative_reference_not_ground_control",
    )


def test_phase_to_los_uses_mintpy_negative_sign() -> None:
    output = phase_to_los_displacement(np.array([4.0 * np.pi]), 0.24)
    assert output[0] == pytest.approx(-0.24)


def test_exact_tree_inversion_validation(tmp_path: Path) -> None:
    stack, phase, wavelength = _write_tree_stack(tmp_path / "stack.h5")
    timeseries = _write_timeseries(
        tmp_path / "timeseries.h5", phase, wavelength
    )
    result = validate_tree_inversion("test", stack, timeseries)
    assert result.passed
    assert result.pixel_count == 30
    assert result.maximum_absolute_error_m < 1e-7


def test_exact_tree_inversion_detects_wrong_sign(tmp_path: Path) -> None:
    stack, phase, wavelength = _write_tree_stack(tmp_path / "stack.h5")
    timeseries = _write_timeseries(
        tmp_path / "timeseries.h5", phase, wavelength
    )
    with h5py.File(timeseries, "r+") as h5:
        h5["timeseries"][1:] *= -1
    result = validate_tree_inversion("wrong", stack, timeseries)
    assert not result.passed
    assert result.maximum_absolute_error_m > 0.01


def test_quality_mask_requires_all_pairs() -> None:
    coherence = np.full((2, 4, 4), 0.8)
    connected = np.ones((2, 4, 4), dtype=np.uint16)
    land = np.ones((4, 4), dtype=bool)
    incidence = np.full((4, 4), 35.0)
    coherence[1, 1, 1] = 0.4
    connected[0, 2, 2] = 0
    output = build_commissioning_quality_mask(
        coherence,
        connected,
        land,
        incidence,
        minimum_pair_coherence=0.5,
    )
    assert output.sum() == 14
    assert not output[1, 1]
    assert not output[2, 2]


def test_uncertainty_propagation_includes_reference_term() -> None:
    uncertainty = np.stack(
        [
            np.full((5, 6), 1.0),
            np.full((5, 6), 2.0),
        ]
    )
    incidence = np.full((5, 6), 30.0)
    lower, conservative = propagate_ionosphere_uncertainty(
        uncertainty,
        0.24,
        incidence,
        reference=_reference(),
    )
    expected_lower_phase = np.sqrt(1.0**2 + 2.0**2)
    expected_conservative_phase = np.sqrt(
        1.0**2 + 1.0**2 + 2.0**2 + 2.0**2
    )
    scale = 0.24 / (4.0 * np.pi) * 1_000.0 / np.cos(np.deg2rad(30.0))
    assert lower[0, 0] == pytest.approx(expected_lower_phase * scale)
    assert conservative[0, 0] == pytest.approx(
        expected_conservative_phase * scale
    )
    assert np.all(conservative > lower)


def test_file_manifest_verification_supports_identical_relocation(
    tmp_path: Path,
) -> None:
    original = tmp_path / "original.dat"
    relocated = tmp_path / "relocated.dat"
    original.write_bytes(b"peatguard-provenance")
    relocated.write_bytes(original.read_bytes())
    manifest = {
        str(original): {
            "bytes": original.stat().st_size,
            "sha256": sha256_file(original),
        }
    }
    result = verify_file_manifest(
        manifest,
        replacements={str(original): relocated},
    )
    assert result[str(original)]["verified"] is True
    assert result[str(original)]["actual_path"] == str(relocated.resolve())


def test_file_manifest_verification_rejects_changed_content(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.dat"
    source.write_bytes(b"original")
    manifest = {
        str(source): {
            "bytes": source.stat().st_size,
            "sha256": sha256_file(source),
        }
    }
    source.write_bytes(b"modified")
    with pytest.raises(ValueError, match="provenance verification failed"):
        verify_file_manifest(manifest)
