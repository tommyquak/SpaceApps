"""Tests for NISAR GUNW correction and readiness gates."""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

import h5py
import numpy as np
import pytest

from peatguard.timeseries.nisar_gunw import (
    ReferenceSelection,
    assess_temporal_readiness,
    prepare_phase_variant,
    select_reference_patch,
)


def _write_stack(
    path: Path,
    phase: np.ndarray,
    *,
    coherence: float = 0.9,
) -> Path:
    pair_count, height, width = phase.shape
    dates = np.array(
        [[b"20260622", b"20260704"], [b"20260704", b"20260716"]][
            :pair_count
        ]
    )
    with h5py.File(path, "w") as h5:
        h5.create_dataset("unwrapPhase", data=phase.astype(np.float32))
        h5.create_dataset(
            "coherence",
            data=np.full(phase.shape, coherence, dtype=np.float32),
        )
        h5.create_dataset(
            "connectComponent",
            data=np.ones(phase.shape, dtype=np.float32),
        )
        h5.create_dataset("date", data=dates)
        h5.create_dataset("bperp", data=np.zeros(pair_count, dtype=np.float32))
        h5.create_dataset(
            "dropIfgram", data=np.ones(pair_count, dtype=bool)
        )
        h5.attrs.update(
            {
                "EPSG": "32750",
                "X_FIRST": "200000",
                "Y_FIRST": "9700000",
                "X_STEP": "80",
                "Y_STEP": "-80",
                "WAVELENGTH": "0.2419632429378531",
            }
        )
    return path


def _write_geometry(path: Path, shape: tuple[int, int]) -> Path:
    with h5py.File(path, "w") as h5:
        h5.create_dataset("waterMask", data=np.ones(shape, dtype=bool))
    return path


def _write_date_network(
    path: Path,
    *,
    acquisition_count: int,
    connection_bandwidth: int,
) -> Path:
    first = datetime(2025, 1, 1)
    dates = [
        (first + timedelta(days=12 * index)).strftime("%Y%m%d")
        for index in range(acquisition_count)
    ]
    pairs = [
        (dates[index], dates[index + lag])
        for index in range(acquisition_count)
        for lag in range(1, connection_bandwidth + 1)
        if index + lag < acquisition_count
    ]
    encoded = np.asarray(pairs, dtype="S8")
    with h5py.File(path, "w") as h5:
        h5.create_dataset("date", data=encoded)
    return path


def test_temporal_readiness_rejects_short_stack_as_velocity(tmp_path: Path) -> None:
    stack = _write_stack(tmp_path / "stack.h5", np.zeros((2, 20, 20)))
    result = assess_temporal_readiness(stack)
    assert result.acquisition_count == 3
    assert result.span_days == 24
    assert result.network_cycle_count == 0
    assert result.bridge_edge_count == 2
    assert not result.redundant_network_ready
    assert not result.tay_network_geometry_ready
    assert result.commissioning_displacement_ready
    assert not result.experimental_trend_ready
    assert not result.annual_velocity_ready
    assert result.release_status == "relative_displacement_commissioning_only"


def test_temporal_readiness_rejects_long_nearest_neighbour_tree(
    tmp_path: Path,
) -> None:
    stack = _write_date_network(
        tmp_path / "nearest.h5",
        acquisition_count=32,
        connection_bandwidth=1,
    )
    result = assess_temporal_readiness(stack)
    assert result.span_days >= 365
    assert result.network_cycle_count == 0
    assert result.bridge_edge_count == 31
    assert not result.redundant_network_ready
    assert not result.experimental_trend_ready
    assert not result.annual_velocity_ready


def test_temporal_readiness_separates_redundancy_from_tay_bandwidth(
    tmp_path: Path,
) -> None:
    stack = _write_date_network(
        tmp_path / "bandwidth_two.h5",
        acquisition_count=32,
        connection_bandwidth=2,
    )
    result = assess_temporal_readiness(stack)
    assert result.redundant_network_ready
    assert result.bridge_edge_count == 0
    assert result.maximum_connection_bandwidth == 2
    assert result.complete_connection_bandwidth == 2
    assert result.annual_velocity_ready
    assert not result.tay_network_geometry_ready


def test_temporal_readiness_rejects_single_long_edge_as_tay_network(
    tmp_path: Path,
) -> None:
    stack = _write_date_network(
        tmp_path / "sparse_long_edge.h5",
        acquisition_count=32,
        connection_bandwidth=2,
    )
    with h5py.File(stack, "r+") as h5:
        dates = h5["date"][:]
        first = dates[0, 0]
        acquisitions = sorted(set(dates.reshape(-1).tolist()))
        del h5["date"]
        h5.create_dataset(
            "date",
            data=np.concatenate(
                [dates, np.asarray([[first, acquisitions[20]]], dtype="S8")],
                axis=0,
            ),
        )

    result = assess_temporal_readiness(stack)
    assert result.redundant_network_ready
    assert result.maximum_connection_bandwidth == 20
    assert result.complete_connection_bandwidth == 2
    assert result.tay_required_pair_fraction < 1.0
    assert result.annual_velocity_ready
    assert not result.tay_network_geometry_ready


def test_temporal_readiness_recognises_tay_bandwidth_geometry(
    tmp_path: Path,
) -> None:
    stack = _write_date_network(
        tmp_path / "bandwidth_twenty.h5",
        acquisition_count=32,
        connection_bandwidth=20,
    )
    result = assess_temporal_readiness(stack)
    assert result.redundant_network_ready
    assert result.maximum_connection_bandwidth == 20
    assert result.complete_connection_bandwidth == 20
    assert result.tay_required_pair_fraction == pytest.approx(1.0)
    assert result.annual_velocity_ready
    assert result.tay_network_geometry_ready
    assert result.release_status == (
        "tay_network_candidate_requires_correction_and_quality_gates"
    )


def test_reference_selection_uses_smooth_high_coherence_interior(
    tmp_path: Path,
) -> None:
    phase = np.zeros((2, 40, 40), dtype=np.float32)
    stack = _write_stack(tmp_path / "stack.h5", phase)
    geometry = _write_geometry(tmp_path / "geometry.h5", phase.shape[1:])
    result = select_reference_patch(
        stack,
        geometry,
        minimum_canal_distance_m=0,
        minimum_peat_distance_m=0,
        minimum_edge_distance_m=400,
    )
    assert 5 <= result.x < 35
    assert 5 <= result.y < 35
    assert result.minimum_pair_coherence == pytest.approx(0.9)
    assert result.reference_status == "relative_reference_not_ground_control"


def test_prepare_phase_variant_subtracts_corrections_and_patch_reference(
    tmp_path: Path,
) -> None:
    rows, columns = np.indices((9, 9))
    raw = np.stack(
        [
            10.0 + rows * 0.1 + columns * 0.2,
            -4.0 + rows * 0.2 - columns * 0.1,
        ]
    ).astype(np.float32)
    correction = np.full_like(raw, 2.0)
    source = _write_stack(tmp_path / "source.h5", raw)
    ionosphere = _write_stack(tmp_path / "ion.h5", correction)
    reference = ReferenceSelection(
        x=4,
        y=4,
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

    output = prepare_phase_variant(
        source,
        tmp_path / "corrected.h5",
        reference,
        correction_stacks={"ionosphere": ionosphere},
    )
    with h5py.File(output, "r") as h5:
        corrected = h5["unwrapPhase"][:]
        assert h5.attrs["PEATGUARD_CORRECTIONS_APPLIED"] == "ionosphere"
        assert h5.attrs["PROCESSOR"] == "nisar"
        for layer in corrected:
            patch = layer[3:6, 3:6]
            assert np.median(patch) == pytest.approx(0.0, abs=1e-6)
    assert np.ptp(corrected[0]) == pytest.approx(np.ptp(raw[0]))
