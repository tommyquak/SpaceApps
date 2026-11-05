from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin
from rasterio.windows import Window

from peatguard.insar.nisar_fine_product import (
    COHERENCE_NODATA,
    COHERENCE_SCALE,
    build_minimum_coherence_cog,
    quantize_coherence,
)
from peatguard.insar.nisar_validation import (
    _reproject_fine_block,
    audit_three_date_phase_closure,
    compare_fine_to_mission,
    wrap_phase,
)


def _write_raster(
    path: Path,
    data: np.ndarray,
    pixel_size: float,
    *,
    reference_date: str | None = "2026-06-22T22:07:11",
    secondary_date: str | None = "2026-07-04T22:07:10",
) -> Path:
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=data.shape[1],
        height=data.shape[0],
        count=1,
        dtype="float32",
        crs="EPSG:32650",
        transform=from_origin(500000.0, 9900000.0, pixel_size, pixel_size),
        nodata=-9999.0,
        tiled=True,
        blockxsize=16,
        blockysize=16,
    ) as dataset:
        dataset.write(data.astype(np.float32), 1)
        if reference_date and secondary_date:
            dataset.update_tags(
                reference_date=reference_date,
                secondary_date=secondary_date,
            )
    return path


def test_wrap_phase_interval() -> None:
    values = np.array([-4 * np.pi, -np.pi, 0, np.pi, 4 * np.pi])
    wrapped = wrap_phase(values)
    assert np.all(wrapped >= -np.pi)
    assert np.all(wrapped < np.pi)
    assert wrapped[2] == 0


def test_compare_selects_opposite_phase_sign(tmp_path: Path) -> None:
    coarse = np.array(
        [
            [0.2, 0.5, 0.8, 1.1],
            [0.3, 0.6, 0.9, 1.2],
            [0.4, 0.7, 1.0, 1.3],
            [0.5, 0.8, 1.1, 1.4],
        ],
        dtype=np.float32,
    )
    fine_phase = np.repeat(np.repeat(coarse, 4, axis=0), 4, axis=1)
    reference_offset = 0.35
    mission_phase = wrap_phase(-coarse + reference_offset).astype(np.float32)
    fine_coherence = np.full(fine_phase.shape, 0.8, dtype=np.float32)
    mission_coherence = np.full(coarse.shape, 0.75, dtype=np.float32)

    report = compare_fine_to_mission(
        _write_raster(tmp_path / "fine_phase.tif", fine_phase, 5.0),
        _write_raster(tmp_path / "fine_coherence.tif", fine_coherence, 5.0),
        _write_raster(tmp_path / "mission_phase.tif", mission_phase, 20.0),
        _write_raster(
            tmp_path / "mission_coherence.tif", mission_coherence, 20.0
        ),
        output_json=tmp_path / "report.json",
    )
    assert report["phase"]["selected_convention"] == "opposite_sign"
    selected = report["phase"]["candidates"]["opposite_sign"]
    assert selected["circular_mae_rad"] < 1e-5
    assert report["coherence"]["bias_fine_minus_mission"] > 0
    assert not report["release_gates"]["supporting_coherence"]["passed"]
    assert not report["release_gates"]["displacement_phase"]["passed"]
    assert (tmp_path / "report.json").exists()
    assert report["aggregation"]["mode"] == "unit_phasor"


def _cross_comparison_inputs(tmp_path: Path, cross_data: np.ndarray) -> tuple[Path, ...]:
    fine_phase = np.angle(cross_data).astype(np.float32)
    fine_phase[~np.isfinite(cross_data) | (np.abs(cross_data) == 0)] = -9999.0
    shape = cross_data.shape
    mission_shape = (shape[0] // 2, shape[1] // 2)
    aggregated = cross_data.reshape(mission_shape[0], 2, mission_shape[1], 2).mean(axis=(1, 3))
    paths = (
        _write_raster(tmp_path / "fine_phase.tif", fine_phase, 5.0),
        _write_raster(tmp_path / "fine_coherence.tif", np.full(shape, 0.8), 5.0),
        _write_raster(tmp_path / "mission_phase.tif", np.angle(aggregated), 10.0),
        _write_raster(tmp_path / "mission_coherence.tif", np.full(mission_shape, 0.8), 10.0),
        tmp_path / "cross.tif",
    )
    common_tags = {
        "reference_date": "2026-06-22T22:07:11",
        "secondary_date": "2026-07-04T22:07:10",
        "reference_granule": "TEST_GSLC_A",
        "secondary_granule": "TEST_GSLC_B",
        "frequency": "frequencyA",
        "polarization": "HH",
        "phase_convention": "reference_times_conjugate_secondary",
    }
    for path in paths[:4]:
        with rasterio.open(path, "r+") as dataset:
            dataset.update_tags(**common_tags)
    with rasterio.open(paths[0]) as phase:
        profile = {**phase.profile, "count": 2, "nodata": np.nan}
    with rasterio.open(paths[4], "w", **profile) as dataset:
        dataset.write(cross_data.real.astype(np.float32), 1)
        dataset.write(cross_data.imag.astype(np.float32), 2)
        dataset.set_band_description(1, "cross_product_real")
        dataset.set_band_description(2, "cross_product_imaginary")
        dataset.update_tags(
            **common_tags,
            complex_representation="real_imaginary_bands",
            amplitude_source="unscaled_reference_times_conjugate_secondary",
            spatial_aggregation="none_native_geocoded_gslc_grid",
        )
    return paths


def test_complex_weighted_aggregation_preserves_nine_plus_i(tmp_path: Path) -> None:
    cross_data = np.array([[9, 1j], [9, 1j]], dtype=np.complex64)
    paths = _cross_comparison_inputs(tmp_path, cross_data)
    with rasterio.open(paths[0]) as phase, rasterio.open(paths[1]) as coherence, rasterio.open(
        paths[2]
    ) as mission, rasterio.open(paths[4]) as cross:
        unit, _, unit_strength = _reproject_fine_block(
            phase, coherence, mission, Window(0, 0, 1, 1)
        )
        weighted, _, weighted_strength = _reproject_fine_block(
            phase, coherence, mission, Window(0, 0, 1, 1), cross
        )
    assert unit[0, 0] == pytest.approx(np.pi / 4)
    assert weighted[0, 0] == pytest.approx(np.arctan2(1, 9))
    assert unit_strength[0, 0] == pytest.approx(np.sqrt(2) / 2)
    assert weighted_strength[0, 0] == pytest.approx(np.sqrt(82) / 10)


def test_complex_weighted_zero_and_nodata_do_not_create_observations(tmp_path: Path) -> None:
    cross_data = np.array([[0, complex(np.nan, np.nan)], [9, 1j]], dtype=np.complex64)
    paths = _cross_comparison_inputs(tmp_path, cross_data)
    with rasterio.open(paths[0]) as phase, rasterio.open(paths[1]) as coherence, rasterio.open(
        paths[2]
    ) as mission, rasterio.open(paths[4]) as cross:
        weighted, _, strength = _reproject_fine_block(
            phase, coherence, mission, Window(0, 0, 1, 1), cross
        )
    assert weighted[0, 0] == pytest.approx(np.arctan2(1, 9))
    assert strength[0, 0] == pytest.approx(np.sqrt(82) / 10)


@pytest.mark.parametrize("first", [-9999 + 0j, 1 + 0j])
def test_signed_cross_components_and_exact_cancellation(tmp_path: Path, first: complex) -> None:
    second = 1j if first.real == -9999 else -1 + 0j
    paths = _cross_comparison_inputs(
        tmp_path, np.array([[first, second], [first, second]], dtype=np.complex64)
    )
    with rasterio.open(paths[0]) as phase, rasterio.open(paths[1]) as coherence, rasterio.open(
        paths[2]
    ) as mission, rasterio.open(paths[4]) as cross:
        weighted, _, strength = _reproject_fine_block(
            phase, coherence, mission, Window(0, 0, 1, 1), cross
        )
    if first.real == -9999:
        # A signed component is not the legacy phase NoData sentinel.
        assert weighted[0, 0] == pytest.approx(np.arctan2(1, -9999))
        assert strength[0, 0] > 0.99
    else:
        # The phase of a zero complex mean is undefined, not zero radians.
        assert np.isnan(weighted[0, 0])
        assert np.isnan(strength[0, 0])


@pytest.mark.parametrize("sign", [1, -1])
def test_complex_weighted_compare_preserves_sign_gates(tmp_path: Path, sign: int) -> None:
    data = np.array(
        [[9, 1j, 1, 1j], [9, 1j, 1, 1j], [1j, -1, 1j, -9], [1j, -1, 1j, -9]],
        dtype=np.complex64,
    )
    paths = _cross_comparison_inputs(tmp_path, data)
    with rasterio.open(paths[2], "r+") as mission:
        mission.write(wrap_phase(sign * mission.read(1) + 0.35).astype(np.float32), 1)
    report = compare_fine_to_mission(
        *paths[:4], aggregation="complex_weighted", fine_cross_product_path=paths[4]
    )
    convention = "same_sign" if sign == 1 else "opposite_sign"
    assert report["phase"]["selected_convention"] == convention
    assert report["phase"]["candidates"][convention]["circular_mae_rad"] < 1e-6
    checks = report["release_gates"]["displacement_phase"]["checks"]
    assert checks["expected_same_sign_selected"] == (sign == 1)
    assert not report["release_gates"]["displacement_phase"]["passed"]
    assert report["aggregation"]["mode"] == "complex_weighted"
    assert not report["aggregation"]["mission_operator_equivalence"]
    assert report["aggregation"]["cross_product_provenance"]["reference_granule"] == "TEST_GSLC_A"
    legacy = compare_fine_to_mission(*paths[:4])
    assert legacy["phase"]["candidates"][convention]["circular_mae_rad"] > 0.1


@pytest.mark.parametrize("fault", ["grid", "date", "amplitude_source", "signal_sign", "bands"])
def test_complex_weighted_rejects_incompatible_cross_product(tmp_path: Path, fault: str) -> None:
    paths = _cross_comparison_inputs(tmp_path, np.full((4, 4), 9 + 1j, dtype=np.complex64))
    with rasterio.open(paths[4], "r+") as cross:
        if fault == "grid":
            cross.transform = from_origin(500005.0, 9900000.0, 5.0, 5.0)
        elif fault == "date":
            cross.update_tags(reference_date="2026-01-01")
        elif fault == "amplitude_source":
            cross.update_tags(amplitude_source="reconstructed_from_coherence")
        elif fault == "signal_sign":
            cross.write(-cross.read(2), 2)
        else:
            cross.set_band_description(2, "unidentified")
    with pytest.raises(ValueError):
        compare_fine_to_mission(
            *paths[:4], aggregation="complex_weighted", fine_cross_product_path=paths[4]
        )


def test_complex_weighted_rejects_wrong_mission_dates(tmp_path: Path) -> None:
    paths = _cross_comparison_inputs(tmp_path, np.full((4, 4), 9 + 1j, dtype=np.complex64))
    for path in paths[2:4]:
        with rasterio.open(path, "r+") as dataset:
            dataset.update_tags(reference_date="2026-07-16T22:07:11")
    with pytest.raises(ValueError, match="ordered acquisition dates"):
        compare_fine_to_mission(
            *paths[:4], aggregation="complex_weighted", fine_cross_product_path=paths[4]
        )


@pytest.mark.parametrize("fault", ["missing", "mismatch", "wrong_pair", "reversed"])
def test_default_mode_requires_valid_matching_pair_dates(tmp_path: Path, fault: str) -> None:
    paths = _cross_comparison_inputs(tmp_path, np.full((4, 4), 9 + 1j, dtype=np.complex64))
    if fault in {"missing", "mismatch"}:
        with rasterio.open(paths[1], "r+") as dataset:
            dataset.update_tags(reference_date="" if fault == "missing" else "2026-01-01")
    else:
        for path in paths[2:4] if fault == "wrong_pair" else paths[:4]:
            with rasterio.open(path, "r+") as dataset:
                dataset.update_tags(reference_date="2026-07-16T22:07:11")
    with pytest.raises(ValueError, match="acquisition-date|ordered acquisition dates"):
        compare_fine_to_mission(*paths[:4])


@pytest.mark.parametrize("index", [1, 3])
def test_comparison_rejects_out_of_range_coherence(tmp_path: Path, index: int) -> None:
    paths = _cross_comparison_inputs(tmp_path, np.full((4, 4), 9 + 1j, dtype=np.complex64))
    with rasterio.open(paths[index], "r+") as dataset:
        dataset.write(np.full(dataset.shape, 1.1, dtype=np.float32), 1)
    with pytest.raises(RuntimeError, match="No overlapping"):
        compare_fine_to_mission(*paths[:4])


def test_comparison_uses_same_coherence_validity_in_both_passes(tmp_path: Path) -> None:
    paths = _cross_comparison_inputs(tmp_path, np.full((4, 4), 9 + 1j, dtype=np.complex64))
    with rasterio.open(paths[3], "r+") as dataset:
        values = dataset.read(1)
        values[0, 0] = 1.1
        dataset.write(values, 1)
    report = compare_fine_to_mission(*paths[:4])
    assert report["phase"]["candidates"]["same_sign"]["valid_pixels"] == 3


@pytest.mark.parametrize("cross_value", [0j, complex(np.nan, np.nan)])
def test_complex_weighted_all_invalid_is_not_a_zero_phase(
    tmp_path: Path, cross_value: complex
) -> None:
    paths = _cross_comparison_inputs(tmp_path, np.full((4, 4), cross_value, dtype=np.complex64))
    with pytest.raises(RuntimeError, match="No overlapping"):
        compare_fine_to_mission(
            *paths[:4], aggregation="complex_weighted", fine_cross_product_path=paths[4]
        )


def test_complex_weighted_requires_explicit_original_complex_input(tmp_path: Path) -> None:
    paths = _cross_comparison_inputs(tmp_path, np.full((4, 4), 9 + 1j, dtype=np.complex64))
    with pytest.raises(ValueError, match="required only"):
        compare_fine_to_mission(*paths[:4], aggregation="complex_weighted")
    with pytest.raises(ValueError, match="required only"):
        compare_fine_to_mission(*paths[:4], fine_cross_product_path=paths[4])


def test_quantize_coherence_reserves_255_for_nodata() -> None:
    values = np.array([[0.0, 0.5, 1.0, np.nan]], dtype=np.float32)
    valid = np.array([[True, True, True, False]])
    encoded = quantize_coherence(values, valid)
    assert encoded.tolist() == [[0, 127, 254, COHERENCE_NODATA]]


def test_build_minimum_coherence_cog_is_scaled_and_ordered(
    tmp_path: Path,
) -> None:
    first = np.linspace(0.0, 1.0, 32 * 32, dtype=np.float32).reshape(32, 32)
    second = np.flipud(first).copy()
    second[0, 0] = -9999.0
    first_path = _write_raster(tmp_path / "first.tif", first, 5.0)
    second_path = _write_raster(tmp_path / "second.tif", second, 5.0)
    output = tmp_path / "minimum.tif"
    report = build_minimum_coherence_cog(
        [first_path, second_path],
        tmp_path / "working.tif",
        output,
    )

    with rasterio.open(output) as dataset:
        encoded = dataset.read(1)
        assert dataset.tags(ns="IMAGE_STRUCTURE")["LAYOUT"] == "COG"
        assert dataset.scales[0] == COHERENCE_SCALE
        assert dataset.nodata == COHERENCE_NODATA
        assert encoded[0, 0] == COHERENCE_NODATA
        decoded = encoded.astype(np.float64) * dataset.scales[0]
        expected = np.minimum(first, second)
        valid = expected != -9999.0
        assert np.max(np.abs(decoded[valid] - expected[valid])) <= (
            COHERENCE_SCALE / 2.0 + 1e-7
        )
    assert report["encoding"]["layout"] == "COG"
    assert report["coverage"]["valid_pixels"] == first.size - 1


def test_three_date_phase_closure_passes_exact_triangle(tmp_path: Path) -> None:
    rng = np.random.default_rng(42)
    phase_a = rng.uniform(-np.pi, np.pi, size=(32, 48))
    phase_b = rng.uniform(-np.pi, np.pi, size=(32, 48))
    phase_c = rng.uniform(-np.pi, np.pi, size=(32, 48))
    coherence = np.full((32, 48), 0.8, dtype=np.float32)
    dates = ("2026-06-22T22:07:11", "2026-07-04T22:07:10", "2026-07-16T22:07:10")

    def pair_rasters(label: str, first: np.ndarray, second: np.ndarray, i: int, j: int):
        phase_path = _write_raster(
            tmp_path / f"phase_{label}.tif",
            wrap_phase(first - second),
            5.0,
            reference_date=dates[i],
            secondary_date=dates[j],
        )
        coherence_path = _write_raster(
            tmp_path / f"coherence_{label}.tif",
            coherence,
            5.0,
            reference_date=dates[i],
            secondary_date=dates[j],
        )
        return phase_path, coherence_path

    phase_ab, coherence_ab = pair_rasters("ab", phase_a, phase_b, 0, 1)
    phase_bc, coherence_bc = pair_rasters("bc", phase_b, phase_c, 1, 2)
    phase_ac, coherence_ac = pair_rasters("ac", phase_a, phase_c, 0, 2)
    report = audit_three_date_phase_closure(
        phase_ab,
        coherence_ab,
        phase_bc,
        coherence_bc,
        phase_ac,
        coherence_ac,
        output_json=tmp_path / "closure.json",
    )

    metrics = report["metrics"]["coherence_masked"]
    assert metrics["valid_pixels"] == phase_a.size
    assert metrics["circular_max_abs_rad"] < 1e-6
    # The production gate also requires enough spatial support.
    assert report["gate"]["checks"]["circular_mae_at_most_0.05_rad"]
    assert not report["gate"]["passed"]
    assert (tmp_path / "closure.json").exists()


def test_three_date_phase_closure_rejects_wrong_date_topology(
    tmp_path: Path,
) -> None:
    phase = np.zeros((16, 16), dtype=np.float32)
    coherence = np.ones((16, 16), dtype=np.float32)
    date_a = "2026-06-22T22:07:11"
    date_b = "2026-07-04T22:07:10"
    date_c = "2026-07-16T22:07:10"

    def write_pair(label: str, reference: str, secondary: str):
        return (
            _write_raster(
                tmp_path / f"phase_{label}.tif",
                phase,
                5.0,
                reference_date=reference,
                secondary_date=secondary,
            ),
            _write_raster(
                tmp_path / f"coherence_{label}.tif",
                coherence,
                5.0,
                reference_date=reference,
                secondary_date=secondary,
            ),
        )

    phase_ab, coherence_ab = write_pair("ab", date_a, date_b)
    phase_bc, coherence_bc = write_pair("bc", date_a, date_c)
    phase_ac, coherence_ac = write_pair("ac", date_a, date_c)

    with np.testing.assert_raises_regex(ValueError, "do not form"):
        audit_three_date_phase_closure(
            phase_ab,
            coherence_ab,
            phase_bc,
            coherence_bc,
            phase_ac,
            coherence_ac,
        )
