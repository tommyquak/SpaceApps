from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
from pyproj import Transformer
import pytest

from peatguard.timeseries.nisar_local_contrast import audit_local_contrasts, decode_gunw_mask


def _fixture(tmp_path, ref="2026-06-22", sec="2026-07-04", *, frame=176, x_shift=0):
    path = tmp_path / f"{ref}_{sec}.h5"
    x = 200000 + (np.arange(31) - 15) * 80 + x_shift
    y = 9700000 - (np.arange(31) - 15) * 80
    yy, xx = np.meshgrid(y, x, indexing="ij")
    target = (xx - 200000)**2 + (yy - 9700000)**2 <= 160**2
    with h5py.File(path, "w") as h:
        ident = h.create_group("science/LSAR/identification")
        for key, value in {
            "granuleId": f"NISAR_L2_PR_TEST_{ref}_{sec}", "compositeReleaseId": "P05023",
            "referenceZeroDopplerStartTime": ref + "T22:00:00",
            "secondaryZeroDopplerStartTime": sec + "T22:00:00",
            "trackNumber": 125, "frameNumber": frame, "orbitPassDirection": "Ascending",
            "productVersion": "1.0.11", "productSpecificationVersion": "1.5.0",
        }.items():
            ident[key] = value
        frequency = h.create_group("science/LSAR/GUNW/grids/frequencyA")
        frequency["centerFrequency"] = 1239000000.0
        unwrap = frequency.create_group("unwrappedInterferogram")
        unwrap["mask"] = np.full((31, 31), 33, dtype=np.uint32)
        group = unwrap.create_group("HH")
        group["xCoordinates"], group["yCoordinates"], group["projection"] = x, y, 32750
        for key, values, unit in (
            ("unwrappedPhase", target.astype(float) * .2, "radians"),
            ("coherenceMagnitude", np.full((31, 31), .9), "1"),
            ("ionospherePhaseScreen", target.astype(float) * .05, "radians"),
            ("connectedComponents", np.ones((31, 31), dtype=np.uint16), "1"),
        ):
            group[key] = values
            group[key].attrs["units"] = unit
        group["connectedComponents"].attrs["_FillValue"] = 65535
        processing = h.create_group("science/LSAR/GUNW/metadata/processingInformation")
        for role, date in (("Reference", ref), ("Secondary", sec)):
            processing[f"inputs/l1{role}SlcGranules"] = np.array([f"RSLC_{date}.h5"], dtype="S")
        for role in ("reference", "secondary", "unwrappedInterferogram"):
            processing[f"parameters/{role}/frequencyA/rangeBandwidth"] = 40000000.0
    lon, lat = Transformer.from_crs(32750, 4326, always_xy=True).transform(200000, 9700000)
    sites = tmp_path / "sites.geojson"
    sites.write_text(json.dumps({"type": "FeatureCollection", "features": [
        {"type": "Feature", "properties": {"site_id": "target"},
         "geometry": {"type": "Point", "coordinates": [lon, lat]}}
    ]}))
    return path, sites


def _run(paths, sites):
    return audit_local_contrasts(paths, sites, target_radius_m=160, background_inner_m=400,
                                background_outer_m=800, min_pixels=2)


def test_uint32_mask_uses_low_byte_and_separates_interpolation():
    values = np.array([33, 133, 30, 3, 255, 33 | (1 << 8), 33 | (1 << 16),
                       33 | (1 << 24), 33 | (1 << 25)], dtype=np.uint32)
    valid, interpolated = decode_gunw_mask(values)
    assert valid.tolist() == [True, False, False, False, False, False, False, True, False]
    assert interpolated.tolist() == [False] * 7 + [True, False]
    with pytest.raises(ValueError, match="schema"):
        decode_gunw_mask(values.astype(np.uint8))


def test_native_contrast_no_reference_and_explicit_corrections(tmp_path):
    path, sites = _fixture(tmp_path)
    report = _run([path], sites)
    assert report["available_sites"] == 1
    assert not report["absolute_datum"]
    assert not report["release_ready"]
    assert not report["annual_velocity_ready"]
    assert report["network"]["cycle_count"] == 0
    row = report["sites"][0]
    variants = row["variants"]
    assert variants["raw"]["contrast_rad"] == pytest.approx([.2])
    assert variants["ionosphere_subtracted"]["contrast_rad"] == pytest.approx([.15])
    assert variants["raw"]["target_count"] == variants["ionosphere_subtracted"]["target_count"]
    assert variants["raw"]["network_solution"]["epoch_covariance_rad2"] is None
    assert len(row["consumed_window_sha256"][0]["native_mask"]) == 64
    json.dumps(report, allow_nan=False)


@pytest.mark.parametrize("bad_mask", [133, 30, 255, 33 | (1 << 8)])
def test_native_validity_rejects_water_fill_and_anomalies(tmp_path, bad_mask):
    path, sites = _fixture(tmp_path)
    with h5py.File(path, "r+") as h:
        h["science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/mask"][:] = bad_mask
    assert _run([path], sites)["available_sites"] == 0


def test_ion_missing_support_applies_to_raw_too(tmp_path):
    path, sites = _fixture(tmp_path)
    with h5py.File(path, "r+") as h:
        h["science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/HH/ionospherePhaseScreen"][:] = np.nan
    row = _run([path], sites)["sites"][0]
    assert not row["variants"]["raw"]["available"]
    assert not row["variants"]["ionosphere_subtracted"]["available"]


@pytest.mark.parametrize("option", ["frame", "grid", "duplicate", "units"])
def test_rejects_incompatible_inputs(tmp_path, option):
    first, sites = _fixture(tmp_path)
    second, _ = _fixture(tmp_path, "2026-07-04", "2026-07-16",
                         frame=177 if option == "frame" else 176,
                         x_shift=5 if option == "grid" else 0)
    if option == "duplicate":
        second = first
    if option == "units":
        with h5py.File(first, "r+") as h:
            h["science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/HH/unwrappedPhase"].attrs["units"] = "m"
    with pytest.raises(ValueError):
        _run([first, second], sites)


def test_outside_site_is_reported_not_silently_dropped(tmp_path):
    path, sites = _fixture(tmp_path)
    content = json.loads(sites.read_text())
    content["features"][0]["geometry"]["coordinates"] = [0, 0]
    sites.write_text(json.dumps(content))
    report = _run([path], sites)
    assert report["total_sites"] == 1
    assert report["available_sites"] == 0
    assert report["sites"][0]["status"] == "outside_full_window"


def test_retains_valid_zero_phase(tmp_path):
    path, sites = _fixture(tmp_path)
    with h5py.File(path, "r+") as h:
        h["science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/HH/unwrappedPhase"][:] = 0
    report = _run([path], sites)
    assert report["available_sites"] == 1
    assert report["sites"][0]["variants"]["raw"]["contrast_rad"] == [0.0]


@pytest.mark.parametrize("field", ["lineage", "bandwidth"])
def test_rejects_changed_source_lineage_or_spectral_support(tmp_path, field):
    first, sites = _fixture(tmp_path)
    second, _ = _fixture(tmp_path, "2026-07-04", "2026-07-16")
    with h5py.File(second, "r+") as h:
        base = "science/LSAR/GUNW/metadata/processingInformation"
        if field == "lineage":
            h[f"{base}/inputs/l1ReferenceSlcGranules"][:] = b"changed"
        else:
            h[f"{base}/parameters/unwrappedInterferogram/frequencyA/rangeBandwidth"][()] = 20000000.0
    with pytest.raises(ValueError, match="lineage|bandwidth"):
        _run([first, second], sites)
