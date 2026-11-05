"""Behavioral tests for retrospective spatial contrast checks and site sampling."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import transform

from peatguard.analysis.contrast_validation import paired_metrics, spatial_group_diagnostics


def script(name):
    path = Path(__file__).parents[1] / "scripts" / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_raster(path, data, nodata=None):
    with rasterio.open(path, "w", driver="GTiff", height=data.shape[0], width=data.shape[1],
                       count=1, dtype="float32", crs="EPSG:32749",
                       transform=from_origin(700000, 9700000, 80, 80), nodata=nodata) as out:
        out.write(data.astype(np.float32), 1)


def test_real_zero_is_valid_and_outside_coordinate_cannot_wrap(tmp_path):
    module = script("site_baseline_table.py")
    path = tmp_path / "zeros.tif"
    write_raster(path, np.zeros((20, 20)))
    lon, lat = transform("EPSG:32749", "EPSG:4326", [700800, 698000], [9699200, 9702000])
    samples = module.sample(path, 1, lon, lat, half_m=80, min_px=1)
    assert samples[0] == (0.0, 0.0, 9)
    assert samples[1] == (None, None, 0)
    assert module.sample(path, 1, lon[:1], lat[:1], half_m=80, min_px=1,
                         zero_is_nodata=True)[0] == (None, None, 0)


def test_sampling_respects_mask_and_shared_support(tmp_path):
    module = script("validate_nisar_site_contrasts.py")
    descending, ascending = tmp_path / "d.tif", tmp_path / "a.tif"
    d, a = np.ones((30, 30)), np.ones((30, 30))*2
    d[14:17, 14] = -9999
    a[14:17, 16] = -9999
    write_raster(descending, d, nodata=-9999)
    write_raster(ascending, a, nodata=-9999)
    lon, lat = transform("EPSG:32749", "EPSG:4326", [701240], [9698760])
    features = [{"geometry": {"type": "Point", "coordinates": [lon[0], lat[0]]},
                 "properties": {"site_id": "A"}}]
    rows, _ = module.sample_common_support(descending, ascending, features,
                                          half_window_m=80, min_pixels=1, min_fraction=0.3)
    assert rows[0]["paired_pixels"] == 3
    assert rows[0]["descending_mm"] == 1
    assert rows[0]["ascending_mm"] == 2
    rows, _ = module.sample_common_support(descending, ascending, features,
                                          half_window_m=80, min_pixels=1, min_fraction=0.5)
    assert rows[0]["status"] == "insufficient_common_support"


def test_many_sites_in_one_group_do_not_claim_many_replicates():
    x = np.arange(50, dtype=float)
    report = spatial_group_diagnostics(x, x, np.zeros((50, 2)), draws=100)
    assert report["groups"] == 1
    assert report["status"] == "insufficient_spatial_groups"
    assert report["pearson_ci95"] is None
    assert report["release_ready"] is False


def test_bootstrap_keeps_pairing_and_never_promotes_release():
    x = np.arange(12, dtype=float)
    xy = np.column_stack((x*6000, np.zeros(12)))
    report = spatial_group_diagnostics(x, x+4, xy, draws=100)
    assert report["groups"] == 12
    assert report["pearson_ci95"] == pytest.approx([1, 1])
    assert report["median_difference_ci95_mm"] == [-4, -4]
    assert report["metrics"]["rmse_difference_mm"] == 4
    assert report["release_ready"] is False


def test_degenerate_and_missing_values_are_explicit():
    report = paired_metrics(np.zeros(4), np.zeros(4))
    assert report["pearson_r"] is None
    assert report["same_sign_fraction_nonzero"] is None
    with pytest.raises(ValueError, match="finite"):
        paired_metrics(np.array([0, np.nan]), np.array([0, 1]))


def test_missing_seasonal_layer_does_not_abort_table(tmp_path, monkeypatch):
    module = script("site_baseline_table.py")
    sites = tmp_path / "sites.json"
    sites.write_text('{"features":[{"geometry":{"coordinates":[114,-2.8]},"properties":{}}]}')
    monkeypatch.setattr(module, "SITES", sites)
    monkeypatch.setattr(module, "OUT", tmp_path / "out")
    monkeypatch.setattr(module, "L", {})
    module.main([])
    assert (tmp_path / "out" / "dam_sites_measured_baseline.csv").exists()


def test_site_table_preserves_site_id_and_rejects_duplicates(tmp_path, monkeypatch):
    import csv
    import json

    module = script("site_baseline_table.py")
    sites = tmp_path / "sites.json"
    feature = {"geometry": {"coordinates": [114, -2.8]},
               "properties": {"site_id": "D-17", "id": "wrong_fallback"}}
    sites.write_text(json.dumps({"features": [feature]}))
    monkeypatch.setattr(module, "L", {})
    module.main(["--sites", str(sites), "--output-dir", str(tmp_path / "out")])
    with (tmp_path / "out/dam_sites_measured_baseline.csv").open() as stream:
        assert next(csv.DictReader(stream))["site"] == "D-17"
    sites.write_text(json.dumps({"features": [feature, feature]}))
    with pytest.raises(ValueError, match="unique"):
        module.main(["--sites", str(sites), "--output-dir", str(tmp_path / "duplicate")])
