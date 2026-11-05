"""Regression tests for ordinary out-of-footprint TiTiler requests."""

from __future__ import annotations

import importlib
from pathlib import Path

import numpy as np
import pytest
import rasterio
from fastapi.testclient import TestClient
from rasterio.warp import transform


@pytest.fixture(autouse=True)
def synthetic_webdata(tmp_path: Path, monkeypatch):
    """Exercise the real reader with tiny fixtures, not ignored local products."""
    dashboard = importlib.import_module("peatguard.dashboard.app")
    webdata = tmp_path / "webdata"
    webdata.mkdir()
    monkeypatch.setattr(dashboard, "_WEBDATA_DIR", webdata)
    coherence = np.full((32, 32), 127, dtype=np.uint8)
    coherence[0, :4] = [255, 0, 127, 254]
    fixtures = (
        ("thermal_stress.tif", np.ones((32, 32), dtype=np.float32), -9999.0, 1.0),
        ("nisar_coherence_fine_5m.tif", coherence, 255, 1.0 / 254.0),
    )
    for filename, values, nodata, scale in fixtures:
        with rasterio.open(
            webdata / f"DEMO_products_{filename}", "w", driver="GTiff",
            width=32, height=32, count=1, dtype=values.dtype,
            crs="EPSG:32750",
            transform=rasterio.transform.from_origin(240000.0, 9700000.0, 5.0, 5.0),
            nodata=nodata, tiled=True, blockxsize=16, blockysize=16,
        ) as dataset:
            dataset.write(values, 1)
            dataset.scales = (scale,)
            dataset.offsets = (0.0,)
    return webdata


def test_outside_point_and_tile_are_not_server_errors(
    tmp_path: Path,
    monkeypatch,
    synthetic_webdata: Path,
) -> None:

    from peatguard.dashboard.app import app

    raster = synthetic_webdata / "DEMO_products_thermal_stress.tif"
    assert raster.exists()
    client = TestClient(app)

    point = client.get(
        "/cog/point/0,0",
        params={"url": str(raster)},
    )
    tile = client.get(
        "/cog/tiles/WebMercatorQuad/1/0/0.png",
        params={"url": str(raster)},
    )

    assert point.status_code == 404
    assert point.json() == {
        "detail": "Point is outside the raster footprint"
    }
    assert tile.status_code == 204
    assert tile.content == b""


def test_cog_reader_rejects_uncatalogued_paths(
    tmp_path: Path,
    monkeypatch,
) -> None:

    from peatguard.dashboard.app import app

    client = TestClient(app)
    local_file = client.get(
        "/cog/info",
        params={"url": "/etc/passwd"},
    )
    internal_url = client.get(
        "/cog/info",
        params={"url": "http://127.0.0.1:8080/health"},
    )

    assert local_file.status_code == 403
    assert internal_url.status_code == 403


@pytest.mark.parametrize("encoded", [0, 127, 254])
def test_fine_coherence_point_applies_uint8_scale(
    tmp_path: Path,
    monkeypatch,
    encoded: int,
) -> None:

    from peatguard.dashboard.app import PRODUCTS, app, resolve_cog_url

    metadata = PRODUCTS["nisar_coherence_fine_5m.tif"]
    assert metadata["point_unscale"] is True
    raster = Path(resolve_cog_url("nisar_coherence_fine_5m.tif"))
    assert raster.exists()

    with rasterio.open(raster) as dataset:
        raw = dataset.read(1, window=((0, 512), (0, 512)))
        rows, columns = np.where(raw == encoded)
        assert rows.size
        row = int(rows[0])
        column = int(columns[0])
        x, y = dataset.xy(row, column)
        longitude, latitude = transform(
            dataset.crs, "EPSG:4326", [x], [y]
        )
        expected = float(raw[row, column]) * dataset.scales[0]

    response = TestClient(app).get(
        f"/cog/point/{longitude[0]},{latitude[0]}",
        params={"url": str(raster), "unscale": "true"},
    )
    assert response.status_code == 200
    assert response.json()["values"][0] == pytest.approx(expected, abs=1e-7)
