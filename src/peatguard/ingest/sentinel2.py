"""Sentinel-2 cloud-free dry-season composite.

Pulls Sentinel-2 L2A surface-reflectance COGs from the Element84 Earth
Search STAC (open AWS bucket, no auth) and builds a cloud-masked median
composite of the bands needed for canal detection:

    green (B03), NIR (B08)  ->  NDWI = (green - NIR) / (green + NIR)

Clouds are screened with the L2A Scene Classification Layer (SCL).
Sentinel-2 is 10 m and free over all of Indonesia, so it can map canals
in areas the InSAR footprint never covered -- which is what lets the
canal-detection step run beyond the original Sebangau AOI.

Output: a green + NIR median composite resampled onto a common EPSG:4326
grid at a requested resolution, returned as numpy arrays + an affine
transform for analysis/canal_detect_s2.py to consume.
"""

from __future__ import annotations

import logging

import numpy as np
import rasterio
import requests
from rasterio.warp import Resampling, reproject
from rasterio.transform import from_bounds

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("sentinel2")

STAC_URL = "https://earth-search.aws.element84.com/v1/search"
# SCL classes that are cloud / shadow / saturated -> masked out.
SCL_BAD = {0, 1, 3, 8, 9, 10, 11}   # nodata, saturated, shadow, clouds, cirrus, snow


def search_scenes(bbox, start, end, max_cloud=30, limit=60) -> list[dict]:
    """STAC search for Sentinel-2 L2A scenes intersecting bbox."""
    body = {
        "collections": ["sentinel-2-l2a"],
        "bbox": list(bbox),
        "datetime": f"{start}T00:00:00Z/{end}T23:59:59Z",
        "query": {"eo:cloud_cover": {"lt": max_cloud}},
        "limit": limit,
    }
    r = requests.post(STAC_URL, json=body, timeout=60)
    r.raise_for_status()
    feats = r.json().get("features", [])
    log.info("Sentinel-2 STAC: %d scenes over bbox, %s..%s, cloud<%d%%",
             len(feats), start, end, max_cloud)
    return feats


def _read_window(href: str, bbox, out_h, out_w):
    """Read a COG band, reprojected/clipped onto the target 4326 grid."""
    dst = np.full((out_h, out_w), np.nan, dtype=np.float32)
    dst_transform = from_bounds(*bbox, out_w, out_h)
    with rasterio.open(f"/vsicurl/{href}") as src:
        reproject(
            source=rasterio.band(src, 1),
            destination=dst,
            dst_transform=dst_transform,
            dst_crs="EPSG:4326",
            dst_nodata=np.nan,
            resampling=Resampling.bilinear,
        )
    return dst


def composite(bbox, start="2024-06-01", end="2024-10-31",
              res_deg=0.00015, max_cloud=30) -> dict:
    """Build a cloud-masked median green + NIR composite over bbox.

    Returns {green, nir, ndwi, transform, bbox, n_scenes}."""
    feats = search_scenes(bbox, start, end, max_cloud=max_cloud)
    if not feats:
        raise RuntimeError("no Sentinel-2 scenes for the requested AOI/window")

    w = int(round((bbox[2] - bbox[0]) / res_deg))
    h = int(round((bbox[3] - bbox[1]) / res_deg))
    log.info("compositing onto %dx%d grid (~%.0f m px)", w, h, res_deg * 111320)

    green_stack, nir_stack = [], []
    used = 0
    for f in feats:
        a = f["assets"]
        try:
            scl = _read_window(a["scl"]["href"], bbox, h, w)
            green = _read_window(a["green"]["href"], bbox, h, w)
            nir = _read_window(a["nir"]["href"], bbox, h, w)
        except Exception as e:  # noqa: BLE001 - skip unreadable scenes
            log.warning("  skip %s: %s", f["id"], e)
            continue
        bad = np.isin(np.nan_to_num(scl, nan=0).astype(np.int16), list(SCL_BAD))
        green[bad] = np.nan
        nir[bad] = np.nan
        if np.isfinite(green).mean() < 0.05:
            continue
        green_stack.append(green)
        nir_stack.append(nir)
        used += 1
    if used == 0:
        raise RuntimeError("no usable (cloud-free enough) Sentinel-2 scenes")
    log.info("composited %d cloud-screened scenes", used)

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        green_med = np.nanmedian(np.stack(green_stack), axis=0).astype(np.float32)
        nir_med = np.nanmedian(np.stack(nir_stack), axis=0).astype(np.float32)

    denom = green_med + nir_med
    ndwi = np.where(denom > 0, (green_med - nir_med) / denom, np.nan).astype(np.float32)
    return {
        "green": green_med, "nir": nir_med, "ndwi": ndwi,
        "transform": from_bounds(*bbox, w, h), "bbox": tuple(bbox),
        "n_scenes": used, "shape": (h, w),
    }
