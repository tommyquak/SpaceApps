"""PlanetScope SuperDove (PSB.SD) 8-band surface-reflectance compositing.

Reads one or more pre-downloaded ortho_analytic_8b_sr GeoTIFFs (Planet Data
API direct asset-download), reprojects + clips each onto a common EPSG:4326
grid at a target resolution, and returns a green + NIR + NDWI median
composite for canal-detection. Pixel size defaults to ~3 m -- the native
SuperDove resolution.

Why a local-files ingester rather than calling Planet's search/activate/
download here: assets take 5-15 min to activate and ~100-1000 MB to fetch.
Downloading inside the analysis run is too brittle and burns Planet API
calls every time the user iterates on detector parameters. Decouple: pull
scenes once with the helper script, then run this against the local files.

SuperDove 8b SR band order (Planet docs):
  1 coastal_blue   485 nm
  2 blue           490 nm
  3 green_i        531 nm
  4 green          565 nm
  5 yellow         610 nm
  6 red            665 nm
  7 rededge        705 nm
  8 nir            865 nm

NDWI = (green - nir) / (green + nir)   uses bands 4 and 8.

SR values are uint16 stored as reflectance * 10000; we keep that scaling
because the detector works on ratios (NDWI) where the 10000x factor cancels.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import rasterio
from rasterio.transform import from_bounds
from rasterio.warp import Resampling, reproject

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("planetscope")

# SuperDove 8b SR band indices (1-based, as rasterio numbers them)
BAND_GREEN = 4
BAND_RED = 6
BAND_REDEDGE = 7
BAND_NIR = 8

# UDM2 band schema (Planet PSScene docs):
#   1 clear  2 snow  3 shadow  4 light_haze  5 heavy_haze  6 cloud
#   7 confidence  8 usable_data
# We treat 3 (shadow) + 5 (heavy_haze) + 6 (cloud) as definitely-bad; light_haze
# is left in (the SR atmospheric correction handles it acceptably for indices).
UDM2_BAD_BANDS = (3, 5, 6)


def _read_band_to_grid(src_path: Path, band: int, bbox, h: int, w: int) -> np.ndarray:
    """Read one band, reproject + clip onto the EPSG:4326 output grid."""
    dst = np.full((h, w), np.nan, dtype=np.float32)
    dst_transform = from_bounds(*bbox, w, h)
    with rasterio.open(src_path) as src:
        # uint16 with 0 = nodata; mask zeros before reprojection so they
        # don't blend into resampled neighbours.
        src_arr = src.read(band, masked=False).astype(np.float32)
        src_arr[src_arr == 0] = np.nan
        reproject(
            source=src_arr,
            destination=dst,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=dst_transform,
            dst_crs="EPSG:4326",
            src_nodata=np.nan,
            dst_nodata=np.nan,
            resampling=Resampling.bilinear,
        )
    return dst


def _read_udm2_bad_mask(udm2_path: Path, bbox, h: int, w: int) -> np.ndarray:
    """Read a UDM2 raster, return a boolean BAD mask on the output grid.

    BAD = any pixel where shadow (3), heavy_haze (5), or cloud (6) is set.
    Reprojected with nearest-neighbour to preserve the mask boundaries.
    """
    bad_total = np.zeros((h, w), dtype=bool)
    dst_transform = from_bounds(*bbox, w, h)
    with rasterio.open(udm2_path) as src:
        for b in UDM2_BAD_BANDS:
            if b > src.count:
                continue
            dst = np.zeros((h, w), dtype=np.uint8)
            reproject(
                source=src.read(b),
                destination=dst,
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=dst_transform,
                dst_crs="EPSG:4326",
                resampling=Resampling.nearest,
            )
            bad_total |= dst > 0
    return bad_total


def _udm2_path_for(scene_path: Path) -> Optional[Path]:
    """Locate the companion UDM2 file for an 8b SR scene.

    Naming convention used by our download scripts:
        <SCENE_ID>_8b_sr.tif  <-> <SCENE_ID>_udm2.tif
    """
    name = scene_path.name
    if "_8b_sr." in name:
        candidate = scene_path.with_name(name.replace("_8b_sr.", "_udm2."))
        return candidate if candidate.exists() else None
    return None


def composite(scene_paths: Iterable[Path], bbox, res_deg: float = 0.00003,
              extended: bool = False,
              ndwi_wet_threshold: float = 0.0,
              nir_water_max: float = 1000.0) -> dict:
    """Build a median green + NIR + NDWI composite from SuperDove 8b SR scenes.

    Args:
        scene_paths: iterable of ortho_analytic_8b_sr GeoTIFFs (any CRS).
        bbox: (W, S, E, N) in EPSG:4326.
        res_deg: output pixel size in degrees. Default 0.00003 ~ 3.3 m near
            the equator, matching the native SuperDove resolution.
        extended: also composite red + rededge medians and a per-pixel
            wet-persistence map (fraction of valid scenes where the pixel was
            wet by NDWI or low-NIR). Used by canal detection v2; off by default
            so the v1 path stays light.
        ndwi_wet_threshold, nir_water_max: the wet criterion for persistence
            (a pixel is "wet" in a scene if NDWI > threshold OR NIR < max).

    Returns:
        dict with keys: green, nir, ndwi (float32), transform, bbox, n_scenes,
        shape (H, W). When extended=True also: red, rededge, persistence
        (float32 in [0,1], nan where never observed). NDWI is nan over
        fully-empty pixels.
    """
    scene_paths = [Path(p) for p in scene_paths]
    if not scene_paths:
        raise RuntimeError("no PlanetScope scenes provided")

    w = int(round((bbox[2] - bbox[0]) / res_deg))
    h = int(round((bbox[3] - bbox[1]) / res_deg))
    log.info("planetscope composite -> %d x %d grid (~%.1f m px)%s",
             w, h, res_deg * 111320, " [extended]" if extended else "")
    log.info("compositing %d scene(s) over bbox=%s", len(scene_paths), bbox)

    green_stack, nir_stack = [], []
    red_stack, rededge_stack = [], []
    # Persistence accumulators (extended only).
    obs_count = np.zeros((h, w), dtype=np.float32)
    wet_count = np.zeros((h, w), dtype=np.float32)
    used = 0
    udm2_used = 0
    for p in scene_paths:
        try:
            g = _read_band_to_grid(p, BAND_GREEN, bbox, h, w)
            n = _read_band_to_grid(p, BAND_NIR, bbox, h, w)
            r = _read_band_to_grid(p, BAND_RED, bbox, h, w) if extended else None
            re = _read_band_to_grid(p, BAND_REDEDGE, bbox, h, w) if extended else None
        except Exception as e:  # noqa: BLE001 -- skip unreadable scenes
            log.warning("  skip %s: %s", p.name, e)
            continue

        # Apply UDM2 cloud / shadow / heavy-haze mask if the companion file
        # exists. Bad pixels go to nan so they fall out of the median composite.
        udm2 = _udm2_path_for(p)
        if udm2 is not None:
            try:
                bad = _read_udm2_bad_mask(udm2, bbox, h, w)
                g[bad] = np.nan
                n[bad] = np.nan
                if extended:
                    r[bad] = np.nan
                    re[bad] = np.nan
                udm2_used += 1
            except Exception as e:  # noqa: BLE001
                log.warning("  UDM2 read failed for %s: %s -- proceeding without", p.name, e)

        # Drop scenes whose AOI coverage is too sparse to be useful (<1%).
        cov = float(np.isfinite(g).mean())
        if cov < 0.01:
            log.info("  skip %s: only %.1f%% AOI coverage", p.name, 100 * cov)
            continue
        udm2_note = " (UDM2-masked)" if udm2 is not None else " (no UDM2)"
        log.info("  used %s: %.1f%% AOI coverage%s", p.name, 100 * cov, udm2_note)
        green_stack.append(g)
        nir_stack.append(n)
        if extended:
            red_stack.append(r)
            rededge_stack.append(re)
            valid = np.isfinite(g) & np.isfinite(n)
            with np.errstate(invalid="ignore"):
                ndwi_s = np.where((g + n) > 0, (g - n) / (g + n), np.nan)
                wet_s = valid & ((ndwi_s > ndwi_wet_threshold) | (n < nir_water_max))
            obs_count += valid.astype(np.float32)
            wet_count += wet_s.astype(np.float32)
        used += 1

    if used == 0:
        raise RuntimeError("no PlanetScope scenes intersect the requested bbox")

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        green_med = np.nanmedian(np.stack(green_stack), axis=0).astype(np.float32)
        nir_med = np.nanmedian(np.stack(nir_stack), axis=0).astype(np.float32)

    denom = green_med + nir_med
    ndwi = np.where(denom > 0, (green_med - nir_med) / denom, np.nan).astype(np.float32)
    log.info("composite finished: %d scenes used, %d with UDM2 cloud-mask applied",
             used, udm2_used)
    out = {
        "green": green_med, "nir": nir_med, "ndwi": ndwi,
        "transform": from_bounds(*bbox, w, h), "bbox": tuple(bbox),
        "n_scenes": used, "n_scenes_udm2_masked": udm2_used,
        "shape": (h, w),
    }
    if extended:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            out["red"] = np.nanmedian(np.stack(red_stack), axis=0).astype(np.float32)
            out["rededge"] = np.nanmedian(np.stack(rededge_stack), axis=0).astype(np.float32)
        out["persistence"] = np.where(
            obs_count > 0, wet_count / np.maximum(obs_count, 1), np.nan
        ).astype(np.float32)
    return out
