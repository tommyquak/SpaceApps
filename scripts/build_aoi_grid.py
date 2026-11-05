"""Assemble the 2x2 AOI grid products.

Two modes:

  mosaic       Mosaic the per-cell Sentinel-2 canal masks into one
               grid-wide canal_mask COG on a common EPSG:4326 grid.

  restoration  From the grid canal_mask + grid fire_history, build the
               grid-wide restoration-priority surface:
                 1. canal-proximity risk -- distance transform of the
                    canal mask, using a project-assumed linear decay over
                    1200 m (the same heuristic as analysis/risk_score.py).
                 2. restoration_zones = 0.62 * proximity + 0.38 * fire
                    (the proximity:fire weight ratio of the project
                    model, renormalised; InSAR subsidence is not
                    available outside the original Sebangau footprint
                    so the grid surface uses proximity + fire, both
                    real grid-wide).

The original InSAR canal_risk stays a separate SBG-C-only layer.

Usage:
    python build_aoi_grid.py mosaic --cells A1=/tmp/canal_A1.tif,... \
        --bbox 114.21,-3.00,114.83,-2.42 --res-deg 0.00018 --out grid_canal.tif
    python build_aoi_grid.py restoration --canal grid_canal.tif \
        --fire grid_fire.tif --out grid_restoration.tif
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.transform import from_bounds
from rasterio.warp import Resampling, reproject
from scipy.ndimage import distance_transform_edt

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("build_aoi_grid")

MAX_INFLUENCE_M = 1200.0   # project sensitivity prior; field calibration pending
W_PROXIMITY = 0.62         # proximity : fire = 0.25 : 0.15, renormalised
W_FIRE = 0.38
COG = dict(driver="COG", compress="deflate")


def _cog_profile(width, height, transform, dtype, nodata):
    return dict(COG, dtype=dtype, count=1, width=width, height=height,
                crs=CRS.from_epsg(4326), transform=transform, nodata=nodata)


def mosaic(cell_paths: list[str], bbox, res_deg: float, out: Path) -> dict:
    """Reproject each cell canal mask onto the common grid; OR-combine."""
    w = int(round((bbox[2] - bbox[0]) / res_deg))
    h = int(round((bbox[3] - bbox[1]) / res_deg))
    transform = from_bounds(*bbox, w, h)
    grid = np.zeros((h, w), dtype=np.uint8)
    for p in cell_paths:
        with rasterio.open(p) as src:
            tmp = np.zeros((h, w), dtype=np.uint8)
            reproject(
                source=src.read(1), destination=tmp,
                src_transform=src.transform, src_crs=src.crs,
                dst_transform=transform, dst_crs="EPSG:4326",
                resampling=Resampling.nearest,
            )
        grid |= (tmp > 0).astype(np.uint8)
    out.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out, "w", **_cog_profile(w, h, transform, "uint8", 0)) as d:
        d.write(grid, 1)
    log.info("mosaic: %dx%d grid, %d canal px (%.2f%%) -> %s",
             w, h, int(grid.sum()), 100.0 * grid.mean(), out)
    return {"width": w, "height": h, "canal_pct": round(100.0 * float(grid.mean()), 3)}


def restoration(canal_path: Path, fire_path: Path, out: Path) -> dict:
    """Build the grid restoration-priority surface from canal proximity
    + fire history."""
    with rasterio.open(canal_path) as r:
        canal = r.read(1)
        transform, crs = r.transform, r.crs
        h, w = canal.shape
        # metres per pixel (latitude-corrected)
        import math
        clat = transform.f + (h / 2) * transform.e
        px_m = abs(transform.a) * 111_320.0 * math.cos(math.radians(clat))

    # Distance transform: metres to the nearest canal pixel.
    dist_px = distance_transform_edt(canal == 0)
    dist_m = dist_px * px_m
    proximity = np.clip(1.0 - dist_m / MAX_INFLUENCE_M, 0.0, 1.0).astype(np.float32)

    # Fire history, reprojected onto the canal grid if necessary.
    with rasterio.open(fire_path) as r:
        if (r.width, r.height) == (w, h):
            fire = r.read(1).astype(np.float32)
        else:
            fire = np.zeros((h, w), dtype=np.float32)
            reproject(
                source=r.read(1), destination=fire,
                src_transform=r.transform, src_crs=r.crs,
                dst_transform=transform, dst_crs=crs,
                resampling=Resampling.bilinear,
            )
    fire = np.clip(np.nan_to_num(fire, nan=0.0), 0.0, 1.0)

    rest = (W_PROXIMITY * proximity + W_FIRE * fire).astype(np.float32)
    rest = np.clip(rest, 0.0, 1.0)

    with rasterio.open(out, "w", **_cog_profile(w, h, transform, "float32", -9999.0)) as d:
        d.write(rest, 1)
    log.info("restoration: mean %.3f, p90 %.3f, max %.3f -> %s",
             float(rest.mean()), float(np.percentile(rest, 90)), float(rest.max()), out)
    return {
        "mean": round(float(rest.mean()), 4),
        "pct_high_gt_0.6": round(100.0 * float((rest > 0.6).mean()), 2),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="mode", required=True)

    m = sub.add_parser("mosaic")
    m.add_argument("--cells", required=True,
                   help="comma-separated list of per-cell canal mask tifs")
    m.add_argument("--bbox", required=True, help="W,S,E,N")
    m.add_argument("--res-deg", type=float, default=0.00018)
    m.add_argument("--out", type=Path, required=True)

    r = sub.add_parser("restoration")
    r.add_argument("--canal", type=Path, required=True)
    r.add_argument("--fire", type=Path, required=True)
    r.add_argument("--out", type=Path, required=True)

    args = ap.parse_args()
    if args.mode == "mosaic":
        bbox = tuple(float(x) for x in args.bbox.split(","))
        cells = [c.strip() for c in args.cells.split(",") if c.strip()]
        mosaic(cells, bbox, args.res_deg, args.out)
    else:
        stats = restoration(args.canal, args.fire, args.out)
        log.info("stats: %s", stats)


if __name__ == "__main__":
    main()
