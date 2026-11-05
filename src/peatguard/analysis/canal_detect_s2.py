"""Finer canal detection from Sentinel-2 optical imagery.

The original SAR canal detector (analysis/canal_detect.py) thresholds a
VV-backscatter percentile. That catches wide canals but misses thin ones:
a 10-20 m drainage canal only partly fills a SAR pixel, so its averaged
backscatter never crosses the threshold.

This module detects canals from Sentinel-2 instead, with two changes that
recover finer canals:

  1. NDWI (Normalised Difference Water Index, green vs NIR) -- drainage
     canals hold water / wet soil and stand out brightly, at 10 m.
  2. A Sato tubeness (ridge) filter -- the key step. Instead of a flat
     NDWI threshold, the ridge filter responds to ELONGATED, locally
     bright structures at canal-width scales. A thin canal whose NDWI is
     diluted by mixed pixels still produces a strong ridge response
     because it is linear and brighter than its immediate surroundings.
     This is what lets the detector pick up fine canals a threshold drops.

  canal = (NDWI > t_ndwi)  OR  (ridge_response > t_ridge)

Output: a binary canal_mask COG on an EPSG:4326 grid. Sentinel-2 is free
and global, so this runs over areas the InSAR footprint never covered.

Usage:
    python -m peatguard.analysis.canal_detect_s2 \
        --bbox 114.21 -2.71 114.52 -2.42 \
        --out  .../canal_mask_s2.tif
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import rasterio
from rasterio.crs import CRS
from scipy.ndimage import label
from skimage.filters import sato

from peatguard.ingest.sentinel2 import composite

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("canal_detect_s2")

# ---- Tunables ---------------------------------------------------------------
RIDGE_SIGMAS = (1.0, 2.0, 3.5)   # canal half-widths in pixels (~15-50 m)
NDWI_THRESHOLD = 0.0             # NDWI above this is wet enough to be canal/water
RIDGE_PERCENTILE = 96.0          # ridge response above this percentile = canal
MIN_SEGMENT_PX = 12              # drop isolated specks smaller than this
# -----------------------------------------------------------------------------


def detect_canals(ndwi: np.ndarray) -> tuple[np.ndarray, dict]:
    """NDWI -> binary canal mask via threshold OR ridge response."""
    finite = np.isfinite(ndwi)
    # Sato needs a filled array; fill gaps with a dry value so they don't
    # masquerade as ridges.
    fill = np.nanpercentile(ndwi[finite], 20) if finite.any() else 0.0
    ndwi_filled = np.where(finite, ndwi, fill).astype(np.float32)

    # Ridge (tubeness) response -- bright elongated structures.
    ridge = sato(ndwi_filled, sigmas=RIDGE_SIGMAS, black_ridges=False)
    ridge = np.nan_to_num(ridge, nan=0.0)
    r_thresh = float(np.percentile(ridge[finite], RIDGE_PERCENTILE)) if finite.any() else 0.0

    wet = finite & (ndwi > NDWI_THRESHOLD)
    ridgey = finite & (ridge > r_thresh)
    canal = wet | ridgey

    # Remove isolated specks (false positives); keeps linear structures.
    lab, n = label(canal)
    if n:
        sizes = np.bincount(lab.ravel())
        small = np.isin(lab, np.where(sizes < MIN_SEGMENT_PX)[0])
        canal &= ~small

    stats = {
        "ndwi_threshold": NDWI_THRESHOLD,
        "ridge_threshold": round(r_thresh, 5),
        "pixels_wet_only": int((wet & ~ridgey).sum()),
        "pixels_ridge_recovered": int((ridgey & ~wet).sum()),
        "canal_pixels": int(canal.sum()),
        "canal_pct": round(100.0 * float(canal.mean()), 3),
    }
    return canal.astype(np.uint8), stats


def build_canal_mask(bbox, out_path: Path, start="2024-06-01",
                     end="2024-10-31", res_deg=0.00015) -> dict:
    """Composite Sentinel-2 over bbox, detect canals, write a COG."""
    comp = composite(bbox, start=start, end=end, res_deg=res_deg)
    canal, stats = detect_canals(comp["ndwi"])
    h, w = comp["shape"]

    profile = {
        "driver": "COG", "dtype": "uint8", "count": 1,
        "width": w, "height": h, "crs": CRS.from_epsg(4326),
        "transform": comp["transform"], "nodata": 0, "compress": "deflate",
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(canal, 1)
    stats["n_scenes"] = comp["n_scenes"]
    stats["output"] = str(out_path)
    stats["grid"] = f"{w}x{h}"
    log.info("wrote %s", out_path)
    return stats


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    # Comma-separated so negative latitudes aren't mistaken for option flags.
    p.add_argument("--bbox", required=True, help="W,S,E,N (comma-separated)")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--res-deg", type=float, default=0.00015)
    p.add_argument("--start", default="2024-06-01")
    p.add_argument("--end", default="2024-10-31")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    bbox = tuple(float(x) for x in args.bbox.split(","))
    stats = build_canal_mask(bbox, args.out,
                             start=args.start, end=args.end, res_deg=args.res_deg)
    log.info("")
    log.info("=== CANAL DETECTION STATS ===")
    for k, v in stats.items():
        log.info("  %s: %s", k, v)


if __name__ == "__main__":
    main()
