"""Fire-history risk layer from MODIS MCD64A1 burned-area granules.

Turns a directory of monthly MCD64A1 v061 HDF granules into a single
fire_history.tif raster (0-1) aligned to the canal_risk grid. The score
combines two physically distinct signals:

  * Recency   -- how recently the pixel last burned, exponential decay
                 (half-life 5 yr). A fresh burn scar is the strongest
                 single indicator that a peatland needs restoration now.
  * Frequency -- how many distinct years the pixel burned in the window.
                 Chronic re-burning marks structurally degraded, drained
                 peat that will keep burning until it is rewetted.

  fire_history = 0.6 * recency + 0.4 * frequency

MCD64A1 is HDF4/HDF-EOS2; this module is read with pyhdf and must run in
an environment that has pyhdf + rasterio (the project .firegdal env).

Usage (with the conda env that has pyhdf):
    .firegdal/bin/python -m peatguard.analysis.fire_history \
        --granules Datasets/fire_mcd64a1 \
        --ref-grid peatguard/.../webdata/DEMO_products_canal_risk.tif \
        --out peatguard/.../webdata/DEMO_products_fire_history.tif
"""

from __future__ import annotations

import argparse
import glob
import logging
import re
from pathlib import Path

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.warp import Resampling, reproject

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("fire_history")

# ---- Tunables ---------------------------------------------------------------
REF_YEAR = 2024          # analysis "present" -- recency measured back from here
RECENCY_HALF_LIFE = 5.0  # years; a burn this old contributes 0.5 recency
FREQ_CAP = 3             # burning in this many distinct years -> frequency 1.0
W_RECENCY = 0.6
W_FREQUENCY = 0.4
# MODIS sinusoidal CRS (the projection MCD64A1 ships in).
MODIS_SINU = CRS.from_proj4(
    "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 "
    "+units=m +no_defs"
)
BURN_DATE_SDS = "Burn Date"
# -----------------------------------------------------------------------------


def _year_of(granule_path: str) -> int:
    """MCD64A1.A2015001.h29v09.061.*.hdf -> 2015."""
    m = re.search(r"\.A(\d{4})\d{3}\.", Path(granule_path).name)
    if not m:
        raise ValueError(f"cannot parse year from {granule_path}")
    return int(m.group(1))


def _tile_bounds_from_structmeta(struct_meta: str) -> tuple[float, float, float, float]:
    """Parse UpperLeftPointMtrs / LowerRightMtrs from HDF-EOS StructMetadata.

    Returns (x_min, y_min, x_max, y_max) in MODIS sinusoidal metres."""
    ul = re.search(r"UpperLeftPointMtrs=\(([-\d.]+),([-\d.]+)\)", struct_meta)
    lr = re.search(r"LowerRightMtrs=\(([-\d.]+),([-\d.]+)\)", struct_meta)
    if not ul or not lr:
        raise ValueError("StructMetadata missing UL/LR corner points")
    x_min, y_max = float(ul.group(1)), float(ul.group(2))
    x_max, y_min = float(lr.group(1)), float(lr.group(2))
    return x_min, y_min, x_max, y_max


def _read_burn_date(hdf_path: str):
    """Read the MCD64A1 'Burn Date' SDS. Returns (array, transform).

    Burn Date values: 0 = unburned, 1-366 = ordinal burn day,
    -1 = unmapped, -2 = water."""
    from pyhdf.SD import SD, SDC

    sd = SD(hdf_path, SDC.READ)
    try:
        arr = sd.select(BURN_DATE_SDS)[:].astype(np.int16)
        struct = sd.attributes()["StructMetadata.0"]
    finally:
        sd.end()
    x_min, y_min, x_max, y_max = _tile_bounds_from_structmeta(struct)
    h, w = arr.shape
    transform = rasterio.transform.from_bounds(x_min, y_min, x_max, y_max, w, h)
    return arr, transform


def build_fire_history(granule_dir: Path, ref_grid: Path, out_path: Path) -> dict:
    """Build fire_history.tif. Returns a dict of validation statistics."""
    granules = sorted(glob.glob(str(granule_dir / "MCD64A1.*.hdf")))
    if not granules:
        raise RuntimeError(f"no MCD64A1 granules in {granule_dir}")
    log.info("found %d MCD64A1 granules", len(granules))

    # Group granules by year.
    by_year: dict[int, list[str]] = {}
    for g in granules:
        by_year.setdefault(_year_of(g), []).append(g)
    years = sorted(by_year)
    log.info("years covered: %d-%d (%d months total)",
             years[0], years[-1], len(granules))

    # Per-pixel: count of distinct years burned, and most recent burn year.
    src_transform = None
    src_shape = None
    burn_count = None
    last_burn_year = None

    for year in years:
        burned_this_year = None
        for g in by_year[year]:
            arr, transform = _read_burn_date(g)
            if src_transform is None:
                src_transform, src_shape = transform, arr.shape
                burn_count = np.zeros(src_shape, dtype=np.int16)
                last_burn_year = np.zeros(src_shape, dtype=np.int16)
            elif arr.shape != src_shape:
                raise RuntimeError(f"tile shape mismatch in {g}")
            month_burned = arr > 0  # 1-366 = burned; <=0 = unburned/water/unmapped
            burned_this_year = (
                month_burned if burned_this_year is None
                else (burned_this_year | month_burned)
            )
        burn_count += burned_this_year.astype(np.int16)
        last_burn_year[burned_this_year] = year
        log.info("  %d: %d burned pixels", year, int(burned_this_year.sum()))

    # ---- Derive the 0-1 fire-history score (in the native MODIS grid) -------
    ever_burned = burn_count > 0
    years_since = np.where(ever_burned, REF_YEAR - last_burn_year, 0).astype(float)
    recency = np.where(ever_burned, np.exp(-years_since / RECENCY_HALF_LIFE), 0.0)
    frequency = np.minimum(1.0, burn_count.astype(float) / FREQ_CAP)
    fire_src = (W_RECENCY * recency + W_FREQUENCY * frequency).astype(np.float32)

    # ---- Reproject onto the canal_risk grid ---------------------------------
    with rasterio.open(ref_grid) as ref:
        dst = np.zeros((ref.height, ref.width), dtype=np.float32)
        reproject(
            source=fire_src,
            destination=dst,
            src_transform=src_transform,
            src_crs=MODIS_SINU,
            dst_transform=ref.transform,
            dst_crs=ref.crs,
            resampling=Resampling.bilinear,
        )
        profile = ref.profile.copy()

    dst = np.clip(dst, 0.0, 1.0)
    # The COG driver auto-tiles internally and rejects the explicit
    # block/tiled/interleave options copied from the reference profile.
    for k in ("blockxsize", "blockysize", "tiled", "interleave"):
        profile.pop(k, None)
    profile.update(
        driver="COG", dtype="float32", count=1, nodata=None,
        compress="deflate",
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out_path, "w", **profile) as out:
        out.write(dst, 1)
    log.info("wrote %s  (%dx%d, EPSG:%s)",
             out_path, dst.shape[1], dst.shape[0], profile["crs"].to_epsg())

    # ---- Validation statistics ----------------------------------------------
    total = ever_burned.size
    burned = int(ever_burned.sum())
    yr_burned = {int(y): int((last_burn_year == y).sum()) for y in years}
    stats = {
        "years": (years[0], years[-1]),
        "modis_tile_pixels": total,
        "pixels_burned_at_least_once": burned,
        "pct_tile_burned": round(100.0 * burned / total, 2),
        "last_burn_year_pixel_counts": yr_burned,
        "max_burn_count": int(burn_count.max()),
        "aoi_fire_history_mean": round(float(dst.mean()), 4),
        "aoi_fire_history_max": round(float(dst.max()), 4),
        "aoi_pct_nonzero": round(100.0 * float((dst > 0.01).mean()), 2),
    }
    return stats


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--granules", type=Path, required=True)
    p.add_argument("--ref-grid", type=Path, required=True,
                   help="raster whose grid the output is aligned to "
                        "(canal_risk.tif)")
    p.add_argument("--out", type=Path, required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    stats = build_fire_history(args.granules, args.ref_grid, args.out)
    log.info("")
    log.info("=== VALIDATION STATS ===")
    for k, v in stats.items():
        log.info("  %s: %s", k, v)


if __name__ == "__main__":
    main()
