"""ECOSTRESS thermal-stress risk layer (diurnal-range method).

Builds thermal_stress.tif (0-1) from a directory of ECOSTRESS
ECO_L2T_LSTE granules. The layer is the DIURNAL LAND-SURFACE-TEMPERATURE
RANGE -- day LST minus night LST -- which is a thermal-inertia proxy for
near-surface moisture.

Physics: water has a very high heat capacity, and wet ground also loses
heat through evapotranspiration. So a wet, intact peat surface barely
warms by day and barely cools by night -> SMALL diurnal range. Drained,
degrading peat has low heat capacity -> it heats hard by day and cools
hard by night -> LARGE diurnal range. The day-night difference is far
less confounded by land cover than single-time LST, because it measures
the AMPLITUDE of the temperature swing, which is governed by subsurface
heat capacity (i.e. moisture), not by how hot the surface is in absolute
terms. This is the classic apparent-thermal-inertia soil-moisture
approach, applied to ECOSTRESS day + night overpasses.

    diurnal       = median(day LST) - median(night LST)     per pixel
    thermal_stress = clip( (diurnal - LO) / (HI - LO), 0, 1 )

ECOSTRESS rides the ISS so overpass local-time drifts; over a multi-month
window each pixel is sampled by a spread of day and night times, and the
median over many overpasses is a stable diurnal-amplitude estimate.

ECO_L2T_LSTE is tiled COG GeoTIFF -> standard rasterio environment.

Usage:
    python -m peatguard.analysis.thermal_stress \
        --granules Datasets/ecostress_lste \
        --ref-grid .../DEMO_products_canal_risk.tif \
        --out      .../DEMO_products_thermal_stress.tif
"""

from __future__ import annotations

import argparse
import glob
import logging
import re
from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import Resampling, reproject

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("thermal_stress")

# ---- Tunables ---------------------------------------------------------------
# Indonesia is UTC+7. Day = local ~08:00-17:00 -> UTC 01:00-10:00.
# Night = local ~20:00-07:00 -> UTC 13:00-24:00 (and 00:00).
DAY_UTC = range(1, 11)
NIGHT_UTC = set(range(13, 24)) | {0}
# Physically plausible tropical-peat LST range (Kelvin); values outside
# are cloud edges / sensor artefacts.
LST_MIN_K, LST_MAX_K = 288.0, 323.0      # 15 C .. 50 C
# Diurnal-range normalisation (Kelvin). Set from the observed dry-season
# diurnal-range distribution over the Sebangau AOI (p5 ~3.4 K, p95 ~9.9 K):
# intact wet peat swamp forest sits near the LOW end, drained / bare
# degraded peat near the HIGH end. A humid tropical peatland has a
# naturally compressed diurnal range (humidity + haze damp the swing),
# and ECOSTRESS day passes drift across 09:00-15:00 local so the median
# does not capture the true peak -- both compress the measured range.
DIURNAL_LO_K = 3.5
DIURNAL_HI_K = 10.0
MIN_GRANULE_PIXELS = 2000
# -----------------------------------------------------------------------------


def _utc_hour(lst_path: str) -> int | None:
    m = re.search(r"_(\d{8})T(\d{2})\d{4}_", Path(lst_path).name)
    return int(m.group(2)) if m else None


def _read_lst(lst_path: str) -> tuple[np.ndarray, dict] | None:
    """Read one LST granule, screen to valid peat-surface temperatures
    in Kelvin. Returns (kelvin, geo) or None if unusable."""
    with rasterio.open(lst_path) as r:
        lst = r.read(1).astype(np.float32)
        geo = {"transform": r.transform, "crs": r.crs}
    valid = np.isfinite(lst) & (lst >= LST_MIN_K) & (lst <= LST_MAX_K)
    water_path = lst_path.replace("_LST.tif", "_water.tif")
    if Path(water_path).exists():
        with rasterio.open(water_path) as r:
            water = r.read(1)
        valid &= water == 0
    if valid.sum() < MIN_GRANULE_PIXELS:
        return None
    return np.where(valid, lst, np.nan).astype(np.float32), geo


def _composite(files: list[str], H: int, W: int, ref_transform, ref_crs) -> tuple[np.ndarray, np.ndarray, int]:
    """Reproject each granule onto the AOI grid; return (median, n_obs, used)."""
    stack: list[np.ndarray] = []
    used = 0
    for f in files:
        res = _read_lst(f)
        if res is None:
            continue
        kelvin, geo = res
        dst = np.full((H, W), np.nan, dtype=np.float32)
        reproject(
            source=kelvin, destination=dst,
            src_transform=geo["transform"], src_crs=geo["crs"],
            dst_transform=ref_transform, dst_crs=ref_crs,
            src_nodata=np.nan, dst_nodata=np.nan,
            resampling=Resampling.bilinear,
        )
        if np.isfinite(dst).sum() > 0:
            stack.append(dst)
            used += 1
    if used == 0:
        return np.full((H, W), np.nan, np.float32), np.zeros((H, W), int), 0
    cube = np.stack(stack, axis=0)
    n_obs = np.isfinite(cube).sum(axis=0)
    import warnings
    with warnings.catch_warnings():
        # All-NaN columns (pixels no granule saw) are expected -> nan, masked later.
        warnings.simplefilter("ignore", category=RuntimeWarning)
        med = np.nanmedian(cube, axis=0).astype(np.float32)
    return med, n_obs, used


def build_thermal_stress(granule_dir: Path, ref_grid: Path, out_path: Path) -> dict:
    """Build thermal_stress.tif. Returns validation statistics."""
    lst_files = sorted(glob.glob(str(granule_dir / "*_LST.tif")))
    if not lst_files:
        raise RuntimeError(f"no *_LST.tif granules in {granule_dir}")
    day = [f for f in lst_files if (_utc_hour(f) or -1) in DAY_UTC]
    night = [f for f in lst_files if (_utc_hour(f) or -1) in NIGHT_UTC]
    log.info("%d LST granules -> %d daytime, %d nighttime",
             len(lst_files), len(day), len(night))
    if not day or not night:
        raise RuntimeError("need both day and night overpasses for diurnal range")

    with rasterio.open(ref_grid) as ref:
        ref_profile = ref.profile.copy()
        H, W = ref.height, ref.width
        ref_transform, ref_crs = ref.transform, ref.crs

    day_med, day_obs, day_used = _composite(day, H, W, ref_transform, ref_crs)
    night_med, night_obs, night_used = _composite(night, H, W, ref_transform, ref_crs)
    log.info("composited %d daytime + %d nighttime granules onto the AOI grid",
             day_used, night_used)

    both = np.isfinite(day_med) & np.isfinite(night_med)
    if both.sum() < 100:
        raise RuntimeError("insufficient day+night overlap for diurnal range")
    diurnal = np.where(both, day_med - night_med, np.nan).astype(np.float32)
    # Diurnal range is physically positive; tiny negatives are noise.
    diurnal = np.where(both, np.clip(diurnal, 0.0, None), np.nan)
    log.info("diurnal range (K): p5 %.2f, p50 %.2f, p95 %.2f",
             *(float(np.nanpercentile(diurnal, p)) for p in (5, 50, 95)))

    thermal = np.clip(
        (diurnal - DIURNAL_LO_K) / (DIURNAL_HI_K - DIURNAL_LO_K), 0.0, 1.0
    ).astype(np.float32)
    thermal[~both] = -9999.0

    for k in ("blockxsize", "blockysize", "tiled", "interleave"):
        ref_profile.pop(k, None)
    ref_profile.update(driver="COG", dtype="float32", count=1,
                       nodata=-9999.0, compress="deflate")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out_path, "w", **ref_profile) as dst:
        dst.write(thermal, 1)
    log.info("wrote %s", out_path)

    valid = thermal != -9999.0
    vals = thermal[valid]
    stats = {
        "day_granules_used": day_used,
        "night_granules_used": night_used,
        "aoi_coverage_pct": round(100.0 * float(valid.mean()), 2),
        "median_diurnal_K": round(float(np.nanmedian(diurnal[both])), 2),
        "thermal_stress_mean": round(float(vals.mean()), 4),
        "thermal_stress_max": round(float(vals.max()), 4),
        "pct_high_gt_0.6": round(100.0 * float((vals > 0.6).mean()), 2),
    }
    return stats


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--granules", type=Path, required=True)
    p.add_argument("--ref-grid", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    stats = build_thermal_stress(args.granules, args.ref_grid, args.out)
    log.info("")
    log.info("=== THERMAL STRESS STATS ===")
    for k, v in stats.items():
        log.info("  %s: %s", k, v)


if __name__ == "__main__":
    main()
