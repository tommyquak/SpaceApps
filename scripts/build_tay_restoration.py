"""Fuse Tay et al. (2025) L-band subsidence + canal proximity + fire into
PeatGuard restoration-priority products for the shifted 2x2 grid.

Subsidence source: Tay, C. et al. (2025), Communications Earth & Environment,
"Satellite radar advances carbon emissions accountability over tropical peat"
-- their published ALOS-2 L-band long-term-rate raster (DR-NTU repository).
This replaces the C-band InSAR velocity, which decorrelated (coherence ~0.17)
and was independently shown to be noise-dominated over the operating area.
The Tay raster remains an independent published input; no fit or calibration
to Tay is applied.

Outputs (EPSG:4326 COGs on the canal-mask grid):
    canal_risk.tif            0.45*proximity + 0.55*subsidence  (proximity-only
                              where Tay has no data; see confidence band)
    canal_risk_confidence.tif 2 = subsidence-backed, 1 = proximity-only
    restoration_zones.tif     0.75*canal_risk + 0.25*fire
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject
from scipy.ndimage import distance_transform_edt

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("tay_restoration")

PROJECT_MAX_INFLUENCE_M = 1200.0  # sensitivity prior; field calibration pending
W_PROX, W_SUBS = 0.45, 0.55
W_RISK, W_FIRE = 0.75, 0.25
NODATA = -9999.0


def _read(path: Path):
    with rasterio.open(path) as r:
        return r.read(1).astype("float64"), r.profile, r.nodata


def _to_grid(src_path: Path, dst_transform, dst_crs, dst_shape):
    """Reproject a raster onto the destination grid."""
    with rasterio.open(src_path) as r:
        out = np.full(dst_shape, np.nan, dtype="float64")
        reproject(
            source=rasterio.band(r, 1),
            destination=out,
            src_transform=r.transform, src_crs=r.crs,
            dst_transform=dst_transform, dst_crs=dst_crs,
            resampling=Resampling.bilinear,
            src_nodata=r.nodata, dst_nodata=np.nan,
        )
        return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--canal", required=True, type=Path)
    p.add_argument("--subsidence", required=True, type=Path, help="Tay L-band rate, mm/yr")
    p.add_argument("--fire", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    args = p.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    # --- canal mask defines the working grid -------------------------------
    canal, cprof, cnd = _read(args.canal)
    canal_bin = (canal > 0) & np.isfinite(canal)
    H, W = canal.shape
    px_m = abs(cprof["transform"].a) * 111_320.0 * np.cos(np.radians(
        cprof["transform"].f + (H / 2) * cprof["transform"].e))
    log.info("grid: %dx%d  pixel ~%.0f m  canal pixels=%d", H, W, px_m, int(canal_bin.sum()))

    # --- canal proximity (project sensitivity prior) ------------------------
    dist_m = distance_transform_edt(~canal_bin) * px_m
    proximity = np.clip(1.0 - dist_m / PROJECT_MAX_INFLUENCE_M, 0.0, 1.0)

    # --- Tay subsidence -> 0..1 risk ---------------------------------------
    subs = _to_grid(args.subsidence, cprof["transform"], cprof["crs"], (H, W))
    subsiding = np.isfinite(subs) & (subs < 0)
    # scale: strongly-subsiding tail (5th pct of subsiding rates) -> risk 1.0
    if subsiding.any():
        sev = abs(np.percentile(subs[subsiding], 5))
        sev = max(sev, 5.0)
    else:
        sev = 15.0
    subs_risk = np.clip(-subs / sev, 0.0, 1.0)   # uplift/0 -> 0 risk
    has_subs = np.isfinite(subs)
    log.info("Tay subsidence: %.1f%% of grid has data, severe-scale=%.1f mm/yr",
             100 * has_subs.mean(), sev)

    # --- canal_risk: blended where Tay exists; NO-DATA elsewhere -----------
    # Restoration priority is only asserted where subsidence is actually
    # measured. Outside Tay coverage we have no subsidence -> explicit
    # no-data, NOT an assumed "low priority". The canal_mask layer still
    # carries the proximity/canal info for those areas.
    canal_risk = np.full((H, W), NODATA, dtype="float64")
    canal_risk[has_subs] = (W_PROX * proximity[has_subs]
                            + W_SUBS * subs_risk[has_subs])
    confidence = np.where(has_subs, 2, 0).astype("float32")  # 2=measured, 0=no data

    # --- fire --------------------------------------------------------------
    fire = _to_grid(args.fire, cprof["transform"], cprof["crs"], (H, W))
    fire = np.where(np.isfinite(fire), fire, 0.0)

    # --- restoration zones (explicit no-data outside subsidence coverage) --
    restoration = np.full((H, W), NODATA, dtype="float64")
    restoration[has_subs] = (W_RISK * canal_risk[has_subs]
                             + W_FIRE * np.clip(fire[has_subs], 0, 1))

    # --- write COGs --------------------------------------------------------
    base = dict(driver="COG", dtype="float32", count=1, height=H, width=W,
                crs=cprof["crs"], transform=cprof["transform"],
                nodata=NODATA, compress="deflate", overview_resampling="nearest")
    for name, arr in [("canal_risk", canal_risk),
                      ("canal_risk_confidence", confidence),
                      ("restoration_zones", restoration)]:
        out = arr.astype("float32")
        path = args.outdir / f"{name}.tif"
        with rasterio.open(path, "w", **base) as dst:
            dst.write(out, 1)
        v = out[np.isfinite(out)]
        log.info("wrote %s  range[%.2f,%.2f] mean=%.3f",
                 path.name, float(v.min()), float(v.max()), float(v.mean()))

    n_sub = int(has_subs.sum())
    log.info("confidence: %.1f%% subsidence-backed, %.1f%% proximity-only",
             100 * n_sub / canal.size, 100 * (canal.size - n_sub) / canal.size)


if __name__ == "__main__":
    main()
