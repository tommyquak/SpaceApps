"""Unified restoration-zones layer.

Fuses the finished PeatGuard risk products into ONE 0-1 surface -- the
single layer an operator reads (beneath canal_mask on the map) to decide
where peatland restoration should be prioritised.

Additive weighted model, with per-pixel renormalisation over the factors
that actually have data at that pixel:

    restoration = sum(Wi * factor_i) / sum(Wi)      for present factors i

The weights trace back to the 4-factor degradation model agreed for
PeatGuard:

    subsidence velocity    0.35  ]
    canal proximity        0.25  ]-> canal_risk already combines these two
    thermal dryness        0.25      (Hooijer-calibrated), so it enters the
    fire history           0.15      fusion as one 0.60-weight term.

canal_risk is the InSAR structural-degradation product and fire_history
(MODIS MCD64A1) is an independent, strongly-validated axis (its burn-year
curve reproduces the documented ENSO fire record). Per-pixel
renormalisation means a pixel still gets a meaningful score where
canal_risk has no InSAR coverage but fire does -- e.g. a burn scar far
from any mapped canal is still flagged for restoration.

thermal_stress (ECOSTRESS diurnal range) is accepted as an OPTIONAL
input but is held OUT of the default fusion: over the Sebangau AOI it
cross-validated only weakly (~0.10 correlation with fire / canal_risk)
and there is no ground truth to anchor it, so it ships as a diagnostic
map layer rather than a headline-risk driver. Pass --thermal-stress
explicitly to fuse it (weight 0.25) once stronger validation exists.

Runs on finished COGs only (no HDF) -> standard rasterio environment.

Usage:
    python -m peatguard.analysis.restoration_zones \
        --canal-risk   .../DEMO_products_canal_risk.tif \
        --fire-history .../DEMO_products_fire_history.tif \
        [--thermal-stress .../DEMO_products_thermal_stress.tif] \
        --out          .../DEMO_products_restoration_zones.tif
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import Resampling, reproject

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("restoration_zones")

# Factor weights in the additive 4-factor model. canal_risk carries the
# subsidence (0.35) + proximity (0.25) share as one combined term.
WEIGHTS = {
    "canal_risk": 0.60,
    "thermal_stress": 0.25,
    "fire_history": 0.15,
}


def _load_aligned(path: Path, ref) -> tuple[np.ndarray, np.ndarray]:
    """Load a raster, reprojecting onto the reference grid if needed.

    Returns (values, valid_mask). valid_mask is False where the source
    raster has nodata / non-finite values."""
    with rasterio.open(path) as src:
        arr = src.read(1).astype(np.float32)
        nod = src.nodata
        if (src.width, src.height, src.transform, src.crs) == (
            ref.width, ref.height, ref.transform, ref.crs
        ):
            data = arr
        else:
            data = np.full((ref.height, ref.width), np.nan, dtype=np.float32)
            reproject(
                source=arr, destination=data,
                src_transform=src.transform, src_crs=src.crs,
                dst_transform=ref.transform, dst_crs=ref.crs,
                src_nodata=nod, dst_nodata=np.nan,
                resampling=Resampling.bilinear,
            )
            nod = None  # nodata folded into NaN
    valid = np.isfinite(data)
    if nod is not None:
        valid &= data != nod
    return data, valid


def build_restoration_zones(
    canal_risk_path: Path,
    fire_history_path: Path,
    out_path: Path,
    thermal_stress_path: Path | None = None,
) -> dict:
    """Build restoration_zones.tif. Returns validation statistics."""
    # canal_risk defines the output grid + profile.
    with rasterio.open(canal_risk_path) as ref:
        ref_profile = ref.profile.copy()
        cr, cr_valid = _load_aligned(canal_risk_path, ref)
        fire, fire_valid = _load_aligned(fire_history_path, ref)
        layers = {"canal_risk": (cr, cr_valid), "fire_history": (fire, fire_valid)}
        if thermal_stress_path is not None and Path(thermal_stress_path).exists():
            th, th_valid = _load_aligned(thermal_stress_path, ref)
            layers["thermal_stress"] = (th, th_valid)

    present = list(layers)
    log.info("fusing %d factor(s): %s", len(present), ", ".join(present))

    h, w = cr.shape
    numerator = np.zeros((h, w), dtype=np.float64)
    denominator = np.zeros((h, w), dtype=np.float64)
    for name, (data, valid) in layers.items():
        wgt = WEIGHTS[name]
        clean = np.where(valid, np.clip(data, 0.0, 1.0), 0.0)
        numerator += wgt * clean * valid
        denominator += wgt * valid
        log.info("  %-15s weight=%.2f  valid=%.1f%%",
                 name, wgt, 100.0 * valid.mean())

    has_any = denominator > 0
    restoration = np.full((h, w), -9999.0, dtype=np.float32)
    restoration[has_any] = (numerator[has_any] / denominator[has_any]).astype(np.float32)
    restoration = np.where(
        restoration == -9999.0, -9999.0, np.clip(restoration, 0.0, 1.0)
    ).astype(np.float32)

    for k in ("blockxsize", "blockysize", "tiled", "interleave"):
        ref_profile.pop(k, None)
    ref_profile.update(driver="COG", dtype="float32", count=1,
                       nodata=-9999.0, compress="deflate")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out_path, "w", **ref_profile) as dst:
        dst.write(restoration, 1)

    valid = restoration != -9999.0
    vals = restoration[valid]
    stats = {
        "factors": present,
        "weights": {k: WEIGHTS[k] for k in present},
        "output": str(out_path),
        "grid": f"{w}x{h}",
        "valid_pct": round(100.0 * float(valid.mean()), 2),
        "mean": round(float(vals.mean()), 4),
        "max": round(float(vals.max()), 4),
        "pct_high_gt_0.6": round(100.0 * float((vals > 0.6).mean()), 2),
        "pct_critical_gt_0.8": round(100.0 * float((vals > 0.8).mean()), 2),
    }
    log.info("wrote %s", out_path)
    return stats


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--canal-risk", type=Path, required=True)
    p.add_argument("--fire-history", type=Path, required=True)
    p.add_argument("--thermal-stress", type=Path, default=None)
    p.add_argument("--out", type=Path, required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    stats = build_restoration_zones(
        args.canal_risk, args.fire_history, args.out, args.thermal_stress,
    )
    log.info("")
    log.info("=== RESTORATION ZONES STATS ===")
    for k, v in stats.items():
        log.info("  %s: %s", k, v)


if __name__ == "__main__":
    main()
