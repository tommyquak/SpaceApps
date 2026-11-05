"""S1-SBAS + ML subsidence with independent Tay validation (the full pipeline).

Trains a Random Forest on the Sentinel-1 SBAS velocity (the measured label),
extrapolates to full coverage from environmental predictors, and validates
against the Tay et al. 2025 (ALOS-2 L-band) raster as an INDEPENDENT benchmark
-- WITHOUT force-calibrating to it (per the no-force-calibration rule). Models
Umarhadi & Siegert 2026.

Run after MintPy produces a geocoded velocity GeoTIFF:
    subsidence_ml_s1.py --velocity /path/to/s1_velocity.tif

LEAKAGE GUARD: predictors are independent environmental layers only. Tay (and
anything derived from it -- restoration_zones, canal_risk, El Nino) is never a
predictor; it is held out purely for validation. Valid because the LABEL here is
S1, a different sensor/time from Tay.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.windows import from_bounds
from rasterio.transform import array_bounds
from rasterio import features
from scipy.ndimage import distance_transform_edt
from shapely.geometry import shape
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.inspection import permutation_importance
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("deliverables/tay_restoration")
CANAL = BASE / "canal_mask.tif"
FIRE = BASE / "fire_history.tif"
TAY = BASE / "subsidence_tay.tif"
WATERWAYS = Path("grid_products/external/osm_waterways_sebangau.geojson")
SETTLE = Path("grid_products/external/osm_settlements_sebangau.geojson")
PEAT = Path("grid_products/external/peatgrids/PEATGRIDS_THICKNESS_CM_MEAN.tif")
OUT = Path("deliverables/subsidence_s1_ml")
OUT.mkdir(parents=True, exist_ok=True)


def warp_val_mask(src_path, transform, crs, shape, resampling=Resampling.bilinear, big=False, bounds=None):
    """Warp a raster to the target grid, returning (values, valid_mask)."""
    with rasterio.open(src_path) as r:
        if big and bounds is not None:
            W, S, E, N = bounds
            win = from_bounds(W - 0.3, S - 0.3, E + 0.3, N + 0.3, r.transform).round_offsets().round_lengths()
            src = r.read(1, window=win, boundless=True, fill_value=(r.nodata or 0)).astype("float32")
            st = r.window_transform(win)
        else:
            src = r.read(1).astype("float32")
            st = r.transform
        snd = r.nodata
        m = np.isfinite(src)
        if snd is not None:
            m &= src != snd
        val = np.zeros((shape[1], shape[0]), "float32")
        reproject(np.where(m, src, 0.0), val, src_transform=st, src_crs=r.crs,
                  dst_transform=transform, dst_crs=crs, resampling=resampling)
        mk = np.zeros((shape[1], shape[0]), "float32")
        reproject(m.astype("float32"), mk, src_transform=st, src_crs=r.crs,
                  dst_transform=transform, dst_crs=crs, resampling=Resampling.nearest)
    return val, mk > 0.5


def rasterize_vec(path, transform, shape, gtypes):
    gj = json.loads(Path(path).read_text())
    geoms = []
    for f in gj.get("features", []):
        g = f.get("geometry")
        if g and g.get("type") in gtypes:
            try:
                geoms.append(shape(g))
            except Exception:
                pass
    if not geoms:
        return np.zeros((shape[1], shape[0]), "uint8")
    return features.rasterize(((g, 1) for g in geoms), out_shape=(shape[1], shape[0]),
                              transform=transform, fill=0, all_touched=True, dtype="uint8")


def edt_m(binary, px_m):
    if binary.sum() == 0:
        return np.full(binary.shape, 1e5, "float32")
    return (distance_transform_edt(~binary.astype(bool)) * px_m).astype("float32")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--velocity", required=True, help="geocoded S1 SBAS velocity GeoTIFF")
    ap.add_argument("--units", default="auto", choices=["auto", "m", "mm"], help="velocity units")
    args = ap.parse_args()

    print("== loading S1 velocity (label) ==")
    with rasterio.open(args.velocity) as r:
        vel = r.read(1).astype("float32")
        transform, crs = r.transform, r.crs
        W, H = r.width, r.height
        nd = r.nodata
    valid = np.isfinite(vel)
    if nd is not None:
        valid &= vel != nd
    valid &= vel != 0.0  # MintPy masks as 0/NaN
    bounds = array_bounds(H, W, transform)
    lat0 = (bounds[1] + bounds[3]) / 2.0
    px_m = abs(transform.a) * 111320.0 * math.cos(math.radians(lat0))
    print(f"  grid {W}x{H} crs={crs} px~{px_m:.0f}m  S1-covered={valid.sum()} ({100*valid.mean():.1f}%)")

    # units -> mm/yr
    med = float(np.median(np.abs(vel[valid]))) if valid.any() else 0.0
    units = args.units if args.units != "auto" else ("m" if med < 0.5 else "mm")
    scale = 1000.0 if units == "m" else 1.0
    vel_mm = vel * scale
    yv = vel_mm[valid]
    lo, hi = np.percentile(yv, [2, 98])
    print(f"  units={units} -> mm/yr; range [{yv.min():.1f},{yv.max():.1f}] clip [{lo:.1f},{hi:.1f}]")
    vel_mm = np.clip(vel_mm, lo, hi)

    print("== predictor stack (independent of subsidence) ==")
    peat_depth, _ = warp_val_mask(PEAT, transform, crs, (W, H), Resampling.bilinear, big=True, bounds=bounds)
    peat_mask = peat_depth > 25.0
    dist_peat_edge = (distance_transform_edt(peat_mask) * px_m).astype("float32")
    canal_density, _ = warp_val_mask(CANAL, transform, crs, (W, H), Resampling.average)
    dist_canal = edt_m(canal_density > 0.02, px_m)
    fire, _ = warp_val_mask(FIRE, transform, crs, (W, H), Resampling.bilinear)
    dist_river = edt_m(rasterize_vec(WATERWAYS, transform, (W, H), {"LineString", "MultiLineString"}).astype(bool), px_m)
    dist_settle = edt_m(rasterize_vec(SETTLE, transform, (W, H), {"Point", "MultiPoint"}).astype(bool), px_m)
    feats = {"peat_depth_cm": peat_depth, "dist_to_peat_edge_m": dist_peat_edge,
             "canal_density": canal_density, "dist_to_canal_m": dist_canal,
             "fire_history": fire, "dist_to_river_m": dist_river, "dist_to_settle_m": dist_settle}
    # multi-temporal NBR (Landsat C2 dry-season composites; the paper's NBR predictor)
    nbr_dir = Path("grid_products/external/nbr")
    for yr in [1990, 1995, 2000, 2005, 2010, 2015, 2020, 2022]:
        p = nbr_dir / f"nbr_{yr}.tif"
        if not p.exists():
            continue
        arr, m = warp_val_mask(p, transform, crs, (W, H), Resampling.bilinear)
        fill = float(np.median(arr[m])) if m.any() else 0.0
        feats[f"nbr_{yr}"] = np.where(m, arr, fill).astype("float32")
    print(f"  NBR epochs loaded: {sum(1 for k in feats if k.startswith('nbr_'))}")
    names = list(feats)
    X_full = np.stack([feats[n].ravel() for n in names], axis=1)

    fv = valid.ravel()
    X, y = X_full[fv], vel_mm.ravel()[fv]
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.30, random_state=42)
    print(f"== RandomForest on {len(ytr)} S1 px (test {len(yte)}) ==")
    rf = RandomForestRegressor(n_estimators=400, max_features="sqrt", min_samples_leaf=4,
                               n_jobs=-1, random_state=42)
    rf.fit(Xtr, ytr)
    pred = rf.predict(Xte)
    r2 = r2_score(yte, pred); mae = mean_absolute_error(yte, pred)
    rmse = float(np.sqrt(mean_squared_error(yte, pred)))
    cv = cross_val_score(rf, X, y, cv=5, scoring="r2", n_jobs=-1)
    print(f"  test R2={r2:.3f} MAE={mae:.2f} RMSE={rmse:.2f} mm/yr | blocked-CV R2 {cv.mean():.3f}+/-{cv.std():.3f}")
    imp = rf.feature_importances_
    perm = permutation_importance(rf, Xte, yte, n_repeats=8, random_state=42, n_jobs=-1)
    fi = [{"feature": names[i], "impurity": float(imp[i]), "permutation": float(perm.importances_mean[i])}
          for i in np.argsort(imp)[::-1]]
    for f in fi:
        print(f"    {f['feature']:18s} {f['impurity']:.3f} | {f['permutation']:.3f}")

    print("== extrapolate to full grid ==")
    full = rf.predict(X_full).reshape(H, W).astype("float32")

    # ---- INDEPENDENT validation vs Tay (no calibration) ----
    print("== validate vs Tay (independent L-band benchmark) ==")
    tay_mm, tay_valid = warp_val_mask(TAY, transform, crs, (W, H), Resampling.bilinear)  # Tay already mm/yr
    # robust sign check on the S1-vs-Tay overlap
    ov = valid & tay_valid
    val_report = {"note": "Tay is an independent L-band InSAR estimate, not ground truth; reported as a benchmark. No calibration applied."}
    if ov.sum() > 50:
        s1o, tyo = vel_mm[ov], tay_mm[ov]
        corr = float(np.corrcoef(s1o, tyo)[0, 1])
        val_report["sign_correlation_s1_vs_tay"] = round(corr, 3)
        val_report["test1_overlap"] = {
            "n": int(ov.sum()),
            "r2": round(float(r2_score(tyo, s1o)), 3),
            "mae_mm_yr": round(float(mean_absolute_error(tyo, s1o)), 2),
            "rmse_mm_yr": round(float(np.sqrt(mean_squared_error(tyo, s1o))), 2),
            "mean_bias_s1_minus_tay": round(float(np.mean(s1o - tyo)), 2),
        }
        print(f"  test1 S1-measured vs Tay (overlap n={ov.sum()}): corr={corr:.2f} "
              f"r2={val_report['test1_overlap']['r2']} rmse={val_report['test1_overlap']['rmse_mm_yr']} "
              f"bias={val_report['test1_overlap']['mean_bias_s1_minus_tay']}")
    # gap test: where S1 has no data but Tay does -> ML extrapolation vs Tay (non-circular)
    gap = (~valid) & tay_valid
    if gap.sum() > 50:
        pg, tg = full[gap], tay_mm[gap]
        val_report["test2_gapfill_vs_tay"] = {
            "n": int(gap.sum()),
            "r2": round(float(r2_score(tg, pg)), 3),
            "mae_mm_yr": round(float(mean_absolute_error(tg, pg)), 2),
            "rmse_mm_yr": round(float(np.sqrt(mean_squared_error(tg, pg))), 2),
            "mean_bias_pred_minus_tay": round(float(np.mean(pg - tg)), 2),
        }
        print(f"  test2 ML-extrapolation vs Tay (S1-gap n={gap.sum()}): "
              f"r2={val_report['test2_gapfill_vs_tay']['r2']} rmse={val_report['test2_gapfill_vs_tay']['rmse_mm_yr']}")
    else:
        print("  test2: too few S1-gap & Tay-covered pixels for the gap-fill test")

    # ---- outputs ----
    op = {"driver": "GTiff", "dtype": "float32", "count": 1, "height": H, "width": W,
          "crs": crs, "transform": transform, "nodata": -9999.0, "compress": "deflate"}
    with rasterio.open(OUT / "s1_ml_subsidence_full.tif", "w", **op) as d:
        d.write(full, 1)
    metrics = {
        "label": "Sentinel-1 C-band SBAS velocity (this run)", "units_in": units,
        "method": "Umarhadi & Siegert 2026 (RF gap-fill)",
        "s1_coverage_pct": round(100 * float(valid.mean()), 2),
        "test_r2": round(float(r2), 3), "test_mae_mm_yr": round(mae, 2), "test_rmse_mm_yr": round(rmse, 2),
        "cv5_r2_mean": round(float(cv.mean()), 3), "feature_importance": fi,
        "tay_validation": val_report,
    }
    (OUT / "metrics.json").write_text(json.dumps(metrics, indent=2))

    fig, ax = plt.subplots(1, 3, figsize=(16, 6))
    vmin, vmax = lo, hi
    ax[0].imshow(np.where(valid, vel_mm, np.nan), cmap="RdYlGn", vmin=vmin, vmax=vmax)
    ax[0].set_title(f"S1 SBAS measured ({100*valid.mean():.0f}%)")
    ax[1].imshow(full, cmap="RdYlGn", vmin=vmin, vmax=vmax)
    ax[1].set_title("S1 + ML full coverage")
    ax[2].imshow(np.where(tay_valid, tay_mm, np.nan), cmap="RdYlGn", vmin=vmin, vmax=vmax)
    ax[2].set_title("Tay L-band (validation)")
    for a in ax:
        a.set_xticks([]); a.set_yticks([])
    fig.tight_layout()
    fig.savefig(OUT / "preview.png", dpi=130, bbox_inches="tight")
    print(f"  wrote {OUT}/s1_ml_subsidence_full.tif, metrics.json, preview.png\nDONE")


if __name__ == "__main__":
    main()
