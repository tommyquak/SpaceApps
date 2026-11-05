"""Prototype: ML gap-filling of sparse InSAR subsidence (Step 1).

Models the Umarhadi & Siegert 2026 approach (Environ. Res. Commun. 8 025025):
train a Random Forest to regress an InSAR-derived subsidence velocity on cheap
environmental predictor maps, then extrapolate it into the un-covered (here,
no-InSAR) gaps to produce a full-coverage subsidence surface.

This prototype uses the EXISTING Tay et al. 2025 (ALOS-2 L-band) subsidence
raster as the sparse label (~1.9% coverage) and predictors PeatGuard already
holds. It processes ZERO Sentinel-1 scenes -- it only proves the predictor
stack + validation harness. Sentinel-1 SBAS would later add an independent,
current, measured label over degraded peat.

LABEL LEAKAGE GUARD: restoration_zones / canal_risk / canal_risk_confidence and
the El Nino drawdown are all derived from the same InSAR velocity, so they are
deliberately EXCLUDED as predictors. Only independent environmental proxies are
used. Peat depth + distance-to-peat-edge (the paper's top predictors) are NOT
available here and are noted as the key missing inputs.
"""
from __future__ import annotations

import json
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
OUT = Path("deliverables/subsidence_ml_proto")
OUT.mkdir(parents=True, exist_ok=True)
LABEL = BASE / "subsidence_tay.tif"
CANAL = BASE / "canal_mask.tif"
FIRE = BASE / "fire_history.tif"
WATERWAYS = Path("grid_products/external/osm_waterways_sebangau.geojson")
SETTLE = Path("grid_products/external/osm_settlements_sebangau.geojson")
PEAT = Path("grid_products/external/peatgrids/PEATGRIDS_THICKNESS_CM_MEAN.tif")

PX_M = 157.0  # approx metres per 0.001416 deg pixel at ~2.5S


def load_label():
    with rasterio.open(LABEL) as r:
        a = r.read(1).astype("float32")
        prof = r.profile
        nd = r.nodata
        valid = np.isfinite(a)
        if nd is not None:
            valid &= a != nd
        return a, valid, prof, r.transform, r.crs, (r.width, r.height)


def warp_to(src_path, transform, crs, shape, resampling):
    with rasterio.open(src_path) as r:
        src = r.read(1).astype("float32")
        src_nd = r.nodata
        if src_nd is not None:
            src = np.where(src == src_nd, 0.0, src)
        dst = np.zeros((shape[1], shape[0]), "float32")
        reproject(src, dst, src_transform=r.transform, src_crs=r.crs,
                  dst_transform=transform, dst_crs=crs, resampling=resampling)
    return dst


def warp_window_to(src_path, transform, crs, shape, resampling, bounds, margin=0.3):
    # For a large (global) source: read only the window over the AOI, then warp.
    W, S, E, N = bounds
    with rasterio.open(src_path) as r:
        win = from_bounds(W - margin, S - margin, E + margin, N + margin, r.transform)
        win = win.round_offsets().round_lengths()
        src = r.read(1, window=win, boundless=True, fill_value=0).astype("float32")
        wt = r.window_transform(win)
        dst = np.zeros((shape[1], shape[0]), "float32")
        reproject(src, dst, src_transform=wt, src_crs=r.crs,
                  dst_transform=transform, dst_crs=crs, resampling=resampling)
    return dst


def rasterize_vec(path, transform, shape, geom_filter=None):
    gj = json.loads(Path(path).read_text())
    geoms = []
    for f in gj.get("features", []):
        g = f.get("geometry")
        if not g:
            continue
        if geom_filter and g.get("type") not in geom_filter:
            continue
        try:
            geoms.append(shape(g))
        except Exception:
            continue
    if not geoms:
        return np.zeros((shape[1], shape[0]), "uint8")
    return features.rasterize(((g, 1) for g in geoms), out_shape=(shape[1], shape[0]),
                              transform=transform, fill=0, all_touched=True, dtype="uint8")


def edt_m(binary):
    # distance (m) to nearest True cell
    if binary.sum() == 0:
        return np.full(binary.shape, 1e5, "float32")
    return (distance_transform_edt(~binary.astype(bool)) * PX_M).astype("float32")


def main():
    print("== loading label (Tay subsidence) ==")
    y_full, valid, prof, transform, crs, (W, H) = load_label()
    yv = y_full[valid]
    lo, hi = np.percentile(yv, [2, 98])
    print(f"  grid {W}x{H} res~{PX_M:.0f}m  covered={valid.sum()} ({100*valid.mean():.2f}%)")
    print(f"  velocity mm/yr: raw range [{yv.min():.1f},{yv.max():.1f}] clip to [{lo:.1f},{hi:.1f}] (p2-p98)")
    y_full = np.clip(y_full, lo, hi)

    print("== building predictor stack on the label grid ==")
    bounds = array_bounds(H, W, transform)  # (W,S,E,N)
    canal_density = warp_to(CANAL, transform, crs, (W, H), Resampling.average)  # 0..1 fraction
    canal_present = canal_density > 0.02
    dist_canal = edt_m(canal_present)
    fire = warp_to(FIRE, transform, crs, (W, H), Resampling.bilinear)
    riv = rasterize_vec(WATERWAYS, transform, (W, H), {"LineString", "MultiLineString"})
    dist_river = edt_m(riv.astype(bool))
    setl = rasterize_vec(SETTLE, transform, (W, H), {"Point", "MultiPoint"})
    dist_settle = edt_m(setl.astype(bool))

    # PEATGRIDS peat thickness (cm) -> peat depth + distance-to-peat-edge
    peat_depth = warp_window_to(PEAT, transform, crs, (W, H), Resampling.bilinear, bounds)
    peat_mask = peat_depth > 25.0  # cm; PEATGRIDS predicts over the peat extent
    # distance INTO the peat from its edge (0 outside peat)
    dist_peat_edge = (distance_transform_edt(peat_mask) * PX_M).astype("float32")
    print(f"  peat depth cm: range [{peat_depth[peat_mask].min():.0f},{peat_depth[peat_mask].max():.0f}] "
          f"peat-covered grid {100*peat_mask.mean():.0f}%; dist-to-edge max {dist_peat_edge.max()/1000:.1f} km")

    feats = {
        "peat_depth_cm": peat_depth,
        "dist_to_peat_edge_m": dist_peat_edge,
        "canal_density": canal_density,
        "dist_to_canal_m": dist_canal,
        "fire_history": fire,
        "dist_to_river_m": dist_river,
        "dist_to_settle_m": dist_settle,
    }
    names = list(feats)
    X_full = np.stack([feats[n].ravel() for n in names], axis=1)
    print("  predictors:", names)

    # training rows = covered label pixels
    flat_valid = valid.ravel()
    X = X_full[flat_valid]
    y = y_full.ravel()[flat_valid]
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.30, random_state=42)

    print(f"== training RandomForestRegressor on {len(ytr)} px (test {len(yte)}) ==")
    rf = RandomForestRegressor(n_estimators=300, max_features="sqrt",
                               min_samples_leaf=4, n_jobs=-1, random_state=42)
    rf.fit(Xtr, ytr)
    pred = rf.predict(Xte)
    r2 = r2_score(yte, pred)
    mae = mean_absolute_error(yte, pred)
    rmse = float(np.sqrt(mean_squared_error(yte, pred)))
    cv = cross_val_score(rf, X, y, cv=5, scoring="r2", n_jobs=-1)
    print(f"  test  R2={r2:.3f}  MAE={mae:.2f}  RMSE={rmse:.2f} mm/yr")
    print(f"  5-fold CV R2: {cv.mean():.3f} +/- {cv.std():.3f}")

    imp = rf.feature_importances_
    perm = permutation_importance(rf, Xte, yte, n_repeats=8, random_state=42, n_jobs=-1)
    print("  feature importance (impurity | permutation):")
    order = np.argsort(imp)[::-1]
    fi = []
    for i in order:
        print(f"    {names[i]:18s} {imp[i]:.3f} | {perm.importances_mean[i]:.3f}")
        fi.append({"feature": names[i], "impurity": float(imp[i]),
                   "permutation": float(perm.importances_mean[i])})

    print("== extrapolating to full grid ==")
    full = rf.predict(X_full).reshape(H, W).astype("float32")

    out_tif = OUT / "subsidence_ml_full.tif"
    op = prof.copy()
    op.update(dtype="float32", count=1, nodata=-9999.0, compress="deflate")
    with rasterio.open(out_tif, "w", **op) as d:
        d.write(full, 1)
    print(f"  wrote {out_tif}  coverage 1.9% -> 100%")

    metrics = {
        "label": "Tay et al. 2025 ALOS-2 L-band subsidence (mm/yr)",
        "method": "Umarhadi & Siegert 2026 ML gap-filling (RandomForest)",
        "label_coverage_pct": round(100 * float(valid.mean()), 2),
        "n_train": int(len(ytr)), "n_test": int(len(yte)),
        "test_r2": round(float(r2), 3), "test_mae_mm_yr": round(float(mae), 2),
        "test_rmse_mm_yr": round(rmse, 2),
        "cv5_r2_mean": round(float(cv.mean()), 3), "cv5_r2_std": round(float(cv.std()), 3),
        "clip_mm_yr": [round(float(lo), 1), round(float(hi), 1)],
        "predictors": names,
        "peat_depth_source": "PEATGRIDS (Widyastuti et al. 2024) ~1km global thickness MEAN, cm",
        "feature_importance": fi,
        "caveats": [
            "Held-out test pixels lie WITHIN Tay coverage -> measures interpolation skill, an optimistic proxy for true no-coverage extrapolation (no independent ground truth in the gaps).",
            "Peat depth is ~1km PEATGRIDS (global DSM), coarser and shallower than LiDAR ground truth; captures dome-scale structure, not fine detail.",
            "Random split ignores spatial autocorrelation -> R2 likely optimistic; a spatially-blocked CV would be stricter.",
            "C-band/L-band subsidence both underestimate; values are relative, not calibrated.",
        ],
    }
    (OUT / "metrics.json").write_text(json.dumps(metrics, indent=2))

    # ---- preview ----
    fig, ax = plt.subplots(1, 3, figsize=(16, 6))
    meas = np.where(valid, y_full, np.nan)
    vmin, vmax = lo, hi
    ax[0].imshow(meas, cmap="RdYlGn", vmin=vmin, vmax=vmax)
    ax[0].set_title(f"Measured Tay L-band ({100*valid.mean():.1f}% coverage)")
    ax[1].imshow(full, cmap="RdYlGn", vmin=vmin, vmax=vmax)
    ax[1].set_title("ML extrapolation (100% coverage)")
    for a in ax[:2]:
        a.set_xticks([]); a.set_yticks([])
    ax[2].scatter(yte, pred, s=3, alpha=0.2, color="#1f4d3a")
    lim = [vmin, vmax]
    ax[2].plot(lim, lim, "--", color="#888", lw=1)
    ax[2].set_xlim(lim); ax[2].set_ylim(lim)
    ax[2].set_xlabel("measured (mm/yr)"); ax[2].set_ylabel("predicted (mm/yr)")
    ax[2].set_title(f"Held-out test: R2={r2:.2f} RMSE={rmse:.1f}")
    fig.suptitle("Subsidence ML gap-filling prototype (Tay L-band label)", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT / "preview.png", dpi=130, bbox_inches="tight")
    print(f"  wrote {OUT/'preview.png'} and {OUT/'metrics.json'}")
    print("DONE")


if __name__ == "__main__":
    main()
