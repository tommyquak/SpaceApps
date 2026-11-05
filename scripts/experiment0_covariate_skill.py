"""Experiment 0: out-of-block covariate signal test on two L-band subsidence labels.

Question: do the environmental covariates PeatGuard holds carry any spatial
information about L-band peat subsidence that survives honest cross-validation,
on (a) the external benchmark (Tay et al. 2025, ALOS-2, 2018 to 2024) and
(b) our own label (ALOS-1 PALSAR mosaic, 2007 to 2009)?

Models per label: M0 mean-only, M1 random forest on x,y only (the interpolation
baseline of Ploton et al. 2020), M2 RF on covariates, M3 RF on covariates + x,y.
Cross-validation per label: random 5-fold (the optimistic number), GroupKFold on
5 km and 20 km spatial blocks (the pessimistic bound), and for PALSAR-1 the
leave-one-orbital-path-out split (paths 422 vs 421 from alos1_frame_id.tif),
an independent-acquisition test. Rows of the PALSAR-1 label are weighted by
1/velocityStd^2. A pre-registered stratum test on the PALSAR-1 label reports
stratum means with a 1.8 km block bootstrap.

Predictors are the seven of the earlier prototype (subsidence_ml_extrapolate.py)
plus a disturbance-year class from the NBR series truncated to 2009 for the
PALSAR-1 label (no leakage from later years). Everything runs on the Tay grid
(~157 m). Output: deliverables/experiment0/metrics.json and a text summary.
Acceptance (from the 2026-09-08 research synthesis): covariate skill exists if
M2 block-CV R2 >= 0.20 and >= M1 block-CV R2 + 0.10; leave-one-path-out R2 >= 0.20.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import Resampling
from scipy.ndimage import distance_transform_edt
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold, KFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
from subsidence_ml_extrapolate import (CANAL, FIRE, LABEL, PEAT, PX_M, SETTLE, WATERWAYS, edt_m,  # noqa: E402
                                       load_label, rasterize_vec, warp_to, warp_window_to)

ALOS1 = Path("grid_products/external/alos1_2007_2009")
NBR = Path("grid_products/external/nbr")
OUT = Path("deliverables/experiment0")
OUT.mkdir(parents=True, exist_ok=True)
RF_KW = dict(n_estimators=300, max_features="sqrt", min_samples_leaf=4, n_jobs=-1, random_state=42)
BLOCK_M = (5000, 20000)
BOOT_BLOCK_PX = 12          # ~1.8 km at 157 m, the BLOCK convention of bw4_regional.py


def r2(y, p, w=None):
    w = np.ones_like(y) if w is None else w
    ybar = np.average(y, weights=w)
    ss_res = np.sum(w * (y - p) ** 2); ss_tot = np.sum(w * (y - ybar) ** 2)
    return float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")


def blocks(rows, cols, size_px):
    return (rows // size_px) * 100000 + (cols // size_px)


def cv_score(X, y, w, groups, n_splits=5, random_split=False):
    """Weighted R2 pooled over folds; groups=None with random_split gives the optimistic number."""
    preds = np.full(y.shape, np.nan)
    if random_split:
        splitter = KFold(n_splits=n_splits, shuffle=True, random_state=42).split(X)
    else:
        splitter = GroupKFold(n_splits=min(n_splits, len(np.unique(groups)))).split(X, y, groups)
    for tr, te in splitter:
        m = RandomForestRegressor(**RF_KW).fit(X[tr], y[tr], sample_weight=w[tr])
        preds[te] = m.predict(X[te])
    return r2(y, preds, w), preds


def block_boot_mean(v, rows, cols, n=400, seed=0):
    bid = blocks(rows, cols, BOOT_BLOCK_PX)
    ub, inv = np.unique(bid, return_inverse=True)
    sums = np.bincount(inv, weights=v); cnt = np.bincount(inv).astype(float)
    rng = np.random.default_rng(seed); means = []
    for _ in range(n):
        pick = rng.integers(0, ub.size, ub.size); means.append(sums[pick].sum() / cnt[pick].sum())
    return float(v.mean()), float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def main():
    print("== grid and predictors (Tay grid) ==")
    y_tay, valid_tay, prof, transform, crs, (W, H) = load_label()
    from rasterio.transform import array_bounds
    bounds = array_bounds(H, W, transform)
    canal_density = warp_to(CANAL, transform, crs, (W, H), Resampling.average)
    dist_canal = edt_m(canal_density > 0.02)
    fire = warp_to(FIRE, transform, crs, (W, H), Resampling.bilinear)
    riv = rasterize_vec(WATERWAYS, transform, (W, H), {"LineString", "MultiLineString"})
    dist_river = edt_m(riv.astype(bool))
    setl = rasterize_vec(SETTLE, transform, (W, H), {"Point", "MultiPoint"})
    dist_settle = edt_m(setl.astype(bool))
    peat_depth = warp_window_to(PEAT, transform, crs, (W, H), Resampling.bilinear, bounds)
    peat_mask = peat_depth > 25.0
    dist_peat_edge = (distance_transform_edt(peat_mask) * PX_M).astype("float32")
    # disturbance year from the NBR series: first epoch below a fixed burn/clearing threshold
    nbr_years = sorted(int(p.stem.split("_")[1]) for p in NBR.glob("nbr_*.tif"))
    dist_year = np.full((H, W), 9999.0, "float32")
    for yr in nbr_years:
        a = warp_to(NBR / f"nbr_{yr}.tif", transform, crs, (W, H), Resampling.bilinear)
        hit = (a < 0.1) & (a != 0) & (dist_year == 9999.0)
        dist_year[hit] = yr
    rows, cols = np.indices((H, W))
    base = {"peat_depth_cm": peat_depth, "dist_to_peat_edge_m": dist_peat_edge, "canal_density": canal_density,
            "dist_to_canal_m": dist_canal, "fire_history": fire, "dist_to_river_m": dist_river, "dist_to_settle_m": dist_settle}
    print("  predictors:", list(base), "+ disturbance_year (NBR epochs", nbr_years, ")")

    # PALSAR-1 label on the same grid (vertical mm/yr, tempCoh >= 0.5, weights 1/std^2, path from frame_id)
    a1 = warp_to(ALOS1 / "alos1_vertical_mmyr.tif", transform, crs, (W, H), Resampling.average)
    a1_tc = warp_to(ALOS1 / "alos1_temporalCoherence.tif", transform, crs, (W, H), Resampling.average)
    a1_sd = warp_to(ALOS1 / "alos1_velocityStd_mmyr.tif", transform, crs, (W, H), Resampling.average)
    a1_fid = warp_to(ALOS1 / "alos1_frame_id.tif", transform, crs, (W, H), Resampling.nearest)
    valid_a1 = (a1 != 0) & np.isfinite(a1) & (a1_tc >= 0.5) & (a1_sd > 0)
    path_a1 = np.where(a1_fid == 4, 421, 422)                      # frame 4 of the mosaic = path 421

    results = {"grid": {"width": W, "height": H, "px_m": PX_M}, "labels": {}}
    for name, y_full, valid, w_full, extra in (
            ("tay_2018_2024", y_tay, valid_tay, None, {}),
            ("palsar1_2007_2009", a1, valid_a1, 1.0 / np.maximum(a1_sd, 1.0) ** 2, {"path": path_a1})):
        print(f"== label {name}: {int(valid.sum())} px ==")
        y = y_full[valid].astype("float64")
        lo, hi = np.percentile(y, [1, 99]); y = np.clip(y, lo, hi)
        w = np.ones_like(y) if w_full is None else w_full[valid].astype("float64")
        w = w / w.mean()
        r, c = rows[valid], cols[valid]
        dy = dist_year[valid].copy()
        if name.startswith("palsar1"):
            dy[dy > 2009] = 9999.0                                # no leakage from post-2009 epochs
        feats = {k: v[valid] for k, v in base.items()}; feats["disturbance_year"] = dy
        Xc = np.stack([feats[k] for k in feats], axis=1).astype("float32")
        Xxy = np.stack([c, r], axis=1).astype("float32")
        Xcxy = np.hstack([Xc, Xxy])
        res = {"n": int(y.size), "label_std_mmyr": float(np.std(y)), "clip": [float(lo), float(hi)], "cv": {}}
        designs = [("random5", None, True)] + [(f"block{bm//1000}km", blocks(r, c, int(round(bm / PX_M))), False) for bm in BLOCK_M]
        if "path" in extra:
            designs.append(("leave_one_path_out", extra["path"][valid], False))
        for dname, groups, rnd in designs:
            row = {}
            for mname, X in (("M1_xy", Xxy), ("M2_cov", Xc), ("M3_cov_xy", Xcxy)):
                sc, _ = cv_score(X, y, w, groups, random_split=rnd)
                row[mname] = round(sc, 3)
            row["M0_mean"] = 0.0
            res["cv"][dname] = row
            print(f"  {dname:20s} " + "  ".join(f"{k}={v:+.3f}" for k, v in row.items()))
        # permutation-free importance from the covariate model fitted on everything (descriptive)
        m = RandomForestRegressor(**RF_KW).fit(Xc, y, sample_weight=w)
        res["importance"] = {k: round(float(v), 3) for k, v in zip(feats, m.feature_importances_)}
        results["labels"][name] = res

    # pre-registered stratum test on the PALSAR-1 label
    print("== PALSAR-1 strata (canal-density tercile x disturbed-by-2009 x burned) ==")
    v = valid_a1
    y = np.clip(a1[v], *np.percentile(a1[v], [1, 99])); r, c = rows[v], cols[v]
    cd = canal_density[v]; t1, t2 = np.percentile(cd, [33.3, 66.7])
    cd_cls = np.where(cd <= t1, 0, np.where(cd <= t2, 1, 2))
    dist_cls = (dist_year[v] <= 2009).astype(int)
    burn_cls = (fire[v] > 0).astype(int)
    all_mean, all_lo, all_hi = block_boot_mean(y, r, c)
    strata = []
    for a in range(3):
        for b in range(2):
            for d in range(2):
                m = (cd_cls == a) & (dist_cls == b) & (burn_cls == d)
                if m.sum() < 500:
                    continue
                mean, lo_, hi_ = block_boot_mean(y[m], r[m], c[m])
                se = (hi_ - lo_) / 3.92
                strata.append({"canal_tercile": a, "disturbed_by_2009": b, "burned": d, "n": int(m.sum()),
                               "mean_mmyr": round(mean, 2), "ci95": [round(lo_, 2), round(hi_, 2)],
                               "distinct_from_all": bool(abs(mean - all_mean) > 2 * se and (lo_ > all_hi or hi_ < all_lo))})
    n_distinct = sum(s["distinct_from_all"] for s in strata)
    results["palsar1_strata"] = {"all_mean_mmyr": round(all_mean, 2), "all_ci95": [round(all_lo, 2), round(all_hi, 2)],
                                 "strata": strata, "n_distinct": n_distinct}
    print(f"  all peat {all_mean:+.2f} [{all_lo:+.2f}, {all_hi:+.2f}] mm/yr; {len(strata)} strata, {n_distinct} distinct at 2 SE")
    for s in strata:
        print(f"    canal{s['canal_tercile']} dist{s['disturbed_by_2009']} burn{s['burned']} n={s['n']:6d} {s['mean_mmyr']:+7.2f} [{s['ci95'][0]:+.2f},{s['ci95'][1]:+.2f}] {'DISTINCT' if s['distinct_from_all'] else ''}")

    # acceptance
    p = results["labels"]["palsar1_2007_2009"]["cv"]
    acc = {"palsar1_block5km_M2_ge_0.20": p["block5km"]["M2_cov"] >= 0.20,
           "palsar1_block5km_M2_beats_xy_by_0.10": p["block5km"]["M2_cov"] >= p["block5km"]["M1_xy"] + 0.10,
           "palsar1_leave_one_path_M2_ge_0.20": p["leave_one_path_out"]["M2_cov"] >= 0.20,
           "strata_ge_3_distinct": n_distinct >= 3}
    results["acceptance"] = acc
    results["verdict"] = "COVARIATE_SKILL" if all(list(acc.values())[:3]) else "STRATUM_MEANS_ONLY"
    (OUT / "metrics.json").write_text(json.dumps(results, indent=1))
    print("acceptance:", acc); print("EXPERIMENT0_VERDICT", results["verdict"])


if __name__ == "__main__":
    main()
