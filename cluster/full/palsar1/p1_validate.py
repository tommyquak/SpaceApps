"""Historical-baseline checks on the ALOS-1 PALSAR 2007-2009 velocity.

Usage: p1_validate.py <mintpy_dir> <tag>

Prints (and writes <mintpy_dir>/validate.json):
  coverage     fraction of the frame and of the dam zone with tempCoh >= 0.6
  rates        median / mean LOS rate (mm/yr, positive = toward the satellite)
               over the dam zone and over the whole peat footprint, with
               block-bootstrap 95 percent CIs (1.8 km blocks), before and
               after re-referencing to the published non-peat reference
               point of Zahro et al. 2025 (Sci Rep 15:23749; -2.213489,
               113.920443, validated against InaCORS GNSS), median of the
               good pixels within 1 km of it, when the frame contains it.
               PeatGrids cannot define mineral ground (381 cm under
               Palangkaraya city) and SRTM heights sit at 45-80 m under
               canopy, so the stable-ground NULL uses ESA WorldCover 2021
               built-up pixels (class 50) instead.
  literature   Zahro et al. 2025 (Sci Rep 15:23749), same ground, same
               era, same processor: mean -1.74 cm/yr, max -15.4 cm/yr,
               76 percent coherent [E]. Hooijer et al. 2012 (Biogeosciences
               9:1053): ~5 cm/yr long-term at 0.7 m water table [E].
  verdict      PASS when (a) dam-zone coverage >= 40 percent at TC_MIN
               (env MINTEMPCOH, default 0.5), (b) the
               re-referenced dam-zone mean is subsidence with a CI clear of
               zero, (c) it lies within 0.5 to 6 cm/yr (the Zahro to Hooijer
               envelope), and (d) stable ground sits within +/- 5 mm/yr.
These are consistency checks against published regional rates, not a
per-pixel validation; there is no independent 2007-2009 raster.
"""
import json
import os
import sys

import h5py
import numpy as np
from osgeo import gdal

ZAHRO_MEAN_CM = -1.74
ZAHRO_MAX_CM = -15.4
ZAHRO_COVER = 0.76
HOOIJER_LONGTERM_CM = -5.0
DAM_ZONE = (114.31, 114.79, -2.985, -2.435)   # lon_min, lon_max, lat_min, lat_max
BLOCK_DEG = 0.0167                             # ~1.8 km
TC_MIN = float(os.environ.get("MINTEMPCOH", "0.5"))
PEAT = os.path.expanduser(os.environ.get("PEAT_TIF", "~/full/peat/PEATGRIDS_THICKNESS_CM_MEAN.tif"))
PEAT_MIN_CM = 50.0                             # peat where PeatGrids says >= 50 cm
ZAHRO_REF = (113.920443, -2.213489)            # Zahro et al. 2025 reference point (non-peat, InaCORS-checked)
ZAHRO_R_DEG = 0.009                            # ~1 km window
WORLDCOVER = os.path.expanduser(os.environ.get("WORLDCOVER_DIR", "~/full/peat/worldcover"))
BUILTUP = 50                                   # ESA WorldCover class: built-up


def geo(h5, key):
    with h5py.File(h5) as f:
        a = f[key][:]
        at = dict(f.attrs)
    x0, dx, y0, dy = (float(at["X_FIRST"]), float(at["X_STEP"]), float(at["Y_FIRST"]), float(at["Y_STEP"]))
    ny, nx = a.shape
    lon = x0 + dx * (np.arange(nx) + 0.5)
    lat = y0 + dy * (np.arange(ny) + 0.5)
    return a, np.meshgrid(lon, lat)


def peat_thickness(lon, lat):
    """PeatGrids mean thickness (cm) sampled at the product's pixel centres."""
    ds = gdal.Open(PEAT)
    a = ds.GetRasterBand(1).ReadAsArray().astype(np.float64)
    nd = ds.GetRasterBand(1).GetNoDataValue()
    if nd is not None:
        a[a == nd] = np.nan
    gt = ds.GetGeoTransform()
    c = np.clip(((lon - gt[0]) / gt[1]).astype(int), 0, a.shape[1] - 1)
    r = np.clip(((lat - gt[3]) / gt[5]).astype(int), 0, a.shape[0] - 1)
    return a[r, c]


def worldcover_class(lon, lat):
    """ESA WorldCover 2021 class sampled at the product's pixel centres (nearest), 0 where no tile."""
    import glob
    out = np.zeros(lon.shape, np.uint8)
    for tif in glob.glob(f"{WORLDCOVER}/ESA_WorldCover_10m_2021_v200_*_Map.tif"):
        ds = gdal.Open(tif); gt = ds.GetGeoTransform()
        nx, ny = ds.RasterXSize, ds.RasterYSize
        inside = (lon >= gt[0]) & (lon < gt[0] + gt[1] * nx) & (lat <= gt[3]) & (lat > gt[3] + gt[5] * ny)
        if not inside.any():
            continue
        # read the tile window covering the product footprint at the product's resolution
        x0 = int((lon[inside].min() - gt[0]) / gt[1]); x1 = int((lon[inside].max() - gt[0]) / gt[1]) + 1
        y0 = int((lat[inside].max() - gt[3]) / gt[5]); y1 = int((lat[inside].min() - gt[3]) / gt[5]) + 1
        x0, y0 = max(x0, 0), max(y0, 0); x1, y1 = min(x1, nx), min(y1, ny)
        a = ds.GetRasterBand(1).ReadAsArray(x0, y0, x1 - x0, y1 - y0)
        c = np.clip(((lon[inside] - gt[0]) / gt[1]).astype(int) - x0, 0, a.shape[1] - 1)
        r = np.clip(((lat[inside] - gt[3]) / gt[5]).astype(int) - y0, 0, a.shape[0] - 1)
        out[inside] = a[r, c]
    return out


def block_boot(v, lon, lat, n=500, seed=0):
    """Mean with a spatial block bootstrap CI (blocks of BLOCK_DEG)."""
    if v.size < 20:
        return float("nan"), (float("nan"), float("nan"))
    bx = np.floor((lon - lon.min()) / BLOCK_DEG).astype(int)
    by = np.floor((lat - lat.min()) / BLOCK_DEG).astype(int)
    bid = bx * 100000 + by
    ub, inv = np.unique(bid, return_inverse=True)
    sums = np.bincount(inv, weights=v)
    cnts = np.bincount(inv).astype(float)
    rng = np.random.default_rng(seed)
    means = []
    for _ in range(n):
        pick = rng.integers(0, ub.size, ub.size)
        means.append(sums[pick].sum() / cnts[pick].sum())
    return float(v.mean()), (float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5)))


def main():
    w, tag = sys.argv[1], sys.argv[2]
    vel, (lon, lat) = geo(f"{w}/geo/geo_velocity.h5", "velocity")
    vel = vel * 1000.0                                   # m/yr -> mm/yr
    tc, _ = geo(f"{w}/geo/geo_temporalCoherence.h5", "temporalCoherence")
    thick = peat_thickness(lon, lat)
    valid = np.isfinite(vel) & np.isfinite(tc) & (tc > 0)
    good = valid & (tc >= TC_MIN)
    dz = (lon >= DAM_ZONE[0]) & (lon <= DAM_ZONE[1]) & (lat >= DAM_ZONE[2]) & (lat <= DAM_ZONE[3])
    peat = np.isfinite(thick) & (thick >= PEAT_MIN_CM)
    wc = worldcover_class(lon, lat)
    stable = good & (wc == BUILTUP) & ~peat                     # built-up, non-peat: the null
    zwin = good & (np.abs(lon - ZAHRO_REF[0]) <= ZAHRO_R_DEG) & (np.abs(lat - ZAHRO_REF[1]) <= ZAHRO_R_DEG)
    # the reference window must sit well inside the frame's valid data: on frame 7130 it
    # falls on the northern azimuth edge, where the phase ramp made the datum 21 mm/yr off
    # the overlap-tied mosaic value. Require 15 km of valid data on every side.
    if zwin.any():
        vy, vx = np.where(valid)
        edge = min(lat[zwin].max() - lat[valid].min(), lat[valid].max() - lat[zwin].min(),
                   lon[zwin].max() - lon[valid].min(), lon[valid].max() - lon[zwin].min())
        if edge < 0.135:            # 15 km: frame 7130 still disagreed with the overlap tie by 21 mm/yr at 6 km
            zwin[:] = False
    out = {"tag": tag, "tc_min": TC_MIN, "n_valid": int(valid.sum()), "n_good": int(good.sum()),
           "n_zahro_window": int(zwin.sum()), "n_builtup": int(stable.sum())}
    out["coverage_frame"] = float(good.sum() / max(valid.sum(), 1))
    out["coverage_damzone"] = float((good & dz).sum() / max((valid & dz).sum(), 1))
    out["coverage_peat"] = float((good & peat).sum() / max((valid & peat).sum(), 1))
    out["coverage_frame_tc06"] = float((valid & (tc >= 0.6)).sum() / max(valid.sum(), 1))
    out["n_stable"] = int(stable.sum())
    out["n_peat_valid"] = int((valid & peat).sum())

    ref, datum = 0.0, "none"
    # a single-frame datum is accepted only when the built-up null exists in the same frame
    # and agrees within 5 mm/yr; frame 7130 reaches the reference point only at its far
    # northern edge (no built-up pixels survive there) and was 25 mm/yr off the overlap tie
    if zwin.sum() >= 20 and stable.sum() >= 50 and abs(float(np.median(vel[stable])) - float(np.median(vel[zwin]))) <= 5.0:
        ref, datum = float(np.median(vel[zwin])), "zahro2025_refpoint"
    elif stable.sum() >= 50:
        ref, datum = float(np.median(vel[stable])), "worldcover_builtup"
    out["datum"] = datum
    out["datum_shift_mmyr_raw"] = ref
    v = vel - ref                                        # re-reference to non-peat ground
    for name, m in (("damzone", good & dz), ("peat", good & peat), ("frame", good), ("stable", stable)):
        if m.sum() < 20:
            out[name] = {"n": int(m.sum())}
            continue
        mean, ci = block_boot(v[m], lon[m], lat[m])
        out[name] = {"n": int(m.sum()), "median_mmyr": float(np.median(v[m])), "mean_mmyr": mean,
                     "ci95": ci, "p10_mmyr": float(np.percentile(v[m], 10)), "p90_mmyr": float(np.percentile(v[m], 90))}
    d = out.get("damzone", {})
    # dam zone minus surrounding peat: the one number that does not depend on the datum
    pz = out.get("peat", {})
    if "mean_mmyr" in d and "mean_mmyr" in pz:
        both = good & (dz | peat)
        contrast = v.copy(); contrast[~dz] = np.nan
        out["damzone_minus_peat_mmyr"] = d["mean_mmyr"] - pz["mean_mmyr"]
        out["damzone_minus_peat_ci95"] = [d["ci95"][0] - pz["ci95"][1], d["ci95"][1] - pz["ci95"][0]]
    no_datum = datum == "none"
    in_frame = (valid & dz).sum() >= 1000
    checks = {
        "coverage_ge_40pct": (out["coverage_damzone"] >= 0.40) if in_frame else "n/a",
        "damzone_subsides_ci_clear": (bool(d.get("ci95") and d["ci95"][1] < 0)) if in_frame else "n/a",
        # absolute rates only mean something once the frame is tied to mineral ground
        "within_zahro_hooijer_envelope": ("n/a" if (no_datum or not in_frame) else
                                          bool(d.get("mean_mmyr") is not None and -60 <= d.get("mean_mmyr", 0) <= -5)),
        # frame 7130 is peat wall to wall on PeatGrids (94 percent >= 50 cm): with no
        # mineral pixels the datum check is not applicable and the mosaic step ties the
        # frame to its neighbour's (Palangkaraya) mineral ground instead
        "stable_within_5mmyr": ("n/a" if out["n_stable"] < 50 else
                                bool(out.get("stable", {}).get("mean_mmyr") is not None and abs(out["stable"]["mean_mmyr"]) <= 5)),
    }
    out["checks"] = checks
    out["verdict"] = "PASS" if all(v is True or v == "n/a" for v in checks.values()) else "FAIL"
    if checks["stable_within_5mmyr"] == "n/a":
        out["verdict"] += "_NO_DATUM"
    out["literature"] = {"zahro2025_mean_cmyr": ZAHRO_MEAN_CM, "zahro2025_max_cmyr": ZAHRO_MAX_CM,
                         "zahro2025_coverage": ZAHRO_COVER, "hooijer2012_longterm_cmyr": HOOIJER_LONGTERM_CM}
    with open(f"{w}/validate.json", "w") as f:
        json.dump(out, f, indent=1)
    print(f"P1_VALIDATE {tag}: tempCoh>={TC_MIN:.2f} coverage frame {out['coverage_frame']:.2f} peat {out['coverage_peat']:.2f} "
          f"damzone {out['coverage_damzone']:.2f} (Zahro 0.76; at 0.6: {out['coverage_frame_tc06']:.2f}); "
          f"datum {datum} (shift {ref:+.1f} mm/yr; Zahro window px {out['n_zahro_window']}, built-up px {out['n_stable']})")
    for k in ("damzone", "peat", "frame", "stable"):
        r = out.get(k, {})
        if "mean_mmyr" in r:
            print(f"  {k:8s} n={r['n']:7d} median {r['median_mmyr']:+7.1f} mean {r['mean_mmyr']:+7.1f} "
                  f"CI [{r['ci95'][0]:+.1f}, {r['ci95'][1]:+.1f}] p10/p90 {r['p10_mmyr']:+.1f}/{r['p90_mmyr']:+.1f} mm/yr")
    if "damzone_minus_peat_mmyr" in out:
        c = out["damzone_minus_peat_ci95"]
        print(f"  dam zone minus peat: {out['damzone_minus_peat_mmyr']:+.1f} mm/yr (conservative CI [{c[0]:+.1f}, {c[1]:+.1f}])")
    print("  checks:", {k: ("n/a" if v == "n/a" else "ok" if v else "FAIL") for k, v in checks.items()})
    print(f"P1_VERDICT {out['verdict']}")


if __name__ == "__main__":
    main()
