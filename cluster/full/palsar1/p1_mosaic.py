"""Mosaic the per-frame PALSAR-1 velocities into one 2007-2009 baseline raster.

Usage: p1_mosaic.py <out_dir> <mintpy_dir_ref> <mintpy_dir_2> [<mintpy_dir_3> ...]

Each MintPy dir must hold geo_velocity.tif, geo_temporalCoherence.tif and
geo_velocityStd.tif (written by p1_mintpy.sbatch). Frames are processed one
at a time (scratch quota), so each has its own reference pixel; the first
directory is the datum and every later frame is shifted by the median
velocity difference over its along-track overlap with the mosaic so far
(needs >= 2000 shared pixels at tempCoh >= 0.6), then merged with
inverse-variance weights where frames overlap. Outputs in <out_dir>:
  alos1_velocity_mmyr.tif, alos1_velocityStd_mmyr.tif,
  alos1_temporalCoherence.tif, alos1_frame_id.tif, mosaic.json
The merged field is re-referenced to the Zahro et al. 2025 non-peat reference
point (-2.213489, 113.920443; median of good pixels within 1 km) when the
mosaic contains it, and the dam-zone rate is reported with a block bootstrap.
The offsets applied are recorded; a shift above 10 mm/yr is flagged because
it means the two frames disagree on their reference ground, not on the peat.
"""
import json
import os
import sys

import numpy as np
from osgeo import gdal

STEP = 0.000833333333333
MIN_OVERLAP = 2000
TC_MIN = float(os.environ.get("MINTEMPCOH", "0.5"))   # same threshold as p1_validate; 0.6 dropped path 421 to 39 percent
ZAHRO_REF = (113.920443, -2.213489)
ZAHRO_R_DEG = 0.009
DAM_ZONE = (114.31, 114.79, -2.985, -2.435)
BLOCK_DEG = 0.0167


def read(path):
    ds = gdal.Open(path)
    a = ds.GetRasterBand(1).ReadAsArray().astype(np.float64)
    nd = ds.GetRasterBand(1).GetNoDataValue()
    if nd is not None:
        a[a == nd] = np.nan
    gt = ds.GetGeoTransform()
    return a, gt


def warp_to(path, x0, y0, nx, ny):
    """Resample a GeoTIFF onto the common lat/lon grid (bilinear)."""
    ds = gdal.Warp("", path, format="MEM", outputBounds=(x0, y0 - ny * STEP, x0 + nx * STEP, y0),
                   xRes=STEP, yRes=STEP, resampleAlg="bilinear", dstNodata=np.nan)
    return ds.GetRasterBand(1).ReadAsArray().astype(np.float64)


def write(path, a, x0, y0, dtype=gdal.GDT_Float32):
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(path, a.shape[1], a.shape[0], 1, dtype, options=["COMPRESS=DEFLATE", "TILED=YES"])
    ds.SetGeoTransform((x0, STEP, 0, y0, 0, -STEP))
    ds.SetProjection('GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]')
    b = ds.GetRasterBand(1)
    b.WriteArray(np.where(np.isfinite(a), a, -9999).astype(np.float32) if dtype == gdal.GDT_Float32 else a)
    b.SetNoDataValue(-9999)
    ds.FlushCache()


def main():
    out, dirs = sys.argv[1], sys.argv[2:]
    os.makedirs(out, exist_ok=True)
    # common grid = union of the frame extents
    xs, ys = [], []
    for d in dirs:
        ds = gdal.Open(f"{d}/geo_velocity.tif"); gt = ds.GetGeoTransform()
        xs += [gt[0], gt[0] + gt[1] * ds.RasterXSize]; ys += [gt[3], gt[3] + gt[5] * ds.RasterYSize]
    x0, y0 = min(xs), max(ys)
    nx, ny = int(np.ceil((max(xs) - x0) / STEP)), int(np.ceil((y0 - min(ys)) / STEP))
    vel = np.full((ny, nx), np.nan); var = np.full((ny, nx), np.nan); tc = np.full((ny, nx), np.nan)
    fid = np.zeros((ny, nx), np.int16)
    log = {"grid": {"x0": x0, "y0": y0, "nx": nx, "ny": ny, "step": STEP}, "frames": []}
    tied = {}                                                    # frame name -> (tied velocity, good mask)
    for i, d in enumerate(dirs, start=1):
        v = warp_to(f"{d}/geo_velocity.tif", x0, y0, nx, ny) * 1000.0      # m/yr -> mm/yr (LOS)
        try:
            inc_ds = gdal.Warp("", f"{d}/geo_incidenceAngle.tif", format="MEM", outputBounds=(x0, y0 - ny * STEP, x0 + nx * STEP, y0),
                               xRes=STEP, yRes=STEP, resampleAlg="bilinear", dstNodata=np.nan)
            inc = inc_ds.GetRasterBand(1).ReadAsArray().astype(np.float64)
        except Exception:
            inc = np.full((ny, nx), 38.7)                                   # ALOS-1 FBD off-nadir 34.3 deg -> ~38.7 deg incidence
        if i == 1:
            incm = inc.copy()
        else:
            incm = np.where(np.isfinite(incm), incm, inc)
        s = warp_to(f"{d}/geo_velocityStd.tif", x0, y0, nx, ny) * 1000.0
        c = warp_to(f"{d}/geo_temporalCoherence.tif", x0, y0, nx, ny)
        good = np.isfinite(v) & np.isfinite(c) & (c >= TC_MIN)
        shift, n_ov = 0.0, 0
        if i > 1:
            ov = good & np.isfinite(vel)
            n_ov = int(ov.sum())
            if n_ov < MIN_OVERLAP:
                print(f"frame {i} ({os.path.basename(d)}): only {n_ov} overlap px, cannot align, SKIPPED")
                log["frames"].append({"dir": d, "skipped": True, "overlap_px": n_ov}); continue
            shift = float(np.median(vel[ov] - v[ov]))
            v = v + shift
        w_new = 1.0 / np.maximum(s, 0.5) ** 2
        w_old = 1.0 / np.maximum(np.sqrt(var), 0.5) ** 2
        both = good & np.isfinite(vel)
        only = good & ~np.isfinite(vel)
        tied[os.path.basename(d)] = (v.copy(), good.copy())
        vel[both] = (vel[both] * w_old[both] + v[both] * w_new[both]) / (w_old[both] + w_new[both])
        var[both] = 1.0 / (w_old[both] + w_new[both])
        vel[only] = v[only]; var[only] = s[only] ** 2; fid[only] = i
        tc[good] = np.fmax(tc[good], c[good])
        log["frames"].append({"dir": d, "shift_mmyr": shift, "overlap_px": n_ov, "good_px": int(good.sum()),
                              "flag": "LARGE_SHIFT" if abs(shift) > 10 else ""})
        print(f"frame {i} ({os.path.basename(d)}): good px {int(good.sum())}, overlap {n_ov}, shift {shift:+.1f} mm/yr")
    # datum: published non-peat reference point, then dam-zone statistics
    lonc = x0 + STEP * (np.arange(nx) + 0.5); latc = y0 - STEP * (np.arange(ny) + 0.5)
    LON, LAT = np.meshgrid(lonc, latc)
    good = np.isfinite(vel)
    zwin = good & (np.abs(LON - ZAHRO_REF[0]) <= ZAHRO_R_DEG) & (np.abs(LAT - ZAHRO_REF[1]) <= ZAHRO_R_DEG)
    # The datum comes from the ONE frame in which the reference window is deep inside the
    # valid data (>= 15 km from every edge), using that frame's tied values only: a frame
    # whose azimuth edge crosses the window carries a ramp there (7130 was 21 mm/yr off).
    best, best_edge = None, 0.0
    for k, (fv, fgood) in tied.items():
        if not (fgood & zwin).any():
            continue
        edge = min(LAT[fgood & zwin].max() - LAT[fgood].min(), LAT[fgood].max() - LAT[fgood & zwin].min(),
                   LON[fgood & zwin].max() - LON[fgood].min(), LON[fgood].max() - LON[fgood & zwin].min())
        if edge > best_edge:
            best, best_edge = k, edge
    if best is not None and best_edge >= 0.135 and (tied[best][1] & zwin).sum() >= 20:
        fv, fgood = tied[best]
        shift = float(np.median(fv[fgood & zwin])); vel = vel - shift
        log["datum"] = {"source": "zahro2025_refpoint", "frame": best, "n": int((fgood & zwin).sum()), "shift_mmyr": shift, "edge_deg": best_edge}
        print(f"datum: Zahro 2025 reference window in {best} ({best_edge*111:.0f} km inside), {int((fgood & zwin).sum())} px, shift {shift:+.1f} mm/yr")
    else:
        log["datum"] = {"source": "none"}; print("datum: NONE (reference point not deep inside any frame)")
    def boot(mask):
        v = vel[mask]; lo = LON[mask]; la = LAT[mask]
        if v.size < 20:
            return None
        bid = np.floor((lo - lo.min()) / BLOCK_DEG).astype(int) * 100000 + np.floor((la - la.min()) / BLOCK_DEG).astype(int)
        ub, inv = np.unique(bid, return_inverse=True); sums = np.bincount(inv, weights=v); cnts = np.bincount(inv).astype(float)
        rng = np.random.default_rng(0); means = []
        for _ in range(500):
            pick = rng.integers(0, ub.size, ub.size); means.append(sums[pick].sum() / cnts[pick].sum())
        return {"n": int(v.size), "median_mmyr": float(np.median(v)), "mean_mmyr": float(v.mean()),
                "ci95": [float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))]}
    dz = good & (LON >= DAM_ZONE[0]) & (LON <= DAM_ZONE[1]) & (LAT >= DAM_ZONE[2]) & (LAT <= DAM_ZONE[3])
    log["damzone"] = boot(dz); log["all"] = boot(good)
    for k in ("damzone", "all"):
        r = log[k]
        if r:
            print(f"  {k:8s} n={r['n']:7d} median {r['median_mmyr']:+.1f} mean {r['mean_mmyr']:+.1f} CI [{r['ci95'][0]:+.1f}, {r['ci95'][1]:+.1f}] mm/yr")
    vert = vel / np.cos(np.deg2rad(np.where(np.isfinite(incm), incm, 38.7)))
    log["vertical_note"] = "vertical = LOS / cos(incidence); positive up, negative = subsidence"
    if log["damzone"]:
        vd = vert[dz]; log["damzone"]["vertical_mean_mmyr"] = float(np.nanmean(vd)); log["damzone"]["vertical_median_mmyr"] = float(np.nanmedian(vd))
        print(f"  damzone vertical mean {log['damzone']['vertical_mean_mmyr']:+.1f} mm/yr (median {log['damzone']['vertical_median_mmyr']:+.1f})")
    if log["all"]:
        va = vert[good]; log["all"]["vertical_mean_mmyr"] = float(np.nanmean(va)); print(f"  all      vertical mean {log['all']['vertical_mean_mmyr']:+.1f} mm/yr")
    write(f"{out}/alos1_vertical_mmyr.tif", vert, x0, y0)
    write(f"{out}/alos1_velocity_mmyr.tif", vel, x0, y0)
    write(f"{out}/alos1_velocityStd_mmyr.tif", np.sqrt(var), x0, y0)
    write(f"{out}/alos1_temporalCoherence.tif", tc, x0, y0)
    write(f"{out}/alos1_frame_id.tif", fid, x0, y0, gdal.GDT_Int16)
    log["mosaic_good_px"] = int(np.isfinite(vel).sum())
    with open(f"{out}/mosaic.json", "w") as f:
        json.dump(log, f, indent=1)
    print(f"P1_MOSAIC_OK {out} good px {log['mosaic_good_px']}")


if __name__ == "__main__":
    main()
