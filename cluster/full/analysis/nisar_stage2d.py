"""NISAR stage 2d: site-level agreement and site-versus-control contrast on the
high-passed fields from stage 2c (descending D and ascending A).

Reads the two GeoTIFFs written by stage 2c, samples 5x5 windows at the 60 dam
sites on each track's own grid, samples 500 random control windows inside
Tay's footprint at least 3 km from any site, and reports: correlation of D and A
across sites (a descriptive consistency diagnostic), site minus control contrast per
track with a bootstrap interval, and the fraction of sites whose D and A agree
in sign. Sign remains the product's own (positive = range increase) [S].

The historical independent-window bootstrap is not a spatial-dependence or
causal restoration-effect test. Product source manifests record the actual dates.
"""
import json
import os
import numpy as np
from osgeo import gdal, osr

gdal.UseExceptions()
OUT = os.path.expanduser(os.environ.get("NISAR_OUT", "~/full/nisar_out"))
HP_KM = float(os.environ.get("HP_KM", "10"))
if not np.isfinite(HP_KM) or HP_KM <= 0:
    raise ValueError("HP_KM must be finite and positive")
TAY = os.path.expanduser("~/full/tay/InSAR_LongTermRate.tif")
SITES = os.path.expanduser("~/full/dam_sites.geojson")
WIN, NCTRL, MIN_KM = 2, 500, 3.0
rng = np.random.default_rng(7)
sites = [(f["properties"].get("site_id"), *f["geometry"]["coordinates"][:2]) for f in json.load(open(SITES))["features"]]


def load(path):
    ds = gdal.Open(path); a = ds.GetRasterBand(1).ReadAsArray().astype(np.float32)
    gt = ds.GetGeoTransform(); srs = osr.SpatialReference(wkt=ds.GetProjection())
    return a, gt, srs, ds


def to_pix(gt, srs, lon, lat):
    src = osr.SpatialReference(); src.ImportFromEPSG(4326); src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    dst = srs.Clone(); dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    X, Y, _ = osr.CoordinateTransformation(src, dst).TransformPoint(float(lon), float(lat))
    return int(round((Y - gt[3]) / gt[5] - 0.5)), int(round((X - gt[0]) / gt[1] - 0.5)), X, Y


def win_mean(a, r, c, w=WIN):
    ny, nx = a.shape
    if r - w < 0 or c - w < 0 or r + w >= ny or c + w >= nx:
        return np.nan
    p = a[r - w:r + w + 1, c - w:c + w + 1]
    return float(np.nanmean(p)) if np.isfinite(p).sum() >= 13 else np.nan


def tay_mask(ds):
    mem = gdal.GetDriverByName("MEM").Create("", ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Float32)
    mem.SetGeoTransform(ds.GetGeoTransform()); mem.SetProjection(ds.GetProjection()); mem.GetRasterBand(1).Fill(float("nan"))
    src = gdal.Open(TAY); gdal.ReprojectImage(src, mem, src.GetProjection(), ds.GetProjection(), gdal.GRA_Average)
    t = mem.GetRasterBand(1).ReadAsArray(); t[t == 0] = np.nan
    return np.isfinite(t)


vals = {}
for tr in ("D", "A"):
    a, gt, srs, ds = load(os.path.join(OUT, f"nisar_{tr}_disp60_hp{HP_KM:g}km_mm.tif"))
    strip = tay_mask(ds)
    sv, sxy = [], []
    for _, lon, lat in sites:
        r, c, X, Y = to_pix(gt, srs, lon, lat); sv.append(win_mean(a, r, c)); sxy.append((X, Y))
    sv = np.array(sv); sxy = np.array(sxy)
    # Controls are footprint windows, not a separately verified peat class.
    rr, cc = np.nonzero(strip & np.isfinite(a))
    pick = rng.choice(len(rr), size=min(20000, len(rr)), replace=False)
    ctrl = []
    for i in pick:
        X = gt[0] + (cc[i] + 0.5) * gt[1]; Y = gt[3] + (rr[i] + 0.5) * gt[5]
        if np.min(np.hypot(sxy[:, 0] - X, sxy[:, 1] - Y)) < MIN_KM * 1000:
            continue
        v = win_mean(a, rr[i], cc[i])
        if np.isfinite(v):
            ctrl.append(v)
        if len(ctrl) >= NCTRL:
            break
    ctrl = np.array(ctrl)
    ins = np.isfinite(sv)
    if ins.sum() < 2 or len(ctrl) < 2:
        raise ValueError(f"[{tr}] insufficient finite sites or controls: {ins.sum()}, {len(ctrl)}")
    d = np.nanmean(sv) - ctrl.mean()
    boots = [np.mean(rng.choice(sv[ins], ins.sum())) - np.mean(rng.choice(ctrl, len(ctrl))) for _ in range(2000)]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    print(f"[{tr}] sites in frame {ins.sum()}/{len(sites)}: mean {np.nanmean(sv):+.1f} mm (sd {np.nanstd(sv):.1f}); "
          f"{len(ctrl)} footprint controls >= {MIN_KM:g} km away: mean {ctrl.mean():+.1f} (sd {ctrl.std():.1f}); "
          f"site minus control {d:+.1f} mm, 95% CI [{lo:+.1f}, {hi:+.1f}]")
    vals[tr] = sv

both = np.isfinite(vals["D"]) & np.isfinite(vals["A"])
if both.sum() >= 10:
    x, y = vals["D"][both], vals["A"][both]
    r = np.corrcoef(x, y)[0, 1]
    boots = [np.corrcoef(*np.array([x, y])[:, rng.integers(0, both.sum(), both.sum())])[0, 1] for _ in range(2000)]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    print(f"\nsite-level two-track agreement, n={both.sum()} sites in both frames: r={r:+.3f} 95% CI [{lo:+.3f}, {hi:+.3f}]; "
          f"same sign {100*np.mean(np.sign(x) == np.sign(y)):.0f}%; mean |D-A| {np.mean(np.abs(x-y)):.1f} mm")
    order = np.argsort(-np.abs((x + y) / 2))[:10]
    ids = np.array([s[0] for s in sites])[both]
    print("largest shared anomalies (site, D, A) mm:", ", ".join(f"{ids[i]} {x[i]:+.0f}/{y[i]:+.0f}" for i in order))
print("STAGE2D_DONE")
