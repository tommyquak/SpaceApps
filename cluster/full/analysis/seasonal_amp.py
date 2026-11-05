"""Component 3: seasonal amplitude of peat surface motion from the 60-day
network series (bw60), the product the short-baseline network exists for.

Per pixel, least-squares fit of d(t) = a + b t + A cos(2 pi t) + B sin(2 pi t)
over the 93 epochs (t in years); amplitude sqrt(A^2 + B^2), phase as the
month of the annual minimum, and the RMS of the residual as the noise level.
Reported for Tay's footprint, the 60 dam sites (5x5 windows), 500 random peat
controls at least 3 km from any site, and stable mineral ground, with a
bootstrap interval on site minus control. The trend b is NOT reported as a
rate: it carries the short-baseline fading bias (see the audit). [E] fit;
[H] amplitude as a rewetting metric.
"""
import datetime as dt
import json
import os
import numpy as np
import h5py
from osgeo import gdal, osr

gdal.UseExceptions()
W = "/mnt/scratch/t/tommy/full/mintpy_bw60"
GEOM = os.path.expanduser("~/full/mintpy_tay/inputs/geometryGeo.h5")
TAY = os.path.expanduser("~/full/tay/InSAR_LongTermRate.tif")
TC20 = os.path.expanduser("~/full/mintpy_tay/temporalCoherence.h5")
SITES = os.path.expanduser("~/full/dam_sites.geojson")
OUT = os.path.expanduser("~/full/analysis_out"); os.makedirs(OUT, exist_ok=True)
SERIES = os.environ.get("SERIES", "timeseries_SET.h5")
W5, NCTRL, MIN_KM = 2, 500, 3.0
rng = np.random.default_rng(3)

with h5py.File(os.path.join(W, SERIES), "r") as f:
    dates = [dt.datetime.strptime(d.decode(), "%Y%m%d") for d in f["date"][:]]
    ts = f["timeseries"][:].astype(np.float32) * 1000.0  # mm, positive toward satellite
    att = dict(f.attrs)
with h5py.File(os.path.join(W, "temporalCoherence.h5"), "r") as f:
    tc = f["temporalCoherence"][:].astype(np.float32)
with h5py.File(TC20, "r") as f:
    stable = f["temporalCoherence"][:] >= 0.7
with h5py.File(GEOM, "r") as f:
    inc = np.deg2rad(f["incidenceAngle"][:].astype(np.float32))
ne, ny, nx = ts.shape
t = np.array([(d - dt.datetime(2000, 1, 1)).days / 365.25 for d in dates], dtype=np.float64)
G = np.stack([np.ones_like(t), t - t.mean(), np.cos(2 * np.pi * t), np.sin(2 * np.pi * t)], 1)
Gp = np.linalg.pinv(G)  # 4 x ne
Y = ts.reshape(ne, -1).astype(np.float64)
coef = Gp @ Y                                  # 4 x npix
resid = Y - G @ coef
amp = np.sqrt(coef[2] ** 2 + coef[3] ** 2).reshape(ny, nx) / np.cos(inc)   # vertical mm
rms = np.sqrt(np.mean(resid ** 2, axis=0)).reshape(ny, nx)
# phase: fraction of year at which the cosine term is minimum (surface lowest)
ph = np.arctan2(coef[3], coef[2]).reshape(ny, nx)          # d ~ R cos(2 pi t - ph)
min_frac = ((ph + np.pi) / (2 * np.pi)) % 1.0              # minimum of cos at 2 pi t - ph = pi
min_month = (min_frac * 12.0 + 1.0)
del Y, resid

x0, dx, y0, dy = float(att["X_FIRST"]), float(att["X_STEP"]), float(att["Y_FIRST"]), float(att["Y_STEP"])
mem = gdal.GetDriverByName("MEM").Create("", nx, ny, 1, gdal.GDT_Float32); mem.SetGeoTransform((x0, dx, 0.0, y0, 0.0, dy))
srs = osr.SpatialReference(); srs.ImportFromEPSG(4326); mem.SetProjection(srs.ExportToWkt()); mem.GetRasterBand(1).Fill(float("nan"))
src = gdal.Open(TAY); gdal.ReprojectImage(src, mem, src.GetProjection(), srs.ExportToWkt(), gdal.GRA_Average)
tay = mem.GetRasterBand(1).ReadAsArray(); tay[tay == 0] = np.nan; strip = np.isfinite(tay)
sites = [(f["properties"].get("site_id"), *f["geometry"]["coordinates"][:2]) for f in json.load(open(SITES))["features"]]
site_rc = [(int(round((lat - y0) / dy)), int(round((lon - x0) / dx))) for _, lon, lat in sites]


def wmean(a, r, c):
    if not (W5 <= r < ny - W5 and W5 <= c < nx - W5):
        return np.nan
    p = a[r - W5:r + W5 + 1, c - W5:c + W5 + 1]; return float(np.nanmean(p)) if np.isfinite(p).sum() >= 13 else np.nan


good = tc >= 0.5
print(f"series {SERIES}, {ne} epochs, fit on {good.sum()} px with bw60 tempCoh >= 0.5 ({100*good.mean():.1f}% of frame)")
for lbl, m in (("Tay footprint", strip & good), ("stable ground", stable & good), ("whole frame", good)):
    print(f"   [{lbl}] n={m.sum():7d}  amplitude median {np.nanmedian(amp[m]):.1f} mm  p90 {np.nanpercentile(amp[m],90):.1f}  | residual RMS median {np.nanmedian(rms[m]):.1f} mm  | month of minimum, median {np.nanmedian(min_month[m]):.1f}")
sa = np.array([wmean(np.where(good, amp, np.nan), r, c) for r, c in site_rc]); sm = np.array([wmean(np.where(good, min_month, np.nan), r, c) for r, c in site_rc])
rr, cc = np.nonzero(strip & good); pick = rng.choice(len(rr), min(20000, len(rr)), replace=False)
sxy = np.array([(x0 + c * dx, y0 + r * dy) for r, c in site_rc]); ctrl = []
for i in pick:
    X, Yy = x0 + cc[i] * dx, y0 + rr[i] * dy
    if np.min(np.hypot((sxy[:, 0] - X) * 111.32 * np.cos(np.deg2rad(Yy)), (sxy[:, 1] - Yy) * 111.32)) < MIN_KM:
        continue
    v = wmean(np.where(good, amp, np.nan), rr[i], cc[i])
    if np.isfinite(v):
        ctrl.append(v)
    if len(ctrl) >= NCTRL:
        break
ctrl = np.array(ctrl); ins = np.isfinite(sa)
d = np.nanmean(sa) - ctrl.mean()
boots = [np.mean(rng.choice(sa[ins], ins.sum())) - np.mean(rng.choice(ctrl, len(ctrl))) for _ in range(2000)]
lo, hi = np.percentile(boots, [2.5, 97.5])
print(f"\n   dam sites (n={ins.sum()}): amplitude mean {np.nanmean(sa):.1f} mm (sd {np.nanstd(sa):.1f}), month of minimum median {np.nanmedian(sm):.1f}")
print(f"   peat controls (n={len(ctrl)}, >= {MIN_KM:g} km from sites): amplitude mean {ctrl.mean():.1f} mm (sd {ctrl.std():.1f})")
print(f"   site minus control: {d:+.1f} mm, 95% CI [{lo:+.1f}, {hi:+.1f}]  (positive = sites breathe more, drained-peat reading [H])")
order = np.argsort(-np.nan_to_num(sa))[:12]
print("   largest site amplitudes: " + ", ".join(f"{sites[k][0]} {sa[k]:.0f}" for k in order if np.isfinite(sa[k])))
for name, a in (("seasonal_amp_mm", amp), ("seasonal_min_month", min_month), ("fit_resid_rms_mm", rms)):
    a = np.where(good, a, np.nan).astype(np.float32)
    ds = gdal.GetDriverByName("GTiff").Create(os.path.join(OUT, f"bw60_{name}.tif"), nx, ny, 1, gdal.GDT_Float32, ["COMPRESS=DEFLATE", "TILED=YES"])
    ds.SetGeoTransform((x0, dx, 0.0, y0, 0.0, dy)); ds.SetProjection(srs.ExportToWkt()); ds.GetRasterBand(1).WriteArray(a); ds.GetRasterBand(1).SetNoDataValue(float("nan")); ds = None
print(f"wrote {OUT}/bw60_seasonal_amp_mm.tif, bw60_seasonal_min_month.tif, bw60_fit_resid_rms_mm.tif")
print("SEASAMP_DONE")
