"""Apply a Zheng bias time series (BIAS) to a MintPy series (BW/timeseries_SET.h5)
to the 60-day-network series and test the corrected rate against Tay's rate
raster. Tay's own remedy for short-baseline rates. Compare, never calibrate.
"""
import datetime as dt
import os
import numpy as np
import h5py
from osgeo import gdal, osr

gdal.UseExceptions()
BW = os.environ.get("BW", "/mnt/scratch/t/tommy/full/mintpy_bw60")
BIAS = os.path.expanduser(os.environ.get("BIAS", "~/full/mintpy_tay/zheng_bw4/timeseriesBiasApprox.h5"))
TAY = os.path.expanduser("~/full/tay/InSAR_LongTermRate.tif")
BLOCK, NBOOT = 20, 300
rng = np.random.default_rng(5)

with h5py.File(os.path.join(BW, "timeseries_SET.h5"), "r") as f:
    ts = f["timeseries"][:].astype(np.float32); dates = [d.decode() for d in f["date"][:]]; att = dict(f.attrs)
with h5py.File(BIAS, "r") as f:
    bias = f["timeseries"][:].astype(np.float32); bdates = [d.decode() for d in f["date"][:]]
assert dates == bdates, "date lists differ between bw60 series and bias estimate"
with h5py.File(os.path.join(BW, "temporalCoherence.h5"), "r") as f:
    tc = f["temporalCoherence"][:].astype(np.float32)
ne, ny, nx = ts.shape
t = np.array([(dt.datetime.strptime(d, "%Y%m%d") - dt.datetime(2000, 1, 1)).days / 365.25 for d in dates]); tm = t - t.mean()


def vel(series):
    s = series.reshape(ne, -1).astype(np.float64); s -= s.mean(axis=0, keepdims=True)
    return ((tm @ s) / (tm @ tm)).reshape(ny, nx) * 1000.0   # mm/yr, positive toward satellite


v0 = vel(ts); v1 = vel(ts - bias); b = vel(bias)
x0, dx, y0, dy = float(att["X_FIRST"]), float(att["X_STEP"]), float(att["Y_FIRST"]), float(att["Y_STEP"])
mem = gdal.GetDriverByName("MEM").Create("", nx, ny, 1, gdal.GDT_Float32); mem.SetGeoTransform((x0, dx, 0.0, y0, 0.0, dy))
srs = osr.SpatialReference(); srs.ImportFromEPSG(4326); mem.SetProjection(srs.ExportToWkt()); mem.GetRasterBand(1).Fill(float("nan"))
src = gdal.Open(TAY); gdal.ReprojectImage(src, mem, src.GetProjection(), srs.ExportToWkt(), gdal.GRA_Average)
tay = mem.GetRasterBand(1).ReadAsArray().astype(np.float64); tay[tay == 0] = np.nan; tay *= 1000.0
strip = np.isfinite(tay)


def boot(a, c, ij):
    bid = (ij[0] // BLOCK) * 100000 + (ij[1] // BLOCK); uniq, inv = np.unique(bid, return_inverse=True)
    order = np.argsort(inv, kind="stable"); starts = np.searchsorted(inv[order], np.arange(len(uniq))); ends = np.r_[starts[1:], len(order)]
    rs = []
    for _ in range(NBOOT):
        pick = rng.integers(0, len(uniq), len(uniq)); idx = np.concatenate([order[starts[p]:ends[p]] for p in pick])
        if idx.size > 30 and a[idx].std() > 0: rs.append(np.corrcoef(a[idx], c[idx])[0, 1])
    return np.percentile(rs, [2.5, 97.5])


print(f"bias velocity (bandwidth-4 quick estimate) over the strip: mean {np.nanmean(b[strip]):+.1f} std {np.nanstd(b[strip]):.1f} mm/yr")
for thr in (0.5, 0.7):
    m = strip & (tc >= thr) & np.isfinite(v0) & np.isfinite(v1) & np.isfinite(tay)
    for lbl, v in (("before Zheng", v0), ("after Zheng bw4", v1)):
        a, c = v[m], tay[m]; lo, hi = boot(a, c, np.nonzero(m))
        print(f"  tempCoh>={thr} n={m.sum():6d} [{lbl:16s}] r={np.corrcoef(a, c)[0,1]:+.3f} CI [{lo:+.3f},{hi:+.3f}]  ours mean {a.mean():+.1f} std {a.std():.1f} | Tay mean {c.mean():+.1f} std {c.std():.1f}")
print("ZAPPLY_DONE")
