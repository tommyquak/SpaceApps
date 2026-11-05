"""Tay et al. 2025 tropospheric correction, implemented as they describe it.

"We reduced the spatially correlated tropospheric noise at each pixel by
utilizing the stacked time series of neighbouring non-deforming pixels":
neighbours with velocity under 1 mm/yr over the whole series, average spatial
coherence above 0.3, within 10 km; the per-epoch mean of at least 5 such
neighbours is subtracted from the target pixel. [E] Tay 2025, Methods.

Implemented as a disc convolution per epoch so it runs in minutes. Reports how
many pixels can be corrected at all, since the criterion needs stable, coherent
neighbours and this frame has few. Writes timeseries_tropoTay.h5 and
velocity_tropoTay.h5 next to the inputs; never modifies them.
"""
import os
import datetime as dt
import numpy as np
import h5py
from scipy.signal import fftconvolve
from osgeo import gdal, osr

W = os.environ.get("W", os.path.expanduser("~/full/mintpy_tay"))
TAY = os.path.expanduser("~/full/tay/InSAR_LongTermRate.tif")
RADIUS_KM, VMAX_MM, COHMIN, MINN = 10.0, 1.0, 0.3, 5


def h5(name, ds=None):
    with h5py.File(os.path.join(W, name), "r") as f:
        if ds is None or ds not in f:
            ds = [k for k in f.keys() if k not in ("date", "bperp")][0]
        return f[ds][:], dict(f.attrs)


for cand in ("timeseries_closureCor.h5", "timeseries_SET.h5", "timeseries.h5"):
    if os.path.exists(os.path.join(W, cand)):
        TS = cand
        break
with h5py.File(os.path.join(W, TS), "r") as f:
    ts = f["timeseries"][:].astype(np.float32)
    dates = [d.decode() for d in f["date"][:]]
    att = dict(f.attrs)
v, _ = h5("velocity.h5", "velocity"); v = v.astype(np.float64) * 1000.0
coh, _ = h5("avgSpatialCoh.h5"); coh = coh.astype(np.float64)
tc, _ = h5("temporalCoherence.h5"); tc = tc.astype(np.float64)
ne, ny, nx = ts.shape
print(f"W={W}  series={TS}  epochs={ne}  grid={ny}x{nx}")

dy_km = abs(float(att["Y_STEP"])) * 111.32
R = int(round(RADIUS_KM / dy_km))
yy, xx = np.mgrid[-R:R + 1, -R:R + 1]
disc = ((yy * yy + xx * xx) <= R * R).astype(np.float32)
print(f"pixel {dy_km*1000:.0f} m, disc radius {R} px, {int(disc.sum())} px in disc")

stable = np.isfinite(v) & (np.abs(v) < VMAX_MM) & np.isfinite(coh) & (coh > COHMIN)
print(f"non-deforming coherent neighbours available: {stable.sum()} px ({100*stable.mean():.2f}% of frame)")
cnt = fftconvolve(stable.astype(np.float32), disc, mode="same")
ok = cnt >= MINN
print(f"pixels with >= {MINN} qualifying neighbours within {RADIUS_KM:g} km: {ok.sum()} ({100*ok.mean():.2f}%)")

ts_c = ts.copy()
for i in range(ne):
    layer = np.where(stable, ts[i], 0.0).astype(np.float32)
    num = fftconvolve(layer, disc, mode="same")
    tropo = np.where(ok, num / np.maximum(cnt, 1.0), 0.0)
    ts_c[i] = ts[i] - tropo

t = np.array([(dt.datetime.strptime(d, "%Y%m%d") - dt.datetime(2000, 1, 1)).days / 365.25 for d in dates])
tm = t - t.mean()


def fit_vel(series):
    s = series.astype(np.float64)
    s -= s.mean(axis=0, keepdims=True)
    return (np.tensordot(tm, s, axes=(0, 0)) / np.sum(tm * tm)) * 1000.0  # mm/yr


v0 = fit_vel(ts); v1 = fit_vel(ts_c)

with h5py.File(os.path.join(W, "timeseries_tropoTay.h5"), "w") as f:
    f.create_dataset("timeseries", data=ts_c, compression="gzip")
    f.create_dataset("date", data=np.array(dates, dtype="S8"))
    f.attrs.update(att); f.attrs["FILE_TYPE"] = "timeseries"
with h5py.File(os.path.join(W, "velocity_tropoTay.h5"), "w") as f:
    f.create_dataset("velocity", data=(v1 / 1000.0).astype(np.float32))
    f.attrs.update(att); f.attrs["FILE_TYPE"] = "velocity"; f.attrs["UNIT"] = "m/year"

# Tay strip on this grid
x0 = float(att["X_FIRST"]); dx = float(att["X_STEP"]); y0 = float(att["Y_FIRST"]); dy = float(att["Y_STEP"])
mem = gdal.GetDriverByName("MEM").Create("", nx, ny, 1, gdal.GDT_Float32)
mem.SetGeoTransform((x0, dx, 0.0, y0, 0.0, dy))
srs = osr.SpatialReference(); srs.ImportFromEPSG(4326); mem.SetProjection(srs.ExportToWkt())
mem.GetRasterBand(1).Fill(float("nan"))
src = gdal.Open(TAY); gdal.ReprojectImage(src, mem, src.GetProjection(), srs.ExportToWkt(), gdal.GRA_Average)
tay = mem.GetRasterBand(1).ReadAsArray().astype(np.float64); tay[tay == 0] = np.nan
strip = np.isfinite(tay)
print(f"\ncorrected pixels inside Tay strip: {(ok & strip).sum()} of {strip.sum()} ({100*(ok & strip).mean()/max(strip.mean(),1e-9):.1f}%)")

for label, m in (("tempCoh>=0.7", tc >= 0.7), ("tempCoh>=0.5", tc >= 0.5), ("Tay strip", strip)):
    mm = m & np.isfinite(v0) & np.isfinite(v1)
    if mm.sum() < 50:
        print(f"[{label}] n={mm.sum()} too few"); continue
    line = f"[{label}] n={mm.sum()}  corrected {100*np.mean(ok[mm]):.1f}%  vel std before {v0[mm].std():.2f} after {v1[mm].std():.2f} mm/yr"
    if label == "Tay strip":
        line += f"  r vs Tay before {np.corrcoef(v0[mm], tay[mm])[0,1]:+.3f} after {np.corrcoef(v1[mm], tay[mm])[0,1]:+.3f}"
    print(line)
print("TROPO_DONE")
