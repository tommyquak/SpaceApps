"""Regional mean rates from the bandwidth-4 series after the full Zheng bias,
with block-bootstrap intervals, and stable ground as the null.

The corrected strip mean landed within 2 mm/yr of Tay's; before this is
claimed, (a) the interval on that mean must be computed with spatial blocks,
(b) stable mineral ground must sit at zero after the same correction, and
(c) the dam zone gets its baseline number with the same interval.
"""
import datetime as dt
import json
import os
import numpy as np
import h5py
from osgeo import gdal, osr

gdal.UseExceptions()
W = "/mnt/scratch/t/tommy/full/mintpy_bw4"
BIAS = os.path.join(W, "zheng", "timeseriesBias.h5")
TAY = os.path.expanduser("~/full/tay/InSAR_LongTermRate.tif")
TC20 = os.path.expanduser("~/full/mintpy_tay/temporalCoherence.h5")
SITES = os.path.expanduser("~/full/dam_sites.geojson")
BLOCK, NBOOT = 20, 500
rng = np.random.default_rng(9)

with h5py.File(os.path.join(W, "timeseries_SET.h5"), "r") as f:
    ts = f["timeseries"][:].astype(np.float32); dates = [d.decode() for d in f["date"][:]]; att = dict(f.attrs)
with h5py.File(BIAS, "r") as f:
    bias = f["timeseries"][:].astype(np.float32)
with h5py.File(os.path.join(W, "temporalCoherence.h5"), "r") as f:
    tc = f["temporalCoherence"][:].astype(np.float32)
with h5py.File(TC20, "r") as f:
    stable = f["temporalCoherence"][:] >= 0.7
ne, ny, nx = ts.shape
t = np.array([(dt.datetime.strptime(d, "%Y%m%d") - dt.datetime(2000, 1, 1)).days / 365.25 for d in dates]); tm = t - t.mean()


def vel(series):
    s = series.reshape(ne, -1).astype(np.float64); s -= s.mean(axis=0, keepdims=True)
    return ((tm @ s) / (tm @ tm)).reshape(ny, nx) * 1000.0


v0, v1 = vel(ts), vel(ts - bias)
x0, dx, y0, dy = float(att["X_FIRST"]), float(att["X_STEP"]), float(att["Y_FIRST"]), float(att["Y_STEP"])
mem = gdal.GetDriverByName("MEM").Create("", nx, ny, 1, gdal.GDT_Float32); mem.SetGeoTransform((x0, dx, 0.0, y0, 0.0, dy))
srs = osr.SpatialReference(); srs.ImportFromEPSG(4326); mem.SetProjection(srs.ExportToWkt()); mem.GetRasterBand(1).Fill(float("nan"))
src = gdal.Open(TAY); gdal.ReprojectImage(src, mem, src.GetProjection(), srs.ExportToWkt(), gdal.GRA_Average)
tay = mem.GetRasterBand(1).ReadAsArray().astype(np.float64); tay[tay == 0] = np.nan; tay *= 1000.0; strip = np.isfinite(tay)
sites = [(f["properties"].get("site_id"), *f["geometry"]["coordinates"][:2]) for f in json.load(open(SITES))["features"]]
lons = np.array([s[1] for s in sites]); lats = np.array([s[2] for s in sites])
dam = np.zeros((ny, nx), bool)
r0, r1 = int((lats.max() - y0) / dy) - 11, int((lats.min() - y0) / dy) + 11; c0, c1 = int((lons.min() - x0) / dx) - 11, int((lons.max() - x0) / dx) + 11
dam[max(r0, 0):r1, max(c0, 0):c1] = True


def block_mean_ci(a, ij):
    bid = (ij[0] // BLOCK) * 100000 + (ij[1] // BLOCK); uniq, inv = np.unique(bid, return_inverse=True)
    bm = np.array([a[inv == k].mean() for k in range(len(uniq))]); w = np.array([(inv == k).sum() for k in range(len(uniq))], float)
    boots = []
    for _ in range(NBOOT):
        pick = rng.integers(0, len(uniq), len(uniq)); boots.append(np.sum(bm[pick] * w[pick]) / np.sum(w[pick]))
    return np.mean(a), np.percentile(boots, [2.5, 97.5]), len(uniq)


print(f"bandwidth-4 series, {ne} epochs; full Zheng bias applied; rates in mm/yr, positive toward satellite; block {BLOCK} px (~1.8 km)")
for thr in (0.5, 0.7):
    print(f"\n--- tempCoh >= {thr} ---")
    for lbl, m in (("Tay strip", strip), ("dam zone box", dam), ("stable mineral ground (null)", stable)):
        mk = m & (tc >= thr) & np.isfinite(v0) & np.isfinite(v1)
        if mk.sum() < 100:
            print(f"  [{lbl}] n={mk.sum()} too few"); continue
        ij = np.nonzero(mk)
        m0, ci0, nb = block_mean_ci(v0[mk], ij); m1, ci1, _ = block_mean_ci(v1[mk], ij)
        line = f"  [{lbl:30s}] n={mk.sum():7d} blocks={nb:4d}  before {m0:+6.1f} [{ci0[0]:+.1f},{ci0[1]:+.1f}]  after Zheng {m1:+6.1f} [{ci1[0]:+.1f},{ci1[1]:+.1f}]"
        if lbl == "Tay strip":
            mt, cit, _ = block_mean_ci(tay[mk], ij); line += f"  | Tay {mt:+.1f} [{cit[0]:+.1f},{cit[1]:+.1f}]"
        print(line)
print("BW4REG_DONE")
