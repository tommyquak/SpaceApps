"""In the region where Tay HAS data, how does coherence decay with temporal baseline?

If coherence is high at short baselines and only collapses at long ones, a
bandwidth-restricted network recovers the Tay region and the project is viable.
If it is already dead at the shortest baseline, no network choice can save it.
"""
import os
import datetime as dt
import numpy as np
import h5py
from osgeo import gdal, osr

W = os.path.expanduser("~/full/mintpy_tay")
TAY = os.path.expanduser("~/full/tay/InSAR_LongTermRate.tif")

with h5py.File(os.path.join(W, "velocity.h5"), "r") as f:
    att = dict(f.attrs)
x0 = float(att["X_FIRST"]); dx = float(att["X_STEP"])
y0 = float(att["Y_FIRST"]); dy = float(att["Y_STEP"])

stack = os.path.join(W, "inputs", "ifgramStack.h5")
with h5py.File(stack, "r") as f:
    print("ifgramStack datasets:", list(f.keys()))
    dates = f["date"][:]
    ny, nx = f["coherence"].shape[1:]
    npair = f["coherence"].shape[0]
    print(f"pairs={npair}  grid={ny}x{nx}")

    # Tay mask on the stack grid
    mem = gdal.GetDriverByName("MEM").Create("", nx, ny, 1, gdal.GDT_Float32)
    mem.SetGeoTransform((x0, dx, 0.0, y0, 0.0, dy))
    srs = osr.SpatialReference(); srs.ImportFromEPSG(4326)
    mem.SetProjection(srs.ExportToWkt())
    mem.GetRasterBand(1).SetNoDataValue(float("nan"))
    mem.GetRasterBand(1).Fill(float("nan"))
    src = gdal.Open(TAY)
    gdal.ReprojectImage(src, mem, src.GetProjection(), srs.ExportToWkt(), gdal.GRA_Average)
    tay = mem.GetRasterBand(1).ReadAsArray().astype(np.float64)
    tay[tay == 0] = np.nan
    tay_ok = np.isfinite(tay)
    print(f"Tay mask on stack grid: {tay_ok.sum()} px")
    if tay_ok.sum() == 0:
        raise SystemExit("stack grid does not match velocity grid")

    tb, coh_tay, coh_all = [], [], []
    for i in range(npair):
        d1 = dt.datetime.strptime(dates[i][0].decode(), "%Y%m%d")
        d2 = dt.datetime.strptime(dates[i][1].decode(), "%Y%m%d")
        c = f["coherence"][i]
        tb.append((d2 - d1).days)
        coh_tay.append(float(np.nanmean(c[tay_ok])))
        coh_all.append(float(np.nanmean(c)))

tb = np.array(tb); coh_tay = np.array(coh_tay); coh_all = np.array(coh_all)
print("\n  tbase(d)   n_pairs   mean coh in TAY region   mean coh whole frame")
edges = [0, 15, 30, 60, 100, 200, 400, 800, 10000]
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (tb > lo) & (tb <= hi)
    if m.sum() == 0:
        continue
    print(f"  {lo:5d}-{hi:<5d} {m.sum():6d}        {np.nanmean(coh_tay[m]):.4f}"
          f"                {np.nanmean(coh_all[m]):.4f}")
print(f"\nshortest baseline present: {tb.min()} d  ({(tb==tb.min()).sum()} pairs)")
print(f"best single-pair coh in Tay region: {np.nanmax(coh_tay):.4f} at tbase {tb[np.nanargmax(coh_tay)]} d")
