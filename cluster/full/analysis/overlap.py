"""Where does Tay have data, and where do we have coherence? Are they the same pixels?"""
import os
import numpy as np
import h5py
from osgeo import gdal, osr

W = os.path.expanduser("~/full/mintpy_tay")
TAY = os.path.expanduser("~/full/tay/InSAR_LongTermRate.tif")

with h5py.File(os.path.join(W, "velocity.h5"), "r") as f:
    att = dict(f.attrs)
    pre = f["velocity"][:].astype(np.float64) * 1000.0
with h5py.File(os.path.join(W, "temporalCoherence.h5"), "r") as f:
    tc = f["temporalCoherence"][:].astype(np.float64)

x0 = float(att["X_FIRST"]); dx = float(att["X_STEP"])
y0 = float(att["Y_FIRST"]); dy = float(att["Y_STEP"])
ny, nx = pre.shape

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
print(f"grid {ny} x {nx} = {tay.size} px at ~92 m")
print(f"Tay has data on          {tay_ok.sum():8d} px  ({100*tay_ok.mean():5.2f}%)")
print(f"Tay std (native units)   {np.nanstd(tay):.5f}   -> {1000*np.nanstd(tay):.2f} mm/yr if m/yr\n")

for thr in (0.3, 0.4, 0.5, 0.6, 0.7):
    ours = tc >= thr
    both = ours & tay_ok
    pct_of_ours = 100.0 * both.sum() / max(ours.sum(), 1)
    print(f"tempCoh>={thr}:  ours {ours.sum():7d} px | overlap with Tay {both.sum():6d} px "
          f"| {pct_of_ours:5.2f}% of our reliable px have a Tay value")

# Is the disjointness geographic? Compare centroids.
def centroid(mask):
    ij = np.nonzero(mask)
    if ij[0].size == 0:
        return None
    lat = y0 + ij[0].mean() * dy
    lon = x0 + ij[1].mean() * dx
    return lat, lon

c_tay = centroid(tay_ok)
c_our = centroid(tc >= 0.7)
print(f"\nTay data centroid          {c_tay}")
print(f"our tempCoh>=0.7 centroid  {c_our}")
