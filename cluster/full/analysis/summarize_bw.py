"""Summary of a MintPy work dir against the Tay strip and Tay's sigma-phi mask."""
import os
import numpy as np
import h5py
from osgeo import gdal, osr

W = os.environ["W"]
TAY = os.path.expanduser("~/full/tay/InSAR_LongTermRate.tif")
MASK = os.path.expanduser("~/full/sigphi/mask_tayrule_nominal.npy")


def h5(name, ds=None):
    with h5py.File(os.path.join(W, name), "r") as f:
        if ds is None or ds not in f:
            ds = list(f.keys())[0]
        return f[ds][:].astype(np.float64), dict(f.attrs)


tc, att = h5("temporalCoherence.h5")
v, _ = h5("velocity.h5", "velocity"); v *= 1000.0
ny, nx = tc.shape
x0 = float(att["X_FIRST"]); dx = float(att["X_STEP"]); y0 = float(att["Y_FIRST"]); dy = float(att["Y_STEP"])
mem = gdal.GetDriverByName("MEM").Create("", nx, ny, 1, gdal.GDT_Float32)
mem.SetGeoTransform((x0, dx, 0.0, y0, 0.0, dy))
srs = osr.SpatialReference(); srs.ImportFromEPSG(4326); mem.SetProjection(srs.ExportToWkt())
mem.GetRasterBand(1).Fill(float("nan"))
src = gdal.Open(TAY)
gdal.ReprojectImage(src, mem, src.GetProjection(), srs.ExportToWkt(), gdal.GRA_Average)
tay = mem.GetRasterBand(1).ReadAsArray().astype(np.float64); tay[tay == 0] = np.nan
tay_ok = np.isfinite(tay)
print(f"W={W}")
print(f"tempCoh median {np.nanmedian(tc):.4f}  p90 {np.nanpercentile(tc,90):.4f}  p99 {np.nanpercentile(tc,99):.4f}")
for t in (0.3, 0.5, 0.7):
    m = tc >= t
    print(f"  tempCoh>={t}: {m.sum():8d} px ({100*m.mean():5.2f}%)  in Tay strip: {(m & tay_ok).sum():6d}")
    if (m & tay_ok).sum() > 100:
        mm = m & tay_ok & np.isfinite(v)
        print(f"     r(ours,Tay)={np.corrcoef(v[mm], tay[mm])[0,1]:+.3f}  our std {v[mm].std():.2f}  Tay std {1000*tay[mm].std():.2f} mm/yr")
if os.path.exists(MASK):
    k = np.load(MASK) & tay_ok & np.isfinite(v)
    print(f"Tay-rule(nominal) AND Tay data: n={k.sum()}  r={np.corrcoef(v[k], tay[k])[0,1]:+.3f}  our std {v[k].std():.2f}")
print("SUMMARY_DONE")
