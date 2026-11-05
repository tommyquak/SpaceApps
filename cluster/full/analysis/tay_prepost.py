"""Does the Zheng closure-bias correction move our velocity toward Tay or away?

Compares corr(ours, Tay) before and after the correction on the SAME pixels,
with a block bootstrap so the CI respects InSAR's spatial correlation.
Tay is an independent benchmark: we compare, we never calibrate to it.
"""
import os
import numpy as np
import h5py
from osgeo import gdal

W = os.path.expanduser(os.environ.get("W", "~/full/mintpy_tay"))   # any MintPy dir with velocity + velocityClosureCor
TAY = os.path.expanduser("~/full/tay/InSAR_LongTermRate.tif")
BLOCK = 20          # pixels per bootstrap block side (~1.8 km)
NBOOT = 400
rng = np.random.default_rng(12345)


def h5(name, ds):
    with h5py.File(os.path.join(W, name), "r") as f:
        if ds not in f:
            ds = list(f.keys())[0]
        return f[ds][:].astype(np.float64), dict(f.attrs)


pre, att = h5("velocity.h5", "velocity")
post, _ = h5("velocityClosureCor.h5", "velocity")
tc, _ = h5("temporalCoherence.h5", "temporalCoherence")
pre *= 1000.0
post *= 1000.0

x0 = float(att["X_FIRST"]); dx = float(att["X_STEP"])
y0 = float(att["Y_FIRST"]); dy = float(att["Y_STEP"])
ny, nx = pre.shape

# Resample Tay onto our exact grid.
out_gt = (x0, dx, 0.0, y0, 0.0, dy)
mem = gdal.GetDriverByName("MEM").Create("", nx, ny, 1, gdal.GDT_Float32)
mem.SetGeoTransform(out_gt)
srs = gdal.osr.SpatialReference(); srs.ImportFromEPSG(4326)
mem.SetProjection(srs.ExportToWkt())
mem.GetRasterBand(1).SetNoDataValue(np.nan)
mem.GetRasterBand(1).Fill(np.nan)
src = gdal.Open(TAY)
gdal.ReprojectImage(src, mem, src.GetProjection(), srs.ExportToWkt(), gdal.GRA_Average)
tay = mem.GetRasterBand(1).ReadAsArray().astype(np.float64)

# Tay's raster declares no NoData, so exact zeros are padding, not data.
tay[tay == 0] = np.nan
print(f"Tay valid pixels on our grid: {int(np.isfinite(tay).sum())} of {tay.size}"
      f"  ({100*np.isfinite(tay).sum()/tay.size:.2f}%)")
print(f"Tay signal std over overlap: {np.nanstd(tay):.2f} mm/yr\n")


def block_boot(a, b, ij, n=NBOOT):
    """Block bootstrap of Pearson r over spatially-correlated pixels."""
    bid = (ij[0] // BLOCK) * 100000 + (ij[1] // BLOCK)
    uniq, inv = np.unique(bid, return_inverse=True)
    order = np.argsort(inv, kind="stable")
    starts = np.searchsorted(inv[order], np.arange(len(uniq)))
    ends = np.r_[starts[1:], len(order)]
    rs = []
    for _ in range(n):
        pick = rng.integers(0, len(uniq), len(uniq))
        idx = np.concatenate([order[starts[p]:ends[p]] for p in pick])
        if idx.size < 30:
            continue
        aa, bb = a[idx], b[idx]
        if aa.std() < 1e-12 or bb.std() < 1e-12:
            continue
        rs.append(np.corrcoef(aa, bb)[0, 1])
    rs = np.array(rs)
    return np.percentile(rs, 2.5), np.percentile(rs, 97.5)


for thr in (0.0, 0.5, 0.7):
    m = np.isfinite(pre) & np.isfinite(post) & np.isfinite(tay)
    if thr > 0:
        m &= tc >= thr
    n = int(m.sum())
    label = "no coh gate" if thr == 0 else f"tempCoh>={thr}"
    if n < 50:
        print(f"[{label}] n={n} too few to compare")
        continue
    ij = np.nonzero(m)
    a, b, t = pre[m], post[m], tay[m]
    r_pre = np.corrcoef(a, t)[0, 1]
    r_post = np.corrcoef(b, t)[0, 1]
    lo1, hi1 = block_boot(a, t, ij)
    lo2, hi2 = block_boot(b, t, ij)
    print(f"[{label}] n={n}")
    print(f"   PRE  Zheng  r = {r_pre:+.3f}   CI [{lo1:+.3f}, {hi1:+.3f}]")
    print(f"   POST Zheng  r = {r_post:+.3f}   CI [{lo2:+.3f}, {hi2:+.3f}]")
    print(f"   change      dr = {r_post - r_pre:+.3f}")
    print()
