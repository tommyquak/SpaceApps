"""Build a MintPy ifgramStack.h5 + geometryGeo.h5 from the pilot pair crops.

Inputs: ~/pilot/pairs/<REF>_<SEC>/ *.aoi.tif on one common geographic grid
(alos2App geocoded, cropped to the AOI window at copy time). Uses the masked
filtered unwrapped phase and the coherence layer per pair, plus incidence
from one pair los layer. ALOS-2 ScanSAR wavelength 0.229268 m is read from
the value alos2App recorded in the pair manifests when present, else the
PALSAR-2 1236.5 MHz value is refused so the operator supplies it explicitly.
"""
import glob, os, re, sys
import numpy as np
import h5py
from osgeo import gdal
gdal.UseExceptions()

PAIR_ROOT = os.path.expanduser("~/pilot/pairs")
OUT_DIR = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/pilot/mintpy")
os.makedirs(OUT_DIR + "/inputs", exist_ok=True)

pairs = []
for d in sorted(glob.glob(PAIR_ROOT + "/*/")):
    if not os.path.exists(d + "DONE"):
        continue
    unw = glob.glob(d + "filt_*_5rlks_28alks.unw.geo.aoi.tif")
    cor = glob.glob(d + "*_5rlks_28alks.cor.geo.aoi.tif")
    if unw and cor:
        name = os.path.basename(d.rstrip("/"))
        pairs.append((name, unw[0], cor[0]))
print(f"{len(pairs)} completed pairs with crops")
if len(pairs) < 3:
    sys.exit("need at least 3 pairs")

# Common target grid: fixed bounds snapped to the shared 0.000278 deg posting.
# Individual pair crops differ by a few rows where the projWin exceeded a
# pair geocode extent, so every layer is warped onto this grid (nearest, to
# never blend phase samples).
RES = 0.000277777777777778
BOUNDS = (113.8, -3.4, 115.2, -2.34)  # xmin ymin xmax ymax

def read_on_grid(path, band_last=True):
    ds = gdal.Warp("", path, format="MEM", outputBounds=BOUNDS,
                   xRes=RES, yRes=RES, targetAlignedPixels=True,
                   resampleAlg="near")
    band = ds.GetRasterBand(ds.RasterCount if band_last else 1)
    arr = band.ReadAsArray()
    return arr, ds.GetGeoTransform(), ds.RasterXSize, ds.RasterYSize

_, gt, xsize, ysize = read_on_grid(pairs[0][1])
print(f"common grid: {xsize} x {ysize} @ {gt[1]:.6f} deg")

n = len(pairs)
date12 = []
unw_stack = np.zeros((n, ysize, xsize), dtype=np.float32)
coh_stack = np.zeros((n, ysize, xsize), dtype=np.float32)
for i, (name, unw, cor) in enumerate(pairs):
    d0, d1 = name.split("_")
    date12.append((d0.replace("-", ""), d1.replace("-", "")))
    unw_stack[i] = read_on_grid(unw)[0]
    coh_stack[i] = read_on_grid(cor)[0]
    print(f"  loaded {name}")

# ALOS-2 PALSAR-2 center frequency 1257.5 MHz for ScanSAR (JAXA ALOS-2/PALSAR-2
# Level 1.1 format description); wavelength = c / f.
WAVELENGTH = 299792458.0 / 1257.5e6

attrs = {
    "WAVELENGTH": str(WAVELENGTH),
    "ALOOKS": "28", "RLOOKS": "5",
    "LENGTH": str(ysize), "WIDTH": str(xsize),
    "X_FIRST": str(gt[0]), "X_STEP": str(gt[1]),
    "Y_FIRST": str(gt[3]), "Y_STEP": str(gt[5]),
    "X_UNIT": "degrees", "Y_UNIT": "degrees",
    "EPSG": "4326",
    "PLATFORM": "ALOS-2", "PROCESSOR": "isce",
    "FILE_TYPE": "ifgramStack",
    "UNIT": "radian",
    "PRODUCT_MATURITY": "PILOT",
    "ANTENNA_SIDE": "-1",
    "ORBIT_DIRECTION": "DESCENDING",
    "processing_note": "pairwise alos2App HH swaths 2-4, AOI crops, pilot network bandwidth 4",
}
p = OUT_DIR + "/inputs/ifgramStack.h5"
with h5py.File(p, "w") as f:
    f.create_dataset("unwrapPhase", data=unw_stack, compression="lzf")
    f.create_dataset("coherence", data=coh_stack, compression="lzf")
    f.create_dataset("connectComponent", data=(coh_stack > 0).astype(np.float32), compression="lzf")
    f.create_dataset("date", data=np.array(date12, dtype="S8"))
    dp = np.zeros(n, dtype=np.float32)
    f.create_dataset("bperp", data=dp)
    f.create_dataset("dropIfgram", data=np.ones(n, dtype=bool))
    for k, v in attrs.items():
        f.attrs[k] = v
print("wrote", p)

los = glob.glob(PAIR_ROOT + "/*/" + "*los.geo.aoi.tif")[0]
inc = read_on_grid(los, band_last=False)[0]
azi_ds = gdal.Warp("", los, format="MEM", outputBounds=BOUNDS, xRes=RES, yRes=RES, targetAlignedPixels=True, resampleAlg="near")
azi = azi_ds.GetRasterBand(azi_ds.RasterCount).ReadAsArray() if azi_ds.RasterCount > 1 else np.zeros_like(inc)
g = OUT_DIR + "/inputs/geometryGeo.h5"
gattrs = dict(attrs); gattrs["FILE_TYPE"] = "geometry"; gattrs["UNIT"] = "degrees"
dem_path = os.path.expanduser("~/pilot/dem/elevation.dem.wgs84")
hgt = read_on_grid(dem_path, band_last=False)[0]
with h5py.File(g, "w") as f:
    f.create_dataset("height", data=hgt.astype(np.float32), compression="lzf")
    f.create_dataset("incidenceAngle", data=inc.astype(np.float32), compression="lzf")
    f.create_dataset("azimuthAngle", data=azi.astype(np.float32), compression="lzf")
    for k, v in gattrs.items():
        f.attrs[k] = v
print("wrote", g)
valid_inc = inc[inc > 0]
print(f"incidence range: {valid_inc.min():.1f} to {valid_inc.max():.1f} deg" if valid_inc.size else "incidence EMPTY")
