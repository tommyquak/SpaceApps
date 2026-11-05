"""Build MintPy ifgramStack.h5 + geometryGeo.h5 for the full baseline.

Sources, unified onto one 3 arcsec (~90 m, Tay-comparable) grid by averaged
warp (proper multilooking of unwrapped phase; zeros excluded as nodata):
Single source: ~/full/pairs/<R>_<S>/*.geo.aoi90.tif (the MAIN-pass output;
pilot and pairwise banks excluded from the corrected stack by decision,
pilot los kept only as a geometry fallback). Perpendicular baselines:
per-pair isce.log when present, else baseline/baseline_center.txt
(refdate secdate Bpar Bperp; reference seeded 0; bperp = b(sec)-b(ref)).
"""
import glob, os, re, sys
import numpy as np
import h5py
from osgeo import gdal
gdal.UseExceptions()

OUT_DIR = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/full/mintpy")
os.makedirs(OUT_DIR + "/inputs", exist_ok=True)

RES = 0.000833333333333333
BOUNDS = (113.8, -3.4, 115.2, -2.34)

# PAIR_ROOT selects which MAIN pass feeds the stack. Defaults to the Tay-faithful
# 17x56 output; override with PAIR_ROOT=... to rebuild from an earlier pass.
PAIR_ROOT = os.environ.get("PAIR_ROOT", os.path.expanduser("~/full/pairs_tay"))

def find_pairs():
    out = {}
    # alosStack finale: discover ONLY the MAIN-pass output root. The pairwise
    # bank (relocated to ~/full/pairs_pairwise_banked) and the pilot pairs are
    # cross-check material and must never enter the corrected stack.
    for root, unw_pat, cor_pat in [
        (PAIR_ROOT, "filt_*.unw.geo.aoi90.tif", "*.cor.geo.aoi90.tif"),
    ]:
        for d in sorted(glob.glob(root + "/*/")):
            if not os.path.exists(d + "DONE"):
                continue
            unw = [u for u in glob.glob(d + unw_pat) if "_msk" not in u and "conncomp" not in u]
            cor = glob.glob(d + cor_pat)
            ccc = glob.glob(d + "filt_*.unw.conncomp.geo.aoi90.tif")
            if unw and cor:
                out[os.path.basename(d.rstrip("/"))] = (unw[0], cor[0], d, ccc[0] if ccc else None)
    return out

pairs = find_pairs()
names = sorted(pairs)
print(f"{len(names)} completed pairs (MAIN pass)")
if len(names) < 100:
    sys.exit("far too few pairs, aborting for review")

def read_on_grid(path, band="last", alg="average", nodata=True):
    kw = dict(srcNodata=0, dstNodata=0) if nodata else {}
    ds = gdal.Warp("", path, format="MEM", outputBounds=BOUNDS,
                   xRes=RES, yRes=RES, targetAlignedPixels=True,
                   resampleAlg=alg, **kw)
    b = ds.GetRasterBand(ds.RasterCount if band == "last" else 1)
    return b.ReadAsArray(), ds.GetGeoTransform(), ds.RasterXSize, ds.RasterYSize

_, gt, xsize, ysize = read_on_grid(pairs[names[0]][0])  # probe grid
print(f"common grid: {xsize} x {ysize} @ {gt[1]:.6f} deg")

BPRE = re.compile(r"perpendicular baseline at center of reference track\s*=\s*(-?[\d.]+)")
_BASELINES = None
def _load_baselines():
    # alosStack bridge (audit r4, format corrected r5): fall back to
    # baseline/baseline_center.txt when no isce.log exists; values relative to
    # the stack reference; pair bperp = b(sec)-b(ref). Row format documented
    # at the parse site below.
    global _BASELINES
    if _BASELINES is not None:
        return _BASELINES
    # Format (compute_baseline.py): 2 header lines, then rows
    # "refdate secdate Bpar Bperp"; the reference date has NO row (seeded 0).
    _BASELINES = {}
    p = os.path.expanduser("~/full/baseline/baseline_center.txt")
    if os.path.exists(p):
        ref = None
        for ln in open(p):
            parts = ln.split()
            if len(parts) == 4 and parts[0].isdigit() and parts[1].isdigit():
                try:
                    _BASELINES[parts[1]] = float(parts[3])
                    ref = parts[0]
                except ValueError:
                    pass
        if ref is not None:
            _BASELINES.setdefault(ref, 0.0)
    return _BASELINES
def bperp_of(pair_dir):
    log = pair_dir + "isce.log"
    if os.path.exists(log):
        m = BPRE.search(open(log, errors="ignore").read())
        if m:
            return float(m.group(1))
    bl = _load_baselines()
    name = os.path.basename(pair_dir.rstrip("/"))
    try:
        a, b = [d.replace("-", "")[2:] for d in name.split("_")]  # YYMMDD
        if a in bl and b in bl:
            return bl[b] - bl[a]
    except ValueError:
        pass
    return None

n = len(names)
date12, bperps, missing_bp, n_realcc = [], [], 0, 0
unw_stack = np.zeros((n, ysize, xsize), dtype=np.float32)
coh_stack = np.zeros((n, ysize, xsize), dtype=np.float32)
cc_stack = np.zeros((n, ysize, xsize), dtype=np.float32)
for i, name in enumerate(names):
    unw, cor, d, ccc = pairs[name]
    d0, d1 = name.split("_")
    date12.append((d0.replace("-", ""), d1.replace("-", "")))
    unw_stack[i] = read_on_grid(unw)[0]
    coh_stack[i] = read_on_grid(cor)[0]
    if ccc is not None:
        # real snaphu components: nearest, labels never blended, 0 is a class
        cc_stack[i] = read_on_grid(ccc, band="first", alg="near", nodata=False)[0]
        n_realcc += 1
    else:
        cc_stack[i] = (coh_stack[i] > 0).astype(np.float32)
    bp = bperp_of(d)
    if bp is None:
        missing_bp += 1; bp = 0.0
    bperps.append(bp)
    if i % 100 == 0:
        print(f"  loaded {i}/{n}")
print(f"bperp parsed for {n - missing_bp}/{n} pairs (missing set to 0)")
print(f"real snaphu conncomp for {n_realcc}/{n} pairs (rest coherence>0 fallback)")

def _looks_from_products(paths):
    """Read the ml2 looks straight out of the alosStack product filenames."""
    found = set()
    for q in paths:
        m = re.search(r"_(\d+)rlks_(\d+)alks", os.path.basename(q))
        if m:
            found.add((int(m.group(1)), int(m.group(2))))
    if len(found) == 1:
        return found.pop()
    if not found:
        raise SystemExit("build_stack2: no _<R>rlks_<A>alks in product names; "
                         "cannot determine looks, refusing to guess")
    raise SystemExit("build_stack2: inconsistent looks across products: %s" % sorted(found))

RLOOKS, ALOOKS = _looks_from_products([pairs[k][0] for k in names])
print("looks parsed from products: RLOOKS=%d ALOOKS=%d" % (RLOOKS, ALOOKS))

WAVELENGTH = 299792458.0 / 1257.5e6
attrs = {
    "WAVELENGTH": str(WAVELENGTH),
    # PySolid (correct_SET) needs the acquisition time of day. Without it
    # solid_earth_tides.py raises KeyError AFTER the inversion completes.
    # 04:02:56.95 UTC start + half the 16.2 s scene duration.
    "CENTER_LINE_UTC": "14585.1",
    # closure_phase_bias (Zheng et al. 2022) needs HEIGHT; utils0.incidence_angle
    # short-circuits on CENTER_INCIDENCE_ANGLE. Both derived from geometryGeo.
    "HEIGHT": "628000.0",
    "CENTER_INCIDENCE_ANGLE": "40.1121",
    "EARTH_RADIUS": "6371000.0",
    # ALOOKS_DERIVED: parsed from the product names so it can never drift from
    # what was actually processed. alosStack encodes ml2 as _<R>rlks_<A>alks.
    "ALOOKS": str(ALOOKS), "RLOOKS": str(RLOOKS),
    "LENGTH": str(ysize), "WIDTH": str(xsize),
    "X_FIRST": str(gt[0]), "X_STEP": str(gt[1]),
    "Y_FIRST": str(gt[3]), "Y_STEP": str(gt[5]),
    "X_UNIT": "degrees", "Y_UNIT": "degrees",
    "EPSG": "4326",
    "PLATFORM": "ALOS-2", "PROCESSOR": "isce",
    "FILE_TYPE": "ifgramStack", "UNIT": "radian",
    "PRODUCT_MATURITY": "BASELINE",
    "ANTENNA_SIDE": "-1", "ORBIT_DIRECTION": "DESCENDING",
    "processing_note": "alosStack HH swaths 2-4, bandwidth-20, stack-level ion correction, 3 arcsec averaged grid, bperp from isce.log or baseline_center fallback",
}
p = OUT_DIR + "/inputs/ifgramStack.h5"
with h5py.File(p, "w") as f:
    f.create_dataset("unwrapPhase", data=unw_stack, compression="lzf")
    f.create_dataset("coherence", data=coh_stack, compression="lzf")
    f.create_dataset("connectComponent", data=cc_stack, compression="lzf")
    f.create_dataset("date", data=np.array(date12, dtype="S8"))
    f.create_dataset("bperp", data=np.array(bperps, dtype=np.float32))
    f.create_dataset("dropIfgram", data=np.ones(n, dtype=bool))
    for k, v in attrs.items():
        f.attrs[k] = v
print("wrote", p)

los_cands = ([u for u in glob.glob(os.path.expanduser("~/full/pairs/*/") + "*los.geo.aoi90.tif")] or
             [u for u in glob.glob(os.path.expanduser("~/pilot/pairs/*/") + "*los.geo.aoi.tif")])
inc = read_on_grid(los_cands[0], band="first", alg="near")[0]
ds = gdal.Warp("", los_cands[0], format="MEM", outputBounds=BOUNDS, xRes=RES, yRes=RES,
               targetAlignedPixels=True, resampleAlg="near", srcNodata=0, dstNodata=0)
azi = ds.GetRasterBand(ds.RasterCount).ReadAsArray() if ds.RasterCount > 1 else np.zeros_like(inc)
hgt = read_on_grid(os.path.expanduser("~/pilot/dem/elevation.dem.wgs84"), band="first")[0]
g = OUT_DIR + "/inputs/geometryGeo.h5"
gattrs = dict(attrs); gattrs["FILE_TYPE"] = "geometry"; gattrs["UNIT"] = "degrees"
# slantRangeDistance from spherical geometry (Re=6371 km, ALOS-2 h=628 km);
# dem_error reads it unconditionally once incidenceAngle exists. A smooth
# approximation suffices: the least-squares DEM-error fit absorbs scale error.
Re, hs = 6371000.0, 628000.0
inc_r = np.deg2rad(np.where(inc > 0, inc, np.nan))
gamma = np.arcsin(Re / (Re + hs) * np.sin(inc_r))
srng = (Re + hs) * np.sin(inc_r - gamma) / np.sin(inc_r)
fill = np.nanmean(srng) if np.isfinite(srng).any() else 850000.0
srng = np.where(np.isfinite(srng), srng, fill).astype(np.float32)
with h5py.File(g, "w") as f:
    f.create_dataset("height", data=hgt.astype(np.float32), compression="lzf")
    f.create_dataset("slantRangeDistance", data=srng, compression="lzf")
    f.create_dataset("incidenceAngle", data=inc.astype(np.float32), compression="lzf")
    f.create_dataset("azimuthAngle", data=azi.astype(np.float32), compression="lzf")
    for k, v in gattrs.items():
        f.attrs[k] = v
print("wrote", g)
v = inc[inc > 0]
print(f"incidence range: {v.min():.1f} to {v.max():.1f} deg" if v.size else "incidence EMPTY")
