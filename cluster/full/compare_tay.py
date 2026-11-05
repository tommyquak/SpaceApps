import os, sys, glob, subprocess, numpy as np, h5py
from workflow_provenance import ProvenanceError, select_velocity
from osgeo import gdal, osr
gdal.UseExceptions()
# Benchmark against Tay et al. 2025 (Comms Earth Environ 6:971).
# COMPARE, NEVER CALIBRATE: nothing is fitted, scaled or offset-corrected.
# This is a SAME-SENSOR benchmark, not independent validation.
W   = os.environ.get("W", os.path.expanduser("~/full/mintpy_tay"))
TAY = os.path.expanduser("~/full/tay/InSAR_LongTermRate.tif")

# guard: never read products a running job may still be writing
# SKIP_RUN_GUARD=1 for products the ALOS-2 MintPy job cannot be writing (PALSAR-1 dirs).
if os.environ.get("SKIP_RUN_GUARD") != "1":
    try:
        queue = subprocess.run(["squeue", "-u", "tommy", "-h", "-n", "full-mintpy,zhapply"],
                               capture_output=True, text=True, check=True)
    except (OSError, subprocess.CalledProcessError):
        sys.exit("cannot verify scheduler state; use SKIP_RUN_GUARD=1 only for a known idle local copy")
    if queue.stdout.strip():
        sys.exit("MintPy or Zheng apply still running; refusing to read partial output")

# Provenance (audit r5 finding 3): velocityClosureCor.h5 is PRESERVED across
# runs, so its mere existence does not mean the current run produced it. A later
# failed or uncorrected run would otherwise be benchmarked using an older job's
# corrected velocity. Require the stamp mintpy2.sbatch writes, and say which job
# made it.
try:
    vf, record = select_velocity(W, allow_unverified=os.environ.get("ALLOW_UNVERIFIED_PRODUCTS") == "1")
except (OSError, ProvenanceError, KeyError) as exc:
    sys.exit(str(exc))
if record.get("status") == "UNVERIFIED_DIAGNOSTIC":
    prov = "UNVERIFIED DIAGNOSTIC: legacy, stale or incomplete provenance"
    print("WARNING:", record["reason"])
else:
    kind = "closure-corrected" if vf.name == "velocityClosureCor.h5" else "UNCORRECTED"
    prov = "%s, recorded %s by job %s" % (kind, record["created_utc"], record["job_id"])
print("our velocity : %s  [%s]" % (os.path.basename(vf), prov))

f = h5py.File(vf,"r"); vel = f["velocity"][:]; a = dict(f.attrs)
unit = a.get("UNIT", "")
unit = unit.decode() if isinstance(unit, bytes) else str(unit)
if unit not in ("m/year", "m/yr"):
    sys.exit("velocity UNIT must be m/year; refusing to infer a conversion")
platform = a.get("PLATFORM", "unknown")
platform = platform.decode() if isinstance(platform, bytes) else str(platform)
if platform.upper().replace("_", "-") == "ALOS-2":
    print("NOTE: same-sensor benchmark (both ALOS-2), NOT independent validation.\n")
else:
    print("NOTE: %s compared with Tay ALOS-2; sensor, epoch and reference differences remain.\n" % platform)
H, Wd = vel.shape
y0,ys = float(a["Y_FIRST"]), float(a["Y_STEP"])
x0,xs = float(a["X_FIRST"]), float(a["X_STEP"])

tc = None
p = os.path.join(W,"temporalCoherence.h5")
if os.path.exists(p):
    with h5py.File(p,"r") as g:
        if "temporalCoherence" not in g:
            sys.exit("temporalCoherence.h5 has no temporalCoherence dataset")
        tc = g["temporalCoherence"][:]
    if tc.shape != vel.shape:
        sys.exit("temporal coherence grid does not match velocity")

# LOS -> vertical, per Tay: divide by cos(local incidence)
inc = None
gp = os.path.join(W,"inputs","geometryGeo.h5")
if os.path.exists(gp):
    gh = h5py.File(gp,"r")
    for k in ("incidenceAngle","incidence_angle"):
        if k in gh: inc = gh[k][:]; break
if inc is None:
    sys.exit("no incidenceAngle; refusing to guess the LOS->vertical projection")
if inc.shape != vel.shape:
    sys.exit("incidence-angle grid does not match velocity")
cosi = np.cos(np.deg2rad(inc)); cosi[~np.isfinite(cosi)|(cosi<0.3)] = np.nan
vert = vel/cosi

# Resample Tay onto our exact grid.
# CRITICAL: their raster declares NO NoData and is ~93% zero-fill. Without
# srcNodata=0 bilinear blends zeros into every valid cell near a boundary and
# biases the comparison toward zero. Warp a validity mask the same way and
# require full validity, so no interpolated cell is contaminated.
def warp(src, alg, srcnd):
    m = gdal.GetDriverByName("MEM").Create("", Wd, H, 1, gdal.GDT_Float32)
    m.SetGeoTransform((x0,xs,0,y0,0,ys))
    srs = osr.SpatialReference(); srs.ImportFromEPSG(4326)
    m.SetProjection(srs.ExportToWkt())
    gdal.Warp(m, src, resampleAlg=alg, srcNodata=srcnd, dstNodata=np.nan)
    return m.GetRasterBand(1).ReadAsArray().astype(np.float64)

tay = warp(TAY, "bilinear", 0)
src = gdal.Open(TAY)
mk = gdal.GetDriverByName("MEM").Create("", src.RasterXSize, src.RasterYSize, 1, gdal.GDT_Float32)
mk.SetGeoTransform(src.GetGeoTransform()); mk.SetProjection(src.GetProjection())
mk.GetRasterBand(1).WriteArray((src.GetRasterBand(1).ReadAsArray()!=0).astype(np.float32))
valid = warp(mk, "average", None)          # 1.0 only where every contributor was valid

m = np.isfinite(vert) & np.isfinite(tay) & (valid > 0.999) & (vert != 0)
print("Tay nonzero cells on our grid : %d (%.2f%%)" % ((valid>0.999).sum(), 100.0*(valid>0.999).mean()))
print("our grid shape                : %d rows x %d columns" % (H, Wd))
print("comparable after both masks   : %d px (%.2f%%)\n" % (m.sum(), 100.0*m.mean()))
if m.sum() < 100: sys.exit("too few comparable pixels")

def report(mask, label):
    o = vert[mask]*1000.0; t = tay[mask]*1000.0
    if o.size < 100: print("  %s: too few px (%d)" % (label,o.size)); return
    r  = np.corrcoef(o,t)[0,1]
    d  = o-t
    dd = (o-o.mean())-(t-t.mean())          # demeaned: datums differ, see below
    print("%s  n=%d" % (label, o.size))
    print("   ours  median %+7.2f  mean %+7.2f  std %6.2f mm/yr" % (np.median(o),o.mean(),o.std()))
    print("   Tay   median %+7.2f  mean %+7.2f  std %6.2f mm/yr" % (np.median(t),t.mean(),t.std()))
    print("   Pearson r        %+.4f" % r)
    print("   bias (ours-Tay)  %+.2f mm/yr   [NOT interpretable: different reference datums]" % d.mean())
    print("   RMSE raw         %7.2f mm/yr   [likewise datum-contaminated]" % np.sqrt((d**2).mean()))
    print("   RMSE demeaned    %7.2f mm/yr   <- the defensible number" % np.sqrt((dd**2).mean()))
    # block bootstrap: InSAR noise is spatially correlated
    ys_,xs_ = np.where(mask); bid=(ys_//40)*10000+(xs_//40); ub=np.unique(bid)
    rng=np.random.default_rng(7); rs=[]
    for _ in range(200):
        pick=rng.choice(ub,size=min(len(ub),150),replace=True)
        idx=np.concatenate([np.where(bid==b)[0] for b in pick])
        if idx.size>500: rs.append(np.corrcoef(o[idx],t[idx])[0,1])
    if rs: print("   block-bootstrap r %+.4f  [%.4f, %.4f] 95%%" % (
        np.mean(rs), np.percentile(rs,2.5), np.percentile(rs,97.5)))
    print()

report(m, "ALL comparable pixels")
if tc is not None:
    for thr in (0.3, 0.5, 0.7):
        report(m & (tc > thr), "temporalCoherence > %.1f" % thr)
print("Reference datums differ (Tay: CPKY GNSS, stable <1 mm/yr; ours: unverified")
print("high-ground pixel). Absolute bias/RMSE are therefore not interpretable;")
print("correlation and demeaned residuals are the defensible measures.")
