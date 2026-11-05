import h5py, numpy as np, os, subprocess, sys
W = os.environ.get("MINTPY_DIR", os.path.expanduser("~/full/mintpy_tay"))
# Guard: never read products a running job may still be writing. This session
# has twice drawn wrong conclusions from partially-written HDF5.
q = subprocess.run(["squeue","-u","tommy","-h","-n","full-mintpy"],
                   capture_output=True, text=True).stdout.strip()
if q:
    sys.exit("mintpy still running; refusing to read partial output")

tc = h5py.File(W+"temporalCoherence.h5","r")
k = [x for x in tc.keys()][0]; T = tc[k][:]
H, Wd = T.shape
g = h5py.File(W+"inputs/geometryGeo.h5","r")
hgt = g["height"][:]
a = dict(h5py.File(W+"inputs/ifgramStack.h5","r").attrs)
y0, ys = float(a["Y_FIRST"]), float(a["Y_STEP"])
x0, xs = float(a["X_FIRST"]), float(a["X_STEP"])

for THR in (0.5, 0.7):
    m = T > THR
    n = int(m.sum())
    print("\n=== temporalCoherence > %.1f : %d px (%.3f%%) ===" % (THR, n, 100.0*n/T.size))
    if n == 0: continue
    ys_, xs_ = np.where(m)
    lat = y0 + ys*ys_; lon = x0 + xs*xs_
    print("  lat %.3f .. %.3f   lon %.3f .. %.3f" % (lat.min(), lat.max(), lon.min(), lon.max()))
    # Is it one clump or spread? Count occupied 20x20-pixel cells (~1.8 km).
    cells = set(zip((ys_//20).tolist(), (xs_//20).tolist()))
    tot_cells = (H//20+1)*(Wd//20+1)
    print("  occupies %d of %d ~1.8km cells (%.2f%% of scene area)" % (
        len(cells), tot_cells, 100.0*len(cells)/tot_cells))
    # Elevation: the peat dome is low ground; the eastern hills are high.
    e = hgt[m]
    allq = np.percentile(hgt[np.isfinite(hgt)], [50, 90])
    print("  elevation of good px: median %.1f m (scene median %.1f, p90 %.1f)" % (
        np.median(e), allq[0], allq[1]))
    frac_hi = float((e > allq[1]).mean())
    print("  fraction on high ground (>scene p90): %.1f%%" % (100*frac_hi))
    verdict = ("CLUSTERED ON HILLS - little value for peat" if frac_hi > 0.5
               else "spread across low ground - usable for peat" if len(cells) > 50
               else "sparse and localised - point-wise use only")
    print("  -> %s" % verdict)

# quadrant spread
m = T > 0.5
if m.sum():
    qy, qx = H//2, Wd//2
    print("\nquadrant counts (>0.5):  NW %d  NE %d  SW %d  SE %d" % (
        int(m[:qy,:qx].sum()), int(m[:qy,qx:].sum()),
        int(m[qy:,:qx].sum()), int(m[qy:,qx:].sum())))
