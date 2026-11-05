import glob, os, re, numpy as np
# Decisive question: is the corruption already present in the PER-PAIR ion
# estimates (so the bug is upstream in subband unwrapping), or does ion_ls
# create it (so the bug is the unweighted least squares)?
B = "/mnt/scratch/t/tommy/full/pairs_ion"
ds = sorted(glob.glob(B + "/*/ion/ion_cal"))
print("pairs with ion_cal:", len(ds))
print("files in first:", sorted(os.listdir(ds[0]))[:12])

def dims(vrt):
    t = open(vrt).read()
    m = re.search(r'rasterXSize="(\d+)"\s+rasterYSize="(\d+)"', t)
    if m: return int(m.group(1)), int(m.group(2))
    m = re.search(r'rasterYSize="(\d+)"\s+rasterXSize="(\d+)"', t)
    return int(m.group(2)), int(m.group(1))

for name in ("filt_ion_80rlks_448alks.ion",):
    rows = []
    for d in ds:
        f = os.path.join(d, name)
        v = f + ".vrt"
        if not (os.path.exists(f) and os.path.exists(v)): continue
        try: W, H = dims(v)
        except Exception: continue
        a = np.fromfile(f, dtype=np.float32)
        if a.size != W * H: continue
        a = a.reshape(H, W)
        g = a[np.isfinite(a)]
        if g.size == 0: continue
        rows.append((float(np.mean(g)), float(g.max() - g.min()),
                     float(np.median(np.abs(np.diff(a, axis=0)))),
                     float(np.median(np.abs(np.diff(a, axis=1))))))
    if not rows: continue
    r = np.array(rows)
    print("\n=== %s : %d pairs, grid %dx%d ===" % (name, len(rows), W, H))
    print("  per-pair mean:     median %10.2f   min %10.2f   max %10.2f rad" % (
        np.median(r[:,0]), r[:,0].min(), r[:,0].max()))
    print("  per-pair ptp:      median %10.2f   max %10.2f rad" % (np.median(r[:,1]), r[:,1].max()))
    print("  |grad| azimuth:    median %10.3f   max %10.3f rad/px" % (np.median(r[:,2]), r[:,2].max()))
    print("  |grad| range:      median %10.3f   max %10.3f rad/px" % (np.median(r[:,3]), r[:,3].max()))
    print("  pairs with |mean| > 1000 rad: %d of %d" % (int((np.abs(r[:,0])>1000).sum()), len(rows)))
    print("  implied TEC (median |mean| / 13.4): %.0f TECU" % (abs(np.median(r[:,0]))/13.4))
