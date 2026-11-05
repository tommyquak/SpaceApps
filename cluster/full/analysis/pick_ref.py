import h5py, numpy as np
W = os.environ.get("MINTPY_DIR", os.environ.get("MINTPY_DIR", "/home/t/tommy/full/mintpy_tay")) + "/"
frac = np.load("/home/t/tommy/ccfrac.npy")
sc = h5py.File(W + "avgSpatialCoh.h5", "r")["coherence"][:]
H, Wd = frac.shape
# A reference must be reliable in as many interferograms as possible AND sit
# away from the scene edge, where support is one-sided and unwrapping is worst.
MARGIN = 60
m = np.zeros_like(frac, dtype=bool)
m[MARGIN:H-MARGIN, MARGIN:Wd-MARGIN] = True
m &= np.isfinite(sc) & (sc > 0)
fr = np.where(m, frac, -1)
print("interior candidates: %d" % int(m.sum()))
print("interior cc-frac: median %.3f  p99 %.3f  max %.3f" % (
    np.median(frac[m]), np.percentile(frac[m], 99), fr.max()))
# among the top cc-coverage pixels, prefer higher coherence, and require the
# 5x5 neighbourhood to be uniformly good so the choice is not a lone pixel
thr = np.percentile(frac[m], 99.9)
cand = np.argwhere(m & (frac >= thr))
print("candidates at/above p99.9 (%.3f): %d" % (thr, len(cand)))
best = None
for y, x in cand:
    nb = frac[y-2:y+3, x-2:x+3]
    if nb.shape != (5, 5): continue
    score = nb.min() * 0.7 + sc[y, x] * 0.3
    if best is None or score > best[0]:
        best = (score, int(y), int(x), float(frac[y, x]), float(nb.min()), float(sc[y, x]))
print("\nSELECTED reference:")
print("  y=%d x=%d   cc>0 in %.1f%% of ifgs (5x5 worst %.1f%%)  avgSpatialCoh %.4f" % (
    best[1], best[2], 100*best[3], 100*best[4], best[5]))
print("  vs MintPy's choice (409,1677): cc 5.0%%, coh 0.5588")
f = h5py.File(W + "inputs/ifgramStack.h5", "r")
a = f.attrs
y0 = float(a["Y_FIRST"]); ys = float(a["Y_STEP"])
x0 = float(a["X_FIRST"]); xs = float(a["X_STEP"])
print("  lat %.5f  lon %.5f" % (y0 + ys*best[1], x0 + xs*best[2]))
print("  dist from edges: N%d S%d W%d E%d px" % (best[1], H-1-best[1], best[2], Wd-1-best[2]))
open("/home/t/tommy/full/ref_yx.txt", "w").write("%d %d\n" % (best[1], best[2]))
print("\nwrote ~/full/ref_yx.txt")
