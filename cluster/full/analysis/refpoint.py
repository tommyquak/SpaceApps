import h5py, numpy as np
W = os.environ.get("MINTPY_DIR", os.environ.get("MINTPY_DIR", "/home/t/tommy/full/mintpy_tay")) + "/"
f = h5py.File(W + "inputs/ifgramStack.h5", "r")
n, H, Wd = f["connectComponent"].shape
print("stack %d pairs, grid %dx%d" % (n, Wd, H))
# per-pixel fraction of interferograms where snaphu produced a valid component.
# The audit's fix: a reference must be reliable in nearly all interferograms,
# not merely high in avgSpatialCoh.
cnt = np.zeros((H, Wd), dtype=np.int32)
step = 40
for i in range(0, n, step):
    cc = f["connectComponent"][i:i+step]
    cnt += (cc != 0).sum(axis=0).astype(np.int32)
frac = cnt / float(n)
ry, rx = 409, 1677
print("\nCHOSEN reference (%d,%d): cc>0 in %.1f%% of interferograms" % (ry, rx, 100*frac[ry, rx]))
sc = h5py.File(W + "avgSpatialCoh.h5", "r")["coherence"][:]
print("   avgSpatialCoh there: %.4f" % sc[ry, rx])
print("\nbest-by-connectComponent pixels:")
flat = np.argsort(frac.ravel())[::-1][:8]
for k in flat:
    y, x = divmod(int(k), Wd)
    print("   (%4d,%4d) cc>0 %.1f%%  avgSpatialCoh %.4f" % (y, x, 100*frac[y, x], sc[y, x]))
print("\nfrac stats: median %.3f  p99 %.3f  max %.3f" % (
    np.median(frac), np.percentile(frac, 99), frac.max()))
# how edge-y is the chosen point
print("\nchosen ref is %d px from E edge, %d from N edge" % (Wd-1-rx, ry))
np.save("/home/t/tommy/ccfrac.npy", frac)
print("saved ccfrac.npy")
