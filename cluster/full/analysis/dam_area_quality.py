import h5py, numpy as np, os, subprocess, sys
W = os.environ.get("MINTPY_DIR", os.path.expanduser("~/full/mintpy_tay"))
if subprocess.run(["squeue","-u","tommy","-h","-n","full-mintpy"],
                  capture_output=True,text=True).stdout.strip():
    print("NOTE: mintpy running; reading temporalCoherence only (written by the")
    print("completed inversion, not touched by the closure step).")
f=h5py.File(W+"temporalCoherence.h5","r"); T=f[[k for k in f.keys()][0]][:]
a=dict(h5py.File(W+"inputs/ifgramStack.h5","r").attrs)
y0,ys=float(a["Y_FIRST"]),float(a["Y_STEP"]); x0,xs=float(a["X_FIRST"]),float(a["X_STEP"])
H,Wd=T.shape
def box(s,n,w,e,label):
    r0=int((n-y0)/ys); r1=int((s-y0)/ys); c0=int((w-x0)/xs); c1=int((e-x0)/xs)
    r0,r1=max(0,min(r0,r1)),min(H,max(r0,r1)); c0,c1=max(0,min(c0,c1)),min(Wd,max(c0,c1))
    sub=T[r0:r1,c0:c1]
    if sub.size==0: print("  %s: empty"%label); return
    g=sub[np.isfinite(sub)]
    print("%s  (%d x %d px)" % (label, sub.shape[1], sub.shape[0]))
    print("   tempCoh median %.4f  p90 %.4f  max %.4f" % (np.median(g),np.percentile(g,90),g.max()))
    for t in (0.3,0.5,0.7):
        print("     > %.1f : %6d px (%.3f%%)" % (t,int((g>t).sum()),100.0*(g>t).mean()))
# the 25 dam sites, padded
box(-2.72,-2.56,114.25,114.39,"DAM SITE AREA (25 candidate dams)")
print()
box(-2.85,-2.45,114.00,114.60,"WIDER DAM/CANAL ZONE")
print()
box(-3.40,-2.90,114.55,115.20,"SE HILLS (where our good pixels are)")
print()
box(-3.40,-2.34,113.80,115.20,"WHOLE AOI")
