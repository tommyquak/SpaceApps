import h5py, numpy as np, os
# Tay et al. 2025 Eq 4:  sigma_phi^2 = (1 - g^2) / (2 N g^2),  N = looks.
# Their exclusion is mean sigma_phi > 0.7 cm. Our 0.35% "reliable" figure used
# MintPy temporal coherence, a DIFFERENT quantity (self-consistency, not formal
# precision), so the two coverage numbers were never comparable.
# N must be the EFFECTIVE number of independent looks, not the nominal 952:
# neighbouring multilooked pixels are correlated, and an earlier audit measured
# L_eff ~10 against 140 nominal at 5x28.
W = os.environ.get("MINTPY_DIR", os.path.expanduser("~/full/mintpy_tay"))
lam = 0.23840354512922465
mm_per_rad = lam/(4*np.pi)*1000.0
NOMINAL = 17*56

f = h5py.File(W+"inputs/ifgramStack.h5","r")
n, H, Wd = f["coherence"].shape
d12 = [(a.decode(), b.decode()) for a,b in f["date"][:]]

# --- effective looks from closure statistics -------------------------------
# For a triplet of independent-noise legs, R = |mean(exp(i*C))| gives
# sigma_tri = sqrt(-2 ln R); sigma_ifg = sigma_tri/sqrt(3). Invert Eq 4 for N.
import itertools, random
dates = sorted(set(x for p in d12 for x in p)); di={d:i for i,d in enumerate(dates)}
idx = {p:i for i,p in enumerate(d12)}
tri=[]
for a,b,c in itertools.combinations(dates,3):
    if (a,b) in idx and (b,c) in idx and (a,c) in idx: tri.append((a,b,c))
random.seed(3); tri = random.sample(tri, min(120,len(tri)))
y0,x0 = H//2-80, Wd//2-80
sig, cohs = [], []
for (a,b,c) in tri:
    try:
        ab=f["unwrapPhase"][idx[(a,b)], y0:y0+160, x0:x0+160].astype(np.float64)
        bc=f["unwrapPhase"][idx[(b,c)], y0:y0+160, x0:x0+160].astype(np.float64)
        ac=f["unwrapPhase"][idx[(a,c)], y0:y0+160, x0:x0+160].astype(np.float64)
        ca=f["coherence"][idx[(a,b)], y0:y0+160, x0:x0+160].astype(np.float64)
    except Exception: continue
    m=(ab!=0)&(bc!=0)&(ac!=0)&np.isfinite(ca)&(ca>0)
    if m.sum()<200: continue
    C=(ab+bc-ac)[m]
    R=abs(np.mean(np.exp(1j*C)))
    if not (0<R<1): continue
    sig.append(np.sqrt(-2*np.log(R))/np.sqrt(3)); cohs.append(np.median(ca[m]))
sig=np.array(sig); cohs=np.array(cohs)
s_med=float(np.median(sig)); g_med=float(np.median(cohs))
L_eff = (1-g_med**2)/(2*s_med**2*g_med**2)
print("EFFECTIVE LOOKS")
print("  triplets used            : %d" % len(sig))
print("  median sigma_ifg         : %.4f rad" % s_med)
print("  median coherence         : %.4f" % g_med)
print("  nominal looks (17x56)    : %d" % NOMINAL)
print("  EFFECTIVE looks L_eff    : %.1f   (ratio %.3f of nominal)" % (L_eff, L_eff/NOMINAL))

# --- sigma_phi map under Tay's Eq 4 ---------------------------------------
sc = h5py.File(W+"avgSpatialCoh.h5","r")
g = sc[[k for k in sc.keys()][0]][:].astype(np.float64)
gm = np.isfinite(g)&(g>0.01)&(g<0.999)
def sigma_cm(N):
    out=np.full(g.shape, np.nan)
    out[gm]=np.sqrt((1-g[gm]**2)/(2*N*g[gm]**2))*mm_per_rad/10.0
    return out
print("\nSIGMA_PHI per interferogram (Tay Eq 4), cm")
for lbl,N in (("nominal 952",NOMINAL),("effective %.0f"%L_eff,max(L_eff,1))):
    s_=sigma_cm(N); v=s_[np.isfinite(s_)]
    print("  N=%-16s median %.3f  p90 %.3f   under 0.7 cm: %.2f%%" % (
        lbl, np.median(v), np.percentile(v,90), 100.0*(v<0.7).mean()))
# time-series sigma propagates as ~ sigma_ifg / sqrt(redundancy per date)
ni = h5py.File(W+"numInvIfgram.h5","r"); nv=ni[[k for k in ni.keys()][0]][:].astype(np.float64)
red = np.nanmedian(nv[nv>0])/len(dates)
print("\n  median ifgs per date (redundancy): %.1f" % red)
for lbl,N in (("nominal",NOMINAL),("effective",max(L_eff,1))):
    s_=sigma_cm(N)/np.sqrt(max(red,1.0)); v=s_[np.isfinite(s_)]
    print("  time-series sigma_phi, N=%-10s median %.3f cm   under 0.7 cm: %.2f%%" % (
        lbl, np.median(v), 100.0*(v<0.7).mean()))
print("\nTay report mean sigma_phi 0.37-0.73 cm at non-forest sites,")
print("0.83 and 0.89 cm at their two FOREST sites (both above their own 0.7 cut).")
