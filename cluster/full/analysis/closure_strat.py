import os, glob, re, itertools, datetime, numpy as np
from osgeo import gdal
gdal.UseExceptions()
ROOT = os.path.expanduser(os.environ.get("ROOT", "~/full/pairs_tay"))
PAT = re.compile(r"^(\d{4}-\d{2}-\d{2})_(\d{4}-\d{2}-\d{2})$")

edges, dates = {}, set()
for d in sorted(glob.glob(ROOT + "/*/")):
    m = PAT.match(os.path.basename(d.rstrip("/")))
    if not m or not os.path.exists(d + "DONE"): continue
    g = [u for u in glob.glob(d + "filt_*.unw.geo.aoi90.tif")
         if "_msk" not in u and "conncomp" not in u]
    if g:
        edges[(m.group(1), m.group(2))] = g[0]; dates.add(m.group(1)); dates.add(m.group(2))
print("completed pairs:", len(edges))
ds = sorted(dates)
def dt(s): return datetime.date(*map(int, s.split("-")))
tri = [(a,b,c) for a,b,c in itertools.combinations(ds,3)
       if (a,b) in edges and (b,c) in edges and (a,c) in edges]
print("closed triplets:", len(tri))

o = gdal.Open(next(iter(edges.values())))
H, W = o.RasterYSize, o.RasterXSize
y0, x0 = H//2-80, W//2-80
cache = {}
def get(k):
    if k not in cache:
        if len(cache) > 260: cache.clear()
        q = gdal.Open(edges[k])
        cache[k] = q.GetRasterBand(q.RasterCount).ReadAsArray(x0, y0, 160, 160).astype(np.float64)
    return cache[k]

# stratify by the triplet's longest leg: does a short-baseline subnetwork close?
bins = [(0,30),(30,60),(60,100),(100,200),(200,400),(400,2000)]
acc = {b: [0,0,[]] for b in bins}
tot_nz = tot_n = 0
import random
random.seed(11); sample = tri if len(tri) <= 1400 else random.sample(tri, 1400)
for (a,b,c) in sample:
    span = (dt(c)-dt(a)).days
    bn = next((x for x in bins if x[0] <= span < x[1]), None)
    if bn is None: continue
    try: ab, bc, ac = get((a,b)), get((b,c)), get((a,c))
    except Exception: continue
    m = (ab!=0)&(bc!=0)&(ac!=0)
    if m.sum() < 100: continue
    C = (ab+bc-ac)[m]; Cw = np.angle(np.exp(1j*C))
    k = np.round((C-Cw)/(2*np.pi))
    nz = int((k!=0).sum()); n = int(k.size)
    acc[bn][0]+=nz; acc[bn][1]+=n; acc[bn][2].append(float(np.angle(np.mean(np.exp(1j*Cw)))))
    tot_nz+=nz; tot_n+=n

lam = 0.2384035; mmrad = lam/(4*np.pi)*1000
print("\n%-14s %8s %10s %14s" % ("max leg (days)","triplets","%% nonzero","wrapped (mm)"))
for bn in bins:
    nz, n, w = acc[bn]
    if not n: continue
    wm = np.angle(np.mean(np.exp(1j*np.array(w))))
    print("%-14s %8d %9.1f%% %13.3f" % ("%d-%d"%bn, len(w), 100.0*nz/n, wm*mmrad))
if tot_n:
    print("\nOVERALL nonzero: %.2f %%   (5x28 run: 96-97 %%; at 372 pairs: 63.6 %%)" % (100.0*tot_nz/tot_n))
