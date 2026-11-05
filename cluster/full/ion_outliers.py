import glob, os, re, numpy as np, json
B = "/mnt/scratch/t/tommy/full/pairs_ion"
N = "filt_ion_80rlks_448alks.ion"
def dims(v):
    t=open(v).read()
    m=re.search(r'rasterXSize="(\d+)"\s+rasterYSize="(\d+)"',t)
    if m: return int(m.group(1)),int(m.group(2))
    m=re.search(r'rasterYSize="(\d+)"\s+rasterXSize="(\d+)"',t); return int(m.group(2)),int(m.group(1))

stats=[]
for d in sorted(glob.glob(B+"/*/ion/ion_cal")):
    pair=d.split("/")[-3]
    f=os.path.join(d,N); v=f+".vrt"
    if not(os.path.exists(f) and os.path.exists(v)): continue
    W,H=dims(v); a=np.fromfile(f,dtype=np.float32)
    if a.size!=W*H: continue
    a=a.reshape(H,W); g=a[np.isfinite(a)]
    if g.size==0: continue
    stats.append(dict(pair=pair, mean=float(g.mean()), ptp=float(g.max()-g.min()),
        gaz=float(np.median(np.abs(np.diff(a,axis=0)))), grg=float(np.median(np.abs(np.diff(a,axis=1)))),
        absmax=float(np.abs(g).max())))
print("pairs measured:", len(stats))

am=np.array([abs(s["mean"]) for s in stats]); pt=np.array([s["ptp"] for s in stats])
for q in (50,75,90,92,95,97,98,99,100):
    print("  p%-3d |mean| %12.2f   ptp %12.2f" % (q,np.percentile(am,q),np.percentile(pt,q)))

# Physical ceiling: a differential L-band ionosphere over this scene cannot
# plausibly exceed a few hundred rad (50 TECU differential is ~670 rad, and
# that is already an extreme storm). Cut generously at 5x that.
BAD = [s for s in stats if abs(s["mean"])>200 or s["ptp"]>1000 or s["gaz"]>1.0 or s["grg"]>1.0]
GOOD= [s for s in stats if s not in BAD]
print("\nBAD pairs: %d   GOOD: %d" % (len(BAD),len(GOOD)))
print("\nworst offenders:")
for s in sorted(BAD,key=lambda s:-abs(s["mean"]))[:12]:
    print("   %s  mean %14.1f  ptp %14.1f  gaz %10.2f" % (s["pair"],s["mean"],s["ptp"],s["gaz"]))

g=np.array([s["mean"] for s in GOOD])
print("\nafter exclusion, GOOD set: mean median %.2f rad, |mean| max %.2f, ptp max %.2f" % (
    np.median(g), np.abs(g).max(), max(s["ptp"] for s in GOOD)))

# network connectivity over the 93 dates using only GOOD pairs
dates=set()
for s in stats: dates.update(s["pair"].split("-"))
dates=sorted(dates)
idx={d:i for i,d in enumerate(dates)}
par=list(range(len(dates)))
def find(x):
    while par[x]!=x: par[x]=par[par[x]]; x=par[x]
    return x
for s in GOOD:
    a,b=s["pair"].split("-"); ra,rb=find(idx[a]),find(idx[b])
    if ra!=rb: par[ra]=rb
comp={}
for d in dates: comp.setdefault(find(idx[d]),[]).append(d)
print("\ndates: %d   connected components after exclusion: %d" % (len(dates),len(comp)))
for k,v in sorted(comp.items(),key=lambda kv:-len(kv[1]))[:4]:
    print("   component size %d  (e.g. %s)" % (len(v), ", ".join(v[:3])))
orphan=[d for d in dates if len([s for s in GOOD if d in s["pair"].split("-")])==0]
print("dates with NO good pair:", orphan if orphan else "none")
json.dump([s["pair"] for s in BAD], open(os.path.expanduser("~/ion_bad_pairs.json"),"w"))
print("\nwrote ~/ion_bad_pairs.json")
