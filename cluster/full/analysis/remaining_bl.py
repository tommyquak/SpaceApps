import os, glob, re, datetime, collections
ROOT = os.path.expanduser("~/full/pairs_tay")
PAT = re.compile(r"^(\d{4}-\d{2}-\d{2})_(\d{4}-\d{2}-\d{2})$")
done = set()
for d in glob.glob(ROOT + "/*/"):
    m = PAT.match(os.path.basename(d.rstrip("/")))
    if m and os.path.exists(d + "DONE"):
        done.add((m.group(1).replace("-", "")[2:], m.group(2).replace("-", "")[2:]))
planned = []
for f in sorted(glob.glob(os.path.expanduser("~/full/main_blocks/*.pairs"))):
    for ln in open(f):
        p = ln.strip()
        if "-" in p: planned.append(tuple(p.split("-")))
def dt(y):
    return datetime.date(2000 + int(y[:2]), int(y[2:4]), int(y[4:6]))
def span(p): return (dt(p[1]) - dt(p[0])).days
missing = [p for p in planned if p not in done]
print("planned %d, done %d, remaining %d" % (len(planned), len(done), len(missing)))
bins = [(0,60),(60,120),(120,240),(240,400),(400,700),(700,3000)]
def hist(ps, label):
    c = collections.Counter()
    for p in ps:
        s = span(p)
        for b in bins:
            if b[0] <= s < b[1]: c[b] += 1; break
    print("\n%s (n=%d)" % (label, len(ps)))
    for b in bins:
        n = c[b]
        print("   %5d-%-5d d : %5d  %s" % (b[0], b[1], n, "#" * int(40*n/max(1,len(ps)))))
hist([p for p in planned if p in done], "COMPLETED")
hist(missing, "REMAINING")
