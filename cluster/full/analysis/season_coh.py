"""Season structure of coherence on the full 1650-pair stack.

The era test showed dry-season pairs (June to October) coherent at 0.32 to
0.42 and wet-season pairs at 0.10 to 0.22 at the same baseline, so the network
should be chosen by season, not by baseline alone. This measures, for every
pair: mean coherence in the Tay strip, at the 60 dam sites (5x5 windows) and
over the frame, with the calendar months of both dates and the baseline. It
reports the month-by-month matrix for short pairs and the inventory of
dry-to-dry annual pairs (both dates June to October, 300 to 430 days apart),
which would carry a rate without the short-baseline fading bias.
"""
import datetime as dt
import json
import os
import numpy as np
import h5py
from osgeo import gdal, osr

gdal.UseExceptions()
W = os.path.expanduser("~/full/mintpy_tay")
TAY = os.path.expanduser("~/full/tay/InSAR_LongTermRate.tif")
SITES = os.path.expanduser("~/full/dam_sites.geojson")
OUT = os.path.expanduser("~/full/analysis_out"); os.makedirs(OUT, exist_ok=True)
DRY = {6, 7, 8, 9, 10}
WIN = 2

with h5py.File(os.path.join(W, "velocity.h5"), "r") as f:
    att = dict(f.attrs)
x0, dx, y0, dy = float(att["X_FIRST"]), float(att["X_STEP"]), float(att["Y_FIRST"]), float(att["Y_STEP"]); ny, nx = int(att["LENGTH"]), int(att["WIDTH"])
mem = gdal.GetDriverByName("MEM").Create("", nx, ny, 1, gdal.GDT_Float32); mem.SetGeoTransform((x0, dx, 0.0, y0, 0.0, dy))
srs = osr.SpatialReference(); srs.ImportFromEPSG(4326); mem.SetProjection(srs.ExportToWkt()); mem.GetRasterBand(1).Fill(float("nan"))
src = gdal.Open(TAY); gdal.ReprojectImage(src, mem, src.GetProjection(), srs.ExportToWkt(), gdal.GRA_Average)
tay = mem.GetRasterBand(1).ReadAsArray(); tay[tay == 0] = np.nan; strip = np.isfinite(tay)
sites = [(f["properties"].get("site_id"), *f["geometry"]["coordinates"][:2]) for f in json.load(open(SITES))["features"]]
site_rc = [(int(round((lat - y0) / dy)), int(round((lon - x0) / dx))) for _, lon, lat in sites]
site_rc = [(r, c) for r, c in site_rc if WIN <= r < ny - WIN and WIN <= c < nx - WIN]

rows = []
with h5py.File(os.path.join(W, "inputs", "ifgramStack.h5"), "r") as f:
    dates = [(a.decode(), b.decode()) for a, b in f["date"][:]]
    n = len(dates)
    for i in range(n):
        c = f["coherence"][i].astype(np.float32); ok = np.isfinite(c) & (c > 0)
        d1 = dt.datetime.strptime(dates[i][0], "%Y%m%d").date(); d2 = dt.datetime.strptime(dates[i][1], "%Y%m%d").date()
        sv = [c[r - WIN:r + WIN + 1, cc - WIN:cc + WIN + 1].mean() for r, cc in site_rc]
        rows.append((dates[i][0], dates[i][1], (d2 - d1).days, d1.month, d2.month, float(c[strip & ok].mean()), float((c[strip & ok] > 0.3).mean()), float(np.mean(sv)), float(c[ok].mean())))
        if i % 200 == 0:
            print(f"  {i}/{n}", flush=True)
rows = np.array(rows, dtype=object)
with open(os.path.join(OUT, "pair_coh_season.csv"), "w") as fh:
    fh.write("date1,date2,tbase,month1,month2,strip_mean,strip_frac_gt03,dams_mean,frame_mean\n")
    for r in rows:
        fh.write(",".join(str(v) for v in r) + "\n")

tb = rows[:, 2].astype(int); m1 = rows[:, 3].astype(int); m2 = rows[:, 4].astype(int)
sc = rows[:, 5].astype(float); sf = rows[:, 6].astype(float); dm = rows[:, 7].astype(float); fr = rows[:, 8].astype(float)
dry1 = np.isin(m1, list(DRY)); dry2 = np.isin(m2, list(DRY))

print(f"\n{n} pairs. Strip coherence by season class (both dates), short pairs <= 60 d:")
short = tb <= 60
for lbl, mk in (("dry-dry", dry1 & dry2), ("dry-wet / wet-dry", dry1 ^ dry2), ("wet-wet", ~dry1 & ~dry2)):
    m = short & mk
    if m.sum():
        print(f"   {lbl:18s} n={m.sum():4d}  strip {sc[m].mean():.3f} (>0.3: {100*sf[m].mean():.0f}%)  dams {dm[m].mean():.3f}  frame {fr[m].mean():.3f}")
print("\nStrip coherence, pairs <= 60 d, by month of date1 (rows) x month of date2 (cols), '.' = none:")
print("     " + " ".join(f"{m:>5d}" for m in range(1, 13)))
for a in range(1, 13):
    line = []
    for b in range(1, 13):
        m = short & (m1 == a) & (m2 == b)
        line.append(f"{sc[m].mean():5.2f}" if m.sum() else "    .")
    print(f"  {a:2d} " + " ".join(line))
print("\nAnnual pairs (300 to 430 d) by season class:")
ann = (tb >= 300) & (tb <= 430)
for lbl, mk in (("dry-dry", dry1 & dry2), ("mixed", dry1 ^ dry2), ("wet-wet", ~dry1 & ~dry2)):
    m = ann & mk
    if m.sum():
        print(f"   {lbl:10s} n={m.sum():4d}  strip {sc[m].mean():.3f} (>0.3: {100*sf[m].mean():.0f}%)  dams {dm[m].mean():.3f}  frame {fr[m].mean():.3f}")
m = ann & dry1 & dry2
order = np.argsort(-sc[m])
print("\ndry-dry annual pairs, best first (date1 date2 tbase strip dams frame):")
for k in order[:30]:
    r = rows[m][k]; print(f"   {r[0]} {r[1]} {int(r[2]):4d} d  {float(r[5]):.3f}  {float(r[7]):.3f}  {float(r[8]):.3f}")
print(f"\nTwo-year dry-dry pairs (660 to 800 d): n={((tb>=660)&(tb<=800)&dry1&dry2).sum()}, strip {sc[(tb>=660)&(tb<=800)&dry1&dry2].mean():.3f}" if ((tb>=660)&(tb<=800)&dry1&dry2).sum() else "\nno two-year dry-dry pairs")
print("SEASON_DONE")
