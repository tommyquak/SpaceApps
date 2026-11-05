"""Era test: is the pre-rewetting ground (2016 to 2017) more coherent than the
same ground in 2019 to 2021 at the SAME temporal baseline?

Reads the geocoded 17x56 coherence of the eight era pairs (coherence-only MAIN,
job 821197) and, for each, every 2019-2021 pair in the finished stack with the
identical baseline in days. Reports mean coherence and the fraction above 0.3
in the Tay strip, at the 60 dam sites (5x5 windows), and over the frame.
Same looks, same grid, same reference geometry, same processing chain; only
the years differ.
"""
import datetime as dt
import glob
import json
import os
import numpy as np
import h5py
from osgeo import gdal, osr

gdal.UseExceptions()
ERA = os.path.expanduser("~/full/pairs_era")
W = os.path.expanduser("~/full/mintpy_tay")
TAY = os.path.expanduser("~/full/tay/InSAR_LongTermRate.tif")
SITES = os.path.expanduser("~/full/dam_sites.geojson")
WIN = 2

with h5py.File(os.path.join(W, "velocity.h5"), "r") as f:
    att = dict(f.attrs)
x0, dx, y0, dy = float(att["X_FIRST"]), float(att["X_STEP"]), float(att["Y_FIRST"]), float(att["Y_STEP"])
ny, nx = int(att["LENGTH"]), int(att["WIDTH"])


def grid_mem(nodata=np.nan):
    m = gdal.GetDriverByName("MEM").Create("", nx, ny, 1, gdal.GDT_Float32); m.SetGeoTransform((x0, dx, 0.0, y0, 0.0, dy))
    s = osr.SpatialReference(); s.ImportFromEPSG(4326); m.SetProjection(s.ExportToWkt()); m.GetRasterBand(1).Fill(nodata); return m


def warp(path_or_ds, alg=gdal.GRA_Average):
    src = gdal.Open(path_or_ds) if isinstance(path_or_ds, str) else path_or_ds
    m = grid_mem(); s = osr.SpatialReference(); s.ImportFromEPSG(4326)
    gdal.ReprojectImage(src, m, src.GetProjection(), s.ExportToWkt(), alg)
    return m.GetRasterBand(1).ReadAsArray().astype(np.float32)


tay = warp(TAY); tay[tay == 0] = np.nan; strip = np.isfinite(tay)
sites = [(f["properties"].get("site_id"), *f["geometry"]["coordinates"][:2]) for f in json.load(open(SITES))["features"]]
site_rc = [(int(round((lat - y0) / dy)), int(round((lon - x0) / dx))) for _, lon, lat in sites]


def region_stats(c):
    ok = np.isfinite(c) & (c > 0)
    out = {}
    out["frame"] = (float(c[ok].mean()), float((c[ok] > 0.3).mean()))
    m = strip & ok; out["strip"] = (float(c[m].mean()), float((c[m] > 0.3).mean()))
    sv = []
    for r, cc in site_rc:
        if WIN <= r < ny - WIN and WIN <= cc < nx - WIN:
            p = c[r - WIN:r + WIN + 1, cc - WIN:cc + WIN + 1]; p = p[np.isfinite(p) & (p > 0)]
            if p.size >= 13: sv.append(p.mean())
    sv = np.array(sv); out["dams"] = (float(sv.mean()), float((sv > 0.3).mean()), len(sv))
    return out


# era pairs
era = []
for d in sorted(glob.glob(os.path.join(ERA, "*"))):
    tifs = glob.glob(os.path.join(d, "*_17rlks_56alks.cor.geo.aoi90.tif"))
    if not tifs or not os.path.exists(os.path.join(d, "DONE_COH")):
        continue
    cal = os.path.basename(d); a, b = cal.split("_")
    tb = (dt.date.fromisoformat(b) - dt.date.fromisoformat(a)).days
    # DEFECT 26: ISCE's .cor is two-band (amplitude, coherence) and the geocoded
    # GeoTIFF keeps both; band 1 is amplitude (values in the thousands).
    ds = gdal.Open(tifs[0]); band = ds.RasterCount
    if ds.RasterYSize != ny or ds.RasterXSize != nx:
        m = grid_mem(); s_ = osr.SpatialReference(); s_.ImportFromEPSG(4326)
        one = gdal.Translate("", ds, format="MEM", bandList=[band])
        gdal.ReprojectImage(one, m, one.GetProjection(), s_.ExportToWkt(), gdal.GRA_Average)
        c = m.GetRasterBand(1).ReadAsArray().astype(np.float32)
    else:
        c = ds.GetRasterBand(band).ReadAsArray().astype(np.float32)
    assert np.nanmax(c) <= 1.01, f"{tifs[0]}: band {band} max {np.nanmax(c)} is not a coherence"
    era.append((a, b, tb, region_stats(c)))
print(f"{len(era)} era pairs")

# 2019-2021 controls at identical baselines
with h5py.File(os.path.join(W, "inputs", "ifgramStack.h5"), "r") as f:
    dates = [(p[0].decode(), p[1].decode()) for p in f["date"][:]]
    need = sorted({tb for _, _, tb, _ in era})
    ctrl = {tb: [] for tb in need}
    for i, (a, b) in enumerate(dates):
        da, db = dt.datetime.strptime(a, "%Y%m%d").date(), dt.datetime.strptime(b, "%Y%m%d").date()
        tb = (db - da).days
        if tb in ctrl and "2019-01-01" <= da.isoformat() <= "2021-12-31":
            ctrl[tb].append((a, b, i))
    print("control pairs per baseline:", {tb: len(v) for tb, v in ctrl.items()})
    ctrl_stats = {}
    for tb, lst in ctrl.items():
        acc = []
        for a, b, i in lst:
            acc.append(region_stats(f["coherence"][i].astype(np.float32)))
        if acc:
            ctrl_stats[tb] = {k: tuple(np.mean([s[k][j] for s in acc]) for j in range(2)) for k in ("frame", "strip", "dams")}
            ctrl_stats[tb]["n"] = len(acc)

print("\n%-24s %5s | %-22s | %-22s | %-22s" % ("era pair", "days", "Tay strip  era / 2019-21", "dam sites  era / 2019-21", "frame      era / 2019-21"))
for a, b, tb, s in era:
    c = ctrl_stats.get(tb)
    fmt = lambda k: f"{s[k][0]:.3f}({100*s[k][1]:3.0f}%) / " + (f"{c[k][0]:.3f}({100*c[k][1]:3.0f}%) n={c['n']}" if c else "no control")
    print("%-24s %5d | %s | %s | %s" % (f"{a}_{b}", tb, fmt("strip"), fmt("dams"), fmt("frame")))
print("\n(values: mean coherence, and percent of pixels above 0.3; dam sites are 5x5 window means at 60 sites)")
print("ERA_COH_DONE")
