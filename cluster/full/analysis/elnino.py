"""2019 El Nino response: our short-baseline series against Tay's published
InSAR_MaxSubsidence-2019ElNino.tif, an episodic signal where a 60-day network
is strongest and the secular short-baseline bias matters least.

Tay's definition (dataset README): lowest displacement between 15 May and
16 Nov 2019 relative to the mean displacement from 12 Oct 2018 to 15 May 2019,
projected to vertical, metres, negative down. Same here from MintPy series
(positive toward the satellite), LOS divided by cos(incidence). Reported for
the bw60 series and, as a null, the bandwidth-20 series, on the same pixels,
with a block bootstrap for the correlation. Tay stays a benchmark.
"""
import datetime as dt
import json
import os
import numpy as np
import h5py
from osgeo import gdal, osr

gdal.UseExceptions()
RUNS = {"bw60": "/mnt/scratch/t/tommy/full/mintpy_bw60", "bw20": os.path.expanduser("~/full/mintpy_tay")}
GEOM = os.path.expanduser("~/full/mintpy_tay/inputs/geometryGeo.h5")
TAYEN = os.path.expanduser("~/full/tay/InSAR_MaxSubsidence-2019ElNino.tif")
SITES = os.path.expanduser("~/full/dam_sites.geojson")
PRE = (dt.date(2018, 10, 12), dt.date(2019, 5, 15)); WIN = (dt.date(2019, 5, 15), dt.date(2019, 11, 16))
BLOCK, NBOOT, W5 = 20, 400, 2
rng = np.random.default_rng(11)

with h5py.File(GEOM, "r") as f:
    inc = np.deg2rad(f["incidenceAngle"][:].astype(np.float64)); att = dict(f.attrs)
x0, dx, y0, dy = float(att["X_FIRST"]), float(att["X_STEP"]), float(att["Y_FIRST"]), float(att["Y_STEP"]); ny, nx = inc.shape
mem = gdal.GetDriverByName("MEM").Create("", nx, ny, 1, gdal.GDT_Float32); mem.SetGeoTransform((x0, dx, 0.0, y0, 0.0, dy))
srs = osr.SpatialReference(); srs.ImportFromEPSG(4326); mem.SetProjection(srs.ExportToWkt()); mem.GetRasterBand(1).Fill(float("nan"))
src = gdal.Open(TAYEN); gdal.ReprojectImage(src, mem, src.GetProjection(), srs.ExportToWkt(), gdal.GRA_Average)
tay = mem.GetRasterBand(1).ReadAsArray().astype(np.float64); tay[tay == 0] = np.nan; tay *= 1000.0
strip = np.isfinite(tay)
print(f"Tay El Nino raster on our grid: {strip.sum()} px; mean {np.nanmean(tay):+.1f} mm, std {np.nanstd(tay):.1f}, p10 {np.nanpercentile(tay,10):+.1f}, p90 {np.nanpercentile(tay,90):+.1f}")
sites = [(f["properties"].get("site_id"), *f["geometry"]["coordinates"][:2]) for f in json.load(open(SITES))["features"]]
site_rc = [(int(round((lat - y0) / dy)), int(round((lon - x0) / dx))) for _, lon, lat in sites]


def block_boot(a, b, ij):
    bid = (ij[0] // BLOCK) * 100000 + (ij[1] // BLOCK)
    uniq, inv = np.unique(bid, return_inverse=True); order = np.argsort(inv, kind="stable")
    starts = np.searchsorted(inv[order], np.arange(len(uniq))); ends = np.r_[starts[1:], len(order)]
    rs = []
    for _ in range(NBOOT):
        pick = rng.integers(0, len(uniq), len(uniq)); idx = np.concatenate([order[starts[p]:ends[p]] for p in pick])
        if idx.size > 30 and a[idx].std() > 0 and b[idx].std() > 0:
            rs.append(np.corrcoef(a[idx], b[idx])[0, 1])
    return np.percentile(rs, [2.5, 97.5]) if rs else (np.nan, np.nan)


def elnino_from(W, series):
    p = os.path.join(W, series)
    if not os.path.exists(p):
        return None
    with h5py.File(p, "r") as f:
        dates = [dt.datetime.strptime(d.decode(), "%Y%m%d").date() for d in f["date"][:]]
        ipre = [i for i, d in enumerate(dates) if PRE[0] <= d <= PRE[1]]; iwin = [i for i, d in enumerate(dates) if WIN[0] <= d <= WIN[1]]
        pre = np.mean(np.stack([f["timeseries"][i] for i in ipre]), axis=0).astype(np.float64)
        win = np.min(np.stack([f["timeseries"][i] for i in iwin]), axis=0).astype(np.float64)
    los = (win - pre) * 1000.0
    return los / np.cos(inc), len(ipre), len(iwin)


for run, W in RUNS.items():
    with h5py.File(os.path.join(W, "temporalCoherence.h5"), "r") as f:
        tc = f["temporalCoherence"][:].astype(np.float64)
    for series in ("timeseries_SET.h5", "timeseries_SET_ramp.h5"):
        r = elnino_from(W, series)
        if r is None:
            print(f"\n[{run} {series}] missing"); continue
        ours, npre, nwin = r
        print(f"\n[{run} {series}] {npre} pre-epochs, {nwin} El Nino epochs; vertical mm, negative down")
        for thr in (0.5, 0.7):
            m = strip & (tc >= thr) & np.isfinite(ours) & np.isfinite(tay)
            if m.sum() < 200:
                print(f"   tempCoh>={thr}: n={m.sum()} too few"); continue
            a, b = ours[m], tay[m]; rr = np.corrcoef(a, b)[0, 1]; lo, hi = block_boot(a, b, np.nonzero(m))
            print(f"   tempCoh>={thr}: n={m.sum():7d}  r={rr:+.3f} CI [{lo:+.3f},{hi:+.3f}]  ours mean {a.mean():+.1f} std {a.std():.1f} | Tay mean {b.mean():+.1f} std {b.std():.1f}")
        sv, tv = [], []
        for rr_, cc in site_rc:
            if W5 <= rr_ < ny - W5 and W5 <= cc < nx - W5:
                po = ours[rr_ - W5:rr_ + W5 + 1, cc - W5:cc + W5 + 1]; pt = tay[rr_ - W5:rr_ + W5 + 1, cc - W5:cc + W5 + 1]
                if np.isfinite(po).sum() >= 13 and np.isfinite(pt).sum() >= 13:
                    sv.append(np.nanmean(po)); tv.append(np.nanmean(pt))
        if len(sv) >= 5:
            print(f"   dam sites with Tay data: {len(sv)}: ours mean {np.mean(sv):+.1f} | Tay {np.mean(tv):+.1f}; r={np.corrcoef(sv, tv)[0,1]:+.3f}")
        else:
            print(f"   dam sites with Tay El Nino data: {len(sv)} (Tay raster does not reach the dam zone)")
print("ELNINO_DONE")
