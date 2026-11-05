"""NISAR GUNW stage 2b: 60-day composites on both tracks with a data-chosen
reference, and their agreement on common ground.

Defect 24 (found by executing stage 2): the ALOS-2 reference pixel is 1.5 km
outside the descending frame and on nodata in the ascending one, so no
composite was produced. The reference is now the 5x5 window of highest
chain-mean coherence inside the ALOS-2 stable-ground mask (tempCoh >= 0.7),
chosen per frame. Sign is reported in the product's own convention: positive
phase is taken as range increase (ISCE lineage), i.e. subsidence, tagged [S].
"""
import glob
import json
import os
import re
import numpy as np
import h5py
from osgeo import gdal, osr
from scipy.ndimage import uniform_filter

gdal.UseExceptions()
D = os.path.expanduser(os.environ.get("GUNW_DIR", "/mnt/scratch/t/tommy/nisar/gunw"))
OUT = os.path.expanduser("~/full/nisar_out"); os.makedirs(OUT, exist_ok=True)
TAY = os.path.expanduser("~/full/tay/InSAR_LongTermRate.tif")
TC = os.path.expanduser("~/full/mintpy_tay/temporalCoherence.h5")
SITES = os.path.expanduser("~/full/dam_sites.geojson")
G = "science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/HH/"
WIN = 2


def find_ds(h, name, prefer):
    hits = []
    h.visititems(lambda n, o: hits.append(n) if isinstance(o, h5py.Dataset) and n.split("/")[-1] == name else None)
    for p in hits:
        if p.startswith(prefer):
            return p
    return hits[0] if hits else None


def read_grid(path):
    with h5py.File(path, "r") as h:
        g = dict(coh=h[G + "coherenceMagnitude"][:].astype(np.float32), unw=h[G + "unwrappedPhase"][:].astype(np.float32),
                 cc=h[G + "connectedComponents"][:], ion=h[G + "ionospherePhaseScreen"][:].astype(np.float32))
        g["x"] = h[find_ds(h, "xCoordinates", G)][:]; g["y"] = h[find_ds(h, "yCoordinates", G)][:]
        v = h[find_ds(h, "epsg", G) or find_ds(h, "projection", G)][()]
        try:
            g["epsg"] = int(v)
        except Exception:
            s = v.decode() if isinstance(v, bytes) else str(v); g["epsg"] = int(re.findall(r"\d{4,5}", s)[-1])
        p = find_ds(h, "centerFrequency", "science"); g["wl"] = 299792458.0 / float(np.ravel(h[p][()])[0]) if p else 0.2420
    return g


def gt(g):
    dx = float(g["x"][1] - g["x"][0]); dy = float(g["y"][1] - g["y"][0])
    return (float(g["x"][0]) - dx / 2, dx, 0.0, float(g["y"][0]) - dy / 2, 0.0, dy)


def mem_like(g, arr=None):
    ny, nx = g["coh"].shape
    m = gdal.GetDriverByName("MEM").Create("", nx, ny, 1, gdal.GDT_Float32); m.SetGeoTransform(gt(g))
    s = osr.SpatialReference(); s.ImportFromEPSG(g["epsg"]); m.SetProjection(s.ExportToWkt())
    if arr is None:
        m.GetRasterBand(1).Fill(float("nan"))
    else:
        m.GetRasterBand(1).WriteArray(np.where(np.isfinite(arr), arr, np.nan).astype(np.float32)); m.GetRasterBand(1).SetNoDataValue(float("nan"))
    return m


def warp_onto(g, src, alg=gdal.GRA_Average):
    m = mem_like(g); s = osr.SpatialReference(); s.ImportFromEPSG(g["epsg"])
    gdal.ReprojectImage(src, m, src.GetProjection(), s.ExportToWkt(), alg)
    return m.GetRasterBand(1).ReadAsArray().astype(np.float32)


def latlon_mem(arr, x0, y0, dx, dy):
    ny, nx = arr.shape
    m = gdal.GetDriverByName("MEM").Create("", nx, ny, 1, gdal.GDT_Float32); m.SetGeoTransform((x0, dx, 0.0, y0, 0.0, dy))
    s = osr.SpatialReference(); s.ImportFromEPSG(4326); m.SetProjection(s.ExportToWkt())
    m.GetRasterBand(1).WriteArray(arr.astype(np.float32)); m.GetRasterBand(1).SetNoDataValue(float("nan")); return m


def to_pix(g, lon, lat):
    src = osr.SpatialReference(); src.ImportFromEPSG(4326); src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    dst = osr.SpatialReference(); dst.ImportFromEPSG(g["epsg"]); dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    X, Y, _ = osr.CoordinateTransformation(src, dst).TransformPoint(float(lon), float(lat))
    dx = float(g["x"][1] - g["x"][0]); dy = float(g["y"][1] - g["y"][0])
    return int(round((Y - float(g["y"][0])) / dy)), int(round((X - float(g["x"][0])) / dx))


def win_mean(a, r, c, w=WIN):
    ny, nx = a.shape
    if r - w < 0 or c - w < 0 or r + w >= ny or c + w >= nx:
        return np.nan
    p = a[r - w:r + w + 1, c - w:c + w + 1]
    return float(np.nanmean(p)) if np.isfinite(p).any() else np.nan


def write_tif(path, g, a):
    d = gdal.GetDriverByName("GTiff").Create(path, a.shape[1], a.shape[0], 1, gdal.GDT_Float32, ["COMPRESS=DEFLATE", "TILED=YES"])
    d.SetGeoTransform(gt(g)); s = osr.SpatialReference(); s.ImportFromEPSG(g["epsg"]); d.SetProjection(s.ExportToWkt())
    d.GetRasterBand(1).WriteArray(a); d.GetRasterBand(1).SetNoDataValue(float("nan")); d.FlushCache()


sites = [(f["properties"].get("site_id"), *f["geometry"]["coordinates"][:2]) for f in json.load(open(SITES))["features"]]
with h5py.File(TC, "r") as f:
    tc = f["temporalCoherence"][:].astype(np.float32); att = dict(f.attrs)
stable_ds = latlon_mem((tc >= 0.7).astype(np.float32), float(att["X_FIRST"]), float(att["Y_FIRST"]), float(att["X_STEP"]), float(att["Y_STEP"]))
tay_ds = gdal.Open(TAY)

files = sorted(glob.glob(os.path.join(D, "*.h5")))
recs = []
for f in files:
    m = re.search(r"_(\d{3})_(\d{3})_([AD])_(\d{3})_\d{3}_\d{4}_SH_(\d{8})T\d{6}_\d{8}T\d{6}_(\d{8})T", os.path.basename(f))
    recs.append(dict(path=f, track=m.group(2), ad=m.group(3), d1=m.group(5), d2=m.group(6)))

results = {}
for ad in ("D", "A"):
    chain = []
    for r in sorted([r for r in recs if r["ad"] == ad], key=lambda r: r["d1"]):
        if not chain or chain[-1]["d2"] == r["d1"]:
            chain.append(r)
    if len(chain) < 2:
        print(f"[{ad}] no chain"); continue
    gs = [read_grid(r["path"]) for r in chain]
    g0 = gs[0]; k = g0["wl"] / (4 * np.pi) * 1000.0
    strip = np.isfinite(warp_onto(g0, tay_ds)) ; tay = warp_onto(g0, tay_ds); tay[tay == 0] = np.nan; strip = np.isfinite(tay)
    stable = warp_onto(g0, stable_ds) > 0.5
    valid = np.ones(g0["coh"].shape, bool)
    for g in gs:
        valid &= np.isfinite(g["unw"]) & (g["cc"] > 0) & np.isfinite(g["coh"])
    chain_coh = np.mean(np.stack([g["coh"] for g in gs]), axis=0)
    # reference: best 5x5 all-valid window inside stable ground
    cand = uniform_filter(np.where(valid, chain_coh, 0.0), 2 * WIN + 1) * (uniform_filter(valid.astype(np.float32), 2 * WIN + 1) > 0.999)
    cand[~stable] = 0.0
    if cand.max() <= 0:
        cand = uniform_filter(np.where(valid, chain_coh, 0.0), 2 * WIN + 1) * (uniform_filter(valid.astype(np.float32), 2 * WIN + 1) > 0.999)
        print(f"[{ad}] no stable-ground window valid through the chain; falling back to the best window anywhere")
    r0, c0 = np.unravel_index(int(np.argmax(cand)), cand.shape)
    print(f"\n=== [{ad}] track {chain[0]['track']} chain {chain[0]['d1']} -> {chain[-1]['d2']} ({len(chain)} pairs), wl {g0['wl']:.4f} m, epsg {g0['epsg']}")
    print(f"   reference window at row {r0} col {c0}, chain-mean coherence {cand[r0, c0]:.3f}, inside stable ground: {bool(stable[r0, c0])}")
    print(f"   valid through all pairs: {100*valid.mean():.1f}% of frame; {100*valid[strip].mean():.1f}% of Tay footprint")
    cum_raw = np.zeros_like(g0["unw"]); cum_ion = np.zeros_like(g0["unw"])
    for g in gs:
        u = g["unw"] - win_mean(g["unw"], r0, c0)
        i = g["ion"] - win_mean(g["ion"], r0, c0)
        cum_raw += u; cum_ion += (u - i)
    cum_raw[~valid] = np.nan; cum_ion[~valid] = np.nan
    dr, di = cum_raw * k, cum_ion * k
    for lbl, a in (("raw", dr), ("ion-corrected", di)):
        sv = np.array([win_mean(a, *to_pix(g0, lon, lat)) for _, lon, lat in sites])
        print(f"   [{lbl}] range-change mm (positive = range increase = subsidence under ISCE lineage [S])")
        print(f"      Tay footprint mean {np.nanmean(a[strip]):+.1f} std {np.nanstd(a[strip]):.1f} | stable ground mean {np.nanmean(a[stable]):+.1f} std {np.nanstd(a[stable]):.1f} (noise floor)")
        print(f"      dam sites n={np.isfinite(sv).sum()}: mean {np.nanmean(sv):+.1f} std {np.nanstd(sv):.1f} min {np.nanmin(sv):+.1f} max {np.nanmax(sv):+.1f}")
    write_tif(os.path.join(OUT, f"nisar_{ad}_disp60_ion_mm.tif"), g0, di)
    write_tif(os.path.join(OUT, f"nisar_{ad}_disp60_raw_mm.tif"), g0, dr)
    write_tif(os.path.join(OUT, f"nisar_{ad}_chain_coh.tif"), g0, chain_coh)
    results[ad] = dict(g=g0, di=di, dr=dr, strip=strip, stable=stable, chain=chain)

if "D" in results and "A" in results:
    gd = results["D"]["g"]
    a_on_d = warp_onto(gd, mem_like(results["A"]["g"], results["A"]["di"]), gdal.GRA_Bilinear)
    d = results["D"]["di"]
    both = np.isfinite(a_on_d) & np.isfinite(d)
    strip = results["D"]["strip"]; stable = results["D"]["stable"]
    print(f"\n=== two-track agreement (ion-corrected, A warped onto D grid), common px {both.sum()}")
    for lbl, mk in (("Tay footprint", strip & both), ("stable ground", stable & both), ("all common", both)):
        if mk.sum() < 100:
            continue
        x, y = d[mk], a_on_d[mk]
        # smooth to 1 km blocks before correlating: per-pixel noise is spatially white, signal is not
        print(f"   [{lbl}] n={mk.sum()}: r={np.corrcoef(x, y)[0,1]:+.3f}; D mean {x.mean():+.1f} std {x.std():.1f}; A mean {y.mean():+.1f} std {y.std():.1f}; mean |D-A| {np.mean(np.abs(x-y)):.1f} mm")
    # 1 km block means
    ny, nx = d.shape; b = 12
    def blk(a):
        m = np.isfinite(a).astype(np.float32); s = np.where(np.isfinite(a), a, 0.0)
        S = uniform_filter(s, b); M = uniform_filter(m, b)
        out = np.where(M > 0.8, S / np.maximum(M, 1e-6), np.nan); return out[b // 2::b, b // 2::b]
    D1, A1, S1 = blk(np.where(both, d, np.nan)), blk(np.where(both, a_on_d, np.nan)), blk(np.where(both & strip, 1.0, np.nan))
    mk = np.isfinite(D1) & np.isfinite(A1)
    if mk.sum() > 30:
        print(f"   [1 km blocks, all common] n={mk.sum()}: r={np.corrcoef(D1[mk], A1[mk])[0,1]:+.3f}")
        mk2 = mk & np.isfinite(S1)
        if mk2.sum() > 30:
            print(f"   [1 km blocks, Tay footprint] n={mk2.sum()}: r={np.corrcoef(D1[mk2], A1[mk2])[0,1]:+.3f}")
print("STAGE2B_DONE")
