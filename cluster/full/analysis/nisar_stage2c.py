"""NISAR GUNW stage 2c (from 2b): contiguous composites on both tracks with a data-chosen
reference, and their agreement on common ground.

Defect 24 (found by executing stage 2): the ALOS-2 reference pixel is 1.5 km
outside the descending frame and on nodata in the ascending one, so no
composite was produced. The reference is now the 5x5 window of highest
chain-mean coherence inside the ALOS-2 stable-ground mask (tempCoh >= 0.7),
chosen per frame. Stage 2b showed the chained phase is dominated by per-track
long-wavelength error (tracks anti-correlate, stable ground moves 16 to 47 mm), so
this stage removes a quadratic ramp per pair fitted on stable ground, masks pixels
whose ionosphere-screen uncertainty is high, and high-passes at 10 km so only
local anomalies remain, which is what site-versus-control needs. Sign is reported in the product's own convention: positive
phase is taken as range increase (ISCE lineage), i.e. subsidence, tagged [S].

NISAR_OUT isolates each run. NISAR_{A,D}_TRACK/FRAME/START/END and NISAR_CRID
can select a single product group and contiguous interval. Historical disp60
filenames are retained for compatibility; the source manifests give actual dates.
"""
import glob
import json
import os
import re
import numpy as np
import h5py
from osgeo import gdal, osr
from scipy.ndimage import uniform_filter
from nisar_safety import deramp, read_pair_identity, select_chain, validate_grids
HP_KM = float(os.environ.get('HP_KM', '10'))
if not np.isfinite(HP_KM) or HP_KM <= 0:
    raise ValueError("HP_KM must be finite and positive")

gdal.UseExceptions()
D = os.path.expanduser(os.environ.get("GUNW_DIR", "/mnt/scratch/t/tommy/nisar/gunw"))
OUT = os.path.expanduser(os.environ.get("NISAR_OUT", "~/full/nisar_out")); os.makedirs(OUT, exist_ok=True)
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
        u = find_ds(h, "ionospherePhaseScreenUncertainty", G)
        g["ionu"] = h[u][:].astype(np.float32) if u else None
        g["x"] = h[find_ds(h, "xCoordinates", G)][:]; g["y"] = h[find_ds(h, "yCoordinates", G)][:]
        v = h[find_ds(h, "epsg", G) or find_ds(h, "projection", G)][()]
        try:
            g["epsg"] = int(v)
        except Exception:
            s = v.decode() if isinstance(v, bytes) else str(v); g["epsg"] = int(re.findall(r"\d{4,5}", s)[-1])
        p = "science/LSAR/GUNW/grids/frequencyA/centerFrequency"
        if p not in h:
            raise ValueError(f"Missing centerFrequency in {path}")
        frequency = float(np.ravel(h[p][()])[0])
        if not np.isfinite(frequency) or frequency <= 0:
            raise ValueError(f"Invalid centerFrequency in {path}")
        g["wl"] = 299792458.0 / frequency
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
recs = [read_pair_identity(path) for path in files]


def highpass(a, km, px_m):
    b = max(1, int(round(km * 1000.0 / px_m)))
    m = np.isfinite(a).astype(np.float32); s = np.where(np.isfinite(a), a, 0.0)
    S = uniform_filter(s, b); M = uniform_filter(m, b)
    local = np.where(M > 0.3, S / np.maximum(M, 1e-6), np.nan)
    return a - local


results = {}
for ad in ("D", "A"):
    chain = select_chain(recs, ad, track=os.environ.get(f"NISAR_{ad}_TRACK"),
                         frame=os.environ.get(f"NISAR_{ad}_FRAME"),
                         crid=os.environ.get("NISAR_CRID"),
                         start=os.environ.get(f"NISAR_{ad}_START"),
                         end=os.environ.get(f"NISAR_{ad}_END"))
    gs = [read_grid(r["path"]) for r in chain]
    validate_grids(gs)
    g0 = gs[0]; k = g0["wl"] / (4 * np.pi) * 1000.0
    px_m = abs(float(g0["x"][1] - g0["x"][0]))
    tay = warp_onto(g0, tay_ds); tay[tay == 0] = np.nan; strip = np.isfinite(tay)
    stable = warp_onto(g0, stable_ds) > 0.5
    valid = np.ones(g0["coh"].shape, bool)
    for g in gs:
        valid &= (np.isfinite(g["unw"]) & (g["cc"] > 0) & (g["cc"] < 65535)
                  & np.isfinite(g["coh"]) & (g["coh"] > 0) & (g["coh"] <= 1)
                  & np.isfinite(g["ion"]))
        if g["ionu"] is not None:
            iu = g["ionu"]; finite = np.isfinite(iu) & (iu >= 0)
            if not finite.any():
                raise ValueError(f"[{ad}] no valid ionosphere uncertainty samples")
            thr = np.percentile(iu[finite], 90)
            valid &= finite & (iu <= thr)
    print(f"\n=== [{ad}] track {chain[0]['track']} chain {chain[0]['d1']} -> {chain[-1]['d2']} ({len(chain)} pairs); ion uncertainty layer: {gs[0]['ionu'] is not None}")
    print(f"   valid (all pairs, ion-uncertainty p90 mask): {100*valid.mean():.1f}% of frame; {100*valid[strip].mean():.1f}% of Tay footprint")
    cum = np.zeros_like(g0["unw"])
    for g in gs:
        u = (g["unw"] - g["ion"]) * k
        u[~valid] = np.nan
        cum += deramp(u, stable)
    cum[~valid] = np.nan
    hp = highpass(cum, HP_KM, px_m)
    for lbl, a in (("deramped", cum), (f"deramped + {HP_KM:g} km high-pass", hp)):
        sv = np.array([win_mean(a, *to_pix(g0, lon, lat)) for _, lon, lat in sites])
        print(f"   [{lbl}] range-change mm; stable ground std = noise floor")
        print(f"      Tay footprint mean {np.nanmean(a[strip]):+.1f} std {np.nanstd(a[strip]):.1f} | stable ground mean {np.nanmean(a[stable]):+.1f} std {np.nanstd(a[stable]):.1f}")
        print(f"      dam sites n={np.isfinite(sv).sum()}: mean {np.nanmean(sv):+.1f} std {np.nanstd(sv):.1f} min {np.nanmin(sv):+.1f} max {np.nanmax(sv):+.1f}")
    write_tif(os.path.join(OUT, f"nisar_{ad}_disp60_hp{HP_KM:g}km_mm.tif"), g0, hp)
    coherence = np.mean([g["coh"] for g in gs], axis=0)
    coherence[~valid] = np.nan
    write_tif(os.path.join(OUT, f"nisar_{ad}_chain_coh.tif"), g0, coherence)
    with open(os.path.join(OUT, f"nisar_{ad}_sources.json"), "w") as stream:
        json.dump(dict(product="PROVISIONAL local-contrast diagnostic", hp_km=HP_KM,
                       wavelength_m=g0["wl"], pairs=chain), stream, indent=2)
    results[ad] = dict(g=g0, hp=hp, cum=cum, strip=strip, stable=stable)

if "D" in results and "A" in results:
    gd = results["D"]["g"]
    for fld in ("cum", "hp"):
        a_on_d = warp_onto(gd, mem_like(results["A"]["g"], results["A"][fld]), gdal.GRA_Bilinear)
        d = results["D"][fld]
        both = np.isfinite(a_on_d) & np.isfinite(d); strip = results["D"]["strip"]
        print(f"\n=== two-track agreement on {fld} (A warped onto D), common px {both.sum()}")
        for lbl, mk in (("Tay footprint", strip & both), ("all common", both)):
            if mk.sum() < 100: continue
            x, y = d[mk], a_on_d[mk]
            print(f"   [{lbl}] n={mk.sum()}: r={np.corrcoef(x, y)[0,1]:+.3f}; D std {x.std():.1f}; A std {y.std():.1f}; mean |D-A| {np.mean(np.abs(x-y)):.1f} mm")
        b = 12
        def blk(a):
            m = np.isfinite(a).astype(np.float32); s = np.where(np.isfinite(a), a, 0.0)
            S = uniform_filter(s, b); M = uniform_filter(m, b)
            return np.where(M > 0.8, S / np.maximum(M, 1e-6), np.nan)[b // 2::b, b // 2::b]
        D1, A1, S1 = blk(np.where(both, d, np.nan)), blk(np.where(both, a_on_d, np.nan)), blk(np.where(both & strip, 1.0, np.nan))
        mk = np.isfinite(D1) & np.isfinite(A1)
        if mk.sum() > 30:
            print(f"   [1 km blocks, all common] n={mk.sum()}: r={np.corrcoef(D1[mk], A1[mk])[0,1]:+.3f}")
            mk2 = mk & np.isfinite(S1)
            if mk2.sum() > 30:
                print(f"   [1 km blocks, Tay footprint] n={mk2.sum()}: r={np.corrcoef(D1[mk2], A1[mk2])[0,1]:+.3f}")
print("STAGE2C_DONE")
