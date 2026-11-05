"""NISAR GUNW stage 2: dam-site coherence, NISAR/ALOS-2 coherence ratio, and a
first 60-day dry-season displacement composite from the descending chain.

Layout confirmed by stage 1 (provisional GUNW):
  science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/HH/{coherenceMagnitude,
  unwrappedPhase, connectedComponents, ionospherePhaseScreen}, UTM grids at 80 m.
The 60-day composite is a magnitude-and-pattern check, not a rate: signal over
64 days is 0.5 to 4 mm of secular subsidence plus any dry-season shrinkage,
against 1 to 2 cm of per-pair noise. Sign convention is reported as found and
must be confirmed against the product specification before any "sinking" claim.
"""
import glob
import json
import os
import re
import numpy as np
import h5py
from osgeo import gdal, osr

gdal.UseExceptions()
D = os.path.expanduser(os.environ.get("GUNW_DIR", "/mnt/scratch/t/tommy/nisar/gunw"))
OUT = os.path.expanduser("~/full/nisar_out"); os.makedirs(OUT, exist_ok=True)
TAY = os.path.expanduser("~/full/tay/InSAR_LongTermRate.tif")
GAMMA = os.path.expanduser("~/full/sigphi/mean_gamma.npy")
TC = os.path.expanduser("~/full/mintpy_tay/temporalCoherence.h5")
SITES = os.path.expanduser("~/full/dam_sites.geojson")
REF_LL = (115.13458333, -3.18458333)   # our mineral-ground reference (lon, lat)
G = "science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/HH/"
WIN = 2  # 5x5 window at 80 m = 400 m


def find_ds(h, name, prefer):
    hits = []
    h.visititems(lambda n, o: hits.append(n) if isinstance(o, h5py.Dataset) and n.split("/")[-1] == name else None)
    for p in hits:
        if p.startswith(prefer):
            return p
    return hits[0] if hits else None


def read_grid(path):
    with h5py.File(path, "r") as h:
        coh = h[G + "coherenceMagnitude"][:].astype(np.float32)
        unw = h[G + "unwrappedPhase"][:].astype(np.float32)
        cc = h[G + "connectedComponents"][:]
        ion = h[G + "ionospherePhaseScreen"][:].astype(np.float32) if (G + "ionospherePhaseScreen") in h else None
        x = h[find_ds(h, "xCoordinates", G)][:]; y = h[find_ds(h, "yCoordinates", G)][:]
        ep = find_ds(h, "epsg", G) or find_ds(h, "projection", G)
        v = h[ep][()]
        try:
            epsg = int(v)
        except Exception:
            s = v.decode() if isinstance(v, bytes) else str(v); epsg = int(re.findall(r"\d{4,5}", s)[-1])
        wl = None
        for cand in ("wavelength", "centerFrequency", "processedCenterFrequency"):
            p = find_ds(h, cand, "science")
            if p is not None:
                val = float(np.ravel(h[p][()])[0]); wl = val if cand == "wavelength" else 299792458.0 / val; break
    return dict(coh=coh, unw=unw, cc=cc, ion=ion, x=x, y=y, epsg=epsg, wl=wl or 0.2384)


def geotransform(g):
    dx = float(g["x"][1] - g["x"][0]); dy = float(g["y"][1] - g["y"][0])
    return (float(g["x"][0]) - dx / 2, dx, 0.0, float(g["y"][0]) - dy / 2, 0.0, dy)


def warp_onto(g, src_ds, alg=gdal.GRA_Average):
    ny, nx = g["coh"].shape
    mem = gdal.GetDriverByName("MEM").Create("", nx, ny, 1, gdal.GDT_Float32)
    mem.SetGeoTransform(geotransform(g))
    srs = osr.SpatialReference(); srs.ImportFromEPSG(g["epsg"]); mem.SetProjection(srs.ExportToWkt())
    mem.GetRasterBand(1).Fill(float("nan"))
    gdal.ReprojectImage(src_ds, mem, src_ds.GetProjection(), srs.ExportToWkt(), alg)
    return mem.GetRasterBand(1).ReadAsArray().astype(np.float32)


def latlon_mem(arr, x0, y0, dx, dy):
    ny, nx = arr.shape
    mem = gdal.GetDriverByName("MEM").Create("", nx, ny, 1, gdal.GDT_Float32)
    mem.SetGeoTransform((x0, dx, 0.0, y0, 0.0, dy))
    srs = osr.SpatialReference(); srs.ImportFromEPSG(4326); mem.SetProjection(srs.ExportToWkt())
    mem.GetRasterBand(1).WriteArray(np.where(np.isfinite(arr), arr, np.nan).astype(np.float32))
    mem.GetRasterBand(1).SetNoDataValue(float("nan"))
    return mem


def to_pix(g, lon, lat):
    src = osr.SpatialReference(); src.ImportFromEPSG(4326); src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    dst = osr.SpatialReference(); dst.ImportFromEPSG(g["epsg"]); dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    t = osr.CoordinateTransformation(src, dst)
    X, Y, _ = t.TransformPoint(float(lon), float(lat))
    dx = float(g["x"][1] - g["x"][0]); dy = float(g["y"][1] - g["y"][0])
    c = int(round((X - float(g["x"][0])) / dx)); r = int(round((Y - float(g["y"][0])) / dy))
    return r, c


def win_mean(a, r, c, w=WIN):
    ny, nx = a.shape
    if r - w < 0 or c - w < 0 or r + w >= ny or c + w >= nx:
        return np.nan
    p = a[r - w:r + w + 1, c - w:c + w + 1]
    return float(np.nanmean(p)) if np.isfinite(p).any() else np.nan


# ---- inputs shared across files
sites = [(f["properties"].get("site_id"), *f["geometry"]["coordinates"][:2]) for f in json.load(open(SITES))["features"]]
gamma = np.load(GAMMA).astype(np.float32); gamma[gamma <= 0.0011] = np.nan
with h5py.File(TC, "r") as f:
    tc = f["temporalCoherence"][:].astype(np.float32); att = dict(f.attrs)
X0, DX, Y0, DY = float(att["X_FIRST"]), float(att["X_STEP"]), float(att["Y_FIRST"]), float(att["Y_STEP"])
gamma_ds = latlon_mem(gamma, X0, Y0, DX, DY)
stable_ds = latlon_mem((tc >= 0.7).astype(np.float32), X0, Y0, DX, DY)
tay_ds = gdal.Open(TAY)

files = sorted(glob.glob(os.path.join(D, "*.h5")))
print(f"{len(files)} files; {len(sites)} dam sites")
per_file = {}
for f in files:
    name = os.path.basename(f)
    m = re.search(r"_(\d{3})_(\d{3})_([AD])_(\d{3})_\d{3}_\d{4}_SH_(\d{8})T\d{6}_\d{8}T\d{6}_(\d{8})T", name)
    track, frame, ad, d1, d2 = m.group(2), m.group(4), m.group(3), m.group(5), m.group(6)
    g = read_grid(f)
    tay = warp_onto(g, tay_ds); tay[tay == 0] = np.nan; strip = np.isfinite(tay)
    gam = warp_onto(g, gamma_ds); stable = warp_onto(g, stable_ds) > 0.5
    ok = np.isfinite(g["coh"]) & (g["coh"] > 0)
    r0, c0 = to_pix(g, *REF_LL)
    ref_ok = 0 <= r0 < g["coh"].shape[0] and 0 <= c0 < g["coh"].shape[1]
    site_coh = np.array([win_mean(g["coh"], *to_pix(g, lon, lat)) for _, lon, lat in sites])
    site_cc = np.array([win_mean((g["cc"] > 0).astype(np.float32), *to_pix(g, lon, lat)) for _, lon, lat in sites])
    both = strip & ok & np.isfinite(gam) & (gam > 0.02)
    ratio = g["coh"][both] / gam[both]
    print(f"\n=== {ad} track {track} frame {frame}  {d1} -> {d2}  epsg {g['epsg']}  wl {g['wl']:.4f} m")
    print(f"   Tay footprint: coh mean {g['coh'][strip & ok].mean():.3f}; ALOS-2 mean gamma on same px {np.nanmean(gam[strip & ok]):.3f}; "
          f"ratio median {np.nanmedian(ratio):.2f} (p10 {np.nanpercentile(ratio,10):.2f}, p90 {np.nanpercentile(ratio,90):.2f})")
    ins = np.isfinite(site_coh)
    print(f"   dam sites inside frame: {ins.sum()}/{len(sites)}; coherence mean {np.nanmean(site_coh):.3f}, "
          f"sites >0.3: {np.nanmean(site_coh[ins] > 0.3)*100:.0f}%, >0.5: {np.nanmean(site_coh[ins] > 0.5)*100:.0f}%; "
          f"unwrapped (conncomp>0) at sites: {np.nanmean(site_cc[ins])*100:.0f}%")
    print(f"   stable-ground (ALOS-2 tempCoh>=0.7) px in frame: {int(stable.sum())}; their NISAR coh mean {np.nanmean(g['coh'][stable & ok]):.3f}")
    print(f"   reference pixel inside frame: {ref_ok}; coh there {win_mean(g['coh'], r0, c0):.3f}" if ref_ok else "   reference pixel outside frame")
    per_file[name] = dict(g=g, strip=strip, stable=stable, ok=ok, ad=ad, track=track, d1=d1, d2=d2, ref=(r0, c0) if ref_ok else None,
                          site_pix=[to_pix(g, lon, lat) for _, lon, lat in sites])

# ---- 60-day composite on the descending chain
desc = sorted([v for v in per_file.values() if v["ad"] == "D"], key=lambda v: v["d1"])
chain = []
for v in desc:
    if not chain or chain[-1]["d2"] == v["d1"]:
        chain.append(v)
if len(chain) >= 2 and all(v["ref"] for v in chain):
    g0 = chain[0]["g"]; wl = g0["wl"]; k = wl / (4 * np.pi) * 1000.0
    cum_raw = np.zeros_like(g0["unw"]); cum_ion = np.zeros_like(g0["unw"]); valid = np.ones(g0["unw"].shape, bool)
    for v in chain:
        g = v["g"]; r0, c0 = v["ref"]
        u = g["unw"].copy(); u[g["cc"] <= 0] = np.nan
        ref_val = win_mean(u, r0, c0)
        u = u - ref_val
        cum_raw += np.nan_to_num(u); valid &= np.isfinite(u)
        if g["ion"] is not None:
            ii = g["ion"] - win_mean(g["ion"], r0, c0)
            cum_ion += np.nan_to_num(u - ii)
        else:
            cum_ion += np.nan_to_num(u)
    cum_raw[~valid] = np.nan; cum_ion[~valid] = np.nan
    dr = cum_raw * k; di = cum_ion * k
    strip = chain[0]["strip"]; stable = chain[0]["stable"]
    print(f"\n=== descending chain {chain[0]['d1']} -> {chain[-1]['d2']} ({len(chain)} pairs), LOS mm, referenced to mineral ground")
    print(f"   valid through all pairs: {100*valid.mean():.1f}% of frame; {100*valid[strip].mean():.1f}% of Tay footprint")
    for lbl, a in (("raw", dr), ("ion-corrected", di)):
        print(f"   [{lbl}] Tay footprint mean {np.nanmean(a[strip]):+.1f}  std {np.nanstd(a[strip]):.1f} | stable ground mean {np.nanmean(a[stable]):+.1f} std {np.nanstd(a[stable]):.1f}")
        sv = np.array([win_mean(a, r, c) for r, c in chain[0]["site_pix"]])
        print(f"          dam sites (n={np.isfinite(sv).sum()}): mean {np.nanmean(sv):+.1f}  std {np.nanstd(sv):.1f}  min {np.nanmin(sv):+.1f} max {np.nanmax(sv):+.1f}")
    drv = gdal.GetDriverByName("GTiff")
    for lbl, a in (("disp60_raw_mm", dr), ("disp60_ion_mm", di)):
        ds = drv.Create(os.path.join(OUT, f"nisar_{lbl}.tif"), a.shape[1], a.shape[0], 1, gdal.GDT_Float32, ["COMPRESS=DEFLATE", "TILED=YES"])
        ds.SetGeoTransform(geotransform(g0)); srs = osr.SpatialReference(); srs.ImportFromEPSG(g0["epsg"]); ds.SetProjection(srs.ExportToWkt())
        ds.GetRasterBand(1).WriteArray(a); ds.GetRasterBand(1).SetNoDataValue(float("nan")); ds = None
    cohmean = np.nanmean(np.stack([v["g"]["coh"] for v in chain]), axis=0)
    ds = drv.Create(os.path.join(OUT, "nisar_coh_mean_desc.tif"), cohmean.shape[1], cohmean.shape[0], 1, gdal.GDT_Float32, ["COMPRESS=DEFLATE", "TILED=YES"])
    ds.SetGeoTransform(geotransform(g0)); ds.SetProjection(srs.ExportToWkt()); ds.GetRasterBand(1).WriteArray(cohmean); ds = None
    print(f"   wrote {OUT}/nisar_disp60_raw_mm.tif, nisar_disp60_ion_mm.tif, nisar_coh_mean_desc.tif")
    print("   SIGN: raw product convention, unconfirmed; see product spec before reading as up/down")
else:
    print("\nno complete descending chain with the reference inside the frame")
print("STAGE2_DONE")
