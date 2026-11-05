"""Peat motion relative to a mineral-ground reference ensemble on native NISAR GUNW chains.

Pre-registration: docs/scientific_audit/NISAR_REFERENCE_ENSEMBLE_PREREG_20260910.md (read it
before changing thresholds). Regional signal preserved: no annulus, no high-pass, no deramp.

Env: TRACK=A|D, GUNW_DIR, REFSEL (frozen candidates json), REFSEL_GRID (label tif),
REFENS_OUT, SITES, PEAT_TIF, TAY_TIF, DEM_VRT.
Corrections are subtracted from the product's own layers (ionosphere at 80 m; hydrostatic, wet
and solid-Earth-tide phases from the metadata cube, interpolated at the DEM height). [E] the
unwrapped phase still contains them (ASF GUNW page; measured 2026-09-10).
Sign [S]: positive unwrapped phase is taken as range increase (ISCE lineage); LOS displacement
toward the sensor = -lambda/(4 pi) * phase; vertical = LOS / cos(incidence), no horizontal motion.
"""
import glob
import hashlib
import json
import os
import re
import numpy as np
import h5py
from osgeo import gdal, osr
from scipy import ndimage

gdal.UseExceptions()
TRACK = os.environ.get("TRACK", "D")
GUNW_DIR = os.path.expanduser(os.environ.get("GUNW_DIR", "/mnt/scratch/t/tommy/nisar/gunw"))
REFSEL = os.path.expanduser(os.environ.get("REFSEL", "~/full/nisar_refsel/refsel_frozen.json"))
REFSEL_GRID = os.path.expanduser(os.environ.get("REFSEL_GRID", "~/full/nisar_refsel/refsel_grid.tif"))
OUT = os.path.expanduser(os.environ.get("REFENS_OUT", "~/full/nisar_refens")); os.makedirs(OUT, exist_ok=True)
SITES = os.path.expanduser(os.environ.get("SITES", "~/full/dam_sites.geojson"))
PEAT = os.path.expanduser(os.environ.get("PEAT_TIF", "~/full/peat/PEATGRIDS_THICKNESS_CM_MEAN.tif"))
TAY = os.path.expanduser(os.environ.get("TAY_TIF", "~/full/tay/InSAR_LongTermRate.tif"))
DEM = os.path.expanduser(os.environ.get("DEM_VRT", "~/pilot/dem/elevation_3s.dem.wgs84.vrt"))
AOI = (113.7, -3.45, 115.3, -2.2)                 # W S E N of the target area
SITE_RADIUS_M = 500.0
MIN_REF_VALID = 0.5
PEAT_MIN_CM = 50.0
PATTERNS = {"A": "*_125_A_176_*_4000_*.h5", "D": "*_061_D_090_*_4000_*.h5"}
G = "science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/"
RG = "science/LSAR/GUNW/metadata/radarGrid/"
C = 299792458.0


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def mem_grid(x, y, epsg, dtype=gdal.GDT_Float32):
    m = gdal.GetDriverByName("MEM").Create("", len(x), len(y), 1, dtype)
    m.SetGeoTransform((float(x[0] - (x[1] - x[0]) / 2), float(x[1] - x[0]), 0.0, float(y[0] - (y[1] - y[0]) / 2), 0.0, float(y[1] - y[0])))
    s = osr.SpatialReference(); s.ImportFromEPSG(int(epsg)); m.SetProjection(s.ExportToWkt())
    return m


def warp_to(src_path, x, y, epsg, alg, dtype=gdal.GDT_Float32, nodata=None):
    src = gdal.Open(src_path); dst = mem_grid(x, y, epsg, dtype)
    if nodata is not None:
        dst.GetRasterBand(1).SetNoDataValue(nodata); dst.GetRasterBand(1).Fill(nodata)
    gdal.ReprojectImage(src, dst, src.GetProjection(), dst.GetProjection(), alg)
    return dst.GetRasterBand(1).ReadAsArray()


def write_tif(path, a, x, y, epsg, nodata=np.nan):
    d = gdal.GetDriverByName("GTiff").Create(path, a.shape[1], a.shape[0], 1, gdal.GDT_Float32, ["COMPRESS=DEFLATE", "TILED=YES"])
    m = mem_grid(x, y, epsg); d.SetGeoTransform(m.GetGeoTransform()); d.SetProjection(m.GetProjection())
    b = d.GetRasterBand(1); b.WriteArray(a.astype(np.float32)); b.SetNoDataValue(float(nodata)); d = None


def decode_mask(mask):
    low = mask & 255
    valid = ((low != 255) & (low < 100) & ((low // 10) > 0) & ((low % 10) > 0)
             & (((mask >> 8) & 65535) == 0) & ((mask >> 25) == 0))
    return valid, ((mask >> 24) & 1).astype(bool)


def relabel(cc, valid):
    """4-connected islands of valid samples within each product component."""
    out = np.zeros(cc.shape, np.int32); nxt = 1
    base = np.where(valid, cc, 0)
    for lab, sl in zip(range(1, int(base.max()) + 1), ndimage.find_objects(base)):
        if sl is None:
            continue
        sub = base[sl] == lab
        l, n = ndimage.label(sub, structure=[[0, 1, 0], [1, 1, 1], [0, 1, 0]])
        out[sl][sub] = (l[sub] + nxt - 1); nxt += n
    return out


def cube_at(h, name, xi, yi, hi):
    cube = h[RG + name][:]
    return ndimage.map_coordinates(cube, [hi, yi, xi], order=1, mode="nearest")


files = sorted(glob.glob(os.path.join(GUNW_DIR, PATTERNS[TRACK])))
if not files:
    raise SystemExit("no GUNW files for track " + TRACK)
pairs = []
for f in files:
    m = re.search(r"_(\d{8})T\d{6}_\d{8}T\d{6}_(\d{8})T\d{6}_", os.path.basename(f))
    pairs.append((m.group(1), m.group(2), f))
pairs.sort()
dates = sorted({d for p in pairs for d in p[:2]})
for i, (r, s, _) in enumerate(pairs):
    if i and r != pairs[i - 1][1]:
        raise SystemExit("pairs do not form a nearest-neighbour chain: " + str([p[:2] for p in pairs]))
print("track", TRACK, "pairs", [(p[0], p[1]) for p in pairs])

with h5py.File(pairs[0][2], "r") as h:
    x = h[G + "xCoordinates"][:]; y = h[G + "yCoordinates"][:]; epsg = int(h[G + "projection"][()])
    wl = C / float(h["science/LSAR/GUNW/grids/frequencyA/centerFrequency"][()])
    cx = h[RG + "xCoordinates"][:]; cy = h[RG + "yCoordinates"][:]; ch = h[RG + "heightAboveEllipsoid"][:]
ny, nx = len(y), len(x)
print("grid", ny, nx, "epsg", epsg, "wavelength %.4f m" % wl)
dem = warp_to(DEM, x, y, epsg, gdal.GRA_Bilinear); dem = np.where(np.isfinite(dem) & (dem > -100), dem, 0.0)
peat = warp_to(PEAT, x, y, epsg, gdal.GRA_Average); peat = np.where(np.isfinite(peat), peat, 0.0)
tay = warp_to(TAY, x, y, epsg, gdal.GRA_Average, nodata=0.0)
reflab = warp_to(REFSEL_GRID, x, y, epsg, gdal.GRA_NearestNeighbour, gdal.GDT_UInt16, 0).astype(np.int32)
ref = json.load(open(REFSEL))
if not ref.get("frozen"):
    raise SystemExit("reference candidates are not frozen (run nisar_refsel_geology.py first)")
cands = {c["rank"]: c for c in ref["candidates"]}
ens_ranks = [r for r, c in cands.items() if c["role"] == "ensemble"]
hold_ranks = [r for r, c in cands.items() if c["role"] == "holdout_control"]
print("references in ensemble", ens_ranks, "holdout", hold_ranks, "pixels per reference",
      {r: int((reflab == r).sum()) for r in cands})
# targets
xx, yy = np.meshgrid(x, y)
s = osr.SpatialReference(); s.ImportFromEPSG(4326); s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
d = osr.SpatialReference(); d.ImportFromEPSG(epsg); d.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
tr = osr.CoordinateTransformation(s, d)
ax0, ay0, _ = tr.TransformPoint(AOI[0], AOI[1]); ax1, ay1, _ = tr.TransformPoint(AOI[2], AOI[3])
aoi = (xx >= min(ax0, ax1)) & (xx <= max(ax0, ax1)) & (yy >= min(ay0, ay1)) & (yy <= max(ay0, ay1))
targets = {"peat_aoi": aoi & (peat >= PEAT_MIN_CM), "tay_footprint": aoi & (tay != 0) & np.isfinite(tay)}
for r in hold_ranks:
    targets["holdout_ref_%d" % r] = reflab == r
sites = json.load(open(SITES))["features"]
site_masks = {}
for i, ft in enumerate(sites):
    lon, lat = ft["geometry"]["coordinates"][:2]; sid = str(ft.get("properties", {}).get("site_id", ft.get("properties", {}).get("id", i + 1)))
    px, py, _ = tr.TransformPoint(lon, lat)
    c0, c1 = np.searchsorted(x, [px - SITE_RADIUS_M, px + SITE_RADIUS_M]); r0, r1 = np.searchsorted(-y, [-py - SITE_RADIUS_M, -py + SITE_RADIUS_M])
    if c0 <= 0 or r0 <= 0 or c1 >= nx or r1 >= ny:
        continue
    sl = (slice(r0, r1), slice(c0, c1)); m = (xx[sl] - px) ** 2 + (yy[sl] - py) ** 2 <= SITE_RADIUS_M ** 2
    site_masks[sid] = (sl, m, lon, lat)
print("targets", {k: int(v.sum()) for k, v in targets.items()}, "sites inside frame", len(site_masks))
crop = ndimage.find_objects(aoi.astype(np.int8))[0]
variants = ("raw", "ion", "all")
ref_mean = {v: {r: [] for r in cands} for v in variants}; ref_valid = {r: [] for r in cands}; ref_label = {r: [] for r in cands}
tgt = {v: {t: {r: [] for r in cands} for t in targets} for v in variants}; tgt_n = {t: {r: [] for r in cands} for t in targets}; tgt_label = {t: [] for t in targets}
site = {v: {sid: {r: [] for r in cands} for sid in site_masks} for v in variants}; site_n = {sid: {r: [] for r in cands} for sid in site_masks}
cum_field = {v: np.zeros(aoi[crop].shape, np.float64) for v in variants}; field_ok = np.zeros(aoi[crop].shape, np.int16)
# amendment A4: mineral-plane datum. Per pair, a plane through the usable ensemble references' mean
# phases at their centroids (needs >= 3); targets are read relative to that plane. Not a peat filter:
# only mineral ground constrains the surface. Reported beside the constant (median) datum.
ref_xy = {}
for r in cands:
    m = reflab == r
    ref_xy[r] = (float(xx[m].mean()), float(yy[m].mean())) if m.any() else (np.nan, np.nan)
plane_tgt = {v: {t: [] for t in targets} for v in variants}; plane_site = {v: {sid: [] for sid in site_masks} for v in variants}; plane_hold = {v: [] for v in variants}; plane_n = []
plane_field = {v: np.zeros(aoi[crop].shape, np.float64) for v in variants}; plane_ok = np.zeros(aoi[crop].shape, np.int16)
diag = []; source_records = []
for ridx, (rd, sd, f) in enumerate(pairs):
    with h5py.File(f, "r") as h:
        unw = h[G + "HH/unwrappedPhase"][:].astype(np.float64); ion = h[G + "HH/ionospherePhaseScreen"][:].astype(np.float64)
        cc = h[G + "HH/connectedComponents"][:]; mask = h[G + "mask"][:]; coh = h[G + "HH/coherenceMagnitude"][:]
        valid, interp = decode_mask(mask)
        valid &= np.isfinite(unw) & np.isfinite(ion) & (cc > 0) & (cc < 65535)
        xi = np.interp(xx.ravel(), cx, np.arange(len(cx))); yi = np.interp(-yy.ravel(), -cy, np.arange(len(cy))); hi = np.interp(dem.ravel(), ch, np.arange(len(ch)))
        hyd = cube_at(h, "hydrostaticTroposphericPhaseScreen", xi, yi, hi).reshape(ny, nx)
        wet = cube_at(h, "wetTroposphericPhaseScreen", xi, yi, hi).reshape(ny, nx)
        set_ = cube_at(h, "slantRangeSolidEarthTidesPhase", xi, yi, hi).reshape(ny, nx)
        if ridx == 0:
            inc = cube_at(h, "incidenceAngle", xi, yi, hi).reshape(ny, nx)
        del xi, yi, hi
    source_records.append(dict(file=os.path.basename(f), bytes=os.path.getsize(f), pair=[rd, sd]))
    lab = relabel(cc, valid)
    ph = {"raw": np.where(valid, unw, np.nan), "ion": np.where(valid, unw - ion, np.nan), "all": np.where(valid, unw - ion - hyd - wet - set_, np.nan)}
    # long-wavelength diagnostic (20 km) on the largest component: does each correction lower the residual?
    big = np.bincount(lab[valid].ravel()).argmax()
    def lp(a, w=250):
        m = np.isfinite(a) & (lab == big); b = np.where(m, a, 0.0); n = ndimage.uniform_filter(m.astype(np.float32), w); v = ndimage.uniform_filter(b.astype(np.float32), w)
        return np.where(n > 0.3, v / np.maximum(n, 1e-6), np.nan)
    lw = {v: lp(ph[v]) for v in variants}
    diag.append(dict(pair=[rd, sd], largest_component_px=int((lab == big).sum()), valid_fraction=float(valid.mean()),
                     ion_interp_fraction_valid=float(interp[valid].mean()),
                     lowpass20km_std_rad={v: float(np.nanstd(lw[v])) for v in variants},
                     tropo_set_std_rad=float(np.nanstd((hyd + wet + set_)[valid]))))
    print("pair %s-%s valid %.3f largest comp %d px; 20 km low-pass std raw %.2f ion %.2f all %.2f rad" % (
        rd, sd, valid.mean(), (lab == big).sum(), diag[-1]["lowpass20km_std_rad"]["raw"], diag[-1]["lowpass20km_std_rad"]["ion"], diag[-1]["lowpass20km_std_rad"]["all"]))
    for r in cands:
        m = (reflab == r) & valid
        frac = m.sum() / max((reflab == r).sum(), 1); ref_valid[r].append(float(frac))
        l = int(np.bincount(lab[m]).argmax()) if m.any() else 0; ref_label[r].append(l)
        for v in variants:
            ref_mean[v][r].append(float(np.mean(ph[v][m & (lab == l)])) if frac >= MIN_REF_VALID and l > 0 else np.nan)
    for t, tm in targets.items():
        mv = tm & valid; tl = int(np.bincount(lab[mv]).argmax()) if mv.any() else 0; tgt_label[t].append(tl)
        for r in cands:
            same = mv & (lab == ref_label[r][-1]) if ref_label[r][-1] > 0 else np.zeros_like(mv)
            tgt_n[t][r].append(int(same.sum()))
            for v in variants:
                tgt[v][t][r].append(float(np.nanmedian(ph[v][same]) - ref_mean[v][r][-1]) if same.sum() >= 10 and np.isfinite(ref_mean[v][r][-1]) else np.nan)
    for sid, (sl, m, lon, lat) in site_masks.items():
        mv = m & valid[sl]; ll = lab[sl]
        for r in cands:
            same = mv & (ll == ref_label[r][-1]) if ref_label[r][-1] > 0 else np.zeros_like(mv)
            site_n[sid][r].append(int(same.sum()))
            for v in variants:
                site[v][sid][r].append(float(np.nanmedian(ph[v][sl][same]) - ref_mean[v][r][-1]) if same.sum() >= 10 and np.isfinite(ref_mean[v][r][-1]) else np.nan)
    # amendment A4: plane datum per pair
    usable = [r for r in ens_ranks if np.isfinite(ref_mean["all"][r][-1])]
    plane_n.append(len(usable))
    for v in variants:
        if len(usable) >= 3:
            Amat = np.array([[ref_xy[r][0], ref_xy[r][1], 1.0] for r in usable]); bvec = np.array([ref_mean[v][r][-1] for r in usable])
            coef = np.linalg.lstsq(Amat, bvec, rcond=None)[0]
            pl = lambda px_, py_: coef[0] * px_ + coef[1] * py_ + coef[2]
            for t, tm in targets.items():
                mv = tm & valid
                plane_tgt[v][t].append(float(np.nanmedian(ph[v][mv] - pl(xx[mv], yy[mv]))) if mv.sum() >= 10 else np.nan)
            for sid, (sl, m, lon, lat) in site_masks.items():
                mv = m & valid[sl]
                plane_site[v][sid].append(float(np.nanmedian(ph[v][sl][mv] - pl(xx[sl][mv], yy[sl][mv]))) if mv.sum() >= 10 else np.nan)
            for r in hold_ranks:
                mh = (reflab == r) & valid
                plane_hold[v].append(float(np.nanmean(ph[v][mh] - pl(xx[mh], yy[mh]))) if mh.sum() >= 10 else np.nan)
            a = ph[v][crop] - pl(xx[crop], yy[crop]); ok = np.isfinite(a); plane_field[v][ok] += a[ok]
            if v == "all":
                plane_ok += ok
            if v == "ion":                       # per-pair plane-relative field, for residual-versus-motion diagnosis
                write_tif(os.path.join(OUT, f"refens_{TRACK}_ion_plane_pair_{rd}_{sd}_los_mm.tif"), np.where(ok, a * (-wl / (4 * np.pi) * 1000.0), np.nan), x[crop[1]], y[crop[0]], epsg)
        else:
            for t in targets: plane_tgt[v][t].append(np.nan)
            for sid in site_masks: plane_site[v][sid].append(np.nan)
            plane_hold[v].append(np.nan)
    # cumulative field relative to the ensemble median reference (all-corrections variant)
    for v in variants:
        rm = np.nanmedian([ref_mean[v][r][-1] for r in ens_ranks])
        a = ph[v][crop] - rm; ok = np.isfinite(a)
        cum_field[v][ok] += a[ok]
        if v == "all":
            field_ok += ok
    del unw, ion, cc, mask, coh, valid, interp, hyd, wet, set_, lab, ph, lw

# chain to epochs (earliest = 0), ensemble statistics, leave-one-out
scale = -wl / (4 * np.pi) * 1000.0              # mm of LOS displacement toward the sensor per radian [S]
def chain(vals):
    out = [0.0]; c = 0.0
    for v in vals:
        c = c + v if np.isfinite(v) and np.isfinite(c) else np.nan; out.append(c)
    return out
def stats(series_by_ref):
    ens = np.array([series_by_ref[r] for r in ens_ranks], float) * scale
    med = np.nanmedian(ens, axis=0); q1, q3 = np.nanpercentile(ens, [25, 75], axis=0)
    loo = [np.nanmedian(np.delete(ens, i, axis=0), axis=0) for i in range(len(ens_ranks))]
    return dict(epochs=dates, ensemble_median_mm=med.tolist(), iqr_mm=(q3 - q1).tolist(), range_mm=(np.nanmax(ens, 0) - np.nanmin(ens, 0)).tolist(),
                per_reference_mm={str(r): (np.array(series_by_ref[r]) * scale).tolist() for r in ens_ranks},
                leave_one_out_final_mm=[float(l[-1]) for l in loo], n_refs_final=int(np.isfinite(ens[:, -1]).sum()))
inc_mean = float(np.nanmean(inc[aoi])); cos_inc = float(np.cos(np.deg2rad(inc_mean)))
report = dict(schema_version=1, track=TRACK, pairs=[[p[0], p[1]] for p in pairs], epochs=dates, wavelength_m=wl, epsg=epsg,
              sign_convention="[S] positive phase = range increase; LOS mm positive toward sensor = -lambda/(4pi)*phase; vertical = LOS/cos(inc), no horizontal motion",
              incidence_mean_aoi_deg=inc_mean, sources=source_records, refsel=dict(path=REFSEL, sha256=sha(REFSEL), candidates_sha256=ref.get("candidates_sha256")),
              references={str(r): dict(role=cands[r]["role"], valid_fraction_per_pair=ref_valid[r], component_per_pair=ref_label[r]) for r in cands},
              reference_consistency={v: stats({r: chain([ref_mean[v][r][i] - np.nanmedian([ref_mean[v][q][i] for q in ens_ranks]) for i in range(len(pairs))]) for r in cands}) for v in variants},
              targets={v: {t: stats({r: chain(tgt[v][t][r]) for r in cands}) for t in targets} for v in variants},
              target_support={t: dict(component_per_pair=tgt_label[t], pixels_shared_with_reference={str(r): tgt_n[t][r] for r in cands}) for t in targets},
              sites={v: {sid: stats({r: chain(site[v][sid][r]) for r in cands}) for sid in site_masks} for v in variants},
              site_support={sid: dict(lon=site_masks[sid][2], lat=site_masks[sid][3], pixels_shared_with_reference={str(r): site_n[sid][r] for r in cands}) for sid in site_masks},
              holdout={v: {t: stats({r: chain(tgt[v][t][r]) for r in cands}) for t in targets if t.startswith("holdout")} for v in variants},
              diagnostics=diag, flags=dict(absolute_datum=False, reference_assumed_stable=True, ground_accuracy_validated=False, annual_velocity_ready=False, release_ready=False, restoration_attribution=False),
              prereg="docs/scientific_audit/NISAR_REFERENCE_ENSEMBLE_PREREG_20260910.md")
# per-reference cumulative values (not relative to the ensemble) for the consistency table
report["reference_cumulative_mm"] = {v: {str(r): (np.array(chain(ref_mean[v][r])) * scale).tolist() for r in cands} for v in variants}
report["plane_datum"] = dict(usable_references_per_pair=plane_n,
                             targets={v: {t: (np.array(chain(plane_tgt[v][t])) * scale).tolist() for t in targets} for v in variants},
                             holdout={v: (np.array(chain(plane_hold[v])) * scale).tolist() for v in variants},
                             sites={v: {sid: (np.array(chain(plane_site[v][sid])) * scale).tolist() for sid in site_masks} for v in variants})
for v in variants:
    print(f"[{v}] plane datum (refs per pair {plane_n}): holdout final {report['plane_datum']['holdout'][v][-1]:+.1f} mm; " + "; ".join(f"{t} {report['plane_datum']['targets'][v][t][-1]:+.1f}" for t in targets))
json.dump(report, open(os.path.join(OUT, f"refens_{TRACK}.json"), "w"), indent=1)
for v in variants:
    a = np.where(field_ok == len(pairs), cum_field[v] * scale, np.nan)
    write_tif(os.path.join(OUT, f"refens_{TRACK}_{v}_los_mm.tif"), a, x[crop[1]], y[crop[0]], epsg)
write_tif(os.path.join(OUT, f"refens_{TRACK}_all_vertical_mm.tif"), np.where(field_ok == len(pairs), cum_field["all"] * scale / cos_inc, np.nan), x[crop[1]], y[crop[0]], epsg)
write_tif(os.path.join(OUT, f"refens_{TRACK}_valid_pairs.tif"), field_ok.astype(np.float32), x[crop[1]], y[crop[0]], epsg)
for v in variants:                                # plane-datum cumulative field, pixels valid in every pair with a plane
    write_tif(os.path.join(OUT, f"refens_{TRACK}_{v}_plane_los_mm.tif"), np.where(plane_ok == sum(1 for n in plane_n if n >= 3), plane_field[v] * scale, np.nan), x[crop[1]], y[crop[0]], epsg)
# console summary against the pre-registered reading
for v in variants:
    rc = report["reference_consistency"][v]
    print(f"[{v}] reference consistency at final epoch: IQR {rc['iqr_mm'][-1]:.1f} mm, range {rc['range_mm'][-1]:.1f} mm, refs {rc['n_refs_final']}")
    for t in targets:
        s_ = report["targets"][v][t]
        print(f"[{v}] {t}: final {s_['ensemble_median_mm'][-1]:+.1f} mm LOS (IQR {s_['iqr_mm'][-1]:.1f}); leave-one-out finals {np.round(s_['leave_one_out_final_mm'], 1).tolist()}")
print("REFENS_OK", TRACK)
