"""Select mineral-ground reference candidates from evidence layers only (no phase read).

Pre-registration: docs/scientific_audit/NISAR_REFERENCE_ENSEMBLE_PREREG_20260910.md.
Evidence: Copernicus 3 arcsec DEM, PeatGrids thickness mean, ESA WorldCover 2021 (where tiles
exist). Geology (Macrostrat) is attached afterwards on a machine with internet access by
nisar_refsel_geology.py; the candidate list is frozen with a SHA-256 after that step.

Outputs (REFSEL_OUT): refsel_grid.tif (uint16 region labels, EPSG:4326, 200 m), refsel.json.
"""
import hashlib
import json
import os
import sys
import glob
import h5py
import numpy as np
from osgeo import gdal, osr
from scipy import ndimage

gdal.UseExceptions()
OUT = os.path.expanduser(os.environ.get("REFSEL_OUT", "~/full/nisar_refsel"))
DEM = os.path.expanduser(os.environ.get("DEM_VRT", "~/pilot/dem/elevation_3s.dem.wgs84.vrt"))
PEAT = os.path.expanduser(os.environ.get("PEAT_TIF", "~/full/peat/PEATGRIDS_THICKNESS_CM_MEAN.tif"))
WC_DIR = os.path.expanduser(os.environ.get("WORLDCOVER_DIR", "~/full/peat/worldcover"))
GUNW_DIR = os.path.expanduser(os.environ.get("GUNW_DIR", "/mnt/scratch/t/tommy/nisar/gunw"))
FOOTPRINT_FILES = {"D": "*_061_D_090_*_4000_*.h5", "A": "*_125_A_176_*_4000_*.h5"}       # amendment A3: coverage only
FOOT_TRACKS = os.environ.get("FOOT_TRACKS", "DA")                                        # amendment A4: per-track lists
BUILTUP_MIN_KM2 = 1.0
ZAHRO_POINT = (113.920443, -2.213489)           # InaCORS-checked stable point, Zahro et al. 2025 [E]
ZAHRO_KM = 3.0
AOI_CENTROID = (114.55, -2.72)
LAND = os.path.expanduser(os.environ.get("LAND_SHP", "~/full/peat/ne_10m_land/ne_10m_land.shp"))   # amendment A2
# common area of descending 61/90 and ascending 125/176, slightly inside the catalog extents
W, S, E, N = 112.35, -4.55, 115.10, -1.62
STEP = 0.0018                                   # about 200 m
MIN_AREA_KM2 = 4.0
MIN_SEP_KM = 10.0                               # amendment A3
MAX_CAND = 12                                   # amendment A3
MAX_AREA_KM2 = 25.0                             # amendment A1: compact datum; 5 km block inside large regions
BLOCK_CELLS = 25                                # 5 km at 200 m
SEED = 7
os.makedirs(OUT, exist_ok=True)
nx = int(round((E - W) / STEP)); ny = int(round((N - S) / STEP))
gt = (W, STEP, 0.0, N, 0.0, -STEP)
srs = osr.SpatialReference(); srs.ImportFromEPSG(4326); wkt = srs.ExportToWkt()


def target(dtype=gdal.GDT_Float32, nodata=None):
    m = gdal.GetDriverByName("MEM").Create("", nx, ny, 1, dtype)
    m.SetGeoTransform(gt); m.SetProjection(wkt)
    if nodata is not None:
        m.GetRasterBand(1).SetNoDataValue(nodata)
    return m


def warp(src_path, alg, dtype=gdal.GDT_Float32, nodata=None):
    src = gdal.Open(src_path)
    dst = target(dtype, nodata)
    gdal.ReprojectImage(src, dst, src.GetProjection(), wkt, alg)
    return dst.GetRasterBand(1).ReadAsArray()


print("grid", nx, "x", ny, "cells at", STEP, "deg")
dem = warp(DEM, gdal.GRA_Average)
dem = np.where(np.isfinite(dem) & (dem > -100), dem, np.nan)
peat_src = gdal.Open(PEAT); pnd = peat_src.GetRasterBand(1).GetNoDataValue()
peat = warp(PEAT, gdal.GRA_Average)
print("PeatGrids nodata", pnd, "grid stats: min %.1f max %.1f, zero fraction %.3f, nan fraction %.3f"
      % (np.nanmin(peat), np.nanmax(peat), np.mean(peat == 0), np.mean(~np.isfinite(peat))))
if pnd is not None:
    peat = np.where(peat == pnd, 0.0, peat)
peat = np.where(np.isfinite(peat), peat, 0.0)
tiles = sorted(f for f in os.listdir(WC_DIR) if f.endswith("_Map.tif"))
vrt = gdal.BuildVRT(os.path.join(OUT, "worldcover.vrt"), [os.path.join(WC_DIR, t) for t in tiles])
vrt = None
wc = warp(os.path.join(OUT, "worldcover.vrt"), gdal.GRA_Mode, gdal.GDT_Byte, 0)
print("WorldCover tiles", tiles, "coverage fraction %.3f" % np.mean(wc > 0))

# rule 1: no peat within 1 km (PeatGrids mean 0 cm); erosion by 5 cells = 1 km
peat0 = peat <= 0.0
peat0_1km = ndimage.binary_erosion(peat0, structure=np.ones((11, 11), bool))
# rule 2: relief: height >= 30 m or >= 20 m above the 5 km minimum
demf = np.where(np.isfinite(dem), dem, 1e4)
local_min = ndimage.minimum_filter(demf, size=51)
relief = np.isfinite(dem) & ((dem >= 30.0) | (dem - local_min >= 20.0))
# rule 3: WorldCover not water / wetland / mangrove where available
wc_bad = np.isin(wc, (80, 90, 95))
# amendment A2: on land, at least 1 km inside the Natural Earth 10 m coastline
land_ds = target(gdal.GDT_Byte, 0)
gdal.Rasterize(land_ds, LAND, burnValues=[1])
land = land_ds.GetRasterBand(1).ReadAsArray() > 0
land_1km = ndimage.binary_erosion(land, structure=np.ones((11, 11), bool))
print("land fraction of grid %.3f" % land.mean())
# amendment A3: inside the valid swath of both tracks (first pair's product mask, coverage only)
def footprint(pattern):
    f = sorted(glob.glob(os.path.join(GUNW_DIR, pattern)))[0]
    with h5py.File(f, "r") as h:
        g = h["science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram"]
        mask = g["mask"][:]; x = g["xCoordinates"][:]; y = g["yCoordinates"][:]; epsg = int(g["projection"][()])
    low = mask & 255
    valid = ((low != 255) & (low < 100) & ((low // 10) > 0) & ((low % 10) > 0) & (((mask >> 8) & 65535) == 0) & ((mask >> 25) == 0))
    src = gdal.GetDriverByName("MEM").Create("", len(x), len(y), 1, gdal.GDT_Byte)
    src.SetGeoTransform((float(x[0] - 40), 80.0, 0.0, float(y[0] + 40), 0.0, -80.0))
    ss = osr.SpatialReference(); ss.ImportFromEPSG(epsg); src.SetProjection(ss.ExportToWkt())
    src.GetRasterBand(1).WriteArray(valid.astype(np.uint8))
    dst = target(gdal.GDT_Float32); gdal.ReprojectImage(src, dst, src.GetProjection(), wkt, gdal.GRA_Average)
    return dst.GetRasterBand(1).ReadAsArray() > 0.5, os.path.basename(f)
foot = {}; foot_src = {}
for trk, pat in FOOTPRINT_FILES.items():
    if trk not in FOOT_TRACKS:
        continue
    foot[trk], foot_src[trk] = footprint(pat)
    print("footprint", trk, foot_src[trk], "fraction of grid %.3f" % foot[trk].mean())
both = ndimage.binary_erosion(np.logical_and.reduce([foot[t] for t in foot]), structure=np.ones((11, 11), bool))
print("both-track footprint after 1 km erosion: %.3f of grid" % both.mean())
qualify = peat0_1km & relief & ~wc_bad & land_1km & both
cell_km2 = (STEP * 111.0) * (STEP * 111.0 * np.cos(np.deg2rad(-3.0)))
yy, xx = np.mgrid[0:ny, 0:nx]
lon_c = W + (xx + 0.5) * STEP; lat_c = N - (yy + 0.5) * STEP
# amendment A3: built-up class (any relief): WorldCover 50 regions >= 1 km2 on PeatGrids 0 or near the Zahro point
bu_lab, bu_n = ndimage.label(wc == 50, structure=[[0, 1, 0], [1, 1, 1], [0, 1, 0]])
zahro = np.hypot((lon_c - ZAHRO_POINT[0]) * np.cos(np.deg2rad(-2.2)), lat_c - ZAHRO_POINT[1]) * 111.0 <= ZAHRO_KM
builtup = np.zeros_like(qualify)
for i, sl in enumerate(ndimage.find_objects(bu_lab), start=1):
    if sl is None:
        continue
    m = np.zeros_like(qualify); m[sl] = bu_lab[sl] == i
    if m.sum() * cell_km2 < BUILTUP_MIN_KM2 or not (m & land_1km & both).mean() > 0:
        continue
    if peat0[m].mean() >= 0.5 or zahro[m].any():
        builtup |= m & land_1km & both
print("built-up candidate cells: %d" % builtup.sum())
qualify |= builtup
print("qualifying cells: %d (%.2f%%)" % (qualify.sum(), 100 * qualify.mean()))

labels, nlab = ndimage.label(qualify, structure=[[0, 1, 0], [1, 1, 1], [0, 1, 0]])
cell_km2 = (STEP * 111.0) * (STEP * 111.0 * np.cos(np.deg2rad(-3.0)))
sizes = ndimage.sum(np.ones_like(labels), labels, index=np.arange(1, nlab + 1))
keep = [i + 1 for i, s in enumerate(sizes) if s * cell_km2 >= MIN_AREA_KM2]
print("regions", nlab, "with area >= %.0f km2: %d" % (MIN_AREA_KM2, len(keep)))
# amendment A1: reduce large regions to the 5 km block with the most qualifying cells (evidence only)
reduced = 0
for lab in keep:
    m = labels == lab
    if m.sum() * cell_km2 <= MAX_AREA_KM2:
        continue
    cnt = ndimage.uniform_filter(m.astype(np.float32), BLOCK_CELLS, mode="constant") * BLOCK_CELLS ** 2
    hgt = ndimage.uniform_filter(np.where(m, np.nan_to_num(dem, nan=0.0), 0.0).astype(np.float32), BLOCK_CELLS, mode="constant")
    score = cnt + 1e-3 * hgt
    r, c = np.unravel_index(int(np.argmax(np.where(m, score, -1))), m.shape)
    h = BLOCK_CELLS // 2
    block = np.zeros_like(m); block[max(r - h, 0):r + h + 1, max(c - h, 0):c + h + 1] = True
    labels[m & ~block] = 0; reduced += 1
print("regions reduced to a 5 km block:", reduced)
yy, xx = np.mgrid[0:ny, 0:nx]
lon_c = W + (xx + 0.5) * STEP; lat_c = N - (yy + 0.5) * STEP
cands = []
for lab in keep:
    m = labels == lab
    lon0, lat0 = float(lon_c[m].mean()), float(lat_c[m].mean())
    hist = {int(k): int(v) for k, v in zip(*np.unique(wc[m], return_counts=True))}
    cands.append(dict(label=int(lab), area_km2=float(m.sum() * cell_km2), lon=lon0, lat=lat0,
                      dem_median_m=float(np.nanmedian(dem[m])), dem_min_m=float(np.nanmin(dem[m])),
                      prominence_m=float(np.nanmedian((dem - local_min)[m])),
                      peat_mean_cm=float(peat[m].mean()), worldcover_hist=hist,
                      worldcover_unavailable_fraction=float(np.mean(wc[m] == 0))))
cands.sort(key=lambda c: -c["area_km2"])


def dist_km(a, b):
    return float(np.hypot((a["lon"] - b["lon"]) * 111.0 * np.cos(np.deg2rad(-3.0)), (a["lat"] - b["lat"]) * 111.0))


for c in cands:                                    # amendment A3: nearest to the peat AOI first
    c["dist_to_aoi_km"] = dist_km(c, {"lon": AOI_CENTROID[0], "lat": AOI_CENTROID[1]})
    m = labels == c["label"]
    c["class"] = "builtup" if (builtup & m).mean() > 0 and (builtup & m).sum() >= 0.5 * m.sum() else "relief"
cands.sort(key=lambda c: c["dist_to_aoi_km"])
chosen = []
for c in cands:
    if len(chosen) >= MAX_CAND:
        break
    if all(dist_km(c, k) >= MIN_SEP_KM for k in chosen):
        chosen.append(c)
rng = np.random.default_rng(SEED)
holdout = int(rng.integers(len(chosen))) if chosen else None
for i, c in enumerate(chosen):
    c["rank"] = i + 1
    c["role"] = "holdout_control" if i == holdout else "ensemble"
    print("cand %d label %d %s area %.1f km2 lon %.4f lat %.4f dem %.0f m prom %.0f m dist %.0f km wc %s role %s" % (
        c["rank"], c["label"], c["class"], c["area_km2"], c["lon"], c["lat"], c["dem_median_m"], c["prominence_m"], c["dist_to_aoi_km"],
        c["worldcover_hist"], c["role"]))

lab_out = np.zeros_like(labels, dtype=np.uint16)
for c in chosen:
    lab_out[labels == c["label"]] = c["rank"]
d = gdal.GetDriverByName("GTiff").Create(os.path.join(OUT, "refsel_grid.tif"), nx, ny, 1, gdal.GDT_UInt16, ["COMPRESS=DEFLATE"])
d.SetGeoTransform(gt); d.SetProjection(wkt); d.GetRasterBand(1).WriteArray(lab_out); d.GetRasterBand(1).SetNoDataValue(0); d = None
report = dict(schema_version=1, stage="evidence_only_no_phase_read", grid=dict(epsg=4326, west=W, south=S, east=E, north=N, step_deg=STEP),
              rules=dict(max_area_km2=MAX_AREA_KM2, block_km=BLOCK_CELLS * STEP * 111.0, peat_mean_cm_max=0.0, peat_buffer_km=1.0, dem_min_m=30.0, prominence_min_m=20.0, prominence_window_km=5.0,
                         worldcover_excluded=[80, 90, 95], land_inside_coast_km=1.0, footprint_both_tracks_eroded_km=1.0, builtup_min_km2=BUILTUP_MIN_KM2, zahro_point=ZAHRO_POINT, zahro_km=ZAHRO_KM, nearest_first=True, footprint_tracks=FOOT_TRACKS, min_area_km2=MIN_AREA_KM2, min_separation_km=MIN_SEP_KM, max_candidates=MAX_CAND,
                         holdout_seed=SEED, coherence_used=False, phase_used=False),
              sources=dict(dem=DEM, peatgrids=PEAT, worldcover_tiles=tiles, land=LAND, footprints=foot_src),
              n_qualifying_regions=len(keep), candidates=chosen, all_regions=cands[:50], geology="pending (nisar_refsel_geology.py)")
with open(os.path.join(OUT, "refsel.json"), "w") as f:
    json.dump(report, f, indent=1)
print("REFSEL_OK", len(chosen), "candidates; holdout rank", (holdout + 1) if holdout is not None else None)
