#!/usr/bin/env python3
"""Suspected degraded peat map and fire-recurrence animation, PeatGuard dam-site AOI (v2).

Run from the repo root:
    .venv/bin/python3 deliverables/degraded_peat_fire_recurrence/build.py [--skip-download]

Inputs
  deliverables/aoi_3d_degradation/            work grid (EPSG:4326, 0.0006 deg, 1000 x 1050), MODIS fire JSONs,
                                              NBR epochs, canals, rivers, Tay 2025 rate, sites, settlements
  grid_products/external/peatgrids/PEATGRIDS_THICKNESS_CM_MEAN.tif   PeatGrids (Widyastuti et al. 2022)
  grid_products/external/firms_viirs/raw/     NASA FIRMS Indonesia yearly CSVs (VIIRS S-NPP 2012 to 2024,
                                              VIIRS NOAA-20 / JPSS-1 2018 to 2024), downloaded on demand
Outputs (this directory)
  map_suspected_degraded_peat.png, map_fire_recurrence_2001_2024.png,
  fires_return_to_degraded_peat_2001_2026.mp4 / .gif, stats.json, README.md, fires_viirs/ (clipped cache)

Sources and tags are in README.md.
"""
import json
import math
import os
import shutil
import subprocess
import sys
import textwrap
import urllib.request

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.spatial import cKDTree
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
from PIL import Image
from skimage.morphology import skeletonize

HERE = os.path.dirname(os.path.abspath(__file__))
# DATA_ROOT is the PeatGuard data root holding deliverables/ and grid_products/
# (the prepared work-grid layers, FIRMS archive and PeatGrids raster are not in
# this repository). Outputs are written next to the inputs.
ROOT = os.environ.get("DATA_ROOT", os.path.abspath(os.path.join(HERE, "..", "..")))
IN = os.path.join(ROOT, "deliverables", "aoi_3d_degradation")
OUT = os.path.join(ROOT, "deliverables", "degraded_peat_fire_recurrence")
FRAMES = os.path.join(OUT, "frames")
VIIRS_CACHE = os.path.join(OUT, "fires_viirs")
VIIRS_RAW = os.path.join(ROOT, "grid_products", "external", "firms_viirs", "raw")
PEAT_TIF = os.path.join(ROOT, "grid_products", "external", "peatgrids", "PEATGRIDS_THICKNESS_CM_MEAN.tif")
FFMPEG = os.environ.get("FFMPEG", shutil.which("ffmpeg") or "ffmpeg")
FIRMS_BASE = "https://firms.modaps.eosdis.nasa.gov/data/country"
FIRMS_SOURCES = {"snpp": ("viirs-snpp", range(2012, 2025)), "jpss1": ("viirs-jpss1", range(2018, 2025))}

# Work grid
NR, NC = 1000, 1050
LON0, LAT0, STEP = 114.2, -2.41, 0.0006
LON1, LAT1 = LON0 + NC * STEP, LAT0 - NR * STEP
EXTENT = [LON0, LON1, LAT1, LAT0]
LAT_MID = 0.5 * (LAT0 + LAT1)
KM_PER_DEG_LAT = 111.2
KM_PER_DEG_LON = 111.32 * math.cos(math.radians(LAT_MID))
CELL_KM = STEP * KM_PER_DEG_LAT
CELL_KM2 = (STEP * KM_PER_DEG_LAT) * (STEP * KM_PER_DEG_LON)

# Fire grids: 500 m (headline, mixed sensor) and 1 km (MODIS-only reference)
FSTEP = 0.0045
FR, FC = math.ceil(NR * STEP / FSTEP), math.ceil(NC * STEP / FSTEP)          # 134 x 140
F_EXTENT = [LON0, LON0 + FC * FSTEP, LAT0 - FR * FSTEP, LAT0]
F_CELL_KM2 = (FSTEP * KM_PER_DEG_LAT) * (FSTEP * KM_PER_DEG_LON)
KSTEP = 0.009
KR, KC = math.ceil(NR * STEP / KSTEP), math.ceil(NC * STEP / KSTEP)          # 67 x 70
K_EXTENT = [LON0, LON0 + KC * KSTEP, LAT0 - KR * KSTEP, LAT0]

# Definitions
PEAT_MIN_CM = 50.0
CANAL_DILATE_CELLS = 15
NBR_LOW = 0.45
NBR_MIN_EPOCHS = 2
NBR_EPOCHS = [1990, 1995, 2000, 2005, 2010, 2015, 2020, 2022, 2023, 2024, 2025]
NBR_EXCLUDED = [2026]
TAY_SINK_MMYR = -5.0
FIRE_YEARS = list(range(2001, 2025))
MODIS_YEARS = list(range(2001, 2012))
VIIRS_YEARS = list(range(2012, 2025))
EL_NINO = [2002, 2006, 2009, 2015, 2019, 2023]
RECUR_N = 3
MODIS_FOOT_KM = 1.0
VIIRS_FOOT_KM = 0.375
DEDUP_KM = 0.375
DEDUP_HOURS = 1.0

# Rendering
SUBFRAMES = 4
FPS = 12
HOLD_FRAMES = 11            # about 0.92 s per year
LAST_HOLD_FRAMES = 36       # 3 s
GIF_FPS = 3
C_RED = "#c8102e"
C_RED_EDGE = "#8a0a1f"
RED_ALPHA = 0.28
C_PEAT = (0.914, 0.863, 0.753)   # #e9dcc0
C_WATER = (0.72, 0.82, 0.92)
C_CANAL = "#5b8fc9"
C_RIVER = "#1f4e8c"
C_SITE = "#1b7f3b"
C_SITE_BURN = "#ff3b3b"
C_INK = "#222222"
C_INK2 = "#555555"
C_CRIT = {"crit_canal": "#b5651d", "crit_burn": "#8c1d18", "crit_sink": "#3b3b8f"}
FIRE_LEVELS = [0.75, 2.0, 4.0]      # in units of one isolated detection's smoothed peak (detection equivalents)
FIRE_AREA_CAP = 3.0                  # drawn glow area <= cap x true footprint area
SIGMA_KM = {"MODIS": 0.9, "VIIRS": 0.5}
FIRE_RGB = [(0.60, 0.05, 0.03), (0.95, 0.40, 0.08), (1.00, 0.72, 0.10)]   # edge, mid, core
FIRE_ALPHA = [0.55, 0.75, 0.90]
CHAR_RGB = (0.122, 0.122, 0.122)     # #1f1f1f
CANAL_BANDS_KM = [0, 0.5, 1, 2, 3, 5, 1e9]

SETTLEMENT_LABELS = ["Pulang Pisau", "Kuala Kapuas", "Marabahan", "Dadahup", "Lamunti Baru",
                     "Sidomulyo", "Tanjung Harapan", "Petak Batuah"]

CAPTION_DEF = ("Suspected degraded peat [S] = PeatGrids thickness >= 50 cm AND (within 1 km of a Dadap 2021 canal OR "
               "Landsat NBR < 0.45 in >= 2 of 11 epochs 1990 to 2025 OR Tay 2025 ALOS-2 rate < -5 mm/yr). "
               "Screening layer, not ground-truthed.")
CAPTION_SRC = ("Sources: PeatGrids (Widyastuti et al. 2022); canals (Dadap et al. 2021); Landsat NBR, thresholds after "
               "Key and Benson 2006; subsidence (Tay et al. 2025); fires NASA FIRMS; hillshade from the AOI DEM. PeatGuard, 2026.")


def load(name):
    return np.load(os.path.join(IN, name + ".npy"))


def jload(name):
    with open(os.path.join(IN, name + ".json")) as fh:
        return json.load(fh)


def disk(r):
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    return (x * x + y * y) <= r * r


def lonlat_to_rc(lon, lat, step):
    return (np.floor((LAT0 - lat) / step).astype(int), np.floor((lon - LON0) / step).astype(int))


def work_to_fire_index(step, nrows, ncols):
    """Flat fire-grid index for every work cell (by cell centre)."""
    rr = np.floor(((np.arange(NR) + 0.5) * STEP) / step).astype(int)
    cc = np.floor(((np.arange(NC) + 0.5) * STEP) / step).astype(int)
    rr = np.clip(rr, 0, nrows - 1)
    cc = np.clip(cc, 0, ncols - 1)
    return (rr[:, None] * ncols + cc[None, :])


def block_count(mask, idx, nrows, ncols):
    return np.bincount(idx.ravel()[mask.ravel()], minlength=nrows * ncols).reshape(nrows, ncols)


def wrap_lines(text, width):
    return textwrap.wrap(text, width=width)


# ---------------------------------------------------------------- peat and criteria
def load_peat_thickness():
    with rasterio.open(PEAT_TIF) as src:
        crs, nodata = src.crs, src.nodata
        dst = np.full((NR, NC), np.nan, np.float32)
        reproject(rasterio.band(src, 1), dst, dst_transform=from_origin(LON0, LAT0, STEP, STEP),
                  dst_crs="EPSG:4326", resampling=Resampling.bilinear, src_nodata=nodata, dst_nodata=np.nan)
    return dst, str(crs), nodata


def build_masks():
    hill = load("hillshade").astype(np.float32)
    canals = load("canals").astype(bool)
    rivers = load("rivers").astype(bool)
    water = load("water_mask").astype(bool)
    tay = load("tay_rate_mmyr").astype(np.float32)
    thick, crs, nodata = load_peat_thickness()
    nbr = np.stack([load(f"nbr_{y}") for y in NBR_EPOCHS]).astype(np.float32)
    valid = ~np.isnan(nbr)
    n_valid = valid.sum(axis=0)
    n_low = ((nbr < NBR_LOW) & valid).sum(axis=0)
    nbr26 = load("nbr_2026").astype(np.float32)
    n_low_incl26 = n_low + ((nbr26 < NBR_LOW) & ~np.isnan(nbr26))
    crit_canal = ndimage.binary_dilation(canals, structure=disk(CANAL_DILATE_CELLS))
    crit_burn = n_low >= NBR_MIN_EPOCHS
    crit_sink = np.where(np.isnan(tay), False, tay < TAY_SINK_MMYR)
    is_water = water | rivers
    peat = (thick >= PEAT_MIN_CM) & ~is_water
    any_crit = crit_canal | crit_burn | crit_sink
    return dict(hill=hill, canals=canals, rivers=rivers, water=is_water, tay=tay, thick=thick, peat_crs=crs,
                peat_nodata=nodata, n_valid=n_valid, n_low=n_low, n_low_incl26=n_low_incl26,
                crit_canal=crit_canal, crit_burn=crit_burn, crit_sink=crit_sink, peat=peat,
                degraded=peat & any_crit, intact=peat & ~any_crit, nonpeat=~peat & ~is_water)


def classify(m, idx, nrows, ncols):
    """Majority class per fire cell among land work cells: 0 degraded, 1 intact, 2 non-peat, -1 no land."""
    stack = np.stack([block_count(m["degraded"], idx, nrows, ncols), block_count(m["intact"], idx, nrows, ncols),
                      block_count(m["nonpeat"], idx, nrows, ncols)])
    cls = np.argmax(stack, axis=0)
    cls[stack.sum(axis=0) == 0] = -1
    return cls


# ---------------------------------------------------------------- fires
def download_firms(skip=False):
    os.makedirs(VIIRS_RAW, exist_ok=True)
    report = {}
    for key, (folder, years) in FIRMS_SOURCES.items():
        for y in years:
            fn = f"{folder}_{y}_Indonesia.csv"
            path = os.path.join(VIIRS_RAW, fn)
            if os.path.exists(path) and os.path.getsize(path) > 1000:
                report[fn] = "cached"
                continue
            if skip:
                report[fn] = "missing (download skipped)"
                continue
            url = f"{FIRMS_BASE}/{folder}/{y}/{fn}"
            try:
                urllib.request.urlretrieve(url, path)
                report[fn] = "downloaded"
            except Exception as exc:   # noqa: BLE001
                report[fn] = f"failed: {exc}"
    return report


def clip_viirs_year(y):
    """Clip the Indonesia CSVs to the AOI, keep vegetation-fire type and nominal/high confidence, dedupe NOAA-20."""
    os.makedirs(VIIRS_CACHE, exist_ok=True)
    cache = os.path.join(VIIRS_CACHE, f"viirs_{y}.json")
    if os.path.exists(cache):
        with open(cache) as fh:
            return json.load(fh)
    parts = {}
    for key, (folder, years) in FIRMS_SOURCES.items():
        if y not in years:
            continue
        path = os.path.join(VIIRS_RAW, f"{folder}_{y}_Indonesia.csv")
        if not os.path.exists(path):
            continue
        cols = ["latitude", "longitude", "acq_date", "acq_time", "confidence", "type"]
        kept, n_raw, n_type, n_lowconf = [], 0, 0, 0
        for chunk in pd.read_csv(path, usecols=cols, chunksize=500_000, dtype={"acq_time": str, "confidence": str}):
            box = chunk[(chunk.longitude >= LON0) & (chunk.longitude < LON1) & (chunk.latitude <= LAT0) & (chunk.latitude > LAT1)]
            n_raw += len(box)
            veg = box[box["type"] == 0]
            n_type += len(box) - len(veg)
            good = veg[veg.confidence.str.lower() != "l"]
            n_lowconf += len(veg) - len(good)
            kept.append(good)
        df = pd.concat(kept) if kept else pd.DataFrame(columns=cols)
        t = pd.to_datetime(df.acq_date + " " + df.acq_time.str.zfill(4), format="%Y-%m-%d %H%M", errors="coerce")
        parts[key] = dict(lon=df.longitude.to_numpy(float), lat=df.latitude.to_numpy(float),
                          hours=(t.astype("int64") / 3.6e12).to_numpy(float),
                          n_raw_in_box=int(n_raw), n_dropped_type=int(n_type), n_dropped_lowconf=int(n_lowconf))
    if "snpp" not in parts:
        return None
    s = parts["snpp"]
    lon, lat, hours = s["lon"], s["lat"], s["hours"]
    sat = np.zeros(len(lon), int)
    n_dup = 0
    if "jpss1" in parts and len(parts["jpss1"]["lon"]):
        j = parts["jpss1"]
        if len(lon):
            tree = cKDTree(np.column_stack([(lon - LON0) * KM_PER_DEG_LON, (lat - LAT0) * KM_PER_DEG_LAT]))
            q = np.column_stack([(j["lon"] - LON0) * KM_PER_DEG_LON, (j["lat"] - LAT0) * KM_PER_DEG_LAT])
            neigh = tree.query_ball_point(q, r=DEDUP_KM)
            keep = np.ones(len(q), bool)
            for i, nb in enumerate(neigh):
                if nb and np.any(np.abs(hours[nb] - j["hours"][i]) <= DEDUP_HOURS):
                    keep[i] = False
        else:
            keep = np.ones(len(j["lon"]), bool)
        n_dup = int((~keep).sum())
        lon = np.concatenate([lon, j["lon"][keep]])
        lat = np.concatenate([lat, j["lat"][keep]])
        hours = np.concatenate([hours, j["hours"][keep]])
        sat = np.concatenate([sat, np.ones(int(keep.sum()), int)])
    rec = dict(year=y, sensor="VIIRS", n=int(len(lon)), n_snpp=int(len(s["lon"])),
               n_jpss1_kept=int((sat == 1).sum()), n_jpss1_duplicates_dropped=n_dup,
               snpp_raw_in_box=s["n_raw_in_box"], snpp_dropped_type=s["n_dropped_type"],
               snpp_dropped_lowconf=s["n_dropped_lowconf"],
               jpss1_raw_in_box=parts.get("jpss1", {}).get("n_raw_in_box", 0),
               jpss1_dropped_type=parts.get("jpss1", {}).get("n_dropped_type", 0),
               jpss1_dropped_lowconf=parts.get("jpss1", {}).get("n_dropped_lowconf", 0),
               points=np.column_stack([lon, lat]).round(5).tolist(), satellite=sat.tolist())
    with open(cache, "w") as fh:
        json.dump(rec, fh)
    return rec


def modis_points(year):
    pts = np.asarray(jload(f"fires_{year}"), dtype=float)
    return pts.reshape(-1, 2)


def load_fire_year(y):
    if y in VIIRS_YEARS:
        rec = clip_viirs_year(y)
        if rec is not None:
            return dict(year=y, sensor="VIIRS", pts=np.asarray(rec["points"], float).reshape(-1, 2), meta=rec)
    pts = modis_points(y)
    return dict(year=y, sensor="MODIS", pts=pts, meta=dict(year=y, sensor="MODIS", n=int(len(pts))))


def week_points():
    w = jload("fires_2026_week")
    cells = np.asarray(w["cells"], dtype=float)
    rng = np.random.default_rng(2026)
    pts = [np.column_stack([rng.uniform(lon, lon + 0.05, int(n)), rng.uniform(lat, lat + 0.05, int(n))])
           for lon, lat, n, _ in cells]
    pts = np.vstack(pts)
    inside = (pts[:, 0] >= LON0) & (pts[:, 0] < LON1) & (pts[:, 1] <= LAT0) & (pts[:, 1] > LAT1)
    return pts[inside], w, cells


def burned_grid(pts, sensor, step, nrows, ncols):
    """Binary 'burned this year' grid. MODIS marks every cell whose centre lies within 500 m of the detection
    (1 km footprint); VIIRS marks the cell containing the detection (375 m footprint)."""
    g = np.zeros((nrows, ncols), bool)
    if len(pts) == 0:
        return g
    if sensor == "MODIS" and step < KSTEP:
        reach = int(math.ceil(0.5 * MODIS_FOOT_KM / (step * KM_PER_DEG_LAT)))
        rc_lat = LAT0 - (np.arange(nrows) + 0.5) * step
        rc_lon = LON0 + (np.arange(ncols) + 0.5) * step
        r0, c0 = lonlat_to_rc(pts[:, 0], pts[:, 1], step)
        for dr in range(-reach, reach + 1):
            for dc in range(-reach, reach + 1):
                r, c = r0 + dr, c0 + dc
                ok = (r >= 0) & (r < nrows) & (c >= 0) & (c < ncols)
                dy = (rc_lat[r[ok]] - pts[ok, 1]) * KM_PER_DEG_LAT
                dx = (rc_lon[c[ok]] - pts[ok, 0]) * KM_PER_DEG_LON
                hit = (dx * dx + dy * dy) <= (0.5 * MODIS_FOOT_KM) ** 2
                g[r[ok][hit], c[ok][hit]] = True
        return g
    r, c = lonlat_to_rc(pts[:, 0], pts[:, 1], step)
    ok = (r >= 0) & (r < nrows) & (c >= 0) & (c < ncols)
    g[r[ok], c[ok]] = True
    return g


def class_stats(fy, cls, cell_km2):
    out = {}
    for name, code in [("suspected_degraded_peat", 0), ("intact_peat", 1), ("non_peat", 2)]:
        sel = cls == code
        n = int(sel.sum())
        v = fy[sel]
        out[name] = dict(cells=n, area_km2=round(n * cell_km2, 1),
                         frac_ge1_fire_years=float((v >= 1).mean()) if n else None,
                         frac_ge3_fire_years=float((v >= 3).mean()) if n else None,
                         frac_ge5_fire_years=float((v >= 5).mean()) if n else None,
                         mean_fire_years_per_cell=float(v.mean()) if n else None,
                         max_fire_years=int(v.max()) if n else None)
    return out


# ---------------------------------------------------------------- rendering helpers
def burn_field(pts, sensor, seed):
    """Organic burn-patch field on the work grid, in detection equivalents (1.0 = the smoothed peak of one
    isolated detection footprint). Returns the field, the per-year sigma factor and the true footprint area (km2)."""
    rng = np.random.default_rng(seed)
    factor = float(rng.uniform(0.8, 1.3))
    sigma = SIGMA_KM[sensor] / CELL_KM * factor
    foot = MODIS_FOOT_KM if sensor == "MODIS" else VIIRS_FOOT_KM
    half = max(1, int(round(0.5 * foot / CELL_KM)))
    if len(pts) == 0:
        return np.zeros((NR, NC), np.float32), factor, 0.0
    jit = pts + rng.normal(0.0, 0.3 / KM_PER_DEG_LAT, pts.shape)
    r, c = lonlat_to_rc(jit[:, 0], jit[:, 1], STEP)
    ok = (r >= 0) & (r < NR) & (c >= 0) & (c < NC)
    ind = np.zeros((NR, NC), bool)
    ind[r[ok], c[ok]] = True
    ind = ndimage.maximum_filter(ind, size=2 * half + 1)
    area_km2 = float(ind.sum() * CELL_KM2)
    field = ndimage.gaussian_filter(ind.astype(np.float32), sigma)
    one = np.zeros((4 * half + 8 * int(sigma) + 3,) * 2, np.float32)
    cc = one.shape[0] // 2
    one[cc - half:cc + half + 1, cc - half:cc + half + 1] = 1.0
    p0 = float(ndimage.gaussian_filter(one, sigma).max())
    field /= p0
    noise = ndimage.gaussian_filter(rng.normal(size=(NR, NC)).astype(np.float32), 6.0)
    noise /= noise.std()
    return (field * np.clip(1.0 + 0.35 * noise, 0.35, 1.8)).astype(np.float32), factor, area_km2


def fire_levels(field, area_km2):
    """Fixed bands, with the outer threshold raised if the drawn glow would exceed FIRE_AREA_CAP x footprint area."""
    edge, mid, core = FIRE_LEVELS
    if area_km2 > 0:
        max_cells = FIRE_AREA_CAP * area_km2 / CELL_KM2
        if (field >= edge).sum() > max_cells:
            vals = np.sort(field[field >= edge])[::-1]
            edge = float(vals[int(max_cells) - 1])
            mid = max(mid, 1.5 * edge)
            core = max(core, 3.0 * edge)
    return [edge, mid, core]


def fire_rgba(field, levels, seed, sub):
    """Banded RGBA fire image with low-frequency flicker (alpha x 0.7 to 1.3, slight core shift)."""
    rng = np.random.default_rng(seed * 100 + sub + 7)
    n1 = ndimage.gaussian_filter(rng.normal(size=(NR, NC)).astype(np.float32), 4.0)
    n1 /= n1.std()
    n2 = ndimage.gaussian_filter(rng.normal(size=(NR, NC)).astype(np.float32), 4.0)
    n2 /= n2.std()
    f = field * (1.0 + 0.12 * n2)
    afac = np.clip(1.0 + 0.3 * n1, 0.7, 1.3)
    out = np.zeros((NR, NC, 4), np.float32)
    for lvl, rgb, a in zip(levels, FIRE_RGB, FIRE_ALPHA):
        sel = f >= lvl
        out[sel, :3] = rgb
        out[sel, 3] = a
    out[..., 3] = np.clip(out[..., 3] * afac, 0, 1)
    return out


def char_layer(cum_counts, idx):
    """Persistent solid charcoal fill on 500 m cells burned in 2 or more earlier years; darker at 4 or more."""
    up = cum_counts.ravel()[idx]
    out = np.zeros((NR, NC, 4), np.float32)
    for sel, a in [((up >= 2) & (up <= 3), 0.50), ((up >= 4), 0.75)]:
        out[sel, :3] = CHAR_RGB
        out[sel, 3] = a
    return out


def canal_distance_stats(m, fy500, idx500):
    """Mean fire years per work cell (2001 to 2024, 500 m grid mapped to the work grid) by distance to the nearest canal."""
    dist = ndimage.distance_transform_edt(~m["canals"], sampling=(STEP * KM_PER_DEG_LAT, STEP * KM_PER_DEG_LON))
    fy = fy500.ravel()[idx500]
    land = ~m["water"]
    bands = []
    for lo, hi in zip(CANAL_BANDS_KM[:-1], CANAL_BANDS_KM[1:]):
        sel = land & (dist >= lo) & (dist < hi)
        n = int(sel.sum())
        v = fy[sel]
        bands.append(dict(band_km=[lo, (None if hi > 1000 else hi)], label=(f"{lo:g} to {hi:g}" if hi < 1000 else f"over {lo:g}"),
                          n_cells=n, area_km2=round(n * CELL_KM2, 1),
                          mean_fire_years=(float(v.mean()) if n else None),
                          frac_ge1=(float((v >= 1).mean()) if n else None),
                          frac_ge3=(float((v >= 3).mean()) if n else None)))
    return dict(bands=bands, max_distance_km=float(dist[land].max()), note="work grid cells (0.0006 deg), water excluded; "
                "fire years from the 500 m mixed-sensor grid; canals.npy is a generous (max-resampled) Dadap 2021 mask")


def draw_canal_bars(ax, cd, dark=False, fontsize=9, annotate=True):
    ink = "#dddddd" if dark else C_INK
    ink2 = "#aaaaaa" if dark else C_INK2
    labels = [b["label"] for b in cd["bands"]]
    vals = [b["mean_fire_years"] or 0 for b in cd["bands"]]
    ns = [b["n_cells"] for b in cd["bands"]]
    x = np.arange(len(vals))
    ax.bar(x, vals, width=0.7, color="#d7301f", edgecolor="none", zorder=3)
    if annotate:
        for xi, v, n in zip(x, vals, ns):
            ax.text(xi, v + 0.03, f"n={n:,}", fontsize=fontsize - 2, ha="center", va="bottom", color=ink2, rotation=90 if dark else 0)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=fontsize, color=ink, rotation=30, ha="right")
    ax.tick_params(axis="y", labelsize=fontsize, colors=ink2, length=0)
    ax.set_ylim(0, max(vals) * (1.45 if dark else 1.25))
    ax.yaxis.grid(True, color="#333333" if dark else "#dddddd", lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(ink2)
    ax.set_xlabel("distance to nearest canal, km", fontsize=fontsize, color=ink)


def make_fig_canal_distance(cd):
    fig = plt.figure(figsize=(8, 5.6), dpi=200)
    ax = fig.add_axes([0.11, 0.22, 0.86, 0.62])
    draw_canal_bars(ax, cd, fontsize=10)
    ax.set_ylabel("mean fire years per cell, 2001 to 2024", fontsize=10, color=C_INK)
    ax.set_title("Fire recurrence against distance to the nearest mapped canal, PeatGuard dam-site AOI",
                 fontsize=11, color=C_INK, loc="left")
    txt = ("Bars: mean number of fire years 2001 to 2024 per work-grid cell (0.0006 deg), fire from the 500 m mixed-sensor grid "
           "(MODIS 1 km to 2011, VIIRS 375 m from 2012), water excluded; n = cells per band. Canal distance from a Euclidean "
           "distance transform of the Dadap et al. 2021 canal mask (max-resampled, so distances are lower bounds). "
           "Canal-fire association in the ex-MRP: Hoscilo et al. 2011, Int. J. Wildland Fire [E]; these band numbers are our "
           "own screening statistics [S]. PeatGuard, 2026.")
    for k, line in enumerate(wrap_lines(txt, 125)):
        fig.text(0.02, 0.085 - k * 0.022, line, fontsize=7.5, color=C_INK2, ha="left", va="top")
    path = os.path.join(OUT, "fire_years_vs_canal_distance.png")
    fig.savefig(path, dpi=200, facecolor="white")
    plt.close(fig)
    return path


def thin_lines(mask, width=1):
    sk = skeletonize(mask)
    if width > 1:
        sk = ndimage.binary_dilation(sk, iterations=width - 1)
    return sk


def rgba_layer(mask, color, alpha):
    out = np.zeros((NR, NC, 4), np.float32)
    out[mask, :3] = color
    out[mask, 3] = alpha
    return out


def compose_base(m):
    """Hillshade grey base, intact peat pale tan, water blue, suspected degraded peat solid red with a dark edge."""
    h = m["hill"]
    h = 0.80 + 0.20 * (h - h.min()) / (h.max() - h.min())
    rgb = np.ones((NR, NC, 3), np.float32) * 0.86
    rgb[m["peat"]] = C_PEAT
    rgb = rgb * h[..., None]
    rgb[m["water"]] = C_WATER
    red = np.array(matplotlib.colors.to_rgb(C_RED), np.float32)
    d = m["degraded"]
    rgb[d] = (1 - RED_ALPHA) * rgb[d] + RED_ALPHA * red
    edge = d & ~ndimage.binary_erosion(d, iterations=1)
    rgb[edge] = matplotlib.colors.to_rgb(C_RED_EDGE)
    return np.clip(rgb, 0, 1)


def draw_lines(ax, m, canals=False, canal_alpha=0.9):
    if canals:
        ax.imshow(rgba_layer(thin_lines(m["canals"]), matplotlib.colors.to_rgb(C_CANAL), canal_alpha),
                  extent=EXTENT, interpolation="nearest", zorder=4)
    ax.imshow(rgba_layer(thin_lines(m["rivers"], 2), matplotlib.colors.to_rgb(C_RIVER), 1.0),
              extent=EXTENT, interpolation="nearest", zorder=4)


def draw_sites_labels(ax, sites, settlements, site_size=70, label_size=11, site_lw=0.8, labels=True):
    burning = np.array([s["burning"] for s in sites])
    lon = np.array([s["lon"] for s in sites])
    lat = np.array([s["lat"] for s in sites])
    ax.scatter(lon[~burning], lat[~burning], marker="^", s=site_size, c=C_SITE, edgecolors="white",
               linewidths=site_lw, zorder=8)
    ax.scatter(lon[burning], lat[burning], marker="^", s=site_size * 1.3, c=C_SITE_BURN, edgecolors="white",
               linewidths=site_lw, zorder=9)
    if not labels:
        return
    halo = [pe.withStroke(linewidth=3, foreground="white")]
    for s in settlements:
        if s["name"] not in SETTLEMENT_LABELS or not (LON0 < s["lon"] < LON1 and LAT1 < s["lat"] < LAT0):
            continue
        dy = 0.008 if s["lat"] < LAT1 + 0.02 else -0.008
        ax.plot(s["lon"], s["lat"], "o", ms=4, color=C_INK, zorder=9)
        ax.text(s["lon"], s["lat"] + dy, s["name"], fontsize=label_size, color=C_INK, ha="center",
                va="bottom" if dy > 0 else "top", zorder=10, path_effects=halo, fontweight="semibold")


def draw_scale_and_north(ax, fontsize=11, lw=3):
    dlon = 10.0 / KM_PER_DEG_LON
    x0, y0 = LON0 + 0.025, LAT1 + 0.022
    halo = [pe.withStroke(linewidth=3, foreground="white")]
    ax.plot([x0, x0 + dlon], [y0, y0], color=C_INK, lw=lw, solid_capstyle="butt", zorder=12,
            path_effects=[pe.withStroke(linewidth=lw + 3, foreground="white")])
    ax.plot([x0, x0 + dlon / 2], [y0, y0], color="white", lw=lw - 1, solid_capstyle="butt", zorder=13)
    for xx, t in [(x0, "0"), (x0 + dlon / 2, "5"), (x0 + dlon, "10 km")]:
        ax.text(xx, y0 + 0.004, t, fontsize=fontsize, ha="center", va="bottom", color=C_INK, zorder=13, path_effects=halo)
    xn, yn = LON1 - 0.03, LAT0 - 0.06
    ax.annotate("", xy=(xn, yn + 0.035), xytext=(xn, yn),
                arrowprops=dict(arrowstyle="-|>", color=C_INK, lw=2, mutation_scale=18), zorder=13)
    ax.text(xn, yn + 0.038, "N", fontsize=fontsize + 3, ha="center", va="bottom", color=C_INK, fontweight="bold",
            zorder=13, path_effects=halo)


def style_map_axes(ax, tick_size=10, ticks=True):
    ax.set_aspect(KM_PER_DEG_LAT / KM_PER_DEG_LON)
    if ticks:
        xt = np.arange(114.3, 114.81, 0.1)
        yt = np.arange(-3.0, -2.45, 0.1)
        ax.set_xticks(xt)
        ax.set_yticks(yt)
        ax.set_xticklabels([f"{v:.1f} E" for v in xt])
        ax.set_yticklabels([f"{-v:.1f} S" for v in yt])
        ax.tick_params(labelsize=tick_size, colors=C_INK2, length=3)
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    ax.set_xlim(LON0, LON1)
    ax.set_ylim(LAT1, LAT0)
    for sp in ax.spines.values():
        sp.set_color(C_INK2)
        sp.set_linewidth(0.8)


def caption_block(fig, x, y_top, lines, fontsize, color, dy):
    for k, line in enumerate(lines):
        fig.text(x, y_top - k * dy, line, fontsize=fontsize, color=color, ha="left", va="top")


# ---------------------------------------------------------------- map 1
def make_map_degraded(m, sites, settlements, areas):
    fig = plt.figure(figsize=(15, 22.4), dpi=200)
    ax = fig.add_axes([0.055, 0.365, 0.925, 0.595])
    ax.imshow(compose_base(m), extent=EXTENT, interpolation="nearest", zorder=1)
    draw_lines(ax, m)
    draw_sites_labels(ax, sites, settlements, site_size=80, label_size=12)
    draw_scale_and_north(ax)
    style_map_axes(ax)
    ax.set_title("Suspected degraded peat, PeatGuard dam-site AOI (ex-Mega Rice Project, Central Kalimantan)",
                 fontsize=16, color=C_INK, loc="left", pad=10)
    a = areas
    handles = [
        Patch(facecolor=C_RED, alpha=RED_ALPHA, edgecolor=C_RED_EDGE,
              label=f"Suspected degraded peat (any criterion): {a['degraded_km2']:.0f} km2, {100*a['degraded_frac_of_box']:.0f}% of the box"),
        Patch(facecolor=C_PEAT, edgecolor=C_INK2, label=f"Intact peat (no criterion met): {a['intact_km2']:.0f} km2"),
        Patch(facecolor="#dddddd", edgecolor=C_INK2, label=f"Non-peat (< 50 cm), hillshade grey: {a['nonpeat_km2']:.0f} km2"),
        Patch(facecolor=C_WATER, edgecolor=C_INK2, label="Water (rivers, open water)"),
        Line2D([0], [0], color=C_RIVER, lw=2.5, label="River (canal lines omitted; see panel a)"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor=C_SITE, markeredgecolor="white", ms=10, label="Dam site (60)"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor=C_SITE_BURN, markeredgecolor="white", ms=11,
               label=f"Dam site burning this week ({a['sites_burning']})"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.055, 0.340), ncol=2, fontsize=11, frameon=False,
               handlelength=2.2, columnspacing=1.5)

    # Secondary panels: which criterion fires
    panels = [("crit_canal", f"(a) Canal-drained: within 1 km of a mapped canal\n{a['crit_canal_km2']:.0f} km2 of peat"),
              ("crit_burn", f"(b) Burn-scarred: Landsat NBR < 0.45 in >= 2 of 11 epochs\n{a['crit_burn_km2']:.0f} km2 of peat"),
              ("crit_sink", f"(c) Measured sinking: Tay 2025 rate < -5 mm/yr\n{a['crit_sink_km2']:.0f} km2 of peat ({a['tay_cover_pct']:.0f}% of box covered)")]
    h = m["hill"]
    h = 0.82 + 0.18 * (h - h.min()) / (h.max() - h.min())
    for k, (key, title) in enumerate(panels):
        px = fig.add_axes([0.055 + k * 0.315, 0.066, 0.295, 0.188])
        rgb = np.ones((NR, NC, 3), np.float32) * 0.94
        rgb[m["peat"]] = C_PEAT
        rgb = rgb * h[..., None]
        rgb[m["water"]] = C_WATER
        sel = m[key] & m["peat"]
        col = np.array(matplotlib.colors.to_rgb(C_CRIT[key]), np.float32)
        rgb[sel] = 0.25 * rgb[sel] + 0.75 * col
        if key == "crit_sink":
            cov = ~np.isnan(m["tay"]) & ~sel
            rgb[cov] = 0.6 * rgb[cov] + 0.4 * np.array([0.75, 0.75, 0.80], np.float32)
        px.imshow(np.clip(rgb, 0, 1), extent=EXTENT, interpolation="nearest")
        if key == "crit_canal":
            px.imshow(rgba_layer(thin_lines(m["canals"]), (0.25, 0.35, 0.6), 0.9), extent=EXTENT, interpolation="nearest")
        px.imshow(rgba_layer(thin_lines(m["rivers"], 2), matplotlib.colors.to_rgb(C_RIVER), 0.8),
                  extent=EXTENT, interpolation="nearest")
        style_map_axes(px, ticks=False)
        px.set_title(title, fontsize=10.5, color=C_INK, loc="left", pad=6)
        if key == "crit_sink":
            px.text(LON0 + 0.01, LAT1 + 0.01, "grey tint = Tay 2025 coverage, not sinking", fontsize=8.5,
                    color=C_INK, ha="left", va="bottom", path_effects=[pe.withStroke(linewidth=2.5, foreground="white")])
    fig.text(0.055, 0.284, "Which criterion fires (secondary panels; criteria overlap, the main map shows their union):",
             fontsize=11, color=C_INK, ha="left", va="top")
    lines = wrap_lines(CAPTION_DEF, 175) + wrap_lines(CAPTION_SRC, 175)
    caption_block(fig, 0.055, 0.050, lines, 9.5, C_INK, 0.0085)
    path = os.path.join(OUT, "map_suspected_degraded_peat.png")
    fig.savefig(path, dpi=200, facecolor="white")
    plt.close(fig)
    return path


# ---------------------------------------------------------------- map 2
def fire_cmap():
    from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
    cmap = LinearSegmentedColormap.from_list("fires", ["#fde3c0", "#fdae61", "#f46d43", "#d7301f", "#a50f15", "#4a0008"])
    bounds = [0.5, 1.5, 2.5, 3.5, 4.5, 6.5, 8.5, 12.5, 24.5]
    return cmap, BoundaryNorm(bounds, cmap.N)


def make_map_recurrence(m, sites, settlements, fy500, counts, stats500, stats1k):
    fig = plt.figure(figsize=(15, 20.6), dpi=200)
    ax = fig.add_axes([0.055, 0.305, 0.925, 0.645])
    base = compose_base(m)
    grey = base.mean(axis=2, keepdims=True) * np.ones((1, 1, 3), np.float32)
    ax.imshow(np.clip(0.5 * base + 0.5 * grey, 0, 1), extent=EXTENT, interpolation="nearest", zorder=1)
    cmap, norm = fire_cmap()
    im = ax.imshow(np.ma.masked_less(fy500, 1).astype(float), extent=F_EXTENT, cmap=cmap, norm=norm,
                   interpolation="nearest", alpha=0.92, zorder=3)
    xs = LON0 + (np.arange(NC) + 0.5) * STEP
    ys = LAT0 - (np.arange(NR) + 0.5) * STEP
    ax.contour(xs, ys, ndimage.gaussian_filter(m["degraded"].astype(np.float32), 2.0), levels=[0.5],
               colors=[C_RED_EDGE], linewidths=1.1, zorder=5)
    draw_lines(ax, m)
    draw_sites_labels(ax, sites, settlements, site_size=80, label_size=12)
    draw_scale_and_north(ax)
    style_map_axes(ax)
    ax.set_title("Fire return count 2001 to 2024, PeatGuard dam-site AOI: years with fire per 500 m cell "
                 "(MODIS 1 km to 2011, VIIRS 375 m from 2012)", fontsize=14.5, color=C_INK, loc="left", pad=10)

    # Right column: legend, colour bar, contrast numbers
    xr = 0.685
    handles = [
        Line2D([0], [0], color=C_RED_EDGE, lw=1.5, label="Suspected degraded peat (outline; fill muted)"),
        Patch(facecolor=C_PEAT, edgecolor=C_INK2, label="Intact peat"),
        Line2D([0], [0], color=C_RIVER, lw=2.5, label="River"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor=C_SITE, markeredgecolor="white", ms=10, label="Dam site"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor=C_SITE_BURN, markeredgecolor="white", ms=11,
               label="Dam site burning this week"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(xr, 0.272), ncol=1, fontsize=10.5, frameon=False)
    cax = fig.add_axes([xr + 0.015, 0.158, 0.27, 0.013])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=[1, 2, 3, 4, 5.5, 7.5, 10.5, 18.5])
    cb.ax.set_xticklabels(["1", "2", "3", "4", "5-6", "7-8", "9-12", "13-24"])
    cb.ax.tick_params(labelsize=10, colors=C_INK)
    cb.set_label("Fire years per 500 m cell, 2001 to 2024", fontsize=10.5, color=C_INK)
    cb.outline.set_linewidth(0.6)
    d, i, n = stats500["suspected_degraded_peat"], stats500["intact_peat"], stats500["non_peat"]
    d1, i1 = stats1k["suspected_degraded_peat"], stats1k["intact_peat"]
    lines = ["Contrast, 500 m cells by majority class (headline):",
             f"degraded n={d['cells']:,}: >=1 yr {100*d['frac_ge1_fire_years']:.0f}%, >=3 yr {100*d['frac_ge3_fire_years']:.0f}%, "
             f">=5 yr {100*d['frac_ge5_fire_years']:.0f}%, mean {d['mean_fire_years_per_cell']:.2f}",
             f"intact n={i['cells']:,}: >=1 yr {100*i['frac_ge1_fire_years']:.0f}%, >=3 yr {100*i['frac_ge3_fire_years']:.0f}%, "
             f">=5 yr {100*i['frac_ge5_fire_years']:.0f}%, mean {i['mean_fire_years_per_cell']:.2f}",
             f"non-peat: {n['cells']} cells (PeatGrids >= 50 cm over the whole box)",
             "MODIS-only 1 km reference (sensor-consistent, 2001 to 2024):",
             f"degraded mean {d1['mean_fire_years_per_cell']:.2f} yr, intact mean {i1['mean_fire_years_per_cell']:.2f} yr; "
             f">=3 yr {100*d1['frac_ge3_fire_years']:.0f}% vs {100*i1['frac_ge3_fire_years']:.0f}%"]
    caption_block(fig, xr, 0.118, lines, 9.2, C_INK, 0.0092)

    # Detections per year (MODIS C6.1, sensor-consistent series), El Nino years marked
    bx = fig.add_axes([0.075, 0.115, 0.565, 0.145])
    yrs = np.array(FIRE_YEARS)
    vals = np.array([counts[str(y)] for y in yrs])
    en = np.isin(yrs, EL_NINO)
    bx.bar(yrs[~en], vals[~en], width=0.72, color="#f46d43", edgecolor="none", zorder=3)
    bx.bar(yrs[en], vals[en], width=0.72, color="#a50f15", edgecolor=C_INK, linewidth=0.8, hatch="//", zorder=3)
    for y, v in zip(yrs[en], vals[en]):
        bx.text(y, v + 60, "El Nino", fontsize=8.5, ha="center", va="bottom", color=C_INK, rotation=90)
    for y, v in zip(yrs, vals):
        if v >= 1800:
            bx.text(y - 0.42, v, f"{v:,}", fontsize=8, ha="right", va="top", color=C_INK)
    bx.set_xticks(yrs)
    bx.set_xticklabels([str(y) for y in yrs], rotation=90, fontsize=9, color=C_INK)
    bx.set_ylim(0, vals.max() * 1.32)
    bx.set_xlim(2000.3, 2024.7)
    bx.tick_params(axis="y", labelsize=9, colors=C_INK2, length=0)
    bx.yaxis.grid(True, color="#dddddd", lw=0.6, zorder=0)
    bx.set_axisbelow(True)
    for sp in ["top", "right", "left"]:
        bx.spines[sp].set_visible(False)
    bx.spines["bottom"].set_color(C_INK2)
    bx.set_ylabel("MODIS detections in AOI", fontsize=10, color=C_INK)
    bx.set_title("NASA FIRMS MODIS C6.1 detections per year inside the AOI, sensor-consistent series (El Nino years hatched, dark)",
                 fontsize=10.5, color=C_INK, loc="left")

    txt = ("Fire return count [E for detections, S for the gridding] = number of years 2001 to 2024 with fire in a 0.0045 deg "
           "(about 500 m) cell. 2001 to 2011: NASA FIRMS MODIS C6.1, 1 km footprint, each detection marks the 500 m cells within "
           "500 m of its centre. 2012 to 2024: VIIRS 375 m (S-NPP, plus NOAA-20 from 2018 deduplicated within 375 m and 1 h), "
           "vegetation-fire type, nominal or high confidence, one cell per detection. VIIRS sees more small fires than MODIS, so "
           "counts are not sensor-consistent across 2011/2012. Detections are thermal anomalies and include agricultural burning.")
    lines = wrap_lines(txt, 200) + wrap_lines(CAPTION_DEF + "  " + CAPTION_SRC, 200)
    caption_block(fig, 0.055, 0.064, lines, 9.0, C_INK, 0.0082)
    path = os.path.join(OUT, "map_fire_recurrence_2001_2024.png")
    fig.savefig(path, dpi=200, facecolor="white")
    plt.close(fig)
    return path


# ---------------------------------------------------------------- animation
def running_counter(cum, cls):
    deg = cls == 0
    n_deg = int(deg.sum())
    return dict(ge_n=int((cum >= RECUR_N).sum()), ge_1=int((cum >= 1).sum()),
                ge_n_degraded=int(((cum >= RECUR_N) & deg).sum()), ge_1_degraded=int(((cum >= 1) & deg).sum()),
                n_degraded=n_deg,
                pct_ge_n_degraded=float(((cum >= RECUR_N) & deg).sum() / max(1, n_deg)),
                pct_ge_1_degraded=float(((cum >= 1) & deg).sum() / max(1, n_deg)))


PANEL_X0 = 0.645
PANEL_W = 0.335


def put_wrapped(fig, x, y, text, fontsize, color, width_frac=PANEL_W, line_h=1.35, **kw):
    """Draw text wrapped to the panel width; returns the y below the block."""
    px_per_char = fontsize * (100.0 / 72.0) * 0.58
    width_chars = max(10, int(width_frac * 1920 / px_per_char))
    lines = textwrap.wrap(text, width=width_chars) or [""]
    dy = fontsize * (100.0 / 72.0) * line_h / 1080.0
    for k, line in enumerate(lines):
        fig.text(x, y - k * dy, line, fontsize=fontsize, color=color, ha="left", va="top", **kw)
    return y - len(lines) * dy


def render_frame(path, label_year, fire_img, char_img, base, m, sites, settlements, annual, running, is_el_nino,
                 sensor_line, count_line, source_lines, highlight_year, cd=None):
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    fig.patch.set_facecolor("#111111")
    ax = fig.add_axes([0.02, 0.05, 0.60, 0.90])
    ax.imshow(base, extent=EXTENT, interpolation="nearest", zorder=1)
    ax.imshow(char_img, extent=EXTENT, interpolation="nearest", zorder=2)
    draw_lines(ax, m)
    ax.imshow(fire_img, extent=EXTENT, interpolation="nearest", zorder=6)
    draw_sites_labels(ax, sites, settlements, site_size=38, label_size=10, site_lw=0.6)
    draw_scale_and_north(ax, fontsize=10, lw=3)
    style_map_axes(ax, tick_size=9)
    ax.tick_params(colors="#cccccc")
    for sp in ax.spines.values():
        sp.set_color("#888888")

    x0 = PANEL_X0
    y = put_wrapped(fig, x0, 0.95, "Fires returning to suspected degraded peat", 17.5, "white", fontweight="semibold")
    y = put_wrapped(fig, x0, y - 0.004, "PeatGuard dam-site AOI, ex-Mega Rice Project, Central Kalimantan", 12.5, "#bbbbbb")
    fig.text(x0, y - 0.005, str(label_year), fontsize=84, color="white", ha="left", va="top", fontweight="bold")
    if is_el_nino:
        fig.text(x0 + 0.205, y - 0.065, "El Nino", fontsize=19, color="#111111", ha="left", va="center", fontweight="bold",
                 bbox=dict(boxstyle="round,pad=0.35", facecolor="#ffb347", edgecolor="none"))
    y = y - 0.135
    y = put_wrapped(fig, x0, y, f"{100 * running['pct_ge_n_degraded']:.0f} percent of the degraded peat has burned in "
                    f"{RECUR_N} or more years so far", 17, "white", fontweight="semibold")
    y = put_wrapped(fig, x0, y - 0.003, f"burned at least once: {100 * running['pct_ge_1_degraded']:.0f} percent "
                    f"(500 m cells, {running['ge_n_degraded']:,} of {running['n_degraded']:,} in {RECUR_N} or more years)",
                    12, "#dddddd")
    y = put_wrapped(fig, x0, y - 0.006, sensor_line, 10.5, "#bbbbbb")
    y = put_wrapped(fig, x0, y - 0.002, count_line, 10.5, "#ffb000")

    bx = fig.add_axes([x0 + 0.008, 0.355, 0.175, 0.135])
    bx.set_facecolor("#111111")
    yrs = np.array(FIRE_YEARS)
    vals = np.array([annual[y] for y in yrs])
    hy = 9999 if highlight_year is None else highlight_year
    cols = ["#ffb000" if y == hy else ("#b0452a" if y < hy else "#5a3a2a") for y in yrs]
    bx.bar(yrs, vals, width=0.72, color=cols, edgecolor="none", zorder=3)
    ymax = vals.max() * 1.25
    for y_ in EL_NINO:
        bx.text(y_, annual[y_] + ymax * 0.02, "EN", fontsize=7.5, ha="center", va="bottom", color="#ffb347")
    bx.axvline(2011.5, color="#888888", lw=0.8, ls="--", zorder=2)
    bx.text(2011.3, ymax * 0.95, "MODIS 1 km", fontsize=8, color="#aaaaaa", ha="right", va="top")
    bx.text(2011.7, ymax * 0.95, "VIIRS 375 m", fontsize=8, color="#aaaaaa", ha="left", va="top")
    bx.set_xticks([2001, 2005, 2010, 2015, 2020, 2024])
    bx.set_xticklabels(["2001", "2005", "2010", "2015", "2020", "2024"], fontsize=9, color="#cccccc")
    bx.tick_params(axis="y", labelsize=8, colors="#aaaaaa", length=0)
    bx.set_xlim(2000.3, 2024.7)
    bx.set_ylim(0, ymax)
    bx.yaxis.grid(True, color="#333333", lw=0.6, zorder=0)
    bx.set_axisbelow(True)
    for sp in ["top", "right", "left"]:
        bx.spines[sp].set_visible(False)
    bx.spines["bottom"].set_color("#666666")
    bx.set_title("500 m cells with fire per year\n(EN = El Nino; sensor change at 2012)", fontsize=8.5, color="#cccccc", loc="left")
    if cd is not None:
        cx = fig.add_axes([x0 + 0.215, 0.355, 0.115, 0.135])
        cx.set_facecolor("#111111")
        draw_canal_bars(cx, cd, dark=True, fontsize=8, annotate=False)
        cx.set_title("mean fire years per cell\nvs distance to canal (km)", fontsize=8.5, color="#cccccc", loc="left")

    handles = [
        Patch(facecolor=FIRE_RGB[2], edgecolor=FIRE_RGB[0], label="Burning this frame (core to edge)"),
        Patch(facecolor=CHAR_RGB, edgecolor="#999999", label="Burned in 2 or more earlier years (charcoal)"),
        Patch(facecolor=C_RED, alpha=RED_ALPHA + 0.2, edgecolor=C_RED_EDGE, label="Suspected degraded peat"),
        Patch(facecolor=C_PEAT, edgecolor="none", label="Intact peat"),
        Line2D([0], [0], color=C_RIVER, lw=2, label="River"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor=C_SITE, markeredgecolor="white", ms=8, label="Dam site"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor=C_SITE_BURN, markeredgecolor="white", ms=9,
               label="Dam site burning this week"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(x0, 0.262), ncol=2, fontsize=9, frameon=False,
               labelcolor="#dddddd", columnspacing=1.0, handlelength=1.6, handletextpad=0.5)
    y = 0.160
    for line in source_lines:
        y = put_wrapped(fig, x0, y, line, 8.5, "#999999", line_h=1.25) - 0.004
    fig.savefig(path, dpi=100, facecolor=fig.get_facecolor())
    plt.close(fig)
    return path


def encode_video(year_frames, mp4_name, size, gif_name=None, gif_size=None):
    """year_frames: list of lists of sub-frame paths (one list per year). Flicker cycles sub-frames at FPS."""
    seq = []
    for k, subs in enumerate(year_frames):
        n = LAST_HOLD_FRAMES if k == len(year_frames) - 1 else HOLD_FRAMES
        seq.extend(subs[i % len(subs)] for i in range(n))
    lst = os.path.join(FRAMES, f"concat_{mp4_name}.txt")
    with open(lst, "w") as fh:
        for p in seq:
            fh.write(f"file '{p}'\nduration {1.0 / FPS:.6f}\n")
        fh.write(f"file '{seq[-1]}'\n")
    mp4 = os.path.join(OUT, mp4_name)
    subprocess.run([FFMPEG, "-y", "-loglevel", "error", "-f", "concat", "-safe", "0", "-i", lst,
                    "-vf", f"fps={FPS},scale={size[0]}:{size[1]}:flags=lanczos,format=yuv420p", "-c:v", "libx264", "-crf", "18",
                    "-profile:v", "high", "-pix_fmt", "yuv420p", "-movflags", "+faststart", mp4], check=True)
    if gif_name is None:
        return mp4, None, len(seq), 0
    gif = os.path.join(OUT, gif_name)
    stride = FPS // GIF_FPS
    gseq = seq[::stride]
    cache = {}
    frames = []
    for p in gseq:
        if p not in cache:
            im = Image.open(p).convert("RGB").resize(gif_size, Image.LANCZOS)
            cache[p] = im.quantize(colors=128, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.FLOYDSTEINBERG)
        frames.append(cache[p])
    frames[0].save(gif, save_all=True, append_images=frames[1:], duration=int(1000 / GIF_FPS), loop=0, optimize=False, disposal=1)
    return mp4, gif, len(seq), len(gseq)


SQ_CROP = 0.015     # deg trimmed from each lon edge so the 0.63 x 0.60 deg box fills a square frame


def render_square_frame(path, label_year, fire_img, char_img, base, m, running, is_el_nino, tag=None):
    """Compact map-only frame for a slide: 1080 x 1080, map edge to edge, large year, headline, scale bar, legend."""
    fig = plt.figure(figsize=(10.8, 10.8), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(base, extent=EXTENT, interpolation="nearest", zorder=1)
    ax.imshow(char_img, extent=EXTENT, interpolation="nearest", zorder=2)
    draw_lines(ax, m)
    ax.imshow(fire_img, extent=EXTENT, interpolation="nearest", zorder=6)
    ax.set_xlim(LON0 + SQ_CROP, LON1 - SQ_CROP)
    ax.set_ylim(LAT1, LAT0)
    ax.set_aspect(KM_PER_DEG_LAT / KM_PER_DEG_LON)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    halo = [pe.withStroke(linewidth=4, foreground="#111111")]
    # Top-left block: year (about 90 px), El Nino tag, headline (about 34 px)
    ax.text(0.03, 0.975, str(label_year), transform=ax.transAxes, fontsize=65, color="white", fontweight="bold",
            ha="left", va="top", zorder=20, path_effects=[pe.withStroke(linewidth=6, foreground="#111111")])
    if is_el_nino:
        ax.text(0.30, 0.935, "El Nino", transform=ax.transAxes, fontsize=24, color="#111111", fontweight="bold",
                ha="left", va="center", zorder=20, bbox=dict(boxstyle="round,pad=0.3", facecolor="#ffb347", edgecolor="none"))
    if tag:
        ax.text(0.30, 0.935, tag, transform=ax.transAxes, fontsize=20, color="white", fontweight="bold",
                ha="left", va="center", zorder=20, bbox=dict(boxstyle="round,pad=0.3", facecolor="#333333", edgecolor="none"))
    ax.text(0.03, 0.875, f"{100 * running['pct_ge_n_degraded']:.0f} percent of degraded peat burned",
            transform=ax.transAxes, fontsize=24.5, color="white", fontweight="semibold", ha="left", va="top", zorder=20,
            path_effects=halo)
    ax.text(0.03, 0.833, f"in {RECUR_N} or more years so far", transform=ax.transAxes, fontsize=24.5, color="white",
            fontweight="semibold", ha="left", va="top", zorder=20, path_effects=halo)
    # Scale bar: 0, 5, 10 km
    dlon = 10.0 / KM_PER_DEG_LON
    x0, y0 = LON1 - SQ_CROP - 0.03 - dlon, LAT1 + 0.03
    ax.plot([x0, x0 + dlon], [y0, y0], color="white", lw=7, solid_capstyle="butt", zorder=20,
            path_effects=[pe.withStroke(linewidth=10, foreground="#111111")])
    ax.plot([x0, x0 + dlon / 2], [y0, y0], color="#111111", lw=5, solid_capstyle="butt", zorder=21)
    for xx, t in [(x0, "0"), (x0 + dlon / 2, "5"), (x0 + dlon, "10 km")]:
        ax.text(xx, y0 + 0.006, t, fontsize=18, ha="center", va="bottom", color="white", zorder=22, path_effects=halo,
                fontweight="semibold")
    # Compact legend, bottom-left, four swatches at about 26 px
    handles = [
        Patch(facecolor=C_RED, alpha=RED_ALPHA + 0.25, edgecolor=C_RED_EDGE, label="Degraded peat"),
        Patch(facecolor=CHAR_RGB, alpha=0.85, edgecolor="#bbbbbb", label="Burned before (2+ years)"),
        Patch(facecolor=FIRE_RGB[2], edgecolor=FIRE_RGB[0], label="Burning now"),
        Line2D([0], [0], color=C_RIVER, lw=3, label="River"),
    ]
    leg = ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.015, 0.015), fontsize=18.7, frameon=True,
                    facecolor="#111111", edgecolor="none", framealpha=0.78, labelcolor="white", handlelength=1.6,
                    handletextpad=0.6, borderpad=0.6, labelspacing=0.45)
    leg.set_zorder(25)
    fig.savefig(path, dpi=100, facecolor="#111111")
    plt.close(fig)
    return path


# ---------------------------------------------------------------- main
def main():
    skip = "--skip-download" in sys.argv
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(FRAMES, exist_ok=True)
    print("firms downloads")
    dl_report = download_firms(skip)
    m = build_masks()
    sites = jload("sites")
    settlements = jload("settlements")
    counts = jload("fires_counts")

    km2 = lambda mask: float(mask.sum() * CELL_KM2)
    areas = dict(
        box_km2=km2(np.ones((NR, NC), bool)), water_km2=km2(m["water"]), peat_km2=km2(m["peat"]),
        nonpeat_km2=km2(m["nonpeat"]), degraded_km2=km2(m["degraded"]), intact_km2=km2(m["intact"]),
        degraded_frac_of_peat=float(m["degraded"].sum() / max(1, m["peat"].sum())),
        degraded_frac_of_box=float(m["degraded"].mean()),
        crit_canal_km2=km2(m["crit_canal"] & m["peat"]), crit_burn_km2=km2(m["crit_burn"] & m["peat"]),
        crit_sink_km2=km2(m["crit_sink"] & m["peat"]),
        crit_canal_only_km2=km2(m["peat"] & m["crit_canal"] & ~m["crit_burn"] & ~m["crit_sink"]),
        crit_burn_only_km2=km2(m["peat"] & ~m["crit_canal"] & m["crit_burn"] & ~m["crit_sink"]),
        crit_sink_only_km2=km2(m["peat"] & ~m["crit_canal"] & ~m["crit_burn"] & m["crit_sink"]),
        crit_all_three_km2=km2(m["peat"] & m["crit_canal"] & m["crit_burn"] & m["crit_sink"]),
        tay_cover_pct=float(100 * (~np.isnan(m["tay"])).mean()),
        sites_burning=int(sum(s["burning"] for s in sites)),
        degraded_incl_nbr2026_km2=km2(m["peat"] & (m["crit_canal"] | (m["n_low_incl26"] >= NBR_MIN_EPOCHS) | m["crit_sink"])),
    )
    thick = m["thick"]
    peat_info = dict(crs=m["peat_crs"], nodata=(None if m["peat_nodata"] is None or np.isnan(m["peat_nodata"]) else m["peat_nodata"]),
                     thickness_cm_min=float(np.nanmin(thick)), thickness_cm_median=float(np.nanmedian(thick)),
                     thickness_cm_max=float(np.nanmax(thick)), nan_fraction=float(np.isnan(thick).mean()),
                     frac_ge_50cm=float((thick >= 50).mean()), frac_ge_100cm=float((thick >= 100).mean()),
                     frac_ge_150cm=float((thick >= 150).mean()), frac_ge_200cm=float((thick >= 200).mean()))

    # Fire data and recurrence grids
    print("fires")
    idx500 = work_to_fire_index(FSTEP, FR, FC)
    idx1k = work_to_fire_index(KSTEP, KR, KC)
    cls500 = classify(m, idx500, FR, FC)
    cls1k = classify(m, idx1k, KR, KC)
    years = {}
    fy500 = np.zeros((FR, FC), np.int32)
    fy1k_modis = np.zeros((KR, KC), np.int32)
    burned500 = {}
    for y in FIRE_YEARS:
        fd = load_fire_year(y)
        b500 = burned_grid(fd["pts"], fd["sensor"], FSTEP, FR, FC)
        b1k = burned_grid(modis_points(y), "MODIS", KSTEP, KR, KC)
        fy500 += b500
        fy1k_modis += b1k
        burned500[y] = b500
        fd["cells500"] = int(b500.sum())
        fd["cells1k_modis"] = int(b1k.sum())
        years[y] = fd
    stats500 = class_stats(fy500, cls500, F_CELL_KM2)
    stats1k = class_stats(fy1k_modis, cls1k, (KSTEP * KM_PER_DEG_LAT) * (KSTEP * KM_PER_DEG_LON))
    fy500_modis_only = sum(burned_grid(modis_points(y), "MODIS", FSTEP, FR, FC) for y in FIRE_YEARS)
    stats500_modis_only = class_stats(fy500_modis_only, cls500, F_CELL_KM2)
    sensitivity = {}
    for cut in [50, 100, 150, 200]:
        mm = dict(m)
        peat = (thick >= cut) & ~m["water"]
        anyc = m["crit_canal"] | m["crit_burn"] | m["crit_sink"]
        mm.update(peat=peat, degraded=peat & anyc, intact=peat & ~anyc, nonpeat=~peat & ~m["water"])
        sensitivity[f"peat_ge_{cut}cm"] = class_stats(fy500, classify(mm, idx500, FR, FC), F_CELL_KM2)
    sr, sc = lonlat_to_rc(np.array([s["lon"] for s in sites]), np.array([s["lat"] for s in sites]), FSTEP)
    site_fire_years = {s["site"]: int(fy500[r, c]) for s, r, c in zip(sites, sr, sc) if 0 <= r < FR and 0 <= c < FC}
    site_class = {s["site"]: ["degraded", "intact", "nonpeat", "water"][cls500[r, c] if cls500[r, c] >= 0 else 3]
                  for s, r, c in zip(sites, sr, sc) if 0 <= r < FR and 0 <= c < FC}

    cd = canal_distance_stats(m, fy500, idx500)
    print("canal distance figure")
    p3 = make_fig_canal_distance(cd)
    print("map 1")
    p1 = make_map_degraded(m, sites, settlements, areas)
    print("map 2")
    p2 = make_map_recurrence(m, sites, settlements, fy500, counts, stats500, stats1k)

    print("frames")
    base = compose_base(m)
    annual = {y: years[y]["cells500"] for y in FIRE_YEARS}
    cum = np.zeros((FR, FC), np.int32)
    year_frames, square_frames, frame_meta = [], [], []
    src_common = ["Peat: PeatGrids (Widyastuti et al. 2022). Canals: Dadap et al. 2021. NBR: Landsat (Key and Benson 2006).",
                  "Subsidence: Tay et al. 2025 (ALOS-2). Suspected degraded peat is a screening layer [S], not ground-truthed.",
                  "Patches: detections rasterised and Gaussian-smoothed (sigma 0.9 km MODIS, 0.5 km VIIRS, seeded per-year "
                  "jitter), drawn where the field exceeds one detection equivalent and capped at 3 x the footprint area; "
                  "glow for visualisation, not burned-area polygons."]
    for y in FIRE_YEARS:
        fd = years[y]
        field, factor, foot_km2 = burn_field(fd["pts"], fd["sensor"], seed=y)
        levels = fire_levels(field, foot_km2)
        cum_now = cum + burned500[y]
        running = running_counter(cum_now, cls500)
        char_img = char_layer(cum, idx500)
        if fd["sensor"] == "MODIS":
            sensor_line = f"NASA FIRMS MODIS C6.1, 1 km footprint, Jan to Dec {y}  (2001 to 2011: MODIS 1 km)"
            src0 = "Fires: NASA FIRMS MODIS C6.1, 1 km footprint; thermal anomalies, includes agricultural burning."
        else:
            sats = "S-NPP + NOAA-20" if y >= 2018 else "S-NPP"
            sensor_line = f"NASA FIRMS VIIRS 375 m ({sats}), Jan to Dec {y}  (2012 onward: VIIRS 375 m)"
            src0 = "Fires: NASA FIRMS VIIRS 375 m, vegetation-fire type, nominal or high confidence; NOAA-20 deduplicated vs S-NPP."
        count_line = f"{len(fd['pts']):,} detections this year, {fd['cells500']:,} cells of 500 m burned"
        subs, sq = [], []
        for k in range(SUBFRAMES):
            fimg = fire_rgba(field, levels, y, k)
            p = os.path.join(FRAMES, f"frame_{y}_{k}.png")
            render_frame(p, y, fimg, char_img, base, m, sites, settlements, annual, running,
                         y in EL_NINO, sensor_line, count_line, [src0] + src_common, y, cd=cd)
            subs.append(p)
            q = os.path.join(FRAMES, f"square_{y}_{k}.png")
            render_square_frame(q, y, fimg, char_img, base, m, running, y in EL_NINO)
            sq.append(q)
        year_frames.append(subs)
        square_frames.append(sq)
        frame_meta.append(dict(year=y, sensor=fd["sensor"], detections=int(len(fd["pts"])), cells500=fd["cells500"],
                               sigma_factor=round(factor, 3), footprint_km2=round(foot_km2, 1),
                               glow_km2=round(float((field >= levels[0]).sum() * CELL_KM2), 1),
                               levels=[round(float(v), 3) for v in levels], running=running))
        cum = cum_now
    wpts, winfo, wcells = week_points()
    field, factor, foot_km2 = burn_field(wpts, "VIIRS", seed=2026)
    levels = fire_levels(field, foot_km2)
    running = running_counter(cum, cls500)
    days = sorted(winfo["byDay"].keys())
    d0, d1 = days[0][8:].lstrip("0"), days[-1][8:].lstrip("0")
    sensor_line = f"This week, VIIRS + MODIS, {d0} to {d1} Oct 2026 (0.05 deg bins, positions approximate)"
    count_line = f"{int(wcells[:, 2].sum()):,} detections this week inside the AOI"
    src0 = "This frame: NASA FIRMS VIIRS S-NPP + NOAA-20 + MODIS, last 7 days, 0.05 deg bins; not comparable with yearly frames."
    subs, sq = [], []
    char_img = char_layer(cum, idx500)
    for k in range(SUBFRAMES):
        fimg = fire_rgba(field, levels, 2026, k)
        p = os.path.join(FRAMES, f"frame_2026_{k}.png")
        render_frame(p, 2026, fimg, char_img, base, m, sites, settlements, annual,
                     running, False, sensor_line, count_line, [src0] + src_common, None, cd=cd)
        subs.append(p)
        q = os.path.join(FRAMES, f"square_2026_{k}.png")
        render_square_frame(q, 2026, fimg, char_img, base, m, running, False, tag=f"this week, {d0} to {d1} Oct")
        sq.append(q)
    year_frames.append(subs)
    square_frames.append(sq)
    frame_meta.append(dict(year=2026, sensor="VIIRS+MODIS week", detections=int(wcells[:, 2].sum()),
                           cells_005deg=int(len(wcells)), sigma_factor=round(factor, 3), footprint_km2=round(foot_km2, 1),
                           glow_km2=round(float((field >= levels[0]).sum() * CELL_KM2), 1), running=running,
                           window=[days[0], days[-1]], fetched_at=winfo.get("fetchedAt")))
    print("encode")
    mp4, _, n_mp4, _ = encode_video(year_frames, "fires_return_to_degraded_peat_2001_2026.mp4", (1920, 1080))
    mp4_sq, gif, n_sq, n_gif = encode_video(square_frames, "fires_return_to_degraded_peat_2001_2026_SLIDE.mp4", (1080, 1080),
                                            "fires_return_to_degraded_peat_2001_2026.gif", (720, 720))

    stats = dict(
        generated_from="deliverables/degraded_peat_fire_recurrence/build.py",
        work_grid=dict(shape=[NR, NC], lon0=LON0, lat0=LAT0, step_deg=STEP, cell_km2=CELL_KM2, extent=EXTENT),
        fire_grid_500m=dict(shape=[FR, FC], step_deg=FSTEP, cell_km2=F_CELL_KM2, extent=F_EXTENT,
                            class_rule="majority of degraded / intact / non-peat work cells by centre; no land = excluded",
                            modis_rule="each MODIS detection marks the 500 m cells whose centres lie within 500 m",
                            viirs_rule="each VIIRS detection marks the 500 m cell containing it"),
        fire_grid_1km_modis=dict(shape=[KR, KC], step_deg=KSTEP, extent=K_EXTENT, note="MODIS-only, all years 2001 to 2024, 1 km"),
        definitions=dict(
            peat=f"PeatGrids mean thickness >= {PEAT_MIN_CM:.0f} cm, minus water_mask and rivers [E dataset, S cutoff]",
            crit_canal=f"binary dilation of canals.npy by a disk of radius {CANAL_DILATE_CELLS} cells (~1 km) [S]",
            crit_burn=f"NBR < {NBR_LOW} in >= {NBR_MIN_EPOCHS} valid epochs of {NBR_EPOCHS}; NaN ignored; {NBR_EXCLUDED} excluded [S]",
            crit_sink=f"Tay et al. 2025 ALOS-2 long-term rate < {TAY_SINK_MMYR} mm/yr where covered [S]",
            suspected_degraded="peat AND (crit_canal OR crit_burn OR crit_sink) [S, screening, not ground-truthed]",
            fire_return_count="years 2001 to 2024 with fire per 500 m cell; MODIS 1 km to 2011, VIIRS 375 m from 2012 [E detections, S gridding]",
            viirs_filter=f"type == 0 (vegetation fire), confidence != low, NOAA-20 dropped if an S-NPP detection within {DEDUP_KM} km and {DEDUP_HOURS} h"),
        peatgrids=peat_info,
        areas_km2=areas,
        fractions_of_box=dict(peat=float(m["peat"].mean()), degraded=float(m["degraded"].mean()), intact=float(m["intact"].mean()),
                              nonpeat=float(m["nonpeat"].mean()), water=float(m["water"].mean()),
                              crit_canal=float((m["crit_canal"] & m["peat"]).mean()), crit_burn=float((m["crit_burn"] & m["peat"]).mean()),
                              crit_sink=float((m["crit_sink"] & m["peat"]).mean())),
        nbr=dict(epochs_used=NBR_EPOCHS, excluded=NBR_EXCLUDED, threshold=NBR_LOW, min_epochs=NBR_MIN_EPOCHS,
                 valid_epochs_per_cell_min=int(m["n_valid"].min()), valid_epochs_per_cell_mean=float(m["n_valid"].mean()),
                 frac_cells_low_ge2=float((m["n_low"] >= 2).mean())),
        class_counts_500m={"degraded": int((cls500 == 0).sum()), "intact": int((cls500 == 1).sum()),
                           "nonpeat": int((cls500 == 2).sum()), "no_land": int((cls500 == -1).sum())},
        recurrence_500m_mixed_sensor=stats500,
        recurrence_500m_modis_only=stats500_modis_only,
        recurrence_1km_modis_only=stats1k,
        recurrence_sensitivity_to_peat_cutoff_500m=sensitivity,
        fire_years_vs_canal_distance=cd,
        fire_years_histogram_500m={str(k): int(v) for k, v in zip(*np.unique(fy500[cls500 >= 0], return_counts=True))},
        per_year={str(y): dict(sensor=years[y]["sensor"], detections=int(len(years[y]["pts"])), cells500=years[y]["cells500"],
                               cells1k_modis=years[y]["cells1k_modis"], modis_detections_reference=counts[str(y)],
                               **{k: v for k, v in years[y]["meta"].items() if k not in ("points", "satellite", "year", "sensor", "n")})
                  for y in FIRE_YEARS},
        firms_downloads=dl_report,
        el_nino_years=EL_NINO,
        running_counter_threshold=RECUR_N,
        site_fire_years_500m=site_fire_years,
        site_class_500m=site_class,
        animation=dict(frames=frame_meta, fps=FPS, subframes_per_year=SUBFRAMES, hold_frames=HOLD_FRAMES,
                       last_hold_frames=LAST_HOLD_FRAMES, mp4_frames=n_mp4, gif_fps=GIF_FPS, gif_frames=n_gif,
                       size_mp4=[1920, 1080], size_mp4_slide=[1080, 1080], size_gif=[720, 720],
                       slide_lon_crop_deg=SQ_CROP,
                       fire_render=f"footprint raster (MODIS 1 km, VIIRS 375 m), 300 m seeded jitter, Gaussian sigma {SIGMA_KM} km x "
                                   f"U(0.8,1.3) per year, normalised to one isolated detection's peak, smoothed multiplicative noise, "
                                   f"bands at {FIRE_LEVELS} detection equivalents, outer level raised if glow area > {FIRE_AREA_CAP} x "
                                   "footprint area; flicker: alpha x clip(1+0.3 n, 0.7, 1.3), core x (1+0.12 n), noise sigma 4 cells, "
                                   "seeded per year and sub-frame"),
        outputs=[os.path.relpath(p, ROOT) for p in [p1, p2, p3, mp4, mp4_sq, gif]],
    )
    with open(os.path.join(OUT, "stats.json"), "w") as fh:
        json.dump(stats, fh, indent=1)
    print(json.dumps(dict(mixed=stats500, modis1k=stats1k), indent=1))
    print(json.dumps({str(y): (years[y]["sensor"], int(len(years[y]["pts"])), years[y]["cells500"]) for y in FIRE_YEARS}))
    return stats


if __name__ == "__main__":
    main()
