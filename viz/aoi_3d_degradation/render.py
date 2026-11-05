#!/usr/bin/env python3
"""Off-screen 3D flyover frames for the dam-site AOI degradation sequence (v4).

Scene: Copernicus DEM (NISAR-supplied, WGS84 ellipsoid heights, relative relief only)
as a textured StructuredGrid in a local km frame with a constant vertical
exaggeration (VEXAG). The texture carries, and only carries:

1. a hillshaded forest-green base (one colour, reads as intact peat forest; no land-cover classes);
2. PRE-EXISTING DEGRADED PEAT as of each year in one muted red, as observed by
   satellites (never inferred from canals): the whole box is peat (PeatGrids
   >= 50 cm) AND Landsat NBR < 0.45 in any composite strictly before the year,
   OR a FIRMS fire footprint in any earlier year, OR the radar-verified sinking
   zone from 2018 (Tay et al. 2025 ALOS-2 rate < VERIFIED_RATE_MMYR, cleaned by
   a binary opening and closing of about 5 cells, blobs under VERIFIED_MIN_CELLS
   dropped); opened with a 3-cell diamond. The red grows through the sequence
   and every flame lands on red (degraded before it burned) or on grey. The
   radar-verified subset is additionally outlined by a dashed cyan line mesh
   lifted above the terrain. A headline under the year gives the share of the
   year's fire that fell on pre-existing degraded peat;
3. black char where the ground burned in an earlier FIRMS year, darker for more
   fire years (build_char), accumulating through the sequence;
4. this year's FIRMS fire: the fine detection density is gated (no flame below
   FLAME_FLOOR_PER_KM2, full flame at FLAME_FULL_PER_KM2) and drives a separate
   translucent flame mesh, z = terrain + FLAME_H * alpha * flicker(tick), coloured
   yellow-orange-red, with a grey-blue smoke mesh above the cores drifting
   north-west; the coarse fused field is baked into the texture only as a low
   ground glow (GROUND_GLOW_STRENGTH);
5. OSM rivers in thin blue; dam sites as white spheres with a dark rim (no
   canals are drawn);
6. for the 2026 segment only, the Planet NICFI September 2026 mosaic replaces
   the green base (hillshade at half strength, red at alpha 0.40), reached by a
   12-frame dissolve from the green 2024 frame.

flicker is a seeded low-frequency noise field that crawls and oscillates 0.7 to
1.3 per frame (tick), so the flames rise, fall and move at 12 fps. From 2018 the
modelled lowering (Tay rate) is applied to z with a further SINK_EXAG (geometry
only, no tint).

Nothing in this scene is a hydrological gradient: the DEM cannot resolve the
peat-dome drainage direction (see memory note project_dem_peat_gradient_unusable).
The canal-proximity (inferred) degradation criterion is deliberately not drawn:
nothing on screen calls an inferred area verified.

Usage (v4 schedule, slide mode, all years, 7 processes):
  .venv/bin/python3 deliverables/aoi_3d_degradation/render.py --all --jobs 7 \
      --minimal --bg "#000000" --outdir deliverables/aoi_3d_degradation/frames_v4

One year of the schedule (moving frames, static hold frames with live flames, and
the 8 pre-roll frames at the previous year's last pose used for the crossfade):
  .venv/bin/python3 deliverables/aoi_3d_degradation/render.py --year 2015 --minimal \
      --bg "#000000" --outdir deliverables/aoi_3d_degradation/frames_v4

Legacy segment mode (any year, explicit camera range, HUD or slide):
  .venv/bin/python3 deliverables/aoi_3d_degradation/render.py --epoch 2015 \
      --t0 0.64 --t1 0.64 --frames 120 --minimal --bg "#000000" --outdir .../proof

Frames are written as YYYY_NNN.png; pre-roll frames as YYYY_pre_NN.png.
assemble.py turns a frames_v4 directory into the MP4.
"""

import argparse
import json
import math
import os
import subprocess
import sys
import time

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage as ndi
from scipy.ndimage import gaussian_filter, map_coordinates

import pyvista as pv

pv.OFF_SCREEN = True

HERE = os.path.dirname(os.path.abspath(__file__))
# Inputs (DEM, NBR, canal, river, fire and lowering grids as .npy/.json, the
# Planet mosaic and the site/settlement lists) and the rendered frames live
# outside this repository. DATA_ROOT is the PeatGuard data root; the
# deliverable folder below it holds the prepared work-grid layers.
DATA_ROOT = os.environ.get("DATA_ROOT", os.path.abspath(os.path.join(HERE, "..", "..")))
DATA_DIR = os.path.join(DATA_ROOT, "deliverables", "aoi_3d_degradation")

# ----------------------------------------------------------------------------
# Constants (all on screen where they matter)
# ----------------------------------------------------------------------------
VEXAG = 40            # vertical exaggeration of the terrain (printed on screen)
SINK_EXAG = 25        # extra exaggeration of the modelled lowering on top of VEXAG
SINK_FROM_EPOCH = 2018
DECIMATE = 2          # grid decimation factor (1000 x 1050 -> 500 x 525)
WINDOW = (1920, 1080)
N_SETTLEMENT_LABELS = 7
BARE = os.environ.get('PG_BARE', '') == '1'   # land only: no sites, labels, north mark, legend or text
LABEL_NAMES = ["Pulang Pisau", "Kuala Kapuas", "Marabahan", "Dadahup",
               "Bina Jaya", "Manuntung", "Lamunti Baru"]
CAM_ELEV_DEG = 35.0
CAM_AZ0_DEG = -40.0
CAM_AZ1_DEG = 40.0
CAM_DIST_KM = 158.0
CAM_VIEW_ANGLE = 22.0   # narrower angle at a longer distance: less near-corner magnification
# pure screen-space shift of the projection (VTK window centre): lifts the box so
# its near corner fits above the bottom edge at the extreme azimuths
CAM_WINDOW_CENTER = (0.0, -0.28)

FIRE_SIGMA_KM = 2.0       # coarse Gaussian that fuses clusters into the ground glow (times a seeded 0.8 to 1.3 per year)
FIRE_SIGMA_FINE_KM = 0.6  # fine Gaussian that keeps detection clusters as distinct flames
FIRE_JITTER_M = 300.0     # seeded positional jitter per detection
FIRE_ALPHA_MAX = 0.85     # flame opacity at a saturated core
FLAME_FLOOR_PER_KM2 = 0.5  # fine density below which there is no flame at all
FLAME_FULL_PER_KM2 = 2.0   # fine density at which the flame has full height and opacity
FLAME_CORE_PER_KM2 = 4.0   # fine density at which the colour reaches the yellow core
FLAME_SINGLE_ALPHA = 0.45  # a lone detection still gets a small distinct flame (quiet years)
FLAME_SINGLE_SIGMA_KM = 0.25  # footprint of that lone-detection flame
FIRE_CORE_PER_KM2 = 3.0   # coarse density that saturates the ground glow
FIRE_FLOOR_PER_KM2 = 0.3  # coarse density below which the ground glow is not drawn
GROUND_GLOW_STRENGTH = 0.25   # the coarse fused field is baked into the texture at a quarter alpha
FIRE_TEXTURE_STRENGTH = 0.5   # the flame field is baked at half strength as that year's ground mark
FLAME_H = 0.60            # flame height in km of the x40 frame at full alpha and flicker 1.0
FLAME_LIFT = 0.05         # constant lift of the fire mesh above the terrain (no z-fighting)
FLAME_MIN_ALPHA = 0.05    # fire mesh is blanked (NaN) where the flame alpha is below this
FLICKER_SIGMA = 3.0       # cells of the work grid
FLICKER_PERIOD = 9        # ticks per oscillation cycle (0.75 s at 12 fps)
FLICKER_DRIFT = (0.8, 1.4)  # cells per tick the noise field crawls (rows, cols)
SMOKE_H = 0.9             # km above the flame base
SMOKE_STRENGTH = 0.35     # smoke alpha = flame alpha times this
SMOKE_DRIFT = 2           # cells per tick the plume drifts north-west
SMOKE_PUFF_SPACING = 4    # cells between puffs along the plume
SMOKE_PUFFS = 8           # puffs per plume (plume length = spacing x puffs cells)
SMOKE_DECAY_CELLS = 18.0  # e-folding length of the puff weight
SITE_LIFT_KM = 0.04
SITE_RADIUS_KM = 0.55
SITE_RIM_KM = 0.15        # the dark rim is a slightly larger back-face-only sphere

VERIFIED_RATE_MMYR = -5.0   # Tay rate below this (sinking) counts as radar-verified degradation
VERIFIED_CLEAN_CELLS = 2    # binary opening and closing with a diamond of this radius (about 5 cells)
VERIFIED_MIN_CELLS = 40     # blobs smaller than this are dropped
DEGRADED_ALPHA = 0.55
DEGRADED_NBR = 0.45          # Landsat NBR below this = vegetation loss in that epoch
DEGRADED_MIN_EPOCHS = 2      # epochs with vegetation loss needed (of the 11 valid composites 1990 to 2025)
DEGRADED_OPEN_CELLS = 1      # binary opening with a diamond of this radius (3 cells across): zones, not speckle
LANDSAT_EPOCHS = [1990, 1995, 2000, 2005, 2010, 2015, 2020, 2022, 2023, 2024, 2025]
OUTLINE_LIFT_KM = 0.08
OUTLINE_DASH_KM = (0.5, 0.3)   # on, off along the line
OUTLINE_WIDTH_PX = 3

ATTRIBUTION = ("Terrain: Copernicus DEM. Degraded peat: Tay et al. 2025 radar sinking, "
               "Landsat burn scar. Fires: NASA FIRMS. PeatGuard")


def hex_rgb(h):
    h = h.lstrip("#")
    return np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)], dtype=np.float32)


# Palette (v4). Fire is the only hot colour on the map; char is black.
PALETTE = {
    "base": hex_rgb("#2f6b45"),       # intact peat forest green; hillshaded to about #244f35 to #3f8458
    "degraded": hex_rgb("#b3261e"),   # degraded peat (verified radar sinking or Landsat burn scar)
    "water": hex_rgb("#3f7fbf"),      # OSM rivers (darker, wider than the canals)
    "char": hex_rgb("#000000"),       # burned in an earlier year
    "fire_core": hex_rgb("#ffe066"),
    "fire_mid": hex_rgb("#ff8c1a"),
    "fire_edge": hex_rgb("#d7301f"),
    "smoke": hex_rgb("#8e9bb0"),
    "verified": "#9fe8ff",            # dashed outline of the radar-verified zone
    "site": "#ffffff",
    "site_rim": "#101418",
}
CHAR_ALPHA1 = 0.15        # black alpha after one earlier fire year
CHAR_ALPHA_MAX = 0.75     # black alpha at CHAR_FULL_YEARS or more fire years
CHAR_FULL_YEARS = 8       # the box has cells with up to 14 FIRMS fire years; 15 percent burned once
RIVER_BLEND = 0.85

LEGEND_PANEL = (10, 14, 20, 184)   # rgba 10,14,20 at 0.72
LEGEND_RADIUS = 10
# the MP4 sits at about half the width of a PowerPoint slide, so the type is
# larger than a full-screen HUD would need
YEAR_PX = 110
LEGEND_LABEL_PX = 26
LEGEND_TITLE_PX = 20
LEGEND_LABEL_MAX_W = 640
EXAG_PX = 24
HEADLINE_PX = 40

FONT_REG = None
FONT_BOLD = None


def _font_paths():
    import matplotlib.font_manager as fm
    return fm.findfont("DejaVu Sans"), fm.findfont("DejaVu Sans:bold")


def font(size, bold=False):
    global FONT_REG, FONT_BOLD
    if FONT_REG is None:
        FONT_REG, FONT_BOLD = _font_paths()
    return ImageFont.truetype(FONT_BOLD if bold else FONT_REG, size)


# ----------------------------------------------------------------------------
# Schedule (v4): years, frames per year, holds, camera chaining, flicker ticks
# ----------------------------------------------------------------------------
FPS = 12
XFADE_FRAMES = 8
XFADE_IMAGERY_FRAMES = 12   # the green 2024 frame dissolves into the imagery 2026 frame over 12
IMAGERY_EPOCH = 2026
IMAGERY_NPY = "planet_nicfi_2026_09"   # Planet NICFI monthly mosaic, Sep 2026, on the work grid
IMAGERY_SHADE = 0.5         # hillshade at half strength on the imagery
IMAGERY_DEGRADED_ALPHA = 0.40
EASE_LINEAR = 0.25          # azimuth easing: this much linear plus the rest smoothstep (see ease)

# Specs. frames: output frames of each year (the camera never stops: a year's
# "hold" is only its data layer persisting while the orbit continues). The last
# xfade_into(next) frames of a year are dissolved with the next year rendered at
# the same global frame indices. prefix: frames before the first year (act 1 of
# the two-act cut). The camera is one function of the global frame index over
# prefix + sum(frames) frames.
SPECS = {
    "v4": {
        "years": [1990, 1995, 2000] + list(range(2001, 2025)) + [2026],
        "frames": lambda y: 36 if y in (1990, 1995, 2000) else (60 if y == 2026 else 26),
        "xfade": lambda nxt: XFADE_IMAGERY_FRAMES if nxt == IMAGERY_EPOCH else XFADE_FRAMES,
        "prefix": 0,
    },
    "twoact": {
        "years": list(range(2006, 2025)) + [2026],
        "frames": lambda y: 51 if y == 2026 else (21 if y == 2024 else (23 if y in (2015, 2019) else 13)),
        "xfade": lambda nxt: XFADE_IMAGERY_FRAMES if nxt == IMAGERY_EPOCH else 4,
        "prefix": 85,
    },
}
SPEC_NAME = "v4"
YEARS = SPECS["v4"]["years"]


def set_spec(name):
    global SPEC_NAME, YEARS
    SPEC_NAME = name
    YEARS = SPECS[name]["years"]


def xfade_frames_for(next_year):
    """Crossfade length into a year (frames of the dissolve)."""
    return SPECS[SPEC_NAME]["xfade"](next_year)


def total_frames():
    sp = SPECS[SPEC_NAME]
    return sp["prefix"] + sum(sp["frames"](y) for y in sp["years"])


def schedule():
    """Per-year dict: n (output frames of the year), tick0 (global frame index of
    its first frame), prev (previous year or None). Year Y occupies global frames
    tick0 .. tick0 + n - 1; its last xfade_frames_for(next) frames are dissolved
    with the next year rendered at those same global indices (its pre-roll)."""
    sp = SPECS[SPEC_NAME]
    out, tick, prev = {}, sp["prefix"], None
    for y in sp["years"]:
        n = sp["frames"](y)
        out[y] = {"n": n, "tick0": tick, "prev": prev}
        tick += n
        prev = y
    return out


# ----------------------------------------------------------------------------
# Grid and geometry
# ----------------------------------------------------------------------------
def load_grid_meta():
    with open(os.path.join(DATA_DIR, "dem.json")) as f:
        m = json.load(f)
    return m["shape"], m["lon0"], m["lat0"], m["step"]


SHAPE, LON0, LAT0, STEP = load_grid_meta()
NROW, NCOL = SHAPE
LON_C = LON0 + STEP * NCOL / 2.0
LAT_C = LAT0 - STEP * NROW / 2.0
KM_PER_DEG_LAT = 110.57
KM_PER_DEG_LON = 111.32 * math.cos(math.radians(LAT_C))
CELL_KM = STEP * KM_PER_DEG_LAT


def lonlat_to_km(lon, lat):
    x = (np.asarray(lon, dtype=np.float64) - LON_C) * KM_PER_DEG_LON
    y = (np.asarray(lat, dtype=np.float64) - LAT_C) * KM_PER_DEG_LAT
    return x, y


def lonlat_to_rowcol(lon, lat):
    col = (np.asarray(lon, dtype=np.float64) - LON0) / STEP
    row = (LAT0 - np.asarray(lat, dtype=np.float64)) / STEP
    return row, col


def rowcol_to_lonlat(row, col):
    lon = LON0 + (np.asarray(col, dtype=np.float64) + 0.5) * STEP
    lat = LAT0 - (np.asarray(row, dtype=np.float64) + 0.5) * STEP
    return lon, lat


def load_npy(name):
    return np.load(os.path.join(DATA_DIR, name + ".npy"))


def build_z_full(epoch):
    """Return z (km, exaggerated) on the full grid plus the lowering field (m)."""
    dem = load_npy("dem").astype(np.float32)
    dem = dem - float(np.nanmin(dem))  # ellipsoid offset removed, relative relief only
    z = dem * VEXAG / 1000.0
    lowering = None
    if epoch >= SINK_FROM_EPOCH:
        p = os.path.join(DATA_DIR, "lowering_%d.npy" % epoch)
        if os.path.exists(p):
            lowering = np.load(p).astype(np.float32)
            low0 = np.nan_to_num(lowering, nan=0.0)
            # lowering is negative where the surface is lowered, so adding it
            # pushes the surface down; scaled by VEXAG * SINK_EXAG (m -> km)
            z = z + low0 * VEXAG * SINK_EXAG / 1000.0
    return z.astype(np.float32), lowering


def work_xy():
    """X, Y (km) of the decimated work grid, shape (nr, nc)."""
    rows = np.arange(0, NROW, DECIMATE)
    cols = np.arange(0, NCOL, DECIMATE)
    lon, _ = rowcol_to_lonlat(np.zeros_like(cols), cols)
    _, lat = rowcol_to_lonlat(rows, np.zeros_like(rows))
    xk, _ = lonlat_to_km(lon, np.full_like(lon, LAT_C))
    _, yk = lonlat_to_km(np.full_like(lat, LON_C), lat)
    X, Y = np.meshgrid(xk, yk)  # X[row, col], Y[row, col]
    return X.astype(np.float32), Y.astype(np.float32)


def build_surface(z_full):
    zs = z_full[::DECIMATE, ::DECIMATE]
    nr, nc = zs.shape
    X, Y = work_xy()
    grid = pv.StructuredGrid(X, Y, zs)
    # texture coordinates: u along columns, v along rows, v = 1 at row 0
    # (north); matches the unflipped north-up texture array (see build_texture).
    u = np.linspace(0.0, 1.0, nc, dtype=np.float32)
    v = np.linspace(1.0, 0.0, nr, dtype=np.float32)
    U, V = np.meshgrid(u, v)
    # StructuredGrid point order is Fortran (x fastest) in pyvista for X[row,col]
    tc = np.column_stack([U.ravel(order="F"), V.ravel(order="F")])
    grid.active_texture_coordinates = tc
    return grid


# ----------------------------------------------------------------------------
# Degraded peat: radar-verified sinking zone and Landsat-observed burn scar
# ----------------------------------------------------------------------------
_MASK_CACHE = {}


def verified_sinking_mask():
    """Radar-verified degraded peat: Tay et al. 2025 ALOS-2 rate below
    VERIFIED_RATE_MMYR (negative = sinking; NaN outside the 10.9 percent
    coverage is never verified), cleaned by a binary opening then closing with a
    diamond of radius VERIFIED_CLEAN_CELLS, blobs under VERIFIED_MIN_CELLS dropped."""
    if "verified" in _MASK_CACHE:
        return _MASK_CACHE["verified"]
    rate = load_npy("tay_rate_mmyr").astype(np.float32)
    m = np.nan_to_num(rate, nan=0.0) < VERIFIED_RATE_MMYR
    st = ndi.iterate_structure(ndi.generate_binary_structure(2, 1), VERIFIED_CLEAN_CELLS)
    m = ndi.binary_closing(ndi.binary_opening(m, structure=st), structure=st)
    lab, n = ndi.label(m)
    if n:
        sizes = ndi.sum(m, lab, index=np.arange(1, n + 1))
        keep = np.where(sizes >= VERIFIED_MIN_CELLS)[0] + 1
        m = np.isin(lab, keep)
    _MASK_CACHE["verified"] = m
    return m


def observed_scar_mask():
    """Landsat-observed cumulative burn scar (NBR < 0.05 on land at any epoch to 2026)."""
    if "scar" not in _MASK_CACHE:
        _MASK_CACHE["scar"] = load_npy("burned_cum_2026").astype(bool)
    return _MASK_CACHE["scar"]


def vegetation_loss_mask():
    """Landsat-observed vegetation loss: NBR below DEGRADED_NBR in at least
    DEGRADED_MIN_EPOCHS of the LANDSAT_EPOCHS composites (NaN gaps, e.g. the 2010
    cloud and SLC-off gaps, do not count), on land (persistent water excluded)."""
    if "vegloss" in _MASK_CACHE:
        return _MASK_CACHE["vegloss"]
    count = np.zeros((NROW, NCOL), dtype=np.int16)
    for e in LANDSAT_EPOCHS:
        nbr = load_npy("nbr_%d" % e).astype(np.float32)
        count += (np.nan_to_num(nbr, nan=1.0) < DEGRADED_NBR)
    water = load_npy("water_mask").astype(bool)
    m = (count >= DEGRADED_MIN_EPOCHS) & ~water
    _MASK_CACHE["vegloss"] = m
    return m


def burned_before_mask(year):
    """Ground with a FIRMS fire footprint in any year before `year` (2001 on)."""
    key = ("burned_before", year)
    if key in _MASK_CACHE:
        return _MASK_CACHE[key]
    m = np.zeros((NROW, NCOL), dtype=bool)
    for y in range(2001, year):
        fp = fire_footprint(y)
        if fp is not None:
            m |= fp
    _MASK_CACHE[key] = m
    return m


def degraded_mask(year):
    """PRE-EXISTING degraded peat as of `year`, as OBSERVED by satellites (never
    inferred from canals); the whole box is peat (PeatGrids >= 50 cm). Land where
    Landsat NBR < DEGRADED_NBR in any composite strictly before `year` (NaN gaps
    do not count, persistent water excluded), OR a FIRMS fire footprint in any
    year before `year`, OR the radar-verified sinking zone from SINK_FROM_EPOCH on
    (the Tay measurement starts in 2018). Cleaned with a binary opening of
    DEGRADED_OPEN_CELLS so it reads as zones. The layer grows through the
    sequence, so every flame lands either on red (degraded before it burned) or
    on grey."""
    key = ("degraded", year)
    if key in _MASK_CACHE:
        return _MASK_CACHE[key]
    m = np.zeros((NROW, NCOL), dtype=bool)
    for e in LANDSAT_EPOCHS:
        if e >= year:
            break
        nbr = load_npy("nbr_%d" % e).astype(np.float32)
        m |= np.nan_to_num(nbr, nan=1.0) < DEGRADED_NBR
    m &= ~load_npy("water_mask").astype(bool)
    m |= burned_before_mask(year)
    if year >= SINK_FROM_EPOCH:
        m |= verified_sinking_mask()
    st = ndi.iterate_structure(ndi.generate_binary_structure(2, 1), DEGRADED_OPEN_CELLS)
    m = ndi.binary_opening(m, structure=st)
    _MASK_CACHE[key] = m
    return m


RECURRENCE_JSON = os.path.join(DATA_ROOT, "deliverables", "degraded_peat_fire_recurrence", "fire_on_preexisting_degraded.json")


def fire_on_degraded_fraction(year):
    """Fraction of this year's fire footprint that lies on peat already degraded
    before it burned: the per-year number from the recurrence deliverable
    (fire_on_preexisting_degraded.json, field frac_pre_degraded) when present,
    otherwise computed the same way here against degraded_mask(year). None when
    the year has no fire record."""
    if os.path.exists(RECURRENCE_JSON):
        with open(RECURRENCE_JSON) as f:
            rows = json.load(f).get("rows", [])
        for r in rows:
            if int(r.get("year", -1)) == year and r.get("frac_pre_degraded") is not None:
                return float(r["frac_pre_degraded"])
    fp = fire_footprint(year)
    if fp is None or not fp.any():
        return None
    return float(degraded_mask(year)[fp].mean())


def mask_contours(mask):
    """Closed contour polylines of a mask as (n, 2) arrays of [row, col]."""
    import contourpy
    cg = contourpy.contour_generator(z=mask.astype(np.float32))
    return [np.column_stack([ln[:, 1], ln[:, 0]]) for ln in cg.lines(0.5) if len(ln) >= 2]


def dashed_line_mesh(mask, z_full, on_km, off_km, lift_km):
    """PolyData of dashes along the contour of a mask, lifted above the terrain.
    Dashes are cut by arc length so they are even on every heading."""
    pts_all, lines = [], []
    for rc in mask_contours(mask):
        lon, lat = rowcol_to_lonlat(rc[:, 0], rc[:, 1])
        x, y = lonlat_to_km(lon, lat)
        seg = np.hypot(np.diff(x), np.diff(y))
        s = np.concatenate([[0.0], np.cumsum(seg)])
        total = s[-1]
        if total <= 0:
            continue
        period = on_km + off_km
        d0 = 0.0
        while d0 < total:
            d1 = min(d0 + on_km, total)
            # sample the dash densely enough to follow the contour
            n = max(2, int((d1 - d0) / (0.5 * CELL_KM)) + 1)
            ss = np.linspace(d0, d1, n)
            xs = np.interp(ss, s, x)
            ys = np.interp(ss, s, y)
            lons, lats = (xs / KM_PER_DEG_LON + LON_C), (ys / KM_PER_DEG_LAT + LAT_C)
            zs = surface_z_at(z_full, lons, lats) + lift_km
            i0 = len(pts_all)
            pts_all.extend(np.column_stack([xs, ys, zs]).tolist())
            lines.append([n] + list(range(i0, i0 + n)))
            d0 += period
    if not pts_all:
        return None
    poly = pv.PolyData(np.asarray(pts_all, dtype=np.float32))
    poly.lines = np.concatenate([np.asarray(ln, dtype=np.int64) for ln in lines])
    return poly


# ----------------------------------------------------------------------------
# Fire fields
# ----------------------------------------------------------------------------
def fire_footprint(year, sigma_km=0.4):
    """Binary 'ground that burned this year' mask: detections rasterised and
    smoothed with a fine Gaussian so each MODIS pixel covers about its own
    footprint, not the fused front used for the glow."""
    pts = load_fires(year)
    if pts is None or len(pts) == 0:
        return None
    row, col = lonlat_to_rowcol(pts[:, 0], pts[:, 1])
    ri = np.clip(np.round(row).astype(int), 0, NROW - 1)
    ci = np.clip(np.round(col).astype(int), 0, NCOL - 1)
    field = np.zeros((NROW, NCOL), dtype=np.float32)
    np.add.at(field, (ri, ci), 1.0)
    sig = sigma_km / CELL_KM
    dens = gaussian_filter(field, sigma=sig, mode="constant")
    peak = 1.0 / (2.0 * math.pi * sig * sig)
    return dens > 0.35 * peak


_CHAR_CACHE = {}


def build_char(year):
    """Persistent char layer: ground that burned in any previous FIRMS year is
    black, more opaque for more fire years: one year CHAR_ALPHA1, CHAR_FULL_YEARS
    or more CHAR_ALPHA_MAX, linear in between. Returns alpha or None."""
    if year in _CHAR_CACHE:
        return _CHAR_CACHE[year]
    count = None
    for y in range(2001, year):
        fp = fire_footprint(y)
        if fp is None:
            continue
        f = fp.astype(np.float32)
        count = f if count is None else count + f
    if count is None:
        _CHAR_CACHE[year] = None
        return None
    k = np.clip((count - 1.0) / (CHAR_FULL_YEARS - 1.0), 0.0, 1.0)   # 0 at one year, 1 at CHAR_FULL_YEARS
    alpha = np.where(count >= 1, CHAR_ALPHA1 + (CHAR_ALPHA_MAX - CHAR_ALPHA1) * k, 0.0).astype(np.float32)
    alpha = gaussian_filter(alpha, sigma=1.0).astype(np.float32)
    _CHAR_CACHE[year] = alpha
    return alpha


def fire_ramp(s):
    """Yellow-orange-red ramp: edge (s = 0) to mid (0.5) to core (1). s is any
    shape; returns (..., 3) float32."""
    s = np.clip(np.asarray(s, dtype=np.float32), 0.0, 1.0)
    lo = np.clip(s * 2.0, 0.0, 1.0)[..., None]
    hi = np.clip(s * 2.0 - 1.0, 0.0, 1.0)[..., None]
    c = PALETTE["fire_edge"] * (1.0 - lo) + PALETTE["fire_mid"] * lo
    c = np.where(s[..., None] > 0.5, PALETTE["fire_mid"] * (1.0 - hi) + PALETTE["fire_core"] * hi, c)
    return c.astype(np.float32)


def load_fires(epoch):
    """Return an (n, 2) array of [lon, lat] detections, or None before 2001.

    For 2026 the input is the live 7-day feed aggregated to 0.05 deg cells
    with a count per cell; each cell is expanded to its count of points spread
    uniformly (seeded) inside the cell so the glow is not a grid of stamps.
    """
    if epoch < 2001:
        return None
    if epoch == 2026:
        p = os.path.join(DATA_DIR, "fires_2026_week.json")
        if os.path.exists(p):
            with open(p) as f:
                w = json.load(f)
            rng = np.random.default_rng(2026)
            out = []
            for c in w["cells"]:
                if not (114.20 <= c[0] <= 114.83 and -3.01 <= c[1] <= -2.41):
                    continue
                n = int(c[2])
                out.append(np.column_stack([c[0] + rng.uniform(-0.025, 0.025, n),
                                            c[1] + rng.uniform(-0.025, 0.025, n)]))
            return np.vstack(out) if out else np.zeros((0, 2))
    p = os.path.join(DATA_DIR, "fires_%d.json" % epoch)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return np.array(json.load(f), dtype=np.float64).reshape(-1, 2)


def build_fire_glow(epoch):
    """Fire fields of a year on the full grid, or None.

    Returns (flame, s, ground):
      flame  in [0, FIRE_ALPHA_MAX]: opacity and height of the 3D flames, from the
             fine detection density (sigma FIRE_SIGMA_FINE_KM): zero below
             FLAME_FLOOR_PER_KM2, full at FLAME_FULL_PER_KM2, so flames stand only
             on actual detection clusters, plus a small FLAME_SINGLE_ALPHA flame
             right at a lone detection;
      s      in [0, 1]: colour ramp position, yellow core at FLAME_CORE_PER_KM2;
      ground in [0, 1]: low ground glow baked into the texture, the coarse fused
             field (sigma FIRE_SIGMA_KM) at GROUND_GLOW_STRENGTH of its alpha plus
             the flame field at FIRE_TEXTURE_STRENGTH.
    Each detection is jittered by FIRE_JITTER_M and the widths vary by a seeded
    0.8 to 1.3 per year so no two years look stamped.
    """
    pts = load_fires(epoch)
    if pts is None or len(pts) == 0:
        return None
    rng = np.random.default_rng(epoch)
    cell_area = CELL_KM * CELL_KM
    jit = FIRE_JITTER_M / (CELL_KM * 1000.0)
    row, col = lonlat_to_rowcol(pts[:, 0], pts[:, 1])
    row = row + rng.normal(0.0, jit, len(row))
    col = col + rng.normal(0.0, jit, len(col))
    ri = np.clip(np.round(row).astype(int), 0, NROW - 1)
    ci = np.clip(np.round(col).astype(int), 0, NCOL - 1)
    field = np.zeros((NROW, NCOL), dtype=np.float32)
    np.add.at(field, (ri, ci), 1.0)
    f_year = rng.uniform(0.8, 1.3)
    sig_c = FIRE_SIGMA_KM / CELL_KM * f_year
    sig_f = FIRE_SIGMA_FINE_KM / CELL_KM * f_year
    dens_c = gaussian_filter(field, sigma=sig_c, mode="constant") / cell_area   # detections per km2
    dens_f = gaussian_filter(field, sigma=sig_f, mode="constant") / cell_area
    noise = gaussian_filter(rng.standard_normal((NROW, NCOL)).astype(np.float32), sigma=6.0)
    noise = noise / max(float(noise.std()), 1e-6)
    shape = np.clip(1.0 + 0.35 * noise, 0.55, 1.25)

    flame = FIRE_ALPHA_MAX * np.clip((dens_f - FLAME_FLOOR_PER_KM2) / (FLAME_FULL_PER_KM2 - FLAME_FLOOR_PER_KM2), 0.0, 1.0) ** 0.8
    # a lone detection: a small flame right at the detection (sharp in density), so
    # quiet years show a few distinct fires instead of nothing
    sig_s = FLAME_SINGLE_SIGMA_KM / CELL_KM
    dens_s = gaussian_filter(field, sigma=sig_s, mode="constant") / cell_area
    single_peak = 1.0 / (2.0 * math.pi * sig_s * sig_s) / cell_area
    a_single = FLAME_SINGLE_ALPHA * np.clip(dens_s / single_peak, 0.0, 1.0) ** 2
    flame = np.clip(np.maximum(flame, a_single) * shape, 0.0, FIRE_ALPHA_MAX)
    flame = np.where(flame < 0.03, 0.0, flame)
    s = np.clip(dens_f / FLAME_CORE_PER_KM2, 0.0, 1.0) ** 0.7

    a_coarse = 0.8 * FIRE_ALPHA_MAX * np.clip(
        (dens_c - FIRE_FLOOR_PER_KM2) / (FIRE_CORE_PER_KM2 - FIRE_FLOOR_PER_KM2), 0.0, 1.0) ** 0.8
    ground = np.maximum(GROUND_GLOW_STRENGTH * a_coarse * shape, FIRE_TEXTURE_STRENGTH * flame)
    return flame.astype(np.float32), s.astype(np.float32), np.clip(ground, 0.0, 1.0).astype(np.float32)


# ----------------------------------------------------------------------------
# Texture
# ----------------------------------------------------------------------------
def build_texture(year, fire, red_field=None, with_char=True):
    """Texture for a year. red_field (float 0 to 1 on the full grid) replaces the
    boolean degraded_mask(year) when given (used to interpolate the red between
    years); with_char=False leaves the char layer off (act 1 of the two-act cut)."""
    hs = load_npy("hillshade").astype(np.float32)
    hmin, hmax = float(hs.min()), float(hs.max())
    shade = 0.75 + 0.50 * (hs - hmin) / max(hmax - hmin, 1e-6)
    imagery = imagery_rgb(year)
    if imagery is None:
        rgb = np.broadcast_to(PALETTE["base"], (NROW, NCOL, 3)).astype(np.float32) * shade[..., None]
        red_alpha = DEGRADED_ALPHA
    else:
        # real imagery for the snapshot year: hillshade at half strength so the
        # imagery keeps its colour, red at a lower alpha so it shows through
        rgb = imagery * (1.0 - IMAGERY_SHADE * (1.0 - shade))[..., None]
        red_alpha = IMAGERY_DEGRADED_ALPHA

    # degraded peat as of this year (pre-existing, grows through the sequence)
    deg = degraded_mask(year).astype(np.float32) if red_field is None else np.clip(red_field, 0.0, 1.0).astype(np.float32)
    a = (red_alpha * deg)[..., None]
    rgb = rgb * (1 - a) + PALETTE["degraded"] * a

    rivers = load_npy("rivers").astype(bool)
    rgb = np.where(rivers[..., None], rgb * (1 - RIVER_BLEND) + PALETTE["water"] * RIVER_BLEND, rgb)

    char = build_char(year) if with_char else None
    if char is not None:
        rgb = rgb * (1 - char[..., None]) + PALETTE["char"] * char[..., None]


    if fire is not None:
        _, s, ground = fire
        rgb = rgb * (1 - ground[..., None]) + fire_ramp(s) * ground[..., None]

    rgb = np.clip(rgb, 0, 255).astype(np.uint8)
    # pv.Texture places array row 0 at v = 1 (verified with a corner-marker
    # test), so the north-up image is passed unflipped with v = 1 at grid row 0
    return pv.Texture(np.ascontiguousarray(rgb))


def imagery_rgb(year):
    """Planet NICFI mosaic (Sep 2026) as float RGB on the work grid for the
    imagery epoch only; None for every other year (the imagery is a 2026
    snapshot, so the historical years keep the green base)."""
    if year != IMAGERY_EPOCH:
        return None
    if "imagery" not in _MASK_CACHE:
        p = os.path.join(DATA_DIR, IMAGERY_NPY + ".npy")
        if not os.path.exists(p):
            raise SystemExit("imagery missing: %s" % p)
        a = np.load(p)
        assert a.shape == (NROW, NCOL, 3), a.shape
        _MASK_CACHE["imagery"] = a.astype(np.float32)
    return _MASK_CACHE["imagery"]


def surface_z_at(z_full, lon, lat):
    row, col = lonlat_to_rowcol(lon, lat)
    row = np.clip(row - 0.5, 0, NROW - 1)
    col = np.clip(col - 0.5, 0, NCOL - 1)
    return map_coordinates(z_full, [row, col], order=1, mode="nearest")


# ----------------------------------------------------------------------------
# Overlays
# ----------------------------------------------------------------------------
def load_sites():
    with open(os.path.join(DATA_DIR, "sites.json")) as f:
        return json.load(f)


def load_settlements():
    with open(os.path.join(DATA_DIR, "settlements.json")) as f:
        s = json.load(f)
    by_name = {x["name"]: x for x in s}
    out = [by_name[n] for n in LABEL_NAMES if n in by_name]
    return out[:N_SETTLEMENT_LABELS]


def points_mesh(lonlat, z_full, lift_km):
    x, y = lonlat_to_km(lonlat[:, 0], lonlat[:, 1])
    z = surface_z_at(z_full, lonlat[:, 0], lonlat[:, 1]) + lift_km
    return pv.PolyData(np.column_stack([x, y, z]).astype(np.float32))


# ----------------------------------------------------------------------------
# Fire animation (3D flame and smoke meshes updated in place per tick)
# ----------------------------------------------------------------------------
def shift2d(a, dr, dc):
    """Shift an array by (dr, dc) cells with zero fill (no wrap-around)."""
    out = np.zeros_like(a)
    nr, nc = a.shape
    r0, r1 = max(dr, 0), min(nr + dr, nr)
    c0, c1 = max(dc, 0), min(nc + dc, nc)
    if r1 <= r0 or c1 <= c0:
        return out
    out[r0:r1, c0:c1] = a[r0 - dr:r1 - dr, c0 - dc:c1 - dc]
    return out


class FireAnimator:
    """Owns the flame and smoke StructuredGrids and rewrites their points and
    rgba scalars for a given tick. Geometry: the work grid (decimated) with
    z = terrain + FLAME_LIFT + FLAME_H * (alpha / FIRE_ALPHA_MAX) * flicker(tick);
    points with alpha below FLAME_MIN_ALPHA get z = NaN (blanked). The flicker
    field is a seeded noise (sigma FLICKER_SIGMA cells) that drifts FLICKER_DRIFT
    cells per tick and swings between two independent fields with period
    FLICKER_PERIOD, clipped to 0.7 to 1.3; it also drives opacity and colour so
    the shimmer reads at 1080p. The same seed and tick drive every year, so the
    crossfade frames of two years share one flicker state."""

    def __init__(self, z_work, X, Y, flame_full, s_full):
        self.nr, self.nc = z_work.shape
        self.z = z_work.astype(np.float32)
        self.Xf = X.ravel(order="F")
        self.Yf = Y.ravel(order="F")
        self.alpha = flame_full[::DECIMATE, ::DECIMATE].astype(np.float32)
        self.s = s_full[::DECIMATE, ::DECIMATE].astype(np.float32)
        rng = np.random.default_rng(7)
        self.na = self._noise(rng)
        self.nb = self._noise(rng)
        self.flame = pv.StructuredGrid(X, Y, self.z + FLAME_LIFT)
        self.smoke = pv.StructuredGrid(X, Y, self.z + FLAME_LIFT + SMOKE_H)
        self.flame.point_data["rgba"] = np.zeros((self.flame.n_points, 4), np.uint8)
        self.smoke.point_data["rgba"] = np.zeros((self.smoke.n_points, 4), np.uint8)
        self.update(0)

    def _noise(self, rng):
        n = gaussian_filter(rng.standard_normal((self.nr, self.nc)).astype(np.float32), sigma=FLICKER_SIGMA)
        return (n / max(float(n.std()), 1e-6)).astype(np.float32)

    def flicker(self, tick):
        th = 2.0 * math.pi * tick / FLICKER_PERIOD
        dr = int(round(FLICKER_DRIFT[0] * tick))
        dc = int(round(FLICKER_DRIFT[1] * tick))
        n = (math.cos(th) * np.roll(self.na, (dr, dc), axis=(0, 1))
             + math.sin(th) * np.roll(self.nb, (-dr, dc), axis=(0, 1)))
        return 1.0 + 0.3 * np.clip(n, -1.0, 1.0)

    def smoke_alpha(self, tick):
        """Plume of SMOKE_PUFFS puffs spaced SMOKE_PUFF_SPACING cells along the
        north-west direction from the flame cores; the whole plume advances
        SMOKE_DRIFT cells per tick and each new puff fades in from zero weight."""
        phase = (SMOKE_DRIFT * tick) % SMOKE_PUFF_SPACING
        out = np.zeros_like(self.alpha)
        for k in range(SMOKE_PUFFS):
            d = phase + SMOKE_PUFF_SPACING * k
            w = min(d / SMOKE_PUFF_SPACING, 1.0) * math.exp(-d / SMOKE_DECAY_CELLS)
            if w <= 0.0:
                continue
            di = int(round(d))
            out = np.maximum(out, w * shift2d(self.alpha, -di, -di))
        return out * SMOKE_STRENGTH

    @staticmethod
    def _rgba_f(rgba):
        # (nr, nc, 4) -> Fortran point order of the StructuredGrid
        return np.ascontiguousarray(rgba.transpose(1, 0, 2).reshape(-1, 4))

    def update(self, tick):
        fl = self.flicker(tick)
        a_norm = self.alpha / FIRE_ALPHA_MAX
        z = self.z + FLAME_LIFT + FLAME_H * a_norm * fl
        keep = self.alpha > FLAME_MIN_ALPHA
        z = np.where(keep, z, np.nan).astype(np.float32)
        opac = np.clip(self.alpha * (0.35 + 0.65 * fl), 0.0, 1.0)
        s_eff = np.clip(self.s + 0.5 * (fl - 1.0), 0.0, 1.0)
        rgba = np.concatenate([fire_ramp(s_eff), (opac * 255.0)[..., None]], axis=-1)
        self.flame.points = np.column_stack([self.Xf, self.Yf, z.ravel(order="F")]).astype(np.float32)
        self.flame.point_data["rgba"] = self._rgba_f(np.clip(rgba, 0, 255).astype(np.uint8))

        sa = self.smoke_alpha(tick)
        zs = self.z + FLAME_LIFT + SMOKE_H + 0.5 * (fl - 1.0)
        zs = np.where(sa > 0.02, zs, np.nan).astype(np.float32)
        srgb = np.broadcast_to(PALETTE["smoke"], sa.shape + (3,))
        srgba = np.concatenate([srgb, (np.clip(sa, 0, 1) * 255.0)[..., None]], axis=-1)
        self.smoke.points = np.column_stack([self.Xf, self.Yf, zs.ravel(order="F")]).astype(np.float32)
        self.smoke.point_data["rgba"] = self._rgba_f(np.clip(srgba, 0, 255).astype(np.uint8))


# ----------------------------------------------------------------------------
# Camera
# ----------------------------------------------------------------------------
def ease(u):
    """Azimuth easing over the whole sequence: EASE_LINEAR of a linear ramp plus
    the rest smoothstep. Pure smoothstep starts at zero velocity, which makes the
    first frame steps grow several-fold frame to frame; the linear share keeps
    every consecutive azimuth step within a few percent of its neighbours while
    the orbit still eases in and out (end velocity EASE_LINEAR of the mean)."""
    u = float(np.clip(u, 0.0, 1.0))
    return EASE_LINEAR * u + (1.0 - EASE_LINEAR) * (3.0 * u * u - 2.0 * u * u * u)


def camera(t):
    """Oblique orbit at parameter t in [0, 1]: azimuth CAM_AZ0 -> CAM_AZ1 (deg,
    0 = from the south), elevation, distance, focal point and up fixed."""
    t = float(np.clip(t, 0.0, 1.0))
    az = math.radians(CAM_AZ0_DEG + (CAM_AZ1_DEG - CAM_AZ0_DEG) * t)
    el = math.radians(CAM_ELEV_DEG)
    focal = (0.0, 0.0, 2.0)
    pos = (focal[0] + CAM_DIST_KM * math.cos(el) * math.sin(az),
           focal[1] - CAM_DIST_KM * math.cos(el) * math.cos(az),
           focal[2] + CAM_DIST_KM * math.sin(el))
    return pos, focal, (0.0, 0.0, 1.0)


def camera_t_at(g):
    """The single camera parameter of global output frame g (0 .. total - 1)."""
    n = total_frames()
    return ease(g / float(max(n - 1, 1)))


def camera_azimuth_deg(g):
    return CAM_AZ0_DEG + (CAM_AZ1_DEG - CAM_AZ0_DEG) * camera_t_at(g)


def check_azimuth_steps(gs):
    """gs: the global frame index of every entry of a concat list, in order.
    Asserts the azimuth is strictly monotonic and that every consecutive step is
    within 10 percent of its neighbours; returns (max_step_deg, max_ratio)."""
    az = np.array([camera_azimuth_deg(g) for g in gs], dtype=np.float64)
    d = np.diff(az)
    assert (d > 0).all(), "azimuth not strictly increasing at list entries %s" % np.where(d <= 0)[0][:10]
    r = d[1:] / d[:-1]
    worst = float(np.max(np.abs(r - 1.0)))
    assert worst <= 0.10, "azimuth step changes by %.1f percent" % (100 * worst)
    return float(d.max()), float(1.0 + worst)


def set_camera(pl, t):
    pos, focal, up = camera(t)
    pl.camera.position = pos
    pl.camera.focal_point = focal
    pl.camera.up = up
    pl.camera.view_angle = CAM_VIEW_ANGLE
    pl.camera.SetWindowCenter(*CAM_WINDOW_CENTER)


# ----------------------------------------------------------------------------
# HUD (PIL compositing on the rendered frame)
# ----------------------------------------------------------------------------
def draw_text(d, xy, text, size, fill=(255, 255, 255), bold=False, shadow=True):
    f = font(size, bold)
    if shadow:
        d.text((xy[0] + 2, xy[1] + 2), text, font=f, fill=(0, 0, 0, 160))
    d.text(xy, text, font=f, fill=fill)
    return d.textbbox(xy, text, font=f)


def _over(base, tint, alpha):
    return tuple(int(round(b * (1 - alpha) + t * alpha)) for b, t in zip(base, tint))


def legend_rows(year, fire_rows=True):
    """(label, swatch) rows; swatch is a colour tuple or a named style.
    fire_rows=False drops the Burned before and Burning rows (act 1 of the cut)."""
    shaded_base = tuple(int(c * 0.88) for c in PALETTE["base"])
    red = _over(shaded_base, PALETTE["degraded"], DEGRADED_ALPHA)
    rows = [
        ("Degraded before this year (Landsat vegetation loss, earlier fire, or radar sinking)", red),
        ("Verified by radar sinking (Tay et al. 2025)", "verified"),
    ]
    if year >= 2002 and fire_rows:
        rows.append(("Burned before (black, darker = more years)", "char"))
    if year >= 2001 and fire_rows:
        rows.append(("Burning now (FIRMS, last 7 days)" if year == 2026 else "Burning this year", "fire"))
    rows.append(("Dam site", "site"))
    return rows, red


def draw_swatch(d, x, y, w, h, swatch, base):
    if swatch == "char":
        for i in range(w):
            k = i / float(w - 1)
            a = CHAR_ALPHA1 + (CHAR_ALPHA_MAX - CHAR_ALPHA1) * k
            d.line([(x + i, y), (x + i, y + h - 1)], fill=_over(base, PALETTE["char"], a) + (255,))
    elif swatch == "fire":
        for i in range(w):
            c = fire_ramp(np.float32(i / float(w - 1)))
            d.line([(x + i, y), (x + i, y + h - 1)], fill=tuple(int(v) for v in c) + (255,))
    elif swatch == "verified":
        # dashed cyan line through the swatch, as drawn on the terrain
        cy = y + h // 2
        col = tuple(int(c) for c in hex_rgb(PALETTE["verified"])) + (255,)
        dash, gap, i = 7, 4, 0
        while i < w:
            d.line([(x + i, cy), (x + min(i + dash, w) - 1, cy)], fill=col, width=3)
            i += dash + gap
    elif swatch == "site":
        cx, cy, r = x + w / 2.0, y + h / 2.0, h / 2.0
        d.ellipse([cx - r, cy - r, cx + r, cy + r], fill=(255, 255, 255, 255), outline=(16, 20, 24, 255), width=2)
    else:
        d.rectangle([x, y, x + w - 1, y + h - 1], fill=tuple(swatch) + (255,))


def draw_legend(d, year, x_left, y_bottom, fire_rows=True):
    """Compact legend panel anchored at its bottom-left corner; long labels wrap."""
    rows, base = legend_rows(year, fire_rows)
    f_title = font(LEGEND_TITLE_PX)
    f_row = font(LEGEND_LABEL_PX)
    title = "Fires: NASA FIRMS; degraded peat: Tay et al. 2025, Landsat"
    foot = "Imagery: Planet NICFI mosaic, Sep 2026" if year == IMAGERY_EPOCH else None
    sw, sh, gap = 28, 18, 12
    pad_x, pad_top, pad_bot = 18, 14, 12
    line_h = LEGEND_LABEL_PX + 6
    row_gap = 10
    wrapped = [(wrap_text(d, lab, f_row, LEGEND_LABEL_MAX_W), sw_) for lab, sw_ in rows]
    lab_w = max(d.textlength(ln, font=f_row) for lines, _ in wrapped for ln in lines)
    lab_w = max(lab_w, d.textlength(title, font=f_title))
    box_w = pad_x + sw + gap + lab_w + pad_x
    title_h = LEGEND_TITLE_PX + 12
    foot_h = LEGEND_TITLE_PX + 10 if foot else 0
    rows_h = sum(len(lines) * line_h + row_gap for lines, _ in wrapped) - row_gap
    box_h = pad_top + title_h + rows_h + foot_h + pad_bot
    x0, y0 = x_left, y_bottom - box_h
    d.rounded_rectangle([x0, y0, x0 + box_w, y0 + box_h], radius=LEGEND_RADIUS, fill=LEGEND_PANEL)
    d.text((x0 + pad_x, y0 + pad_top), title, font=f_title, fill=(255, 255, 255, 178))
    y = y0 + pad_top + title_h
    for lines, swatch in wrapped:
        draw_swatch(d, x0 + pad_x, y + (line_h - sh) // 2, sw, sh, swatch, base)
        for ln in lines:
            d.text((x0 + pad_x + sw + gap, y), ln, font=f_row, fill=(255, 255, 255, 255))
            y += line_h
        y += row_gap
    if foot:
        d.text((x0 + pad_x, y + 2), foot, font=f_title, fill=(255, 255, 255, 178))
    return box_w, box_h


def wrap_text(d, text, f, max_w):
    words = text.split()
    lines, cur = [], ""
    for w in words:
        trial = (cur + " " + w).strip()
        if d.textlength(trial, font=f) <= max_w or not cur:
            cur = trial
        else:
            lines.append(cur)
            cur = w
    if cur:
        lines.append(cur)
    return lines


def hud_alpha_for(k, n):
    """HUD opacity for crossfade frame k of n: the text fades out over the two
    frames before the midpoint and back in over the two after it, so the year,
    headline and legend never double-expose; the year switches at the midpoint."""
    mid = n // 2
    d = k - mid            # -1, -2 before the switch; 0, 1 after
    if d in (-2, 1):
        return 0.66
    if d in (-1, 0):
        return 0.33
    return 1.0


def compose_minimal(frame, epoch, degraded_pct=None, hud_alpha=1.0):
    """Slide mode: terrain, the year top-left, the legend bottom-left and the
    exaggeration note bottom-right. degraded_pct (act 1 of the two-act cut)
    replaces the fire headline with "Peat degraded: N percent of the box" and
    drops the fire rows from the legend. hud_alpha scales the whole overlay
    (crossfade text fade)."""
    im = Image.fromarray(np.asarray(frame)).convert("RGBA")
    W, H = im.size
    if BARE:
        # land only, except the year so the viewer can follow time
        ov = Image.new("RGBA", im.size, (0, 0, 0, 0))
        d = ImageDraw.Draw(ov)
        draw_text(d, (40, 24), str(epoch), 96, bold=True)
        if hud_alpha < 1.0:
            a = ov.split()[3].point(lambda v: int(v * hud_alpha)); ov.putalpha(a)
        return Image.alpha_composite(im, ov).convert("RGB")
    ov = Image.new("RGBA", im.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(ov)
    draw_text(d, (40, 20), str(epoch), YEAR_PX, bold=True)
    if degraded_pct is None:
        draw_headline(d, epoch, 44, 20 + YEAR_PX + 14, W - 88)
    else:
        draw_text(d, (44, 20 + YEAR_PX + 14), "Peat degraded: %d percent of the box" % int(round(degraded_pct)), HEADLINE_PX)
    f = font(EXAG_PX)
    txt = "terrain height x%d" % VEXAG + ("; sinking x%d" % SINK_EXAG if epoch >= SINK_FROM_EPOCH else "")
    tw = d.textlength(txt, font=f)
    draw_text(d, (W - 40 - tw, H - 40 - EXAG_PX - 6), txt, EXAG_PX)
    draw_legend(d, epoch, 40, H - 40, fire_rows=degraded_pct is None)
    if hud_alpha < 1.0:
        a = ov.getchannel("A").point(lambda v: int(v * hud_alpha))
        ov.putalpha(a)
    return Image.alpha_composite(im, ov).convert("RGB")


def headline_text(epoch):
    """One line under the year from 2001: the share of this year's fire that fell
    on peat already degraded before it burned (None before the FIRMS record)."""
    if epoch < 2001:
        return None
    frac = fire_on_degraded_fraction(epoch)
    if frac is None:
        return None
    pct = int(round(100.0 * frac))
    if epoch == 2026:
        # two balanced lines: the 7-day qualifier does not fit on one
        return ["%d percent of the fire in the FIRMS last 7 days" % pct,
                "is on peat already degraded before it burned"]
    return ["%d percent of this year's fire is on peat already degraded before it burned" % pct]


def draw_headline(d, epoch, x, y, max_w):
    lines = headline_text(epoch)
    if not lines:
        return
    f = font(HEADLINE_PX)
    for txt in lines:
        for ln in wrap_text(d, txt, f, max_w):
            draw_text(d, (x, y), ln, HEADLINE_PX)
            y += HEADLINE_PX + 8


def compose_hud(frame, epoch, caption, has_sinking):
    im = Image.fromarray(frame).convert("RGBA")
    W, H = im.size
    ov = Image.new("RGBA", im.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(ov)
    right_col_w = 420
    cap_max_w = W - 48 - right_col_w - 60 - 48

    # top-left: year, caption (wrapped to two lines)
    draw_text(d, (44, 22), str(epoch), YEAR_PX, bold=True)
    y = 150
    if caption:
        f28 = font(28)
        for ln in wrap_text(d, caption, f28, cap_max_w)[:2]:
            draw_text(d, (48, y), ln, 28)
            y += 38

    # top-right: exaggeration statement (always on screen)
    xr = W - 48
    f26 = font(26, bold=True)
    txt = "terrain height x%d" % VEXAG
    draw_text(d, (xr - d.textlength(txt, font=f26), 30), txt, 26, bold=True)
    if has_sinking:
        txt2 = "sinking x%d on top of terrain x%d (modelled from the Tay et al. 2025 rate)" % (SINK_EXAG, VEXAG)
        draw_text(d, (xr - d.textlength(txt2, font=font(22, True)), 70), txt2, 22, bold=True)

    # bottom-left legend, bottom-right attribution
    draw_legend(d, epoch, 40, H - 40)
    f19 = font(19)
    draw_text(d, (xr - d.textlength(ATTRIBUTION, font=f19), H - 44), ATTRIBUTION, 19,
              fill=(225, 225, 225), shadow=True)
    return Image.alpha_composite(im, ov).convert("RGB")


# ----------------------------------------------------------------------------
# Scene
# ----------------------------------------------------------------------------
BG_OVERRIDE = ""


def build_scene(epoch, with_fire=True):
    """with_fire=False: no flame or smoke mesh and no char or glow in the texture."""
    z_full, lowering = build_z_full(epoch)
    grid = build_surface(z_full)
    fire = build_fire_glow(epoch) if with_fire else None
    tex = build_texture(epoch, fire, with_char=with_fire)

    pl = pv.Plotter(off_screen=True, window_size=WINDOW, lighting="none")
    if BG_OVERRIDE:
        pl.set_background(BG_OVERRIDE)
    else:
        pl.set_background("#0a0f1c", top="#1a2a44")
    pl.enable_anti_aliasing("fxaa")
    pl.enable_depth_peeling(number_of_peels=6, occlusion_ratio=0.0)

    pl.add_mesh(grid, texture=tex, smooth_shading=True, show_edges=False,
                ambient=0.50, diffuse=0.62, specular=0.03, specular_power=12,
                name="terrain")

    # soft key light from the north-west, high; a weak fill from the camera
    key = pv.Light(position=(-60, 80, 90), focal_point=(0, 0, 0), light_type="scene light",
                   intensity=0.95, color="white")
    fill = pv.Light(light_type="headlight", intensity=0.25)
    pl.add_light(key)
    pl.add_light(fill)

    # radar-verified zone outline: dashed cyan line lifted above the terrain
    outline = dashed_line_mesh(verified_sinking_mask(), z_full, OUTLINE_DASH_KM[0], OUTLINE_DASH_KM[1],
                               OUTLINE_LIFT_KM)
    if outline is not None:
        pl.add_mesh(outline, color=PALETTE["verified"], line_width=OUTLINE_WIDTH_PX, lighting=False,
                    name="verified_outline")

    # dam sites: white spheres with a thin dark rim (a larger sphere drawn
    # back-faces only, so only its silhouette shows around the white one)
    sites = [] if BARE else load_sites()
    if sites:
        ll = np.array([[s["lon"], s["lat"]] for s in sites], dtype=np.float64)
        pts = points_mesh(ll, z_full, SITE_LIFT_KM)
        rim = pts.glyph(geom=pv.Sphere(radius=SITE_RADIUS_KM + SITE_RIM_KM, theta_resolution=16, phi_resolution=16),
                        scale=False, orient=False)
        pl.add_mesh(rim, color=PALETTE["site_rim"], lighting=False, culling="front", name="site_rim")
        core = pts.glyph(geom=pv.Sphere(radius=SITE_RADIUS_KM, theta_resolution=16, phi_resolution=16),
                         scale=False, orient=False)
        pl.add_mesh(core, color=PALETTE["site"], ambient=0.75, diffuse=0.35, specular=0.1,
                    smooth_shading=True, name="site_core")

    # settlement labels
    setts = [] if BARE else load_settlements()
    if setts:
        ll = np.array([[s["lon"], s["lat"]] for s in setts], dtype=np.float64)
        pm = points_mesh(ll, z_full, 0.25)
        pl.add_point_labels(pm, [s["name"] for s in setts], font_size=16, text_color="white",
                            shape=None, shadow=False, show_points=True, point_color="white",
                            point_size=5, always_visible=True, bold=False,
                            font_family="arial")

    # north marker at the middle of the northern edge
    xn, yn = lonlat_to_km(LON_C, LAT0)
    zn = float(surface_z_at(z_full, np.array([LON_C]), np.array([LAT0 - STEP]))[0])
    if not BARE:
      pl.add_point_labels(pv.PolyData(np.array([[xn, yn + 1.5, zn + 0.4]], dtype=np.float32)),
                        ["N"], font_size=20, text_color="#e6e6e6", shape=None, shadow=False,
                        show_points=False, always_visible=True, bold=True, font_family="arial")

    # animated fire: flame mesh and smoke mesh, updated in place per tick
    anim = None
    if fire is not None and float(fire[0].max()) > FLAME_MIN_ALPHA:
        X, Y = work_xy()
        anim = FireAnimator(z_full[::DECIMATE, ::DECIMATE], X, Y, fire[0], fire[1])
        pl.add_mesh(anim.flame, scalars="rgba", rgba=True, lighting=False, name="flame")
        pl.add_mesh(anim.smoke, scalars="rgba", rgba=True, lighting=False, name="smoke")

    pl.camera.view_angle = CAM_VIEW_ANGLE
    return pl, anim, lowering is not None


def render_frame(pl, anim, t, tick):
    """One frame at camera parameter t (the legacy path) and flicker tick."""
    set_camera(pl, t)
    if anim is not None:
        anim.update(tick)
    pl.render()
    return pl.screenshot(return_img=True, window_size=WINDOW)


def render_global_frame(pl, anim, g):
    """One frame at global output index g: camera and flicker both from g."""
    return render_frame(pl, anim, camera_t_at(g), g)


def compose(frame, epoch, args, has_sinking):
    if args.minimal:
        return compose_minimal(frame, epoch)
    return compose_hud(frame, epoch, args.caption, has_sinking)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def render_year(year, args):
    """One year of the active spec: its pre-roll (the year rendered at the global
    frame indices of the previous year's last xfade frames, same camera, same
    flicker, used only inside the dissolve), then its own n frames. Every camera
    pose comes from camera_t_at(g) with g the global frame index."""
    sch = schedule()
    s = sch[year]
    t_build = time.time()
    pl, anim, has_sinking = build_scene(year)
    print("%d: scene built in %.1f s" % (year, time.time() - t_build), file=sys.stderr)
    written = []

    def emit(path, g, raw=False):
        """raw=True also writes the scene without the HUD (<name>.raw.png) for
        the crossfade frames: assemble.py blends the raw scenes and composes the
        HUD once per frame, so the text never double-exposes."""
        t_f = time.time()
        frame = render_global_frame(pl, anim, g)
        if raw:
            Image.fromarray(frame).save(raw_path(path), compress_level=3)
        compose(frame, year, args, has_sinking).save(path, compress_level=3)
        written.append(path)
        print("%s  g=%d az=%.3f  %.2f s" % (os.path.basename(path), g, camera_azimuth_deg(g), time.time() - t_f),
              file=sys.stderr)

    n_x = xfade_frames_for(year) if s["prev"] is not None else 0
    if s["prev"] is not None:
        p = sch[s["prev"]]
        for k in range(n_x):
            emit(os.path.join(args.outdir, "%d_pre_%02d.png" % (year, k)), p["tick0"] + p["n"] - n_x + k, raw=True)
    elif SPECS[SPEC_NAME]["prefix"] > 0:
        # first year after the act 1 prefix: pre-roll over the prefix's last frames
        n_x = xfade_frames_for(year)
        for k in range(n_x):
            emit(os.path.join(args.outdir, "%d_pre_%02d.png" % (year, k)), s["tick0"] - n_x + k, raw=True)
    n_out = xfade_frames_for(YEARS[YEARS.index(year) + 1]) if year != YEARS[-1] else 0
    for i in range(s["n"]):
        emit(os.path.join(args.outdir, "%d_%03d.png" % (year, i)), s["tick0"] + i, raw=i >= s["n"] - n_out)
    pl.close()
    return written


def raw_path(path):
    return path[:-4] + ".raw.png" if path.endswith(".png") else path + ".raw.png"


def render_segment(args):
    """Legacy mode: explicit camera range for one epoch."""
    t_build = time.time()
    pl, anim, has_sinking = build_scene(args.epoch)
    print("scene built in %.1f s" % (time.time() - t_build), file=sys.stderr)
    written = []
    for i in range(args.frames):
        t = args.t0 + (args.t1 - args.t0) * (i / float(args.frames))
        t_f = time.time()
        frame = render_frame(pl, anim, t, args.start_index + i)
        path = os.path.join(args.outdir, "%d_%03d.png" % (args.epoch, args.start_index + i))
        compose(frame, args.epoch, args, has_sinking).save(path, compress_level=3)
        written.append(path)
        print("%s  t=%.4f  %.2f s" % (path, t, time.time() - t_f), file=sys.stderr)
    pl.close()
    return written


def run_all(args):
    """Render every year of the schedule, args.jobs processes at a time."""
    base = [sys.executable, os.path.abspath(__file__), "--outdir", args.outdir, "--spec", SPEC_NAME]
    if args.minimal:
        base.append("--minimal")
    if args.bg:
        base += ["--bg", args.bg]
    pending = list(YEARS)
    running = []
    failed = []
    while pending or running:
        while pending and len(running) < args.jobs:
            y = pending.pop(0)
            log = open(os.path.join(args.outdir, "render_%d.log" % y), "w")
            running.append((y, subprocess.Popen(base + ["--year", str(y)], stdout=log, stderr=subprocess.STDOUT), log))
        time.sleep(0.5)
        still = []
        for y, p, log in running:
            if p.poll() is None:
                still.append((y, p, log))
            else:
                log.close()
                print("%d finished with code %d" % (y, p.returncode), file=sys.stderr)
                if p.returncode != 0:
                    failed.append(y)
        running = still
    return failed


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--year", type=int, help="render one year of the v4 schedule")
    ap.add_argument("--all", action="store_true", help="render every year of the v4 schedule")
    ap.add_argument("--jobs", type=int, default=7, help="parallel processes for --all")
    ap.add_argument("--spec", default="v4", choices=sorted(SPECS), help="frame schedule: v4 (full) or twoact (act 2 of the cut)")
    ap.add_argument("--epoch", type=int, help="legacy segment mode: year to render")
    ap.add_argument("--t0", type=float, default=0.0)
    ap.add_argument("--t1", type=float, default=1.0)
    ap.add_argument("--frames", type=int, default=24)
    ap.add_argument("--caption", type=str, default="")
    ap.add_argument("--outdir", type=str, default=os.path.join(DATA_DIR, "frames_v4"))
    ap.add_argument("--start-index", type=int, default=0, help="first frame number (legacy mode)")
    ap.add_argument("--minimal", action="store_true", help="slide mode: year, legend and exaggeration note only")
    ap.add_argument("--bg", type=str, default="", help="flat background colour for slide mode, e.g. #000000")
    ap.add_argument("--no-canals", action="store_true", help="accepted for compatibility; canals are never drawn in v4")
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)
    global BG_OVERRIDE
    BG_OVERRIDE = args.bg
    set_spec(args.spec)

    if args.all:
        failed = run_all(args)
        print(json.dumps({"failed": failed}))
        sys.exit(1 if failed else 0)
    if args.year is not None:
        written = render_year(args.year, args)
        print(json.dumps({"year": args.year, "frames": written}))
        return
    if args.epoch is None:
        ap.error("one of --year, --all or --epoch is required")
    written = render_segment(args)
    print(json.dumps({"epoch": args.epoch, "frames": written}))


if __name__ == "__main__":
    main()
