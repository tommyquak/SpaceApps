"""Two-class hydro layer from PlanetScope: drainage canals vs rivers/water.

The standard canal detector (canal_detect_planet) only outputs a single
binary canal class. For the SBG-CE deployed view, the user wants to also
see the major rivers and natural water bodies on the same map -- separately
colour-coded so restoration targets (canals) stay visually distinct from
contextual hydrography (rivers).

Classification rule:
    1. Compute NDWI from green + NIR; mark all pixels with NDWI > NDWI_T as
       "wet". This catches every water surface: rivers, lakes, ponds, and
       canals all hit the same threshold.
    2. Apply morphological closing (small radius) to fill speckle gaps in
       genuinely wet features without merging unrelated ones.
    3. Connected-component label the wet mask.
    4. For each component, measure "thickness" by counting how many
       iterations of binary erosion it survives. Components surviving
       >= WIDE_EROSION_DEPTH iterations have a core wider than
       2 * WIDE_EROSION_DEPTH pixels (in any direction) -- these are
       rivers / large water bodies.
    5. Components surviving fewer iterations are narrow -- candidate
       canals. The same aspect-ratio shape filter as canal_detect_planet
       drops near-circular blobs (forest gaps mis-firing on the ridge).

Output: a uint8 COG with values:
    0 = nothing / dry
    1 = canal (narrow linear, restoration target)
    2 = river / large water body (context only)

Usage:
    python -m peatguard.analysis.canal_and_water_planet \\
        --bbox 114.21,-3.0,114.83,-2.42 \\
        --scenes 'Datasets/planet/big_sebangau/*_8b_sr.tif' \\
        --out deliverables/canal_planet_big_aoi/canal_and_water.tif
"""

from __future__ import annotations

import argparse
import glob
import logging
from pathlib import Path

import numpy as np
import rasterio
from rasterio.crs import CRS
from scipy.ndimage import binary_closing, binary_erosion, label
from skimage.filters import sato
from skimage.measure import regionprops

from peatguard.ingest.planetscope import composite

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("canal_and_water_planet")

# ---- Tunables --------------------------------------------------------------
# Wet detection: standard NDWI fails on tropical peat rivers because their
# tannin / sediment load suppresses the green-NIR contrast. Diagnosed on
# the Sebangau system (2026-05-19): NDWI > 0 catches 0.14% of AOI;
# NDWI > -0.1 OR NIR < 1500 catches 2.6% and traces the actual river
# system. NIR is the discriminating band -- tannic water still absorbs
# NIR strongly even when its visible-band signature is muted.
NDWI_THRESHOLD = -0.1          # loosened for tannic tropical rivers
NIR_WATER_MAX = 1000           # uint16 SR scale (10% reflectance);
                               # tightened from 1500 (which over-caught
                               # burnt-peat and wet plantation patches).
                               # NIR < 1000 isolates the actual river
                               # water surface in tropical tannin systems.
# Ridge detection: narrow dark linear features (dry-bed drainage canals
# in plantations) -- not water, but elongated and locally bright in NDWI.
RIDGE_SIGMAS = (1.0, 1.5, 2.5, 4.0)
RIDGE_PERCENTILE = 97.0
CLOSING_RADIUS = 1             # fill 1-pixel gaps inside features
# Erosion-survival threshold: a component that survives N erosions has a
# minimum local thickness of about 2*N+1 pixels. At ~5.5 m output:
#   N = 3 -> width >= ~35 m -> rivers and large drainage channels
#   N = 2 -> width >= ~25 m
# Sebangau river is ~80-200 m wide -> safely caught; Kapuas tributaries 50+ m.
WIDE_EROSION_DEPTH = 3
# Canal aspect-ratio (same as canal_detect_planet) -- drop near-circular
# components that are not rivers either (forest gaps, ponds).
CANAL_MIN_ASPECT = 3.0
CANAL_MIN_PX = 60              # ~1800 m^2 at 5.5 m output
# -----------------------------------------------------------------------------


def classify(ndwi: np.ndarray, nir: np.ndarray) -> tuple[np.ndarray, dict]:
    """NDWI + NIR -> 2-class hydro raster.

    Args:
        ndwi: NDWI raster (green-NIR/green+NIR), nan over invalid.
        nir: NIR band reflectance (uint16 SR scale, i.e. reflectance * 10000).

    Returns:
        out: uint8 with 0=dry, 1=canal, 2=river.
        stats: dict with pixel counts + percentages.
    """
    finite = np.isfinite(ndwi)

    # 1. Wet pixels: loosened NDWI threshold OR low NIR. The OR catches
    #    tropical tannin rivers that fail NDWI > 0 but still absorb NIR.
    wet_ndwi = finite & (ndwi > NDWI_THRESHOLD)
    wet_nir = finite & np.isfinite(nir) & (nir < NIR_WATER_MAX)
    wet = wet_ndwi | wet_nir
    if CLOSING_RADIUS > 0:
        wet_closed = binary_closing(wet, iterations=CLOSING_RADIUS)
    else:
        wet_closed = wet
    log.info("wet pixels: NDWI=%d  NIR=%d  union=%d",
             int(wet_ndwi.sum()), int(wet_nir.sum()), int(wet.sum()))

    # 2. Sato ridge response catches narrow dark linear features (dry
    #    drainage canals in plantations, which appear as darker lines but
    #    not necessarily as water surfaces).
    fill = float(np.nanpercentile(ndwi[finite], 20)) if finite.any() else 0.0
    ndwi_filled = np.where(finite, ndwi, fill).astype(np.float32)
    ridge = sato(ndwi_filled, sigmas=RIDGE_SIGMAS, black_ridges=False)
    ridge = np.nan_to_num(ridge, nan=0.0)
    r_thresh = float(np.percentile(ridge[finite], RIDGE_PERCENTILE)) if finite.any() else 0.0
    ridgey = finite & (ridge > r_thresh)

    # 3. RIVER classification operates on the WET mask only (real surface
    #    water), not on the ridge response. A real river is wide AND wet.
    #    Sato ridge thickens narrow canals to ~3-5 px in the response;
    #    using it for river classification would mis-promote canals to
    #    rivers. NDWI-wet is the physical criterion for surface water.
    wet_lab, n_wet = label(wet_closed)
    eroded = binary_erosion(wet_closed, iterations=WIDE_EROSION_DEPTH)
    river_labels = np.unique(wet_lab[eroded])
    river_labels = river_labels[river_labels > 0]
    is_river = np.isin(wet_lab, river_labels)
    log.info("wet components: %d  ->  rivers (wide+wet): %d", n_wet, len(river_labels))

    # 4. CANAL candidates = everything we detected (wet OR ridge) that is
    #    NOT classified as a river. Includes narrow wet canals and dry
    #    plantation drainage lines.
    hydro = wet_closed | ridgey
    is_canal_candidate = hydro & ~is_river

    # Shape-filter the canal class: drop tiny specks and near-circular blobs.
    lab_c, n_c = label(is_canal_candidate)
    keep = np.zeros(n_c + 1, dtype=bool)
    dropped_small = 0
    dropped_round = 0
    for r in regionprops(lab_c):
        if r.area < CANAL_MIN_PX:
            dropped_small += r.area
            continue
        minor = r.axis_minor_length
        if minor <= 0:
            keep[r.label] = True
            continue
        if r.axis_major_length / minor >= CANAL_MIN_ASPECT:
            keep[r.label] = True
        else:
            dropped_round += r.area
    is_canal = keep[lab_c]

    out = np.zeros(ndwi.shape, dtype=np.uint8)
    out[is_canal] = 1
    out[is_river] = 2

    stats = {
        "ndwi_threshold": NDWI_THRESHOLD,
        "nir_water_max": NIR_WATER_MAX,
        "ridge_sigmas": list(RIDGE_SIGMAS),
        "ridge_percentile": RIDGE_PERCENTILE,
        "ridge_threshold_value": round(r_thresh, 5),
        "closing_radius": CLOSING_RADIUS,
        "wide_erosion_depth": WIDE_EROSION_DEPTH,
        "canal_min_aspect": CANAL_MIN_ASPECT,
        "canal_min_px": CANAL_MIN_PX,
        "wet_pixels": int(wet.sum()),
        "ridge_pixels": int(ridgey.sum()),
        "wet_components_total": int(n_wet),
        "river_components": int(len(river_labels)),
        "canal_pixels": int(is_canal.sum()),
        "river_pixels": int(is_river.sum()),
        "dropped_small_canal_px": int(dropped_small),
        "dropped_round_canal_px": int(dropped_round),
        "canal_pct_of_aoi": round(100.0 * float(is_canal.mean()), 3),
        "river_pct_of_aoi": round(100.0 * float(is_river.mean()), 3),
        "valid_aoi_coverage_pct": round(100.0 * float(finite.mean()), 1),
    }
    return out, stats


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bbox", required=True, help="W,S,E,N (comma-separated)")
    p.add_argument("--scenes", required=True,
                   help="glob of ortho_analytic_8b_sr GeoTIFFs (quote it)")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--res-deg", type=float, default=0.00005,
                   help="output pixel size in deg (default ~5.5 m)")
    args = p.parse_args()

    bbox = tuple(float(x) for x in args.bbox.split(","))
    paths = sorted(glob.glob(args.scenes))
    if not paths:
        raise SystemExit(f"no scenes match: {args.scenes}")
    log.info("scenes (%d):", len(paths))
    for pth in paths:
        log.info("  %s", pth)

    comp = composite(paths, bbox, res_deg=args.res_deg)
    out, stats = classify(comp["ndwi"], comp["nir"])
    h, w = comp["shape"]

    profile = {
        "driver": "COG", "dtype": "uint8", "count": 1,
        "width": w, "height": h, "crs": CRS.from_epsg(4326),
        "transform": comp["transform"], "nodata": 0, "compress": "deflate",
        "overview_resampling": "nearest",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(args.out, "w", **profile) as dst:
        dst.write(out, 1)
    stats["n_scenes"] = comp["n_scenes"]
    stats["output"] = str(args.out)
    stats["grid"] = f"{w}x{h}"
    stats["pixel_size_m_approx"] = round(args.res_deg * 111320, 2)
    log.info("wrote %s", args.out)
    log.info("")
    log.info("=== HYDRO CLASSIFICATION STATS (PlanetScope) ===")
    for k, v in stats.items():
        log.info("  %s: %s", k, v)


if __name__ == "__main__":
    main()
