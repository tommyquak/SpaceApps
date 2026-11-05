"""Canal detection from PlanetScope SuperDove 8-band surface reflectance.

Same NDWI-threshold OR Sato-ridge structure as analysis/canal_detect_s2.py,
re-tuned for the 3 m pixel size (Sentinel-2 is 10 m). Detector reads any
number of pre-downloaded ortho_analytic_8b_sr scenes (via the helper in
ingest/planetscope.py) and writes a binary canal mask COG on an EPSG:4326
grid.

Why a separate module rather than reusing canal_detect_s2.py:
  * The Sato sigmas scale with canal width measured in pixels. A 5-15 m
    canal is 0.5-1.5 px at 10 m S2 (sub-pixel, why S2 underperforms) but
    1.7-5 px at 3 m PlanetScope. Different sigma stack is needed.
  * The minimum-segment-size in pixels scales as 1/res^2: a 12-px cleanup
    at 10 m is ~1200 m2; the equivalent at 3 m is ~133 px.

The detector itself stays simple (no CNN). Tay et al. (2025) used a CNN on
PlanetScope and got cleaner results; a CNN can be layered on later. Goal
of this module is the 3 m resolution win over S2, not a method change.

Usage:
    python -m peatguard.analysis.canal_detect_planet \\
        --bbox 113.94,-3.12,114.09,-2.99 \\
        --scenes Datasets/planet/tay_sw_aoi/*.tif \\
        --out  deliverables/canal_planet_sw_aoi/canal_mask.tif
"""

from __future__ import annotations

import argparse
import glob
import logging
from pathlib import Path

import numpy as np
import rasterio
from rasterio.crs import CRS
from scipy.ndimage import binary_dilation, binary_erosion, label
from skimage.filters import sato
from skimage.measure import regionprops

from peatguard.ingest.planetscope import composite

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("canal_detect_planet")

# ---- Tunables (3 m PlanetScope) --------------------------------------------
# Canal half-widths in pixels for 5-15 m wide canals at 3 m px:
#   half_width = (5..15) / 2 / 3 = 0.83..2.5 px  -> sigmas span 1..3
# Plus an extra wider sigma for larger drainage canals (~15-25 m wide):
RIDGE_SIGMAS = (1.0, 1.5, 2.5, 4.0)
NDWI_THRESHOLD = 0.0             # standard NDWI water/wet threshold
RIDGE_PERCENTILE = 97.0          # ridge response above this percentile = canal
MIN_SEGMENT_PX = 110             # ~1000 m^2 at 3 m (matches 12 px @ 10 m S2)
# Shape filter: a real canal is elongated (length >> width). Forest clearings
# and bright vegetation patches are roughly round and produce a Sato ridge
# response but fail this test. Aspect ratio is major-axis / minor-axis from
# the fitted ellipse on each connected component. A 3:1 threshold keeps real
# canals (typically 5:1 or higher) while dropping near-circular blobs.
MIN_ASPECT_RATIO = 3.0
# River-exclusion: the Sato ridge filter fires on the bright/dark edges of
# wide water bodies (rivers, lakes), producing spurious "canal" lines that
# trace the banks. Detect wide-and-wet components internally and zero out
# any canal pixel within a buffer around them. The output is still canal-
# only -- river detection is used as a mask, not as an output class.
RIVER_NDWI_THRESHOLD = -0.1     # loosened for tannic tropical rivers
RIVER_NIR_MAX = 1000             # uint16 SR (10% reflectance)
RIVER_WIDE_EROSION_DEPTH = 2     # components surviving N erosions are wide
                                 # (lowered from 3 -> 2 so narrower
                                 # tributaries also classify as rivers and
                                 # their banks are excluded)
RIVER_BUFFER_PX = 25             # exclude canals within ~140 m at 5.5 m
                                 # output (widened from 8 -> 25; the Sato
                                 # ridge response along riverbanks extends
                                 # well past the visible bank, especially
                                 # for the wider Sebangau / Kahayan rivers).
                                 # Conservative bias: better to drop a few
                                 # real canals near rivers than to keep
                                 # un-dammable bank-edge false positives.
# -----------------------------------------------------------------------------


def _build_river_exclusion(ndwi: np.ndarray, nir: np.ndarray | None,
                           finite: np.ndarray) -> tuple[np.ndarray, int]:
    """Identify wide water bodies and return a buffered exclusion mask.

    A "river" pixel here is anything that is (a) wet enough -- by either
    loosened NDWI threshold or low NIR -- and (b) part of a connected
    component wide enough to survive RIVER_WIDE_EROSION_DEPTH erosions.
    The mask is then dilated by RIVER_BUFFER_PX to cover bank edges where
    Sato ridge would otherwise spuriously detect canals.

    Returns the exclusion mask and the count of river components.
    """
    wet_ndwi = finite & (ndwi > RIVER_NDWI_THRESHOLD)
    if nir is not None:
        wet_nir = finite & np.isfinite(nir) & (nir < RIVER_NIR_MAX)
        wet = wet_ndwi | wet_nir
    else:
        wet = wet_ndwi
    if not wet.any():
        return np.zeros_like(finite), 0
    wet_lab, _ = label(wet)
    eroded = binary_erosion(wet, iterations=RIVER_WIDE_EROSION_DEPTH)
    river_labels = np.unique(wet_lab[eroded])
    river_labels = river_labels[river_labels > 0]
    is_river = np.isin(wet_lab, river_labels)
    if RIVER_BUFFER_PX > 0 and is_river.any():
        is_river = binary_dilation(is_river, iterations=RIVER_BUFFER_PX)
    return is_river, int(len(river_labels))


def detect_canals(ndwi: np.ndarray, nir: np.ndarray | None = None
                  ) -> tuple[np.ndarray, dict]:
    """NDWI (+ optional NIR) -> binary canal mask.

    Detection: NDWI threshold OR Sato ridge response.
    Cleanup: minimum-size + aspect-ratio shape filter + river-buffer
    exclusion (the Sato ridge fires on river banks; if NIR is supplied,
    we mask those out internally).
    """
    finite = np.isfinite(ndwi)
    # Sato needs a filled array; fill nans with a dry value so they don't
    # masquerade as bright ridges at the AOI boundary.
    fill = float(np.nanpercentile(ndwi[finite], 20)) if finite.any() else 0.0
    ndwi_filled = np.where(finite, ndwi, fill).astype(np.float32)

    # Tubeness response -- bright elongated structures (canals are wetter
    # so brighter in NDWI; black_ridges=False finds bright lines).
    ridge = sato(ndwi_filled, sigmas=RIDGE_SIGMAS, black_ridges=False)
    ridge = np.nan_to_num(ridge, nan=0.0)
    r_thresh = float(np.percentile(ridge[finite], RIDGE_PERCENTILE)) if finite.any() else 0.0

    wet = finite & (ndwi > NDWI_THRESHOLD)
    ridgey = finite & (ridge > r_thresh)
    canal = wet | ridgey

    # Stage 0 cleanup: zero out canal pixels inside the river-buffer
    # exclusion (riverbank edges that Sato ridge spuriously picks up).
    river_mask, n_rivers = _build_river_exclusion(ndwi, nir, finite)
    river_buffer_dropped = int((canal & river_mask).sum())
    canal &= ~river_mask

    # Stage 1 cleanup: drop isolated specks below the minimum area.
    lab, n = label(canal)
    pre_shape_pixels = int(canal.sum())
    if n:
        sizes = np.bincount(lab.ravel())
        small = np.isin(lab, np.where(sizes < MIN_SEGMENT_PX)[0])
        canal &= ~small

    # Stage 2 cleanup: drop near-circular components (forest clearings,
    # bright vegetation patches). A canal is elongated; aspect_ratio is
    # major-axis / minor-axis from the fitted ellipse via skimage
    # regionprops. Below MIN_ASPECT_RATIO -> not canal-shaped, drop.
    lab2, n2 = label(canal)
    round_dropped = 0
    if n2:
        keep_labels = np.zeros(n2 + 1, dtype=bool)
        for r in regionprops(lab2):
            minor = r.axis_minor_length
            if minor <= 0:
                # degenerate (e.g. single-pixel-wide) -> treat as elongated
                keep_labels[r.label] = True
                continue
            ratio = r.axis_major_length / minor
            if ratio >= MIN_ASPECT_RATIO:
                keep_labels[r.label] = True
            else:
                round_dropped += int(r.area)
        canal = keep_labels[lab2]

    stats = {
        "ndwi_threshold": NDWI_THRESHOLD,
        "ridge_threshold": round(r_thresh, 5),
        "ridge_sigmas": list(RIDGE_SIGMAS),
        "min_segment_px": MIN_SEGMENT_PX,
        "min_aspect_ratio": MIN_ASPECT_RATIO,
        "river_buffer_px": RIVER_BUFFER_PX,
        "river_components_excluded": int(n_rivers),
        "pixels_dropped_river_buffer": river_buffer_dropped,
        "pixels_wet_only": int((wet & ~ridgey & ~river_mask).sum()),
        "pixels_ridge_recovered": int((ridgey & ~wet & ~river_mask).sum()),
        "pixels_pre_shape_filter": pre_shape_pixels,
        "pixels_dropped_round": round_dropped,
        "canal_pixels": int(canal.sum()),
        "canal_pct_of_aoi": round(100.0 * float(canal.mean()), 3),
        "valid_aoi_coverage_pct": round(100.0 * float(finite.mean()), 1),
    }
    return canal.astype(np.uint8), stats


def build_canal_mask(bbox, scene_paths, out_path: Path,
                     res_deg: float = 0.00003) -> dict:
    """Composite PlanetScope scenes over bbox, detect canals, write a COG."""
    comp = composite(scene_paths, bbox, res_deg=res_deg)
    canal, stats = detect_canals(comp["ndwi"], comp.get("nir"))
    h, w = comp["shape"]

    profile = {
        "driver": "COG", "dtype": "uint8", "count": 1,
        "width": w, "height": h, "crs": CRS.from_epsg(4326),
        "transform": comp["transform"], "nodata": 0, "compress": "deflate",
        "overview_resampling": "nearest",
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(canal, 1)
    stats["n_scenes"] = comp["n_scenes"]
    stats["output"] = str(out_path)
    stats["grid"] = f"{w}x{h}"
    stats["pixel_size_deg"] = res_deg
    stats["pixel_size_m_approx"] = round(res_deg * 111320, 2)
    log.info("wrote %s", out_path)
    return stats


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bbox", required=True, help="W,S,E,N (comma-separated)")
    p.add_argument("--scenes", required=True,
                   help="glob of ortho_analytic_8b_sr GeoTIFFs (quote it)")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--res-deg", type=float, default=0.00003,
                   help="output pixel size in deg (default ~3.3 m)")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    bbox = tuple(float(x) for x in args.bbox.split(","))
    paths = sorted(glob.glob(args.scenes))
    if not paths:
        raise SystemExit(f"no scenes match glob: {args.scenes}")
    log.info("scenes:")
    for p in paths:
        log.info("  %s", p)
    stats = build_canal_mask(bbox, paths, args.out, res_deg=args.res_deg)
    log.info("")
    log.info("=== CANAL DETECTION STATS (PlanetScope 3 m) ===")
    for k, v in stats.items():
        log.info("  %s: %s", k, v)


if __name__ == "__main__":
    main()
