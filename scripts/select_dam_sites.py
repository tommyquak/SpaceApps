"""Automated canal-block dam site selection.

Given the canal_mask (binary canal detection in UTM 32750) and canal_risk
(combined subsidence + proximity risk in EPSG:4326), pick a set of dam
sites that:

    1. Lie on detected canal centerlines (skeleton of canal_mask)
    2. Are spaced at least MIN_SPACING_M apart (typical BRG / Hooijer
       guidance is one earthen dam every 100-200 m on flat peatdomes;
       wider spacing on lower-slope canals)
    3. Score highly on canal_risk in a small buffer

Output: GeoJSON FeatureCollection of Point geometries with properties
{site_id, risk, lat, lon}.  Also emits a Python literal list ready for
pasting into a dam-site seed list.

Usage:
    python3 scripts/select_dam_sites.py \
        --mask data/sample/canal_mask.tif \
        --risk data/sample/canal_risk.tif \
        --out  /tmp/dam_sites.geojson
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.transform import xy as rio_xy
from scipy.ndimage import distance_transform_edt
from skimage.morphology import skeletonize

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("dam_sites")


# ---- Tunables ---------------------------------------------------------------
MIN_SPACING_M = 250    # min distance between dams (centerline-walk + final NMS)
SAMPLE_STEP_M = 50     # candidate centerline sample interval
RISK_BUFFER_M = 60     # buffer around each candidate to average canal_risk over
RISK_FLOOR = 0.55      # only consider candidates whose buffered risk >= this
TARGET_N = 40          # cap on number of dam sites returned

# Width filter: at each centerline candidate, look up the half-width to the
# nearest non-canal pixel (= the distance transform on canal_mask).  Canals
# wider than this are open water / rivers / large clearings -- not actual
# drainage canals you can earthen-dam.  Drainage canals targeted by BRG /
# Hooijer guidance are typically 2-15 m wide; we keep the upper bound at
# ~30 m total width (15 m half-width) which still excludes the Sebangau
# river and the wide canals around plantation clearings.
MAX_DAM_HALF_WIDTH_M = 50.0
# -----------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--mask", required=True, type=Path,
                   help="canal_mask COG (binary)")
    p.add_argument("--risk", required=True, type=Path,
                   help="restoration-priority COG (EPSG:4326, 0-1 float)")
    p.add_argument("--out", required=True, type=Path,
                   help="output GeoJSON path")
    p.add_argument("--target", type=int, default=TARGET_N,
                   help="number of dam sites to return")
    p.add_argument("--bbox", default=None,
                   help="optional W,S,E,N -- restrict sites to this cell")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    # ---- Load canal_mask + skeletonize ------------------------------------
    with rasterio.open(args.mask) as ms:
        mask = (ms.read(1) > 0).astype(np.uint8)
        mask_transform = ms.transform
        mask_crs = ms.crs
        # Pixel size in METRES. For a projected (UTM) mask transform.a is
        # already metres; for a geographic (EPSG:4326) mask it is degrees,
        # so convert with the latitude-corrected degree length.
        if mask_crs.is_geographic:
            import math
            clat = mask_transform.f + (mask.shape[0] / 2) * mask_transform.e
            px_m = abs(mask_transform.a) * 111_320.0 * math.cos(math.radians(clat))
        else:
            px_m = abs(mask_transform.a)
    log.info("canal_mask: %s shape=%s pixel=%.1fm CRS=%s canal-pixels=%d",
             args.mask.name, mask.shape, px_m, mask_crs, int(mask.sum()))

    skel = skeletonize(mask.astype(bool)).astype(np.uint8)
    log.info("skeleton:    centerline-pixels=%d", int(skel.sum()))

    # ---- Distance transform: drop centerline pixels in wide water -------
    # `dist[r,c]` = pixels to the nearest non-canal cell.  At centerlines
    # this approximates the canal half-width.  Anything > MAX_DAM_HALF_WIDTH
    # is open water (Sebangau river, lakes), not a dammable drainage canal.
    dist_px = distance_transform_edt(mask)
    half_width_max_px = MAX_DAM_HALF_WIDTH_M / px_m
    narrow_skel = skel & (dist_px <= half_width_max_px)
    log.info("narrow:      half-width <= %.1fm (%.2fpx) -> %d centerline pixels"
             " (%.0f%% of skeleton)",
             MAX_DAM_HALF_WIDTH_M, half_width_max_px,
             int(narrow_skel.sum()),
             100.0 * narrow_skel.sum() / max(1, skel.sum()))

    # ---- Sample candidates along narrow centerlines -----------------------
    rows, cols = np.where(narrow_skel > 0)
    step_px = max(1, int(round(SAMPLE_STEP_M / px_m)))
    rows = rows[::step_px]
    cols = cols[::step_px]
    log.info("sampling:    every %d px (~%dm) -> %d candidates",
             step_px, step_px * px_m, len(rows))

    # ---- Project candidate UTM pixels to lat/lon -------------------------
    utm_xs, utm_ys = rio_xy(mask_transform, rows.tolist(), cols.tolist())
    utm_xs = np.asarray(utm_xs)
    utm_ys = np.asarray(utm_ys)
    to_4326 = Transformer.from_crs(mask_crs, "EPSG:4326", always_xy=True)
    lons, lats = to_4326.transform(utm_xs, utm_ys)

    # ---- Sample canal_risk via buffered mean -----------------------------
    with rasterio.open(args.risk) as rs:
        risk = rs.read(1)
        risk_transform = rs.transform
        risk_nodata = rs.nodata if rs.nodata is not None else -9999
        # canal_risk pixel size in degrees
        deg_per_px_lat = abs(risk_transform.e)
        deg_per_px_lon = abs(risk_transform.a)

    # Convert RISK_BUFFER_M to a degree-radius (rough; latitude ~constant)
    deg_per_m_lat = 1.0 / 111_320
    buf_deg = RISK_BUFFER_M * deg_per_m_lat
    buf_rows = max(1, int(round(buf_deg / deg_per_px_lat)))
    buf_cols = max(1, int(round(buf_deg / deg_per_px_lon)))

    def buffered_risk(lat: float, lon: float) -> float:
        # convert lat/lon -> risk pixel
        col = int((lon - risk_transform.c) / risk_transform.a)
        row = int((lat - risk_transform.f) / risk_transform.e)
        r0 = max(0, row - buf_rows)
        r1 = min(risk.shape[0], row + buf_rows + 1)
        c0 = max(0, col - buf_cols)
        c1 = min(risk.shape[1], col + buf_cols + 1)
        if r1 <= r0 or c1 <= c0:
            return -1.0
        win = risk[r0:r1, c0:c1]
        valid = (win != risk_nodata) & np.isfinite(win)
        if not valid.any():
            return -1.0
        return float(win[valid].mean())

    # Optional bbox restriction (W,S,E,N) -- used to select dams per AOI cell.
    cell_bbox = None
    if args.bbox:
        cell_bbox = tuple(float(x) for x in args.bbox.split(","))

    candidates = []
    for r, c, lat, lon in zip(rows, cols, lats, lons):
        if cell_bbox is not None:
            w, s, e, n = cell_bbox
            if not (w <= lon <= e and s <= lat <= n):
                continue
        rscore = buffered_risk(lat, lon)
        if rscore < RISK_FLOOR:
            continue
        candidates.append((float(lat), float(lon), rscore, int(r), int(c)))
    log.info("candidates:  %d above risk floor %.2f%s",
             len(candidates), RISK_FLOOR,
             f" within bbox {cell_bbox}" if cell_bbox else "")

    # ---- Greedy NMS by min spacing in metres ------------------------------
    # Sort desc by score; walk down the list keeping each candidate that's
    # >= MIN_SPACING_M from every already-kept point.
    candidates.sort(key=lambda x: x[2], reverse=True)
    kept = []
    for cand in candidates:
        lat, lon, score, _r, _c = cand
        ok = True
        for klat, klon, _ks, _kr, _kc in kept:
            # great-circle approximation good enough at this scale
            dlat_m = (lat - klat) * 111_320
            dlon_m = (lon - klon) * 111_320 * np.cos(np.deg2rad(lat))
            if (dlat_m * dlat_m + dlon_m * dlon_m) ** 0.5 < MIN_SPACING_M:
                ok = False
                break
        if ok:
            kept.append(cand)
        if len(kept) >= args.target:
            break
    log.info("dam sites:   %d (after NMS, min spacing %dm)",
             len(kept), MIN_SPACING_M)

    # ---- Emit GeoJSON ----------------------------------------------------
    features = []
    for i, (lat, lon, score, _r, _c) in enumerate(kept, 1):
        features.append({
            "type": "Feature",
            "geometry": {"type": "Point", "coordinates": [lon, lat]},
            "properties": {
                "site_id": f"D-{i:02d}",
                "risk": round(score, 3),
                "lat": lat,
                "lon": lon,
            },
        })
    gj = {"type": "FeatureCollection", "features": features}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(gj, indent=2))
    log.info("wrote %s (%d features)", args.out, len(features))

    # ---- Emit Python literal for seed.py ---------------------------------
    print("\n# ----- dam-site seed list -----")
    print("DAM_SITES = [")
    for f in features:
        p = f["properties"]
        print(f"    ({p['site_id']!r}, {p['lat']:.5f}, {p['lon']:.5f}, {p['risk']:.3f}),")
    print("]")


if __name__ == "__main__":
    main()
