"""Graph-based canal-block dam siting (v2) -- "GraphChainage".

Consumes the v2 canal centerline GRAPH (canal_detect_planet_v2 ->
canal_centerlines.geojson) plus its paired canal_mask.tif, instead of
skeletonizing a raster like the original select_dam_sites.py. Every step now
uses real reach topology + attributes:

  - reaches loaded as a vertex graph (UTM 50S); confluences = degree>=3 nodes,
    outlets = degree-1 tips.
  - dammability from width (recomputed from the canal mask as 2*EDT, because
    the v2 geojson mean_width_m can be deflated by skeleton drift), sinuosity,
    and reach length -- all as SOFT ramps, not hard gates.
  - candidates placed at confluences (highest hydrological leverage -- one
    block backs water up several tributaries), at outlets (built first so
    works don't drain), and filled along long reaches by chainage.
  - score = restoration-priority x Tay-subsidence x dammability x leverage,
    with restoration_zones (Dadap-lineage canal/proximity/fire) and Tay 2025
    L-band subsidence sampled INDEPENDENTLY and never fused (user mandate).
  - spacing enforced as GRAPH-GEODESIC distance (~500 m, clamp 300-700 m),
    not Euclidean -- tropical peat dams influence ~400-600 m upstream
    (Urzainki/Cobb 2023) and the Sumatra 4,800 ha trial used ~500 m intervals
    / ~2 dams/km and halved subsidence (Sci Reports 2024). The old 250 m
    Euclidean NMS was temperate-bog logic and over-densified.

Citations: Ritzema et al. 2014 (canal blocking, build downstream-first);
Jaenicke et al. 2010 (canal-slope / head-drop method); Hooijer et al. 2012
(drainage, subsidence, and carbon-process interpretation, not the project
weights); Pavelsky &
Smith 2008 (channel width from distance transform); Sumatra 2024 trial
(Sci Reports 10.1038/s41598-024-60462-3). Tay et al. 2025 (L-band InSAR
subsidence) is kept methodologically separate from Dadap et al. 2021 /
PlanetScope optical canal geometry.

Usage:
    python3 select_dam_sites_v2.py \
        --centerlines deliverables/canal_planet_v2_sbgce/canal_centerlines.geojson \
        --mask        deliverables/canal_planet_v2_sbgce/canal_mask.tif \
        --restoration deliverables/tay_restoration/restoration_zones.tif \
        --subsidence  deliverables/tay_restoration/subsidence_tay.tif \
        --bbox 114.6108,-2.7441,114.7692,-2.5859 --aoi SBG-CE \
        --target 15 --spacing-m 500 --out /tmp/dam_sites_sbgce.geojson
"""

from __future__ import annotations

import argparse
import json
import logging
import math
from pathlib import Path

import numpy as np
import networkx as nx
import rasterio
from pyproj import Transformer
from scipy.ndimage import distance_transform_edt
from skimage.morphology import skeletonize
from shapely.geometry import shape

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("dam_sites_v2")

# ---- Tunables ---------------------------------------------------------------
TARGET_CRS = "EPSG:32750"          # UTM 50S -- correct for ~114.7 E (NOT 49S)
SNAP_M = 3.0                       # endpoint/vertex snap grid for graph nodes
# Dammability soft ramps (full channel width, metres)
W_FULL_LO, W_FULL_HI = 3.0, 12.0   # full weight band (village-buildable)
W_TAPER_HI = 15.0                  # taper end of village-dam range
# Wide collector canals (15-30 m) stay SELECTABLE at reduced weight -- they
# carry real drainage and are surfaced as dam sites with an "engineered dam"
# note in the per-site remark (not a separate tier). Beyond this they are
# rivers (dammability 0). Lets dense wide-canal grids (e.g. the NE plantation
# block) get sites instead of being silently excluded.
W_ENGINEERED_MAX = 30.0
SIN_OK, SIN_ZERO = 1.4, 1.8        # sinuosity ramp (canals straight ~1.0)
LEN_ZERO, LEN_FULL = 60.0, 150.0   # reach-length ramp (median reach ~80 m)
# Scoring
REST_FLOOR_W = 0.45                # score = (0.45 + 0.55*R_prox) * ...
SUB_FLOOR_W = 0.30                 #         (0.30 + 0.70*R_sub) * ...
SUB_NORM_MM = 50.0                 # subsidence normalisation (mm/yr)
# Audit fix: restoration_zones already contains the project risk composite
# 0.45*proximity + 0.55*clip(|v|/40) (risk_score.py), so multiplying it with
# the standalone subsidence factor double-counted Tay where measured. We
# algebraically recover the proximity share before scoring:
#   r_prox = clip((r_rest - 0.55*clip(|v_pt|/40,0,1)) / 0.45, 0, 1)
# where v_pt is the point-sampled velocity; where unmeasured the composite is
# proximity-backed by construction, so r_prox = r_rest. Subsidence then
# enters the score exactly once (the upstream-integrated factor).
COMPOSITE_SUB_W = 0.55
COMPOSITE_PROX_W = 0.45
COMPOSITE_SEVERE_MM = 40.0
REST_BUFFER_M = 45.0               # restoration sampling window
SUB_BUFFER_M = 200.0               # subsidence window (~157 m Tay pixel)
# Audit fix (2026-06-13): was 1.6 -- the "one dam backs up several branches"
# story does not hold (a dam is one wall across ONE channel; a junction needs
# >=2 structures to impound -- Ritzema et al. 2014 block each canal), and 3 m
# skeleton junctions include artifacts. Kept as a mild candidate-richness
# bonus only; flow accumulation now carries the leverage idea.
LEV_CONFLUENCE = 1.2
LEV_OUTLET = 1.3                   # TRUE outlets only (tip near a river)
LEV_CANAL_END = 1.1                # mid-peat dead end: mild bonus
LEV_EDGE_ARTIFACT = 0.9            # canal truncated by the AOI bbox: penalty
TIP_RIVER_NEAR_M = 150.0           # tip within this of an OSM river = outlet
TIP_EDGE_NEAR_M = 120.0            # tip within this of the bbox = artifact
# Flow-accumulation leverage: how much canal network drains through the
# candidate (water routed toward true outlets along the graph). Aligned with
# placement-by-hydrological-impact (Urzainki et al., canal-blocking
# optimisation); greedy selection by impact is near-optimal there.
ACC_FLOOR = 0.5                    # score *= (0.5 + 0.5*acc_norm)
ACC_NORM_M = 2500.0                # upstream canal length that saturates acc
UPSTREAM_SAMPLE_M = 500.0          # integrate subsidence this far upstream
# El Nino 2019 max-drawdown factor (Tay et al. 2025 episodic motion):
# peat that draws down hardest under drought gains the most from rewetting.
# Mild multiplier because coverage is sparse (~2% of pixels measured).
EN_FLOOR = 0.7                     # score *= (0.7 + 0.3*r_en)
EN_NORM_CM = 8.0                   # drawdown (cm) that saturates the factor
EN_BUFFER_M = 200.0
# Quality floor: better to return fewer good sites than pad to target.
# Relative to the best candidate (the absolute scale shifts whenever a new
# multiplicative factor joins the score), with a small absolute backstop.
MIN_SCORE_FRACTION = 0.35
MIN_SITE_SCORE_ABS = 0.08
# Navigation-consultation proxy: a dammable point within this graph distance
# of a river outlet sits on a corridor villagers may use for boat transport.
CONSULT_RIVER_M = 1500.0
# Access/cost factor (audit feature B): sites a crew can actually reach with
# materials score higher. Two honest proxies from real data -- straight-line
# distance to the nearest OSM settlement (crew origin) and to navigable water
# (OSM river; boats carry sandbags/timber). A soft multiplier, never a gate:
#   access = ACC_MIN + (1-ACC_MIN) * 0.5*(reach_settlement + reach_river)
# 0.72 floor: a maximally-remote site keeps 72% of score -- a clear demotion
# (feature B intent) without gating valid high-subsidence sites out. Tuned
# 2026-06-13 after a 0.6 floor over-compressed scores and the relative quality
# floor then cut too many remote-but-valid sites (49->34); 0.72 -> ~44.
ACCESS_MIN = 0.72                  # remotest sites keep 72% of score
SETTLEMENT_FULL_M = 1500.0         # <=this from a settlement = full marks
SETTLEMENT_ZERO_M = 8000.0         # >=this = no settlement credit
RIVER_FULL_M = 800.0               # <=this from navigable water = full marks
RIVER_ZERO_M = 6000.0
# Consultation upgrade (audit feature C): flag a candidate when it lies on a
# plausible village<->river boat route -- i.e. its graph distance to a
# settlement-snapped node PLUS to a river outlet is close to the direct
# settlement->river graph distance (the reach is "on the way").
CONSULT_DETOUR_FRAC = 1.35         # on-route if d_via <= 1.35 * d_direct
CONSULT_SETTLEMENT_SNAP_M = 400.0  # settlement must be within this of a canal node
# Spacing
SPACING_M = 500.0                  # default; clamped to [300, 700]
SPACING_MIN, SPACING_MAX = 300.0, 700.0
# Width calibration (audit fix B-2, 2026-06-13): our detector's mask paints
# canals as thicker bands than they are, so 2*EDT over-reads channel width by
# ~1.6x. Validated against the independent Dadap et al. 2021 canal map across
# all four AOIs (width ratio ours/Dadap: CE 1.41, NE 1.41, NW 1.58, S 2.00;
# mean 1.60, std 0.24). Divide the EDT width by this factor so the dammability
# gate and dam-type classes operate on calibrated widths -- recovering genuine
# ~10-13 m canals that were wrongly rejected as >15 m. [Dadap et al. 2021]
WIDTH_CALIBRATION = 1.6
# -----------------------------------------------------------------------------


def _ramp(x, lo, hi):
    """0 below lo, 1 above hi, linear between."""
    if hi <= lo:
        return 1.0 if x >= hi else 0.0
    return float(np.clip((x - lo) / (hi - lo), 0.0, 1.0))


def dammability(width_m, sinuosity, length_m):
    # width weight: full across [LO,HI]; taper 1.0->0.55 across [HI,TAPER_HI]
    # (still village-buildable); reduced 0.55->0.15 across the engineered band
    # [TAPER_HI, ENGINEERED_MAX] (wide collectors -- selectable but lower
    # priority, flagged "engineered" in the remark); 0 below LO / above MAX.
    if width_m < W_FULL_LO:
        w = 0.0
    elif width_m <= W_FULL_HI:
        w = 1.0
    elif width_m <= W_TAPER_HI:
        w = 1.0 - 0.45 * _ramp(width_m, W_FULL_HI, W_TAPER_HI)
    elif width_m <= W_ENGINEERED_MAX:
        w = 0.55 - 0.40 * _ramp(width_m, W_TAPER_HI, W_ENGINEERED_MAX)
    else:
        w = 0.0
    s = 1.0 - _ramp(sinuosity, SIN_OK, SIN_ZERO)
    l = _ramp(length_m, LEN_ZERO, LEN_FULL)
    return w * s * l


def dam_type(width_m):
    if width_m <= 8.0:
        return "compacted-peat"
    if width_m <= 15.0:
        return "composite"
    return "engineered"


def reaches_from_mask(mask_path, bbox, to_utm):
    """FALLBACK for AOIs without a v2 centerline geojson (SBG-NW/NE/S use the
    deployed 5.6 m PlanetScope mask): skeletonize the mask, build a reach graph
    with sknw, emit the same GeoJSON-feature shape build_graph expects
    (geometry.coordinates in lon/lat + properties length_m, sinuosity)."""
    import sknw
    with rasterio.open(mask_path) as ms:
        mask = ms.read(1) > 0
        tr = ms.transform
    skel = skeletonize(mask)
    g = sknw.build_sknw(skel, multi=False)
    feats = []
    for s, e in g.edges():
        pts = g.get_edge_data(s, e)["pts"]  # (row, col)
        if len(pts) < 2:
            continue
        lons, lats = rasterio.transform.xy(tr, pts[:, 0], pts[:, 1])
        coords = list(zip(np.atleast_1d(lons), np.atleast_1d(lats)))
        # length + sinuosity in UTM metres
        xy = np.array([to_utm.transform(lo, la) for lo, la in coords])
        seglen = np.hypot(*np.diff(xy, axis=0).T).sum()
        d_end = np.hypot(*(xy[-1] - xy[0]))
        sinu = float(seglen / d_end) if d_end > 1e-6 else 1.0
        # bbox filter (drop reaches whose midpoint is outside the AOI cell)
        if bbox:
            w, s_, e_, n_ = bbox
            mlon, mlat = coords[len(coords) // 2]
            if not (w <= mlon <= e_ and s_ <= mlat <= n_):
                continue
        feats.append({"type": "Feature",
                      "geometry": {"type": "LineString",
                                   "coordinates": [[float(x), float(y)] for x, y in coords]},
                      "properties": {"length_m": float(seglen), "sinuosity": sinu}})
    log.info("mask fallback: skeletonized %s -> %d reaches", mask_path.name, len(feats))
    return feats


def build_graph(feats, mask_path, to_utm):
    """Vertex graph in UTM from a list of reach features; each reach's vertices
    are nodes, consecutive vertices are edges weighted by segment length.
    Per-reach width recomputed from the canal mask (2*EDT). Returns graph G
    with node attrs (lon, lat, x, y) and edge attrs (width_m, sinuosity,
    reach_len_m)."""

    # Recompute per-reach width from the mask: 2 * EDT(mask) at reach pixels.
    with rasterio.open(mask_path) as ms:
        mask = ms.read(1) > 0
        mtr = ms.transform
        clat = mtr.f + (mask.shape[0] / 2) * mtr.e
        px_m = abs(mtr.a) * 111320.0 * math.cos(math.radians(clat))
        # 2*EDT over-reads thick-painted canals ~1.6x vs Dadap ground truth;
        # calibrate so dammability + dam-type see true channel width.
        fullwidth = (2.0 / WIDTH_CALIBRATION) * distance_transform_edt(mask).astype(np.float32) * px_m

    def lonlat_to_px(lon, lat):
        col = int((lon - mtr.c) / mtr.a)
        row = int((lat - mtr.f) / mtr.e)
        return row, col

    def snap(x, y):
        return (round(x / SNAP_M), round(y / SNAP_M))

    G = nx.Graph()
    nonzero_w = 0
    for rid, f in enumerate(feats):
        coords = f["geometry"]["coordinates"]
        if len(coords) < 2:
            continue
        props = f["properties"]
        length_m = props.get("length_m", 0.0)
        sinu = props.get("sinuosity", 1.0)
        # recompute width: 75th pct of fullwidth along the reach pixels
        ws = []
        for lon, lat in coords:
            r, c = lonlat_to_px(lon, lat)
            if 0 <= r < mask.shape[0] and 0 <= c < mask.shape[1]:
                ws.append(fullwidth[r, c])
        width_m = float(np.nanpercentile(ws, 75)) if ws else 0.0
        if width_m > 0:
            nonzero_w += 1
        # add vertices as nodes, consecutive as edges
        prev = None
        for lon, lat in coords:
            x, y = to_utm.transform(lon, lat)
            key = snap(x, y)
            if key not in G:
                G.add_node(key, lon=lon, lat=lat, x=x, y=y)
            if prev is not None and prev != key:
                seg = math.hypot(G.nodes[key]["x"] - G.nodes[prev]["x"],
                                 G.nodes[key]["y"] - G.nodes[prev]["y"])
                G.add_edge(prev, key, weight=seg, reach_id=rid,
                           width_m=width_m, sinuosity=sinu, reach_len_m=length_m)
            prev = key
    log.info("graph: %d nodes, %d edges; %d/%d reaches got width>0 (recomputed from mask)",
             G.number_of_nodes(), G.number_of_edges(), nonzero_w, len(feats))
    return G


def classify_tips(G, bbox, osm_path):
    """Split degree-1 tips into river outlets / bbox artifacts / canal ends.

    The old selector bonused EVERY dead end as an "outlet"; most tips are
    either detection truncations at the AOI boundary (the canal continues
    outside) or canals petering out mid-peat -- neither is the
    dam-downstream-first outlet of Ritzema et al. (2014)."""
    tips = [n for n in G.nodes() if G.degree(n) == 1]
    edge_tips, river_tips = set(), set()

    if bbox:
        w, s, e, nn = bbox
        for n in tips:
            lon, lat = G.nodes[n]["lon"], G.nodes[n]["lat"]
            coslat = math.cos(math.radians(lat))
            dx = min(lon - w, e - lon) * 111320.0 * coslat
            dy = min(lat - s, nn - lat) * 111320.0
            if min(dx, dy) < TIP_EDGE_NEAR_M:
                edge_tips.add(n)

    if osm_path and Path(osm_path).exists():
        from shapely.geometry import Point, shape
        from shapely.strtree import STRtree
        fc = json.load(open(osm_path))
        geoms = []
        for f in fc.get("features", []):
            if f.get("properties", {}).get("waterway") not in ("river", "stream"):
                continue
            try:
                geoms.append(shape(f["geometry"]))
            except Exception:  # noqa: BLE001
                continue
        if geoms:
            tree = STRtree(geoms)
            for n in tips:
                if n in edge_tips:
                    continue
                lat = G.nodes[n]["lat"]
                deg = TIP_RIVER_NEAR_M / (111320.0 * math.cos(math.radians(lat)))
                pt = Point(G.nodes[n]["lon"], lat)
                try:  # shapely 2.x
                    near = len(tree.query_nearest(pt, max_distance=deg)) > 0
                except (AttributeError, TypeError):  # shapely 1.8
                    g = tree.nearest(pt)
                    near = g is not None and pt.distance(g) <= deg
                if near:
                    river_tips.add(n)
    canal_ends = set(tips) - edge_tips - river_tips
    log.info("tips: %d river outlets, %d bbox artifacts, %d canal ends",
             len(river_tips), len(edge_tips), len(canal_ends))
    return river_tips, edge_tips, canal_ends


def flow_accumulation(G, roots):
    """Upstream drained canal length per node, water routed along the graph
    toward the nearest root (true outlet). Each edge's length is charged to
    its upstream endpoint and routed downstream D8-style (single downstream
    neighbour), so reconverging paths are not double-counted."""
    if not roots:
        return {}, {}, {}
    d = nx.multi_source_dijkstra_path_length(G, roots, weight="weight")
    dn = {}
    for n in d:
        best, bestd = None, d[n]
        for nbr in G.neighbors(n):
            if d.get(nbr, math.inf) < bestd:
                best, bestd = nbr, d[nbr]
        if best is not None:
            dn[n] = best
    acc = {n: 0.0 for n in d}
    for u, v, data in G.edges(data=True):
        du, dv = d.get(u, math.inf), d.get(v, math.inf)
        if math.isinf(du) and math.isinf(dv):
            continue
        up = u if du >= dv else v
        if up in acc:
            acc[up] += data.get("weight", 0.0)
    for n in sorted(d, key=d.get, reverse=True):
        if n in dn:
            acc[dn[n]] = acc.get(dn[n], 0.0) + acc.get(n, 0.0)
    return acc, d, dn


def upstream_nodes(G, dn, start, max_len_m):
    """Walk upstream (against the downstream pointers) from ``start``,
    following the largest-accumulation branch, until ``max_len_m``."""
    ups_of = {}
    for child, parent in dn.items():
        ups_of.setdefault(parent, []).append(child)
    out, cur, walked = [start], start, 0.0
    while walked < max_len_m:
        ups = ups_of.get(cur, [])
        if not ups:
            break
        nxt = ups[0] if len(ups) == 1 else max(
            ups, key=lambda m: G.degree(m))
        if not G.has_edge(cur, nxt):
            break
        walked += G.edges[cur, nxt].get("weight", 0.0)
        out.append(nxt)
        cur = nxt
    return out


def _ramp_down(d, full_m, zero_m):
    """1 at/below full_m, 0 at/above zero_m, linear between (nearness score)."""
    if d <= full_m:
        return 1.0
    if d >= zero_m:
        return 0.0
    return float((zero_m - d) / (zero_m - full_m))


def load_settlements(path, to_utm):
    """OSM settlement points -> list of UTM (x, y). Used for access + consent."""
    if not path or not Path(path).exists():
        return []
    fc = json.load(open(path))
    pts = []
    for f in fc.get("features", []):
        try:
            lon, lat = f["geometry"]["coordinates"]
            pts.append(to_utm.transform(lon, lat))
        except Exception:  # noqa: BLE001
            continue
    return pts


def river_points_utm(osm_path, to_utm, bbox=None):
    """Vertices of OSM river/stream lines as UTM points (for nearest-water
    distance). Restricted to a padded bbox so we don't load the whole planet
    extract if the file is large."""
    if not osm_path or not Path(osm_path).exists():
        return []
    fc = json.load(open(osm_path))
    pts = []
    for f in fc.get("features", []):
        if f.get("properties", {}).get("waterway") not in ("river", "stream"):
            continue
        geom = f.get("geometry", {})
        coords = geom.get("coordinates", [])
        rings = coords if geom.get("type") == "MultiLineString" else [coords]
        for ring in rings:
            for c in ring[::3]:  # decimate; water-distance is coarse anyway
                lon, lat = c[0], c[1]
                if bbox and not (bbox[0]-0.05 <= lon <= bbox[2]+0.05 and bbox[1]-0.05 <= lat <= bbox[3]+0.05):
                    continue
                pts.append(to_utm.transform(lon, lat))
    return pts


def access_factor(x, y, settle_tree, river_tree):
    """Soft 0.6-1.0 multiplier from nearness to a settlement (crew origin) and
    to navigable water (material transport). Both from real OSM data."""
    rs = rr = 0.0
    if settle_tree is not None:
        d, _ = settle_tree.query([x, y])
        rs = _ramp_down(float(d), SETTLEMENT_FULL_M, SETTLEMENT_ZERO_M)
    if river_tree is not None:
        d, _ = river_tree.query([x, y])
        rr = _ramp_down(float(d), RIVER_FULL_M, RIVER_ZERO_M)
    near = 0.5 * (rs + rr) if (settle_tree is not None and river_tree is not None) else max(rs, rr)
    return ACCESS_MIN + (1 - ACCESS_MIN) * near, rs, rr


def node_reach_attrs(G, node):
    """Best (max-width) incident reach attributes at a node."""
    best = None
    for nbr in G.neighbors(node):
        d = G.edges[node, nbr]
        if best is None or d["width_m"] > best["width_m"]:
            best = d
    return best or {"width_m": 0.0, "sinuosity": 1.0, "reach_len_m": 0.0}


def sample_buffer(ds_arr, transform, nodata, lon, lat, buffer_m, signed_neg=False):
    """Buffered-mean sample, half-pixel-centred. If signed_neg, use
    |min(0, value)| (subsidence: negative = sinking)."""
    a, c, e, f = transform.a, transform.c, transform.e, transform.f
    col = (lon - c) / a - 0.5
    row = (lat - f) / e - 0.5
    col, row = int(round(col)), int(round(row))
    deg_lat = buffer_m / 111320.0
    br = max(1, int(round(deg_lat / abs(e))))
    bc = max(1, int(round(deg_lat / abs(a))))
    r0, r1 = max(0, row - br), min(ds_arr.shape[0], row + br + 1)
    c0, c1 = max(0, col - bc), min(ds_arr.shape[1], col + bc + 1)
    if r1 <= r0 or c1 <= c0:
        return None
    win = ds_arr[r0:r1, c0:c1]
    valid = np.isfinite(win) & (win != nodata)
    if not valid.any():
        return None
    v = win[valid]
    if signed_neg:
        v = np.abs(np.minimum(0.0, v))
    return float(v.mean())


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--centerlines", type=Path, default=None,
                   help="v2 centerline geojson (preferred). If omitted, the "
                        "--mask is skeletonized as a fallback (5.6 m AOIs).")
    p.add_argument("--mask", type=Path, required=True)
    p.add_argument("--restoration", type=Path, required=True)
    p.add_argument("--subsidence", type=Path, default=None)
    p.add_argument("--bbox", default=None, help="W,S,E,N restrict to AOI cell")
    p.add_argument("--aoi", default="SBG-CE", help="AOI code for seed tuples")
    p.add_argument("--target", type=int, default=15)
    p.add_argument("--spacing-m", type=float, default=SPACING_M)
    p.add_argument("--osm-rivers", type=Path, default=None,
                   help="OSM waterways geojson: classifies tips into true "
                        "river outlets vs bbox artifacts vs canal ends")
    p.add_argument("--elnino", type=Path, default=None,
                   help="Tay 2025 El Nino max-drawdown raster (cm): "
                        "drought-vulnerability score factor")
    p.add_argument("--settlements", type=Path, default=None,
                   help="OSM settlement points geojson: access factor + "
                        "village<->river consultation routing")
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()

    spacing = float(np.clip(args.spacing_m, SPACING_MIN, SPACING_MAX))
    to_utm = Transformer.from_crs("EPSG:4326", TARGET_CRS, always_xy=True)
    bbox = tuple(float(x) for x in args.bbox.split(",")) if args.bbox else None

    if args.centerlines and args.centerlines.exists():
        feats = json.load(open(args.centerlines))["features"]
        log.info("source: v2 centerline graph (%s), %d reaches", args.centerlines.name, len(feats))
    else:
        feats = reaches_from_mask(args.mask, bbox, to_utm)
        log.info("source: mask-skeleton fallback (%s)", args.mask.name)
    G = build_graph(feats, args.mask, to_utm)
    if G.number_of_nodes() == 0:
        raise SystemExit("empty graph")

    confluences = [n for n in G.nodes() if G.degree(n) >= 3]
    outlets = [n for n in G.nodes() if G.degree(n) == 1]
    log.info("topology: %d confluences (deg>=3), %d tips (deg1)", len(confluences), len(outlets))

    # ---- tip classification + flow accumulation ---------------------------
    river_tips, edge_tips, canal_ends = classify_tips(G, bbox, args.osm_rivers)
    # Roots per connected component: prefer that component's river outlets,
    # else its mid-peat canal ends, else any tip -- so EVERY component gets
    # routed (a single global root set leaves most components unrouted and
    # blanket-penalised).
    # Audit fix: track which components actually contain a RIVER outlet --
    # only there is the routed direction real drainage. Everywhere else the
    # roots are arbitrary canal ends, and accumulation is a network-centrality
    # proxy, which the score and remarks must say honestly.
    roots = set()
    river_directed: set = set()
    comp_of: dict = {}
    comp_best_d: dict = {}
    for ci, comp in enumerate(nx.connected_components(G)):
        comp = set(comp)
        for n in comp:
            comp_of[n] = ci
        if river_tips & comp:
            river_directed |= comp
        for pool in (river_tips & comp, canal_ends & comp,
                     {n for n in comp if G.degree(n) == 1}, comp):
            if pool:
                roots |= pool
                break
    acc, dist_to_root, dn = flow_accumulation(G, roots)
    if acc:
        top = max(acc.values())
        log.info("flow accumulation: %d routed nodes, max upstream length %.0f m", len(acc), top)
    # Navigation-consultation proxy: graph distance to a RIVER outlet (not to
    # arbitrary per-component roots) -- canals plugged into the river network
    # may be boat corridors and need community sign-off before blocking.
    dist_to_river = (
        nx.multi_source_dijkstra_path_length(G, river_tips, weight="weight")
        if river_tips else {}
    )

    # ---- access (feature B) + consultation routing (feature C) -----------
    from scipy.spatial import cKDTree
    settle_utm = load_settlements(args.settlements, to_utm)
    river_utm = river_points_utm(args.osm_rivers, to_utm, bbox)
    settle_tree = cKDTree(settle_utm) if settle_utm else None
    river_tree = cKDTree(river_utm) if river_utm else None
    log.info("access data: %d settlements, %d river vertices",
             len(settle_utm), len(river_utm))

    # Consent route-set: graph distance from settlement-snapped nodes, and the
    # direct settlement->river graph distance, so a candidate is "on a boat
    # route" if reaching it from a settlement then continuing to the river is
    # not much longer than going settlement->river directly.
    node_xy = np.array([[G.nodes[n]["x"], G.nodes[n]["y"]] for n in G.nodes()])
    node_list = list(G.nodes())
    settle_nodes = set()
    if settle_utm and len(node_xy):
        ntree = cKDTree(node_xy)
        for sx, sy in settle_utm:
            d, idx = ntree.query([sx, sy])
            if d <= CONSULT_SETTLEMENT_SNAP_M:
                settle_nodes.add(node_list[idx])
    dist_to_settle = (
        nx.multi_source_dijkstra_path_length(G, settle_nodes, weight="weight")
        if settle_nodes else {}
    )
    # direct settlement->river graph distance per component (min over routes)
    direct_sr = math.inf
    for sn in settle_nodes:
        dr = dist_to_river.get(sn)
        if dr is not None:
            direct_sr = min(direct_sr, dr)
    log.info("consultation routing: %d settlement-snapped canal nodes", len(settle_nodes))

    # ---- candidate nodes: confluences + outlets + chainage fill -----------
    cands = set(confluences) | set(outlets)
    # chainage fill: walk each connected component, sample nodes ~spacing apart
    for comp in nx.connected_components(G):
        sub = G.subgraph(comp)
        # pick the lowest-leverage spanning walk: BFS from an outlet if any
        seeds = [n for n in comp if G.degree(n) == 1] or [next(iter(comp))]
        seen = {}
        for seed in seeds:
            d = nx.single_source_dijkstra_path_length(sub, seed, weight="weight")
            for n, dist in d.items():
                if n not in seen or dist < seen[n]:
                    seen[n] = dist
        # bucket nodes by chainage / spacing, keep one per bucket
        buckets = {}
        for n, dist in seen.items():
            b = int(dist // spacing)
            if b not in buckets:
                buckets[b] = n
        cands |= set(buckets.values())

    # bbox filter
    if args.bbox:
        w, s, e, nn = (float(x) for x in args.bbox.split(","))
        cands = {n for n in cands
                 if w <= G.nodes[n]["lon"] <= e and s <= G.nodes[n]["lat"] <= nn}
    log.info("candidates: %d (confluences + outlets + chainage)", len(cands))

    # ---- score each candidate --------------------------------------------
    with rasterio.open(args.restoration) as rs:
        rest = rs.read(1).astype(np.float32)
        rest_tr, rest_nd = rs.transform, (rs.nodata if rs.nodata is not None else -9999)
    sub_arr = sub_tr = sub_nd = None
    if args.subsidence:
        with rasterio.open(args.subsidence) as ss:
            sub_arr = ss.read(1).astype(np.float32)
            sub_tr, sub_nd = ss.transform, (ss.nodata if ss.nodata is not None else -9999)
    en_arr = en_tr = en_nd = None
    if args.elnino and args.elnino.exists():
        with rasterio.open(args.elnino) as es:
            en_arr = es.read(1).astype(np.float32)
            en_tr, en_nd = es.transform, (es.nodata if es.nodata is not None else -9999)

    scored = []
    for n in cands:
        lon, lat = G.nodes[n]["lon"], G.nodes[n]["lat"]
        attrs = node_reach_attrs(G, n)
        dam = dammability(attrs["width_m"], attrs["sinuosity"], attrs["reach_len_m"])
        if dam <= 0:
            continue
        r_rest = sample_buffer(rest, rest_tr, rest_nd, lon, lat, REST_BUFFER_M)
        if r_rest is None:
            r_rest = 0.0
        # De-duplicate subsidence out of the composite (see header note):
        # point-sample the velocity at the candidate with the same support as
        # the composite, invert the known 0.45/0.55 blend.
        r_prox = r_rest
        if sub_arr is not None:
            v_pt = sample_buffer(sub_arr, sub_tr, sub_nd, lon, lat,
                                 REST_BUFFER_M, signed_neg=True)
            if v_pt is not None:
                sub_share = float(np.clip(v_pt / COMPOSITE_SEVERE_MM, 0, 1))
                r_prox = float(np.clip(
                    (r_rest - COMPOSITE_SUB_W * sub_share) / COMPOSITE_PROX_W, 0, 1))
        # Subsidence integrated along the reach the dam would rewet: mean of
        # buffered samples walked up to UPSTREAM_SAMPLE_M upstream, not just
        # the point under the dam.
        r_sub = 0.0
        sub_mm = None
        if sub_arr is not None:
            pts = upstream_nodes(G, dn, n, UPSTREAM_SAMPLE_M) if dn else [n]
            vals = []
            for m in pts[:5]:
                v = sample_buffer(sub_arr, sub_tr, sub_nd,
                                  G.nodes[m]["lon"], G.nodes[m]["lat"],
                                  SUB_BUFFER_M, signed_neg=True)
                if v is not None:
                    vals.append(v)
            if vals:
                rs_val = float(np.mean(vals))
                r_sub = float(np.clip(rs_val / SUB_NORM_MM, 0, 1))
                sub_mm = -rs_val  # report signed (negative = sinking)
        # El Nino drought-drawdown factor (sparse coverage -> gentle floor).
        r_en = 0.0
        en_cm = None
        if en_arr is not None:
            ev = sample_buffer(en_arr, en_tr, en_nd, lon, lat, EN_BUFFER_M, signed_neg=True)
            if ev is not None:
                r_en = float(np.clip(ev / EN_NORM_CM, 0, 1))
                en_cm = -ev
        is_conf = G.degree(n) >= 3
        tip = G.degree(n) == 1
        if is_conf:
            lev, tip_kind = LEV_CONFLUENCE, None
        elif tip and n in river_tips:
            lev, tip_kind = LEV_OUTLET, "river-outlet"
        elif tip and n in edge_tips:
            lev, tip_kind = LEV_EDGE_ARTIFACT, "bbox-artifact"
        elif tip:
            lev, tip_kind = LEV_CANAL_END, "canal-end"
        else:
            lev, tip_kind = 1.0, None
        # Flow-accumulation leverage: upstream canal length routed through n.
        acc_m = acc.get(n, 0.0)
        acc_norm = float(np.clip(acc_m / ACC_NORM_M, 0, 1)) if n in acc else 0.3
        # Access (feature B): reachability for a crew with materials.
        acc_fac, reach_settle, reach_river = access_factor(
            G.nodes[n]["x"], G.nodes[n]["y"], settle_tree, river_tree)
        # Consultation (feature C): on a plausible village<->river boat route?
        d_set = dist_to_settle.get(n, math.inf)
        d_riv = dist_to_river.get(n, math.inf)
        on_route = (math.isfinite(d_set) and math.isfinite(d_riv)
                    and math.isfinite(direct_sr)
                    and (d_set + d_riv) <= CONSULT_DETOUR_FRAC * direct_sr)
        consult = bool(on_route or (d_riv < CONSULT_RIVER_M and math.isfinite(d_set)
                                    and d_set < 3000))
        score = (REST_FLOOR_W + (1 - REST_FLOOR_W) * r_prox) * \
                (SUB_FLOOR_W + (1 - SUB_FLOOR_W) * r_sub) * \
                (ACC_FLOOR + (1 - ACC_FLOOR) * acc_norm) * \
                (EN_FLOOR + (1 - EN_FLOOR) * r_en) * dam * lev * acc_fac
        scored.append({
            "node": n, "lon": lon, "lat": lat, "score": score,
            "width_m": round(attrs["width_m"], 1), "sinuosity": round(attrs["sinuosity"], 3),
            "reach_len_m": round(attrs["reach_len_m"], 0),
            "r_rest": round(r_rest, 3), "r_prox": round(r_prox, 3),
            "sub_mm_yr": round(sub_mm, 1) if sub_mm is not None else None,
            "elnino_cm": round(en_cm, 1) if en_cm is not None else None,
            "upstream_m": round(acc_m, 0),
            "drainage_directed": n in river_directed,
            "comp": comp_of.get(n, -1),
            "d_to_root": round(dist_to_root.get(n, -1), 0),
            "is_confluence": is_conf, "is_outlet": tip and n in river_tips,
            "tip_kind": tip_kind,
            "access": round(acc_fac, 3),
            "reach_settle_m": round(float(settle_tree.query([G.nodes[n]["x"], G.nodes[n]["y"]])[0]), 0) if settle_tree is not None else None,
            "consultation_required": consult,
            "dam_type": dam_type(attrs["width_m"]),
        })
    scored.sort(key=lambda d: d["score"], reverse=True)
    log.info("scored candidates above dammability>0: %d", len(scored))

    # ---- graph-geodesic spacing (greedy) ----------------------------------
    kept = []
    floor = max(MIN_SITE_SCORE_ABS,
                MIN_SCORE_FRACTION * scored[0]["score"]) if scored else 0.0
    for cand in scored:
        # Quality floor: return fewer, better sites rather than pad to target.
        if cand["score"] < floor:
            log.info("quality floor %.3f reached -- stopping at %d sites", floor, len(kept))
            break
        n = cand["node"]
        ok = True
        if kept:
            dist = nx.single_source_dijkstra_path_length(G, n, cutoff=spacing, weight="weight")
            for k in kept:
                if k["node"] in dist:
                    ok = False
                    break
        if ok:
            kept.append(cand)
        if len(kept) >= args.target:
            break
    log.info("dam sites: %d (graph-geodesic spacing %.0f m)", len(kept), spacing)

    def build_why(c, rank, total):
        """Plain-language 'why this site' remark composed from the actual
        scoring factors -- deterministic, every clause traces to a number the
        selector computed (no generative text, nothing to hallucinate)."""
        parts = []
        if c["sub_mm_yr"] is not None and c["sub_mm_yr"] <= -10:
            parts.append(f"peat here is sinking {abs(c['sub_mm_yr']):.0f} mm/yr "
                         f"(Tay et al. 2025 InSAR)")
        elif c["sub_mm_yr"] is not None:
            parts.append(f"measured subsidence {abs(c['sub_mm_yr']):.0f} mm/yr (Tay et al. 2025)")
        else:
            parts.append("no direct InSAR coverage; ranked on canal-network position")
        if c["elnino_cm"] is not None and c["elnino_cm"] <= -2:
            parts.append(f"surface dropped {abs(c['elnino_cm']):.0f} cm in the 2019 "
                         f"El Nino drought (drought-vulnerable)")
        # Honest drainage wording: "upstream" is only meaningful where the
        # component demonstrably reaches a river; elsewhere it is connected
        # network length, stated as such.
        # Above ~25 km the routed figure reflects a whole grid component
        # funnelled through one node (a single-flow-path model artifact), so
        # state the situation qualitatively rather than quote a silly number.
        if c["upstream_m"] >= 25_000:
            parts.append("a large share of this canal network drains toward this point")
        elif c["upstream_m"] >= 400 and c["drainage_directed"]:
            parts.append(f"blocks ~{c['upstream_m'] / 1000:.1f} km of upstream drainage")
        elif c["upstream_m"] >= 400:
            parts.append(f"sits on ~{c['upstream_m'] / 1000:.1f} km of connected canal network")
        if c["is_confluence"]:
            parts.append("sits at a canal junction (several branches within reach)")
        if c["tip_kind"] == "river-outlet":
            parts.append("true river outlet -- built first so upstream works hold water")
        if c.get("reach_settle_m") is not None and c["reach_settle_m"] <= 2000:
            parts.append(f"reachable: ~{c['reach_settle_m']/1000:.1f} km from the nearest village")
        elif c.get("reach_settle_m") is not None and c["reach_settle_m"] >= 5000:
            parts.append(f"remote: ~{c['reach_settle_m']/1000:.0f} km to the nearest village (plan logistics)")
        if c["dam_type"] == "engineered":
            parts.append(f"about {c['width_m']:.0f} m wide (confirm on site): a wide collector "
                         f"canal -- needs an ENGINEERED dam (heavier materials / machinery), "
                         f"not a village-built structure")
        else:
            parts.append(f"about {c['width_m']:.0f} m wide (confirm on site): a {c['dam_type']} "
                         f"dam a village crew can build")
        why = "Why here: " + "; ".join(parts) + f". Rank {rank}/{total} in this zone."
        if c["consultation_required"]:
            why += " River-connected canal: community consultation required before blocking."
        return why

    # ---- emit GeoJSON + seed literal --------------------------------------
    # Build order, audit fix: REAL downstream-first (Ritzema et al. 2014 --
    # downstream structures go in first so upstream works hold water).
    # Components are ranked by their best site's score; within a component
    # whose drainage direction is known (contains a river outlet), sites
    # build in ascending distance-to-outlet. Where direction is unknown
    # (arbitrary roots), fall back to score order rather than pretend.
    comp_best = {}
    for c in kept:
        comp_best[c["comp"]] = max(comp_best.get(c["comp"], 0.0), c["score"])
    kept_sorted = sorted(kept, key=lambda c: (
        -comp_best[c["comp"]],
        c["comp"],
        c["d_to_root"] if c["drainage_directed"] else -c["score"],
    ))
    feats = []
    for i, c in enumerate(kept_sorted, 1):
        feats.append({
            "type": "Feature",
            "geometry": {"type": "Point", "coordinates": [c["lon"], c["lat"]]},
            "properties": {
                "site_id": f"D-{i:02d}", "score": round(c["score"], 3),
                "lat": c["lat"], "lon": c["lon"], "width_m": c["width_m"],
                "sinuosity": c["sinuosity"], "reach_len_m": c["reach_len_m"],
                "r_rest": c["r_rest"], "r_prox": c["r_prox"], "sub_mm_yr": c["sub_mm_yr"],
                "elnino_cm": c["elnino_cm"], "upstream_m": c["upstream_m"],
                "drainage_directed": c["drainage_directed"],
                "access": c["access"], "reach_settle_m": c.get("reach_settle_m"),
                "is_confluence": c["is_confluence"], "is_outlet": c["is_outlet"],
                "tip_kind": c["tip_kind"],
                "consultation_required": c["consultation_required"],
                "dam_type": c["dam_type"], "build_order": i, "aoi": args.aoi,
                "why": build_why(c, i, len(kept_sorted)),
            },
        })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"type": "FeatureCollection", "features": feats}, indent=2))
    log.info("wrote %s (%d sites)", args.out, len(feats))
    log.info("score: min %.3f max %.3f; confluences %d outlets %d; widths %.1f-%.1f m",
             min(f["properties"]["score"] for f in feats) if feats else 0,
             max(f["properties"]["score"] for f in feats) if feats else 0,
             sum(f["properties"]["is_confluence"] for f in feats),
             sum(f["properties"]["is_outlet"] for f in feats),
             min(f["properties"]["width_m"] for f in feats) if feats else 0,
             max(f["properties"]["width_m"] for f in feats) if feats else 0)


if __name__ == "__main__":
    main()
