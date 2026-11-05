"""Canal detection v2 -- continuous dammable centerlines + wide-river rejection.

Rebuilds the PlanetScope canal detector around two goals the v1 detector
(canal_detect_planet.py) misses:

  (1) CONTINUITY: emit continuous canal centerlines instead of the fragments
      v1 produces. v1 hard-thresholds a Sato ridge at one percentile and runs
      the area filter BEFORE any linking, so every sub-threshold gap (bridge,
      planting strip, cloud-shadow) splits a line and small real fragments are
      deleted before they can reconnect.
  (2) RIVER REJECTION: remove wide un-dammable rivers (Sebangau ~80-200 m,
      Kahayan tributaries) and their banks, which v1's 140 m isotropic buffer
      handles crudely and which also truncates canals tapping off the river.

Pipeline (plan: canal_detection_v2_plan.md):
  composite (extended: red, rededge, persistence)
    -> features (NDWI, NDRE)
    -> response front-end: RivaMap Multiscale Singularity Index (edge-
       suppressing) if available, else max(frangi, meijering) fallback
    -> persistence weighting (soft)
    -> hysteresis double-threshold linking (recover dim-but-connected runs)
    -> river rejection: distance-transform width + morphological
       open-and-reclaim (never deletes whole components, so canals tapping
       off the river survive), with OSM rivers as soft negatives and the
       Dadap 2021 canal raster as positive protection
    -> component cleanup AFTER linking (size + aspect)
    -> skeleton -> graph: prune spurs, bridge near-collinear endpoints
    -> outputs: canal_mask.tif (uint8 COG) + canal_centerlines.geojson

Citations (see canal_detection_v2_plan.md for the full list with confidence
tags): Isikdogan et al. 2017 (RivaMap MSI); Frangi 1998; Meijering 2004;
Canny 1986 (hysteresis); Pavelsky & Smith 2008 (distance-transform width);
Dadap et al. 2021 (CNN canal precedent + positive-protection reference);
Steger 1998 (centerline graph); Mueller 1968 / Leopold 1964 (sinuosity).
"""

from __future__ import annotations

import argparse
import glob
import json
import logging
from pathlib import Path

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.features import rasterize
from rasterio.warp import Resampling, reproject
from scipy.ndimage import (
    binary_closing,
    binary_dilation,
    binary_opening,
    distance_transform_edt,
    label,
)
from skimage.filters import frangi, meijering
from skimage.measure import regionprops
from skimage.morphology import skeletonize

from peatguard.ingest.planetscope import composite

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("canal_detect_v2")

# Optional response front-end: RivaMap Multiscale Singularity Index suppresses
# the riverbank step-edge response that Frangi/Sato fire on. Falls back to
# Frangi+Meijering if rivamap is unavailable.
try:
    from rivamap.singularity_index import SingularityIndexFilters, applyMMSI
    _HAVE_RIVAMAP = True
except Exception:  # noqa: BLE001
    _HAVE_RIVAMAP = False

# Optional skeleton->graph helper.
try:
    import sknw
    _HAVE_SKNW = True
except Exception:  # noqa: BLE001
    _HAVE_SKNW = False

# ---- Tunables (set for the OUTPUT resolution, default 3 m) ------------------
DEFAULT_RES_DEG = 0.00003          # ~3.3 m native SuperDove
NDWI_THRESHOLD = 0.0
NIR_WATER_MAX = 1000.0
# Response: ridge half-widths for 3-15 m canals at ~3 m (half-width 0.5-2.5 px)
FRANGI_SIGMAS = (1.0, 1.5, 2.5)
RIVAMAP_MIN_SCALE = 1.0
RIVAMAP_NR_SCALES = 5
# Hysteresis double-threshold (percentiles of the response over finite pixels)
HYST_HIGH_PCT = 97.0               # confident seed
HYST_LOW_PCT = 90.0                # grow, kept only if connected to a seed
# Persistence weighting OFF by default: multiplying the response by wet-
# persistence collapsed recall ~10x in the few-scene regime (only ~4 scenes
# cover SBG-CE, so most dry-season canals score low persistence and get
# suppressed). Empirically (2026-05-30) persistence ON gave 0.19% canal
# density vs 1.94% OFF, with the width-based river rejection keeping rivers
# out either way. Persistence is still computed and carried as a per-reach
# confidence attribute; it just no longer gates detection. Re-enable only
# with a dense multi-season stack (and prefer an additive, not multiplicative,
# boost).
PERSISTENCE_WEIGHT = False
# River rejection
RIVER_WIDTH_M = 13.0               # half-width above which a wet pixel is river
RIVER_BANK_BUFFER_M = 25.0         # extra buffer around river to kill banks
OSM_RIVER_BUFFER_M = 40.0          # soft-negative buffer around OSM rivers
DADAP_PROTECT = True               # exempt Dadap-overlapping canals from cut
# Component cleanup (AFTER linking)
MIN_SEGMENT_M2 = 1000.0            # ~ drop specks below this area
MIN_ASPECT_RATIO = 3.0
# Continuity: directional morphological closing bridges gaps ALONG the canal's
# own local orientation before skeletonization (a cheap path-closing that
# exploits the straight plantation-grid prior without vendoring a path
# operator). Gaps up to ~this length along any tested orientation are healed.
# 22.5-degree steps: the previous 45-degree set missed canals at intermediate
# headings, leaving their gaps unhealed.
DIRECTIONAL_CLOSE_LEN_M = 45.0
DIRECTIONAL_CLOSE_ANGLES = (0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5)
# Skeleton graph
SPUR_PRUNE_M = 30.0                # drop dangling spurs shorter than this
BRIDGE_MAX_M = 220.0               # bridge endpoint gaps shorter than this
BRIDGE_MAX_TURN_DEG = 38.0         # only if near-collinear
BRIDGE_PASSES = 3                  # iterate prune+bridge (each pass extends reach)
# Geodesic escalation (plan step 13, Cohen & Kimmel 1997): endpoint pairs the
# straight-line bridge rejects are routed through the inverse-response cost
# surface with skimage.graph.route_through_array; accepted only if the path
# follows elevated response (a dim canal trace), so it cannot invent canals
# across blank peat or reinforce sinuous river traces.
GEO_BRIDGE_MAX_M = 300.0           # search radius for geodesic bridging
GEO_BRIDGE_MIN_R = 0.06            # mean response along path must exceed this
GEO_BRIDGE_MARGIN_PX = 24          # window margin around each endpoint pair
RIVER_SINUOSITY = 1.4             # a reach is river only if wide AND sinuous
# -----------------------------------------------------------------------------


def _pixel_m(res_deg: float) -> float:
    return res_deg * 111320.0


def compute_response(ndwi: np.ndarray, ndre: np.ndarray | None,
                     finite: np.ndarray) -> np.ndarray:
    """Curvilinear response R in [0,1]. RivaMap MSI if available (suppresses
    riverbank edges), else max(frangi, meijering) over NDWI (+NDRE)."""
    fill = float(np.nanpercentile(ndwi[finite], 20)) if finite.any() else 0.0
    ndwi_f = np.where(finite, ndwi, fill).astype(np.float32)

    if _HAVE_RIVAMAP:
        filters = SingularityIndexFilters(minScale=RIVAMAP_MIN_SCALE,
                                          nrScales=RIVAMAP_NR_SCALES)
        psi, _widths, _orient = applyMMSI(ndwi_f, filters)
        R = np.nan_to_num(psi, nan=0.0).astype(np.float32)
        log.info("response: RivaMap MSI (edge-suppressing)")
    else:
        r1 = frangi(ndwi_f, sigmas=FRANGI_SIGMAS, black_ridges=False)
        r = np.nan_to_num(r1, nan=0.0)
        if ndre is not None:
            ndre_f = np.where(np.isfinite(ndre), ndre,
                              float(np.nanpercentile(ndre[np.isfinite(ndre)], 20))
                              if np.isfinite(ndre).any() else 0.0).astype(np.float32)
            r = np.maximum(r, np.nan_to_num(meijering(ndwi_f, sigmas=FRANGI_SIGMAS,
                                                       black_ridges=False), nan=0.0))
        R = r.astype(np.float32)
        log.info("response: Frangi+Meijering fallback (rivamap unavailable)")
    # normalise to [0,1]
    mx = float(R[finite].max()) if finite.any() and R[finite].max() > 0 else 1.0
    return (R / mx).astype(np.float32)


def hysteresis_link(R: np.ndarray, finite: np.ndarray,
                    hi_pct: float, lo_pct: float) -> tuple[np.ndarray, float]:
    """Double-threshold: keep low-threshold pixels only if 8-connected to a
    high-threshold seed. Recovers dim-but-connected canal runs.

    The curvilinear response (RivaMap MSI / Frangi) is near-zero off canals,
    so percentiles must be taken over the POSITIVE response -- a percentile
    over all pixels collapses to 0 (the majority) and the hysteresis becomes
    a no-op that keeps everything.
    """
    pos = finite & (R > 0)
    if pos.sum() < 50:
        return np.zeros_like(finite), 0.0
    hi = float(np.percentile(R[pos], hi_pct))
    lo = float(np.percentile(R[pos], lo_pct))
    if hi <= 0:
        hi = float(R[pos].max()) * 0.5
    if lo <= 0:
        lo = hi * 0.3
    seed = finite & (R >= hi)
    grow = finite & (R >= lo)
    lab, n = label(grow)
    if n == 0:
        return np.zeros_like(finite), lo
    keep_labels = np.unique(lab[seed])
    keep_labels = keep_labels[keep_labels > 0]
    out = np.isin(lab, keep_labels)
    log.info("hysteresis: hi=%.5f lo=%.5f (over %d positive px) -> %d seeded components",
             hi, lo, int(pos.sum()), len(keep_labels))
    return out, lo


def build_wet_mask(ndwi: np.ndarray, nir: np.ndarray, finite: np.ndarray) -> np.ndarray:
    wet = finite & ((ndwi > NDWI_THRESHOLD) | (np.isfinite(nir) & (nir < NIR_WATER_MAX)))
    return binary_closing(wet, iterations=1)


def river_mask_from_width(wet: np.ndarray, pixel_m: float) -> np.ndarray:
    """Open-and-reclaim: morphologically open the wet mask by ~half the
    river/canal width gap. Wide river core + banks survive the opening;
    narrow canals vanish and are reclaimed (kept as canal). Returns the
    river mask (wide water), buffered to cover banks."""
    r_open = max(1, int(round(RIVER_WIDTH_M / pixel_m)))
    struct = _disk(r_open)
    river = binary_opening(wet, structure=struct)
    # buffer to kill the step-edge banks the ridge filter fires on
    buf = max(1, int(round(RIVER_BANK_BUFFER_M / pixel_m)))
    river = binary_dilation(river, iterations=buf)
    return river


def _disk(radius: int) -> np.ndarray:
    y, x = np.ogrid[-radius:radius + 1, -radius:radius + 1]
    return (x * x + y * y <= radius * radius)


def _line_se(length_px: int, angle_deg: float) -> np.ndarray:
    """A 1-px-thick line structuring element of the given length + angle."""
    rad = np.radians(angle_deg)
    dx, dy = np.cos(rad), np.sin(rad)
    half = max(1, length_px // 2)
    pts = set()
    for t in np.linspace(-half, half, length_px * 2 + 1):
        pts.add((int(round(t * dy)), int(round(t * dx))))
    ys = [p[0] for p in pts]; xs = [p[1] for p in pts]
    miny, minx = min(ys), min(xs)
    se = np.zeros((max(ys) - miny + 1, max(xs) - minx + 1), dtype=bool)
    for y, x in pts:
        se[y - miny, x - minx] = True
    return se


def directional_close(mask: np.ndarray, pixel_m: float,
                      length_m: float | None = None,
                      angles: tuple = DIRECTIONAL_CLOSE_ANGLES) -> np.ndarray:
    """Close gaps ALONG the canal orientation: union of binary closings by a
    line structuring element at several angles. A closing with a length-L line
    bridges a gap shorter than L only when the two pieces are collinear with
    the line, so parallel canals are NOT merged (the gap there is
    perpendicular). This is a cheap directional path-closing that exploits the
    straight plantation-grid prior [standard morphology; cf. path closings,
    Heijmans/Talbot 2005]."""
    L = max(3, int(round((length_m or DIRECTIONAL_CLOSE_LEN_M) / pixel_m)))
    out = mask.copy()
    for ang in angles:
        se = _line_se(L, ang)
        out |= binary_closing(mask, structure=se)
    return out


def rasterize_osm_rivers(geojson_path: Path, transform, shape,
                         pixel_m: float) -> np.ndarray:
    """Rasterize OSM waterway=river/stream lines, buffered, as a soft river
    negative mask. Returns a boolean grid (True = near an OSM river)."""
    try:
        fc = json.load(open(geojson_path))
    except Exception as e:  # noqa: BLE001
        log.warning("OSM rivers unreadable (%s) -- skipping", e)
        return np.zeros(shape, dtype=bool)
    from shapely.geometry import shape as shp_shape
    buf_deg = OSM_RIVER_BUFFER_M / 111320.0
    geoms = []
    for f in fc.get("features", []):
        if f.get("properties", {}).get("waterway") not in ("river", "stream"):
            continue
        try:
            geoms.append(shp_shape(f["geometry"]).buffer(buf_deg))
        except Exception:  # noqa: BLE001
            continue
    if not geoms:
        log.info("OSM: no river/stream features in AOI")
        return np.zeros(shape, dtype=bool)
    mask = rasterize([(g, 1) for g in geoms], out_shape=shape, transform=transform,
                     fill=0, dtype="uint8")
    log.info("OSM: %d river/stream features rasterized as soft negatives", len(geoms))
    return mask.astype(bool)


def load_dadap_canals(dadap_path: Path, transform, shape, dst_crs="EPSG:4326") -> np.ndarray:
    """Reproject the Dadap 2021 canal raster onto the grid. Returns a boolean
    canal-presence mask (True where Dadap mapped a canal)."""
    try:
        with rasterio.open(dadap_path) as src:
            dst = np.zeros(shape, dtype=np.float32)
            reproject(source=rasterio.band(src, 1), destination=dst,
                      src_transform=src.transform, src_crs=src.crs,
                      dst_transform=transform, dst_crs=dst_crs,
                      resampling=Resampling.max, src_nodata=src.nodata, dst_nodata=0)
        m = dst > 0
        log.info("Dadap: %d canal pixels overlap the AOI (positive protection)", int(m.sum()))
        return m
    except Exception as e:  # noqa: BLE001
        log.warning("Dadap raster unreadable (%s) -- skipping positive protection", e)
        return np.zeros(shape, dtype=bool)


def shape_cleanup(mask: np.ndarray, pixel_m: float,
                  linked: np.ndarray | None = None,
                  R: np.ndarray | None = None,
                  r_floor: float | None = None) -> np.ndarray:
    """Drop specks below MIN_SEGMENT_M2 and near-circular blobs (aspect < 3).

    When ``linked`` (the directionally-closed mask) is given, size/aspect are
    judged on the LINKED connected component each pixel belongs to, and the
    verdict is applied back to the un-thickened ``mask``. A 40 m real canal
    fragment that joins a kilometre-long linked run therefore survives, while
    an isolated speck still dies -- this is the true "cleanup AFTER linking"
    the v2 plan called for (previously cleanup ran before the linkers and
    deleted fragments they could have reconnected)."""
    min_px = int(MIN_SEGMENT_M2 / (pixel_m * pixel_m))
    judge = linked if linked is not None else mask
    lab, n = label(judge)
    if n == 0:
        return mask
    # Mean ridge response per linked component: separates real canal NETWORKS
    # (crossing grid canals merge into aspect~1 patches but ride strong
    # response) from speckle fields the close glued together (barely above the
    # hysteresis floor). Shape alone cannot make this distinction.
    if R is not None and r_floor is not None and n > 0:
        from scipy.ndimage import mean as nd_mean
        comp_mean_r = nd_mean(np.nan_to_num(R, nan=0.0), labels=lab, index=np.arange(1, n + 1))
        r_gate = 1.6 * max(r_floor, 1e-6)
    else:
        comp_mean_r = None
        r_gate = None
    keep = np.zeros(n + 1, dtype=bool)
    for r in regionprops(lab):
        if r.area < min_px:
            continue
        # Aspect only filters SMALL isolated blobs. A large linked component
        # may be a canal NETWORK whose whole-component major/minor ratio
        # approaches 1 -- admit it only when its mean response clears the
        # strength gate (no response info -> fall back to the aspect rule).
        if r.area >= 10 * min_px and comp_mean_r is not None:
            if comp_mean_r[r.label - 1] >= r_gate:
                keep[r.label] = True
                continue
        minor = r.axis_minor_length
        if minor <= 0 or (r.axis_major_length / minor) >= MIN_ASPECT_RATIO:
            keep[r.label] = True
    return mask & keep[lab]


def _geodesic_bridge(graph, R: np.ndarray, pixel_m: float) -> int:
    """Escalation bridging (Cohen & Kimmel 1997 minimal path): for dangling
    endpoint pairs the straight-line bridge rejected, route through the
    inverse-response cost surface and keep the path only when it rides
    elevated response (a dim but real canal trace). Straightness-gated so it
    cannot reinforce sinuous river skeletons."""
    from scipy.spatial import cKDTree
    from skimage.graph import route_through_array

    max_px = GEO_BRIDGE_MAX_M / pixel_m
    ends = [n for n in graph.nodes() if graph.degree(n) == 1]
    if len(ends) < 2:
        return 0
    pos = np.array([graph.nodes[n]["o"] for n in ends], dtype=float)
    tree = cKDTree(pos)
    pairs = sorted(tree.query_pairs(r=max_px),
                   key=lambda ij: np.hypot(*(pos[ij[0]] - pos[ij[1]])))
    cost = (1.0 / (np.maximum(R, 0.0) + 0.05)).astype(np.float32)
    cos_thresh = np.cos(np.radians(BRIDGE_MAX_TURN_DEG + 12))
    H, W = R.shape
    m = GEO_BRIDGE_MARGIN_PX
    added = 0
    for i, j in pairs:
        n1, n2 = ends[i], ends[j]
        # Re-check degree: an earlier bridge this pass may have consumed one.
        if graph.has_edge(n1, n2) or graph.degree(n1) != 1 or graph.degree(n2) != 1:
            continue
        d1 = _endpoint_outward_dir(graph, n1)
        d2 = _endpoint_outward_dir(graph, n2)
        if d1 is None or d2 is None:
            continue
        c = pos[j] - pos[i]
        cn = np.hypot(*c)
        if cn < 1e-6:
            continue
        cu = c / cn
        if np.dot(d1, cu) < cos_thresh or np.dot(d2, -cu) < cos_thresh:
            continue
        r0 = max(0, int(min(pos[i][0], pos[j][0])) - m)
        r1 = min(H, int(max(pos[i][0], pos[j][0])) + m + 1)
        c0 = max(0, int(min(pos[i][1], pos[j][1])) - m)
        c1 = min(W, int(max(pos[i][1], pos[j][1])) + m + 1)
        start = (int(pos[i][0]) - r0, int(pos[i][1]) - c0)
        end = (int(pos[j][0]) - r0, int(pos[j][1]) - c0)
        try:
            path, _tc = route_through_array(cost[r0:r1, c0:c1], start, end,
                                            fully_connected=True, geometric=True)
        except Exception:  # noqa: BLE001
            continue
        path = np.asarray(path, dtype=int)
        mean_r = float(R[path[:, 0] + r0, path[:, 1] + c0].mean())
        if mean_r < GEO_BRIDGE_MIN_R:
            continue
        pts = path + np.array([r0, c0])
        graph.add_edge(n1, n2, pts=pts, weight=float(cn))
        added += 1
    return added


def centerline_graph(mask: np.ndarray, transform, pixel_m: float,
                     width_m: np.ndarray, R: np.ndarray | None = None):
    """Skeletonize -> graph -> prune spurs + bridge collinear endpoints
    (straight first, then geodesic through the response surface).
    Returns a list of GeoJSON LineString features (lon/lat) with per-reach
    length_m, mean_width_m, sinuosity."""
    skel = skeletonize(mask)
    if not _HAVE_SKNW:
        log.warning("sknw unavailable -- emitting skeleton pixels without graph linking")
        return [], skel
    graph = sknw.build_sknw(skel, multi=False)

    def to_lonlat(rc):
        # sknw points are (row, col)
        xs, ys = rasterio.transform.xy(transform, rc[:, 0], rc[:, 1])
        return list(zip(np.atleast_1d(xs), np.atleast_1d(ys)))

    spur_px = SPUR_PRUNE_M / pixel_m
    # Iterate prune-then-bridge: each pass that bridges a gap creates new
    # through-paths whose remaining dangling ends the next pass can extend.
    total_bridged = 0
    for _ in range(BRIDGE_PASSES):
        deg = dict(graph.degree())
        to_drop = []
        for s, e in list(graph.edges()):
            pts = graph.get_edge_data(s, e)["pts"]
            if (deg.get(s, 0) == 1 or deg.get(e, 0) == 1) and _polyline_len(pts) < spur_px:
                to_drop.append((s, e))
        graph.remove_edges_from(to_drop)
        n_bridged = _bridge_endpoints(graph, pixel_m)
        total_bridged += n_bridged
        if n_bridged == 0:
            break
    if total_bridged:
        log.info("centerline graph: bridged %d collinear endpoint gaps", total_bridged)

    if R is not None:
        n_geo = _geodesic_bridge(graph, R, pixel_m)
        if n_geo:
            log.info("centerline graph: geodesic-bridged %d response-following gaps", n_geo)

    feats = []
    for s, e in graph.edges():
        data = graph.get_edge_data(s, e)
        pts = data["pts"]
        if len(pts) < 2:
            continue
        coords = to_lonlat(pts)
        length_px = _polyline_len(pts)
        length_m = length_px * pixel_m
        # endpoint straight-line distance for sinuosity
        d_end = np.hypot(pts[0, 0] - pts[-1, 0], pts[0, 1] - pts[-1, 1]) * pixel_m
        sinuosity = (length_m / d_end) if d_end > 1e-6 else 1.0
        # width along the reach -- 75th percentile is robust to skeleton drift
        # (which deflates a mean) and to confluence spikes (which inflate a max).
        ww = width_m[pts[:, 0].astype(int), pts[:, 1].astype(int)]
        mean_w = float(np.nanpercentile(ww, 75)) if ww.size else 0.0
        feats.append({
            "type": "Feature",
            "properties": {
                "length_m": round(length_m, 1),
                "mean_width_m": round(mean_w, 1),
                "sinuosity": round(float(sinuosity), 3),
            },
            "geometry": {"type": "LineString",
                         "coordinates": [[float(x), float(y)] for x, y in coords]},
        })
    log.info("centerline graph: %d reaches after spur prune", len(feats))
    return feats, skel


def _polyline_len(pts: np.ndarray) -> float:
    d = np.diff(pts.astype(float), axis=0)
    return float(np.hypot(d[:, 0], d[:, 1]).sum())


def _endpoint_outward_dir(graph, node) -> np.ndarray | None:
    """Unit (row,col) direction pointing OUT of the line at a degree-1 node
    (from a few px inside the incident edge toward the node)."""
    nbrs = list(graph.neighbors(node))
    if not nbrs:
        return None
    pts = np.asarray(graph.get_edge_data(node, nbrs[0])["pts"], dtype=float)
    o = np.asarray(graph.nodes[node]["o"], dtype=float)
    # orient pts so the END is at the node
    if np.hypot(*(pts[0] - o)) < np.hypot(*(pts[-1] - o)):
        pts = pts[::-1]
    k = min(5, len(pts) - 1)
    if k < 1:
        return None
    d = o - pts[-1 - k]
    nrm = np.hypot(*d)
    return d / nrm if nrm > 1e-6 else None


def _bridge_endpoints(graph, pixel_m: float) -> int:
    """Add straight bridge edges between near-collinear dangling endpoints
    whose gap is shorter than BRIDGE_MAX_M (Steger 1998 centerline linking)."""
    from scipy.spatial import cKDTree
    bridge_px = BRIDGE_MAX_M / pixel_m
    ends = [n for n in graph.nodes() if graph.degree(n) == 1]
    if len(ends) < 2:
        return 0
    pos = np.array([graph.nodes[n]["o"] for n in ends], dtype=float)
    tree = cKDTree(pos)
    pairs = tree.query_pairs(r=bridge_px)
    added = 0
    cos_thresh = np.cos(np.radians(BRIDGE_MAX_TURN_DEG))
    for i, j in pairs:
        n1, n2 = ends[i], ends[j]
        if n1 == n2 or graph.has_edge(n1, n2):
            continue
        d1 = _endpoint_outward_dir(graph, n1)
        d2 = _endpoint_outward_dir(graph, n2)
        if d1 is None or d2 is None:
            continue
        c = pos[j] - pos[i]
        cn = np.hypot(*c)
        if cn < 1e-6:
            continue
        c = c / cn
        # n1's outward dir should align with the connector toward n2, and
        # n2's outward dir with the connector toward n1 (-c).
        if np.dot(d1, c) >= cos_thresh and np.dot(d2, -c) >= cos_thresh:
            steps = max(2, int(cn))
            bridge_pts = np.linspace(pos[i], pos[j], steps).round().astype(int)
            graph.add_edge(n1, n2, pts=bridge_pts, weight=float(cn))
            added += 1
    return added


def detect(scene_paths, bbox, res_deg=DEFAULT_RES_DEG,
           osm_rivers: Path | None = None, dadap: Path | None = None) -> dict:
    """Run the full v2 pipeline. Returns dict with mask, centerlines (GeoJSON
    FeatureCollection), transform, shape, stats."""
    pixel_m = _pixel_m(res_deg)
    comp = composite(scene_paths, bbox, res_deg=res_deg, extended=True,
                     ndwi_wet_threshold=NDWI_THRESHOLD, nir_water_max=NIR_WATER_MAX)
    ndwi = comp["ndwi"]; nir = comp["nir"]
    rededge = comp.get("rededge"); persistence = comp.get("persistence")
    finite = np.isfinite(ndwi)
    h, w = comp["shape"]
    transform = comp["transform"]

    # NDRE for the fallback response (drained-bank vegetation stress channel)
    ndre = None
    if rededge is not None:
        denom = nir + rededge
        ndre = np.where(denom > 0, (nir - rededge) / denom, np.nan).astype(np.float32)

    # 1-2. response front-end
    R = compute_response(ndwi, ndre, finite)
    # 3. persistence weighting (soft)
    if PERSISTENCE_WEIGHT and persistence is not None:
        R = R * np.nan_to_num(persistence, nan=0.0).astype(np.float32)
    # 4. hysteresis linking
    canal, r_floor = hysteresis_link(R, finite, HYST_HIGH_PCT, HYST_LOW_PCT)
    # also union with explicit wet (water-filled canals)
    wet = build_wet_mask(ndwi, nir, finite)
    canal = canal | (wet & ~_wide_only(wet, pixel_m))

    # 5. river rejection: open-and-reclaim + OSM soft + Dadap protect
    river = river_mask_from_width(wet, pixel_m)
    if osm_rivers is not None:
        river = river | rasterize_osm_rivers(osm_rivers, transform, (h, w), pixel_m)
    pre_river_px = int(canal.sum())
    canal_excl = canal & ~river
    # Dadap positive protection: re-add canal pixels the river cut removed that
    # overlap a Dadap-mapped canal (a real canal tapping off the river).
    if DADAP_PROTECT and dadap is not None:
        dmask = load_dadap_canals(dadap, transform, (h, w))
        reclaimed = canal & river & dmask
        canal_excl = canal_excl | reclaimed
        log.info("Dadap protection reclaimed %d near-river canal pixels", int(reclaimed.sum()))
    river_dropped = pre_river_px - int(canal_excl.sum())
    canal = canal_excl

    # 6. component cleanup on RAW components (strict size + aspect). Ablation
    # (2026-06-10, SBG-CE 4-scene stack) showed this filter is load-bearing
    # for precision: judging on a linked mask -- even a conservative 18 m
    # close, even response-gated -- let speckle fields glued by the close
    # slip through and inflated canal density 1.9% -> 4.3% of the AOI.
    # Continuity is therefore won in the LINKING stages below (45 m 8-angle
    # skeleton close, 220 m collinear bridges, geodesic escalation), not by
    # loosening cleanup.
    canal = shape_cleanup(canal, pixel_m)

    # 6b. directional close is applied ONLY to a skeletonization copy, to heal
    # short along-orientation gaps so centerlines come out continuous. It is
    # NOT written to the output mask, because closing thickens/merges canals
    # and would corrupt the EDT width estimate.
    canal_for_skel = directional_close(canal, pixel_m)
    log.info("directional close (skeleton only): +%d px bridged",
             int(canal_for_skel.sum()) - int(canal.sum()))

    # 7. width map from the CLEAN canal mask (full width = 2 x EDT; measured on
    # the canal mask, not the wet mask -- dry ridge-detected canals are absent
    # from the wet mask, which made the old wet-EDT width read ~0).
    width_m = (2.0 * distance_transform_edt(canal) * pixel_m).astype(np.float32)
    feats, skel = centerline_graph(canal_for_skel, transform, pixel_m, width_m, R=R)

    stats = {
        "pixel_m": round(pixel_m, 2),
        "response": "rivamap_msi" if _HAVE_RIVAMAP else "frangi_meijering",
        "canal_pixels": int(canal.sum()),
        "canal_pct_of_aoi": round(100.0 * float(canal.mean()), 3),
        "river_pixels_dropped": river_dropped,
        "centerline_reaches": len(feats),
        "valid_aoi_coverage_pct": round(100.0 * float(finite.mean()), 1),
        "n_scenes": comp["n_scenes"],
    }
    return {
        "mask": canal.astype(np.uint8), "skeleton": skel.astype(np.uint8),
        "centerlines": {"type": "FeatureCollection", "features": feats},
        "transform": transform, "shape": (h, w), "stats": stats,
    }


def _wide_only(wet: np.ndarray, pixel_m: float) -> np.ndarray:
    """Pixels belonging to wide (river-scale) water -- used to avoid pulling
    the whole river into the canal class via the explicit-wet union."""
    r_open = max(1, int(round(RIVER_WIDTH_M / pixel_m)))
    return binary_opening(wet, structure=_disk(r_open))


def write_outputs(result: dict, out_mask: Path, out_lines: Path) -> None:
    h, w = result["shape"]
    profile = {"driver": "COG", "dtype": "uint8", "count": 1, "width": w, "height": h,
               "crs": CRS.from_epsg(4326), "transform": result["transform"],
               "nodata": 0, "compress": "deflate", "overview_resampling": "nearest"}
    out_mask.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out_mask, "w", **profile) as dst:
        dst.write(result["mask"], 1)
    json.dump(result["centerlines"], open(out_lines, "w"))
    log.info("wrote %s and %s", out_mask, out_lines)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bbox", required=True, help="W,S,E,N")
    p.add_argument("--scenes", required=True, help="glob of 8b SR scenes (quote it)")
    p.add_argument("--out", type=Path, required=True, help="output canal_mask.tif")
    p.add_argument("--out-lines", type=Path, default=None, help="centerlines geojson")
    p.add_argument("--res-deg", type=float, default=DEFAULT_RES_DEG)
    p.add_argument("--osm-rivers", type=Path, default=None)
    p.add_argument("--dadap", type=Path, default=None)
    args = p.parse_args()

    bbox = tuple(float(x) for x in args.bbox.split(","))
    paths = sorted(glob.glob(args.scenes))
    if not paths:
        raise SystemExit(f"no scenes match: {args.scenes}")
    out_lines = args.out_lines or args.out.with_suffix(".geojson")
    res = detect(paths, bbox, res_deg=args.res_deg,
                 osm_rivers=args.osm_rivers, dadap=args.dadap)
    write_outputs(res, args.out, out_lines)
    log.info("")
    log.info("=== CANAL DETECTION v2 STATS ===")
    for k, v in res["stats"].items():
        log.info("  %s: %s", k, v)


if __name__ == "__main__":
    main()
