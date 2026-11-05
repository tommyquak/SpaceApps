"""Combine per-AOI dam-site GeoJSONs from select_dam_sites_v2.py into a single
globally-ranked set, emit the DAM_SITE_TASKS python literal for seed.py, a
combined GeoJSON for the map, and a preview PNG.

Cosmetic demo fields (payout, deadline, village, deliverables) follow the exact
cycling scheme already in seed.py; the substantive fields (site_id, aoi, the
verbatim selector rationale, polygon, score) come straight from the selector.
"""
from __future__ import annotations

import json
import sys
from datetime import date, timedelta
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

INPUTS = {
    "SBG-CE": "/tmp/dam_ce.geojson",
    "SBG-NE": "/tmp/dam_ne.geojson",
    "SBG-NW": "/tmp/dam_nw.geojson",
    "SBG-S": "/tmp/dam_s.geojson",
}
CENTERLINES = {
    "SBG-CE": "deliverables/canal_planet_v2_sbgce/canal_centerlines.geojson",
    "SBG-NE": "deliverables/canal_planet_v2_sbg-ne/canal_centerlines.geojson",
    "SBG-NW": "deliverables/canal_planet_v2_sbg-nw/canal_centerlines.geojson",
    "SBG-S": "deliverables/canal_planet_v2_sbg-s/canal_centerlines.geojson",
}
PAYOUTS = [200000, 220000, 240000, 260000, 280000]
DELIVERABLES = [
    "3 photos: before, during, after canal block",
    "GPS track of dam edge (>= 6 m linear)",
    "Timestamp within working window 06:00-17:00 WIB",
    "Voice note: rationale for blockage",
]
VILLAGE = "vil-hampangen"
START = date(2026, 7, 1)


def load_all() -> list[dict]:
    feats = []
    for aoi, path in INPUTS.items():
        d = json.loads(Path(path).read_text())
        for f in d["features"]:
            p = dict(f["properties"])
            p["aoi"] = aoi
            feats.append(p)
    feats.sort(key=lambda p: p["score"], reverse=True)
    return feats


def square(lon: float, lat: float, half: float = 0.0005) -> dict:
    lon0, lon1 = round(lon - half, 5), round(lon + half, 5)
    lat0, lat1 = round(lat - half, 5), round(lat + half, 5)
    return {
        "type": "Polygon",
        "coordinates": [[[lon0, lat0], [lon1, lat0], [lon1, lat1], [lon0, lat1], [lon0, lat0]]],
    }


def main() -> None:
    feats = load_all()
    n = len(feats)
    print(f"combined sites: {n}")
    from collections import Counter
    print("per-AOI:", dict(Counter(p["aoi"] for p in feats)))
    print("dam_type:", dict(Counter(p["dam_type"] for p in feats)))

    tuples_src = []
    combined = {"type": "FeatureCollection", "features": []}
    for i, p in enumerate(feats):
        tid = f"D-{i + 1:02d}"
        lon, lat = round(p["lon"], 5), round(p["lat"], 5)
        payout = PAYOUTS[i % len(PAYOUTS)]
        deadline = (START + timedelta(days=i % 21)).isoformat() + "T17:00"
        poly = square(lon, lat)
        why = p["why"]
        tuples_src.append(
            "    ({tid!r}, {aoi!r}, {vil!r}, 'canal_block',\n"
            "     {title!r},\n"
            "     {why!r},\n"
            "     {payout}, {deadline!r}, 'available', None, 'DANA',\n"
            "     {deliv!r},\n"
            "     {poly},\n"
            "     {lat}, {lon}, {score}),".format(
                tid=tid, aoi=p["aoi"], vil=VILLAGE, title=f"Block canal · {tid}",
                why=why, payout=payout, deadline=deadline, deliv=DELIVERABLES,
                poly=json.dumps(poly), lat=lat, lon=lon, score=round(p["score"], 3),
            )
        )
        cf = {"type": "Feature", "geometry": {"type": "Point", "coordinates": [lon, lat]},
              "properties": {**p, "site_id": tid}}
        combined["features"].append(cf)

    Path("/tmp/dam_all.geojson").write_text(json.dumps(combined, indent=2))
    literal = "DAM_SITE_TASKS = [\n" + "\n".join(tuples_src) + "\n]\n"
    Path("/tmp/dam_seed_tuples.py").write_text(literal)
    print("wrote /tmp/dam_all.geojson and /tmp/dam_seed_tuples.py")

    # ---- preview image ----
    fig, ax = plt.subplots(figsize=(11, 11))
    for aoi, clp in CENTERLINES.items():
        try:
            cl = json.loads(Path(clp).read_text())
            for f in cl["features"]:
                g = f["geometry"]
                lines = g["coordinates"] if g["type"] == "MultiLineString" else [g["coordinates"]]
                for ln in lines:
                    xs = [c[0] for c in ln]
                    ys = [c[1] for c in ln]
                    ax.plot(xs, ys, color="#9fb8c8", lw=0.3, alpha=0.5, zorder=1)
        except Exception as e:
            print(f"  (centerlines {aoi} skipped: {e})")

    sc = ax.scatter(
        [p["lon"] for p in feats], [p["lat"] for p in feats],
        c=[p["score"] for p in feats], cmap="YlOrRd",
        s=[max(18, p["width_m"] * 6) for p in feats],
        edgecolors="#222", linewidths=0.5, zorder=3,
    )
    plt.colorbar(sc, ax=ax, shrink=0.5, label="dam-site score")
    for i, p in enumerate(feats):
        ax.annotate(f"{i+1}", (p["lon"], p["lat"]), fontsize=5, ha="center", va="center", zorder=4)
    ax.set_title(f"PeatGuard dam sites v6 (calibrated + engineered band): {n} sites across 4 AOIs", fontsize=11)
    ax.set_xlabel("lon"); ax.set_ylabel("lat")
    ax.set_aspect("equal")
    legend = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#fdae61", markersize=8, label="marker size = canal width"),
    ]
    ax.legend(handles=legend, loc="lower left", fontsize=8)
    out = "deliverables/dam_sites_v6_preview.png"
    fig.tight_layout()
    fig.savefig(out, dpi=130, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
