#!/usr/bin/env python3
"""Recheck provisional NISAR contrasts on common spatial support at dam sites.

Reads existing, filtered, millimetre contrast rasters. Does not apply further
deramping, fit either track to the other, or annualize the observations.
Missing acquisition/correction metadata remains an explicit release blocker.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.warp import transform
from rasterio.windows import Window

from peatguard.analysis.contrast_validation import spatial_group_diagnostics


def fingerprint(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"path": str(Path(path).resolve()), "sha256": digest.hexdigest()}


def sample_common_support(descending, ascending, features, *, half_window_m=500,
                          min_pixels=10, min_fraction=0.5):
    if half_window_m <= 0 or min_pixels < 1 or not 0 < min_fraction <= 1:
        raise ValueError("invalid window coverage parameters")
    with rasterio.open(descending) as d, rasterio.open(ascending) as a:
        if not d.crs or not a.crs or not d.crs.is_projected:
            raise ValueError("rasters need a CRS; descending grid must be projected")
        if d.crs.linear_units_factor[1] != 1 or d.transform.b or d.transform.d:
            raise ValueError("descending grid must be unrotated and use metres")
        hx, hy = (int(math.ceil(half_window_m / abs(v)))
                  for v in (d.transform.a, d.transform.e))
        width, height = 2 * hx + 1, 2 * hy + 1
        rows = []
        with WarpedVRT(a, crs=d.crs, transform=d.transform, width=d.width, height=d.height,
                       resampling=Resampling.nearest, nodata=np.nan, dtype="float64") as av:
            for index, feature in enumerate(features):
                if feature["geometry"]["type"] != "Point":
                    raise ValueError("site geometries must be points")
                lon, lat = feature["geometry"]["coordinates"][:2]
                if not np.isfinite([lon, lat]).all() or not -180 <= lon <= 180 or not -90 <= lat <= 90:
                    raise ValueError("site coordinate outside longitude/latitude bounds")
                xs, ys = transform("EPSG:4326", d.crs, [lon], [lat])
                row, col = d.index(xs[0], ys[0])
                rec = {"site": feature["properties"].get("site_id", feature["properties"].get("id", index+1)),
                       "lon": lon, "lat": lat, "xy_metres": [xs[0], ys[0]],
                       "status": "outside_full_window", "paired_pixels": 0,
                       "window_pixels": width * height}
                rows.append(rec)
                if row-hy < 0 or col-hx < 0 or row+hy >= d.height or col+hx >= d.width:
                    continue
                window = Window(col-hx, row-hy, width, height)
                dv, aa = d.read(1, window=window, masked=True), av.read(1, window=window, masked=True)
                valid = (~np.ma.getmaskarray(dv) & ~np.ma.getmaskarray(aa)
                         & np.isfinite(dv.data) & np.isfinite(aa.data))
                rec["paired_pixels"] = int(valid.sum())
                rec["paired_fraction"] = float(valid.mean())
                rec["status"] = "insufficient_common_support"
                if valid.sum() < min_pixels or valid.mean() < min_fraction:
                    continue
                rec.update(status="paired_diagnostic", descending_mm=float(np.median(dv.data[valid])),
                           ascending_mm=float(np.median(aa.data[valid])))
        return rows, {"crs": str(d.crs), "window_width_m": width*abs(d.transform.a),
                      "window_height_m": height*abs(d.transform.e), "resampling": "nearest"}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--descending", type=Path, required=True)
    parser.add_argument("--ascending", type=Path, required=True)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--block-km", type=float, nargs="+", default=[3, 5, 10])
    parser.add_argument("--half-window-m", type=float, default=500)
    parser.add_argument("--draws", type=int, default=2000)
    args = parser.parse_args(argv)
    if args.output.exists() or args.output.resolve() in {p.resolve() for p in (args.sites, args.ascending, args.descending)}:
        parser.error("output must be a new path, distinct from every input")
    if any(not math.isfinite(k) or k <= 0 for k in args.block_km):
        parser.error("block sizes must be positive and finite")
    features = json.loads(args.sites.read_text())["features"]
    rows, grid = sample_common_support(args.descending, args.ascending, features,
                                      half_window_m=args.half_window_m)
    valid = [r for r in rows if r["status"] == "paired_diagnostic"]
    trials = []
    if len(valid) >= 2:
        for km in args.block_km:
            trials.append(spatial_group_diagnostics(
                np.array([r["descending_mm"] for r in valid]),
                np.array([r["ascending_mm"] for r in valid]),
                np.array([r["xy_metres"] for r in valid]),
                block_metres=km*1000, draws=args.draws))
    report = {
        "schema_version": 1, "observable": "provisional_filtered_local_contrast_mm",
        "release_ready": False, "annual_velocity_ready": False,
        "inputs": {key: fingerprint(getattr(args, key)) for key in ("ascending", "descending", "sites")},
        "spatial_support": grid, "site_count": len(rows), "paired_sites": len(valid),
        "trials": trials, "sites": rows,
        "limitations": [
            "Source rasters lack verified acquisition lists, correction lineage, and per-pair component linkage.",
            "Filtered contrasts remove long spatial wavelengths and cannot recover frame-wide displacement.",
            "Cross-geometry agreement is consistency evidence; LOS geometry and time windows may differ.",
            "Sites were not reserved before filter selection; this is retrospective, not held-out validation.",
            "Spatial groups reduce pseudoreplication but do not establish independence; inspect block-size sensitivity.",
            "Requires a subsequent acquisition and independent withheld stable controls before operational use.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(f"Paired sites: {len(valid)}/{len(rows)}; diagnostic report: {args.output}")
    for trial in trials:
        print(f"{trial['block_metres']/1000:g} km: groups={trial['groups']}, "
              f"r={trial['metrics']['pearson_r']}, CI={trial['pearson_ci95']}, {trial['status']}")
    return 0 if trials else 2


if __name__ == "__main__":
    raise SystemExit(main())
