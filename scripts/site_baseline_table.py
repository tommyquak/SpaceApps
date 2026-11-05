"""Measured per-site baseline table for the 60 ranked dam sites.

Each column preserves its sensor and era. Windows are approximately 1 km,
with their odd pixel counts set by native posting. Spread is within-window
variability, not a measurement-error or datum-uncertainty estimate.

Columns:
  palsar1_2007_2009_vertical_mmyr   ALOS-1 historical mosaic using a published reference location
  tay_2018_2024_mmyr                 Tay et al. 2025 published benchmark rate raster
  (the ALOS-2 1.5 km block rate is deliberately NOT a column: r 0.12 vs Tay and an
   unreferenced datum make it misleading at site scale; it stays a coarse map layer)
  alos2_seasonal_amp_mm              exploratory fitted amplitude, residual screen is not validation
  nisar_D_60d_hp_mm, nisar_A_60d_hp_mm   NISAR provisional 60-day high-passed local contrast, Jun-Sep 2026
                                     (retrospective consistency evidence only)
Output: deliverables/site_baseline/dam_sites_measured_baseline.csv
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import transform as rtransform

SITES = Path("deliverables/tay_full_grid/dam_sites.geojson")
OUT = Path("deliverables/site_baseline")
L = {
    "palsar1_2007_2009_vertical_mmyr": ("grid_products/external/alos1_2007_2009/alos1_vertical_mmyr.tif", 1.0),
    "palsar1_tempcoh": ("grid_products/external/alos1_2007_2009/alos1_temporalCoherence.tif", 1.0),
    "tay_2018_2024_mmyr": ("deliverables/tay_restoration/subsidence_tay.tif", 1.0),
    "alos2_seasonal_amp_mm": ("grid_products/external/alos2_block_2018_2024/bw60_seasonal_amp_mm.tif", 1.0),
    "alos2_seasonal_resid_rms_mm": ("grid_products/external/alos2_block_2018_2024/bw60_fit_resid_rms_mm.tif", 1.0),
    "nisar_D_60d_hp_mm": ("grid_products/external/nisar/20260905/nisar_D_disp60_hp10km_mm.tif", 1.0),
    "nisar_A_60d_hp_mm": ("grid_products/external/nisar/20260905/nisar_A_disp60_hp10km_mm.tif", 1.0),
}


def sample(path, scale, lon, lat, half_m=500.0, min_px=10, *, zero_is_nodata=False):
    """Sample in-bounds windows, respecting the raster mask and valid zeros.

    Pixel spread describes within-window variability, not measurement accuracy.
    Zero padding must be explicitly declared for products that lack nodata metadata.
    """
    if (len(lon) != len(lat) or not np.isfinite([half_m, scale]).all()
            or half_m <= 0 or min_px < 1):
        raise ValueError("invalid coordinates or sampling parameters")
    with rasterio.open(path) as r:
        if r.crs is None or r.transform.b != 0 or r.transform.d != 0:
            raise ValueError("sampling requires a CRS and an unrotated raster")
        if r.crs.to_epsg() != 4326:
            xs, ys = rtransform("EPSG:4326", r.crs, lon, lat)
        else:
            xs, ys = lon, lat
        px = abs(r.transform.a); py = abs(r.transform.e)
        out = []
        for x, y, latitude in zip(xs, ys, lat):
            if not np.isfinite([x, y, latitude]).all():
                out.append((None, None, 0)); continue
            if r.crs.is_projected:
                metres_per_unit = r.crs.linear_units_factor[1]
                hx = int(round(half_m / (px * metres_per_unit)))
                hy = int(round(half_m / (py * metres_per_unit)))
            elif r.crs.to_epsg() == 4326:
                hx = int(round(half_m / (px * 111320 * np.cos(np.deg2rad(latitude)))))
                hy = int(round(half_m / (py * 111320)))
            else:
                raise ValueError("unsupported geographic CRS")
            rr, c = r.index(x, y)
            if not (0 <= rr < r.height and 0 <= c < r.width):
                out.append((None, None, 0)); continue
            top, left = max(rr - hy, 0), max(c - hx, 0)
            bottom, right = min(rr + hy + 1, r.height), min(c + hx + 1, r.width)
            w = r.read(1, window=rasterio.windows.Window(left, top, right-left, bottom-top),
                       masked=True).astype("float64")
            ok = ~np.ma.getmaskarray(w) & np.isfinite(w.data)
            if zero_is_nodata:
                ok &= w.data != 0
            if ok.sum() < min_px:
                out.append((None, None, int(ok.sum()))); continue
            v = w.data[ok] * scale
            out.append((float(np.median(v)), float(np.std(v)), int(ok.sum())))
        return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", type=Path, default=SITES)
    parser.add_argument("--output-dir", type=Path, default=OUT)
    args = parser.parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    feats = json.loads(args.sites.read_text())["features"]
    if not feats:
        raise ValueError("site collection is empty")
    lon = [f["geometry"]["coordinates"][0] for f in feats]; lat = [f["geometry"]["coordinates"][1] for f in feats]
    rows = [{"site": f["properties"].get("site_id", f["properties"].get("id", i + 1)), "lon": round(x, 5), "lat": round(y, 5), "risk": f["properties"].get("risk", "")}
            for i, (f, x, y) in enumerate(zip(feats, lon, lat))]
    identifiers = [str(row["site"]) for row in rows]
    if any(row["site"] is None or not str(row["site"]).strip() for row in rows):
        raise ValueError("site identifiers must be nonempty")
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("site identifiers must be unique")
    for name, (path, scale) in L.items():
        if not Path(path).exists():
            print("missing", path); continue
        vals = sample(path, scale, lon, lat)
        for row, (med, sd, n) in zip(rows, vals):
            row[name] = "" if med is None else round(med, 1)
            row[name + "_npx"] = n
            if name in ("palsar1_2007_2009_vertical_mmyr", "tay_2018_2024_mmyr", "nisar_D_60d_hp_mm"):
                row[name + "_spread"] = "" if sd is None else round(sd, 1)
    for row in rows:
        rr = row.get("alos2_seasonal_resid_rms_mm")
        # This is the legacy exploratory residual screen, not a validation test.
        row["alos2_seasonal_valid"] = "" if rr in (None, "") else int(float(rr) < 134)
        row["alos2_seasonal_screen_status"] = "exploratory_residual_screen_only"
        row["nisar_status"] = "provisional_filtered_contrast_not_rate"
        row["palsar1_status"] = "historical_relative_rate_reference_stability_not_reverified"
        got = [k for k in ("palsar1_2007_2009_vertical_mmyr", "tay_2018_2024_mmyr", "nisar_D_60d_hp_mm") if row.get(k, "") != ""]
        row["measured_layers"] = len(got)
    keys = list(rows[0].keys())
    with open(args.output_dir / "dam_sites_measured_baseline.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    def cnt(k): return sum(1 for r in rows if r.get(k, "") != "")
    print("sites with: palsar1", cnt("palsar1_2007_2009_vertical_mmyr"), "| tay", cnt("tay_2018_2024_mmyr"),
          "| seasonal", cnt("alos2_seasonal_amp_mm"), "| nisar D", cnt("nisar_D_60d_hp_mm"), "| nisar A", cnt("nisar_A_60d_hp_mm"))
    import statistics as st
    for k in ("palsar1_2007_2009_vertical_mmyr", "tay_2018_2024_mmyr", "alos2_seasonal_amp_mm", "nisar_D_60d_hp_mm"):
        v = [float(r[k]) for r in rows if r.get(k, "") != ""]
        if v: print(f"  {k:36s} n={len(v):2d} median {st.median(v):+7.1f} min {min(v):+7.1f} max {max(v):+7.1f}")
    print("SITE_TABLE_OK", args.output_dir / "dam_sites_measured_baseline.csv")


if __name__ == "__main__":
    main()
