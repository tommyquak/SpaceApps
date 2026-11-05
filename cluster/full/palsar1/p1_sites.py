"""Per-site 2007-2009 rates at the ranked dam sites from the PALSAR-1 mosaic.

Usage: p1_sites.py <mosaic_dir> <dam_sites.geojson> <out_csv>

For every site: LOS and vertical rate as the median inside a 500 m window
(11 x 11 pixels at 90 m), the pixel spread inside the window (a local noise
proxy, not the formal uncertainty), the temporal coherence, the number of good
pixels, and the frame the value came from. Sites with fewer than 20 good pixels
are reported as no-data rather than filled.
"""
import csv
import json
import sys

import numpy as np
from osgeo import gdal

HALF = 5  # pixels: 11 x 11 window at 90 m posting ~ 1 km, matches the block scale of the product


def read(path):
    ds = gdal.Open(path)
    a = ds.GetRasterBand(1).ReadAsArray().astype(np.float64)
    nd = ds.GetRasterBand(1).GetNoDataValue()
    if nd is not None:
        a[a == nd] = np.nan
    return a, ds.GetGeoTransform()


def main():
    mdir, sites, out = sys.argv[1], sys.argv[2], sys.argv[3]
    los, gt = read(f"{mdir}/alos1_velocity_mmyr.tif")
    vert, _ = read(f"{mdir}/alos1_vertical_mmyr.tif")
    tc, _ = read(f"{mdir}/alos1_temporalCoherence.tif")
    fid, _ = read(f"{mdir}/alos1_frame_id.tif")
    feats = json.load(open(sites))["features"]
    rows = []
    for k, f in enumerate(feats):
        lon, lat = f["geometry"]["coordinates"][:2]
        p = f.get("properties", {})
        c = int((lon - gt[0]) / gt[1]); r = int((lat - gt[3]) / gt[5])
        r0, r1 = max(r - HALF, 0), min(r + HALF + 1, los.shape[0]); c0, c1 = max(c - HALF, 0), min(c + HALF + 1, los.shape[1])
        w = los[r0:r1, c0:c1]; wv = vert[r0:r1, c0:c1]; wt = tc[r0:r1, c0:c1]; wf = fid[r0:r1, c0:c1]
        ok = np.isfinite(w)
        n = int(ok.sum())
        row = {"site": p.get("id", p.get("name", k + 1)), "lon": round(lon, 5), "lat": round(lat, 5), "n_good_px": n}
        if n >= 20:
            row.update({"los_mmyr": round(float(np.median(w[ok])), 1), "vertical_mmyr": round(float(np.nanmedian(wv[ok])), 1),
                        "window_spread_mmyr": round(float(np.std(w[ok])), 1), "tempcoh": round(float(np.nanmedian(wt[ok])), 2),
                        "frame_id": int(np.nanmedian(wf[ok]))})
        else:
            row.update({"los_mmyr": "", "vertical_mmyr": "", "window_spread_mmyr": "", "tempcoh": "", "frame_id": ""})
        for key in ("rank", "score", "priority", "risk"):
            if key in p:
                row[key] = p[key]
        rows.append(row)
    keys = list(rows[0].keys())
    with open(out, "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=keys); wtr.writeheader(); wtr.writerows(rows)
    have = [r for r in rows if r["los_mmyr"] != ""]
    v = np.array([r["vertical_mmyr"] for r in have], float)
    print(f"P1_SITES {len(have)}/{len(rows)} sites with data; vertical rate median {np.median(v):+.1f} mm/yr, "
          f"p10/p90 {np.percentile(v, 10):+.1f}/{np.percentile(v, 90):+.1f}; {int((v < 0).sum())} subsiding, {int((v >= 0).sum())} not")


if __name__ == "__main__":
    main()
