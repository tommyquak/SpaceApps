"""Stage 1: strip window in radar coordinates and the radar-coordinate geometry.

The 80x448 latitude/longitude grids written by the ionosphere stage (163 x 120) are
interpolated to the 17x56 grid (770 x 960); the window is the bounding box of the
pixels inside the strip; incidence, azimuth, height and slant range are sampled from
the existing 90 m geocoded geometry at each radar pixel's lat/lon (nearest cell).
"""
import json
import os

import h5py
import numpy as np
from scipy.ndimage import map_coordinates

from strip_common import GEOGEO, ION, L_FULL, OUT, REFSTACK, STRIP, W_FULL

os.makedirs(OUT, exist_ok=True)
lat8 = np.fromfile(f"{ION}/lat_80rlks_448alks.lat", dtype=np.float64).reshape(120, 163)
lon8 = np.fromfile(f"{ION}/lon_80rlks_448alks.lon", dtype=np.float64).reshape(120, 163)
valid8 = (lat8 != 0) & (lon8 != 0)
# fill the zero (no-data) cells by nearest valid value so the interpolation does not bleed zeros
from scipy.ndimage import distance_transform_edt
idx = distance_transform_edt(~valid8, return_distances=False, return_indices=True)
lat8f, lon8f = lat8[tuple(idx)], lon8[tuple(idx)]
rows = (np.arange(L_FULL) + 0.5) * 56 / 448 - 0.5
cols = (np.arange(W_FULL) + 0.5) * 17 / 80 - 0.5
rr, cc = np.meshgrid(rows, cols, indexing="ij")
lat = map_coordinates(lat8f, [rr, cc], order=1, mode="nearest")
lon = map_coordinates(lon8f, [rr, cc], order=1, mode="nearest")
valid = map_coordinates(valid8.astype(np.float32), [rr, cc], order=0, mode="nearest") > 0.5
inside = valid & (lat >= STRIP["lat_min"]) & (lat <= STRIP["lat_max"]) & (lon >= STRIP["lon_min"]) & (lon <= STRIP["lon_max"])
r, c = np.where(inside)
r0, r1, c0, c1 = int(r.min()), int(r.max()) + 1, int(c.min()), int(c.max()) + 1
win = dict(r0=r0, r1=r1, c0=c0, c1=c1, length=r1 - r0, width=c1 - c0, strip_px=int(inside.sum()))
json.dump(win, open(f"{OUT}/window.json", "w"), indent=1)
print("window", win)

with h5py.File(GEOGEO) as g:
    a = dict(g.attrs)
    x0, dx, y0, dy = float(a["X_FIRST"]), float(a["X_STEP"]), float(a["Y_FIRST"]), float(a["Y_STEP"])
    geo = {k: g[k][:] for k in ("incidenceAngle", "azimuthAngle", "height", "slantRangeDistance")}
gi = np.clip(((lat[r0:r1, c0:c1] - y0) / dy).astype(int), 0, geo["incidenceAngle"].shape[0] - 1)
gj = np.clip(((lon[r0:r1, c0:c1] - x0) / dx).astype(int), 0, geo["incidenceAngle"].shape[1] - 1)
with h5py.File(REFSTACK) as f:
    attrs = {k: v for k, v in f.attrs.items()}
for k in ("X_FIRST", "Y_FIRST", "X_STEP", "Y_STEP", "X_UNIT", "Y_UNIT", "REF_Y", "REF_X", "REF_LAT", "REF_LON", "EPSG", "FILE_TYPE", "UNIT"):
    attrs.pop(k, None)
attrs.update(LENGTH=str(r1 - r0), WIDTH=str(c1 - c0), FILE_TYPE="geometry", processing_note="strip experiment: radar coords, geometry sampled from the 90 m geocoded geometry at interpolated 80x448 lat/lon")
with h5py.File(f"{OUT}/geometryRadar.h5", "w") as g:
    g.create_dataset("latitude", data=lat[r0:r1, c0:c1].astype(np.float32))
    g.create_dataset("longitude", data=lon[r0:r1, c0:c1].astype(np.float32))
    for k in geo:
        v = geo[k][gi, gj].astype(np.float32); v[~np.isfinite(v)] = 0
        g.create_dataset(k, data=v)
    g.create_dataset("strip_mask", data=inside[r0:r1, c0:c1].astype(np.uint8))
    for k, v in attrs.items():
        g.attrs[k] = v
inc = geo["incidenceAngle"][gi, gj]
print("geometry written: incidence %.1f..%.1f, lat %.3f..%.3f, lon %.3f..%.3f" % (np.nanmin(inc), np.nanmax(inc), lat[r0:r1, c0:c1].min(), lat[r0:r1, c0:c1].max(), lon[r0:r1, c0:c1].min(), lon[r0:r1, c0:c1].max()))
print("STRIP_PREP_OK")
