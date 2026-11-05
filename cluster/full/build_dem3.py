"""Build elevation_3s.dem.wgs84 (3 arcsec) from the 1-arcsec elevation.dem.wgs84.

alosStack geocode.py takes the output grid posting from the DEM it is given, so
the geocoding DEM must be 3 arcsec (~90 m, matching 5x28 looks); with the
1-arcsec DEM every geocoded product would come out 9x larger. The earlier
interactive attempt failed inside ISCE's stitcher ("the attribute corresponding
to the key reference is not present in Image": a DEM xml loaded as a plain
Image), so this builds the raster directly and renders the xml through
isceobj.createDemImage, which does carry the reference property.

The 1-arcsec grid is pixel-is-point (startingvalue 111.0, delta 1/3600, size
21601 = exactly 6 degrees inclusive), so the 3-arcsec grid keeps every third
point and averages the 3x3 neighbourhood around it: the kept points stay on the
exact same geodetic locations, with aliasing suppressed.

Run on a compute node (the source is 2.2 GB):
  srun -p normal -t 30 --mem=16G -c 4 python3 build_dem3.py <in> <out>
"""
import os
import sys

import numpy as np

STEP = 3


def main(src, dst):
    import xml.etree.ElementTree as ET

    import isce  # noqa: F401
    import isceobj

    # read the source metadata straight from the xml: ISCE's own parser returns
    # a flat property map that does not carry the coordinate components
    root = ET.parse(src + ".xml").getroot()
    prop = {}
    def walk(node, path=""):
        for p in node.findall("property"):
            v = p.find("value")
            prop[path + p.get("name").lower()] = v.text.strip() if v is not None else None
        for c in node.findall("component"):
            walk(c, path + c.get("name").lower() + ".")
    walk(root)
    w, l = int(prop["width"]), int(prop["length"])
    x0, dx = float(prop["coordinate1.startingvalue"]), float(prop["coordinate1.delta"])
    y0, dy = float(prop["coordinate2.startingvalue"]), float(prop["coordinate2.delta"])
    dt = str(prop["data_type"]).upper()
    if dt != "FLOAT":
        raise SystemExit(f"build_dem3: source data_type is {dt}, expected FLOAT")
    if os.path.getsize(src) != w * l * 4:
        raise SystemExit(f"build_dem3: {src} is not {w}x{l} float32")
    # keep every third point: indices 0, 3, ... <= w-1
    w3, l3 = (w - 1) // STEP + 1, (l - 1) // STEP + 1
    print(f"source {w}x{l} at {dx:.9f} deg, origin {x0}/{y0} -> {w3}x{l3} at {dx * STEP:.9f} deg")
    a = np.memmap(src, dtype="float32", mode="r", shape=(l, w))
    out = np.empty((l3, w3), dtype="float32")
    # 3x3 mean centred on each kept sample, clipped at the edges, in row blocks
    for i3 in range(l3):
        i = i3 * STEP
        r0, r1 = max(i - 1, 0), min(i + 2, l)
        block = np.asarray(a[r0:r1], dtype="float64")
        acc = np.zeros(w3, dtype="float64")
        n = np.zeros(w3, dtype="float64")
        cols = np.arange(w3) * STEP
        for off in (-1, 0, 1):
            c = cols + off
            ok = (c >= 0) & (c < w)
            acc[ok] += block[:, c[ok]].sum(axis=0)
            n[ok] += block.shape[0]
        out[i3] = (acc / n).astype("float32")
        if i3 % 1000 == 0:
            print(f"  row {i3}/{l3}")
    del a
    tmp = dst + ".tmp"
    out.tofile(tmp)
    os.replace(tmp, dst)
    print(f"wrote {dst}: {out.shape[1]}x{out.shape[0]}, "
          f"{out.min():.1f}..{out.max():.1f} m, {os.path.getsize(dst)} B")

    d, base = os.path.split(os.path.abspath(dst))
    img = isceobj.createDemImage()
    img.setFilename(base)
    img.setWidth(w3)
    img.setLength(l3)
    img.setAccessMode("READ")
    img.dataType = "FLOAT"
    img.scheme = str(prop["scheme"]).upper()
    img.bands = 1
    img.reference = prop.get("reference") or "WGS84"
    img.firstLatitude = y0
    img.firstLongitude = x0
    img.deltaLatitude = dy * STEP
    img.deltaLongitude = dx * STEP
    img.coord1.coordStart = x0
    img.coord1.coordDelta = dx * STEP
    img.coord1.coordSize = w3
    img.coord2.coordStart = y0
    img.coord2.coordDelta = dy * STEP
    img.coord2.coordSize = l3
    cwd = os.getcwd()
    try:
        os.chdir(d)
        img.renderHdr()
    finally:
        os.chdir(cwd)
    print("rendered", dst + ".xml")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2])
