"""Stage 3 (variant A'): pad one unwrapped crop back to the full 770 x 960 frame and
geocode it with alosStack geocode.py exactly as main.sbatch does (unw, cor, conncomp),
then translate to GeoTIFFs under the finale names build_stack2 globs for."""
import os
import subprocess
import sys

import numpy as np
import isce  # noqa: F401
import isceobj

from strip_common import L_FULL, ROOT, W_FULL, cal_to_p, window

cal = sys.argv[1]; track = sys.argv[2]; dem3 = sys.argv[3]; bbox = sys.argv[4]
w = window(); p = cal_to_p(cal); ML = "_17rlks_56alks"
od = f"{ROOT}/pairsA/{cal}"
if os.path.exists(f"{od}/DONE"):
    print("SKIP", cal); sys.exit(0)
os.makedirs(od, exist_ok=True); os.chdir(od)
unw = np.zeros((L_FULL, 2, W_FULL), np.float32); coh = np.zeros((L_FULL, 2, W_FULL), np.float32); cc = np.zeros((L_FULL, W_FULL), np.uint8)
u = np.load(f"{ROOT}/pairs/{cal}/unw.npy"); c = np.load(f"{ROOT}/pairs/{cal}/coh.npy"); k = np.load(f"{ROOT}/pairs/{cal}/cc.npy")
unw[w["r0"]:w["r1"], 1, w["c0"]:w["c1"]] = u; unw[w["r0"]:w["r1"], 0, w["c0"]:w["c1"]] = 1.0
coh[w["r0"]:w["r1"], 1, w["c0"]:w["c1"]] = c; coh[w["r0"]:w["r1"], 0, w["c0"]:w["c1"]] = 1.0
cc[w["r0"]:w["r1"], w["c0"]:w["c1"]] = k
names = {f"filt_{p}{ML}.unw": (unw, 2, "FLOAT", "unw"), f"{p}{ML}.cor": (coh, 2, "FLOAT", "cor"), f"filt_{p}{ML}.unw.conncomp": (cc, 1, "BYTE", "byte")}
for fn, (arr, nb, dt, kind) in names.items():
    arr.tofile(fn)
    img = isceobj.createImage(); img.setFilename(fn); img.setWidth(W_FULL); img.setLength(L_FULL); img.setBands(nb); img.setDataType(dt); img.setInterleavedScheme("BIL"); img.setAccessMode("READ")
    if kind == "unw": img.setImageType("unw")
    img.renderHdr(); img.renderVRT()
    r = subprocess.run(["geocode.py", "-ref_date_stack_track", track, "-dem", dem3, "-input", fn, "-nrlks", "17", "-nalks", "56", "-bbox", " " + bbox], capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(fn + ".geo.vrt"):
        print("GEOCODE_FAIL", cal, fn, r.stderr[-300:]); sys.exit(1)
conv = {f"filt_{p}{ML}.unw.geo": f"filt_{cal}{ML}.unw.geo.aoi90.tif", f"{p}{ML}.cor.geo": f"{cal}{ML}.cor.geo.aoi90.tif", f"filt_{p}{ML}.unw.conncomp.geo": f"filt_{cal}{ML}.unw.conncomp.geo.aoi90.tif"}
for src, dst in conv.items():
    subprocess.run(["gdal_translate", "-q", "-of", "GTiff", "-co", "COMPRESS=DEFLATE", src + ".vrt", dst], check=True)
for f in os.listdir("."):
    if not f.endswith(".aoi90.tif"): os.remove(f)
open("DONE", "w").write(""); print("GEOA_OK", cal)
