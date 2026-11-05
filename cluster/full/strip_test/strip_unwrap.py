"""Stage 2: crop one filtered interferogram to the strip window and unwrap it (snaphu)."""
import os
import subprocess
import sys

import numpy as np

from strip_common import ROOT, cal_to_p, phase_consistency, read_int, window

cal = sys.argv[1]
w = window()
od = f"{ROOT}/pairs/{cal}"
if os.path.exists(f"{od}/DONE"):
    print("SKIP", cal); sys.exit(0)
os.makedirs(od, exist_ok=True)
z = read_int(cal)[w["r0"]:w["r1"], w["c0"]:w["c1"]].astype(np.complex64)
coh = phase_consistency(z)
nl, nw = z.shape
z.tofile(f"{od}/crop.int"); coh.tofile(f"{od}/crop.coh")
conf = f"""INFILE {od}/crop.int
LINELENGTH {nw}
INFILEFORMAT COMPLEX_DATA
OUTFILE {od}/crop.unw
OUTFILEFORMAT ALT_LINE_DATA
CORRFILE {od}/crop.coh
CORRFILEFORMAT FLOAT_DATA
STATCOSTMODE SMOOTH
INITMETHOD MCF
NCORRLOOKS 41
CONNCOMPFILE {od}/crop.conncomp
CONNCOMPOUTTYPE UCHAR
"""
open(f"{od}/snaphu.conf", "w").write(conf)
r = subprocess.run(["snaphu", "-f", f"{od}/snaphu.conf"], capture_output=True, text=True)
if r.returncode != 0 or not os.path.exists(f"{od}/crop.unw"):
    print("SNAPHU_FAIL", cal, r.stderr[-400:]); sys.exit(1)
unw = np.fromfile(f"{od}/crop.unw", dtype=np.float32).reshape(nl, 2, nw)[:, 1, :]
cc = np.fromfile(f"{od}/crop.conncomp", dtype=np.uint8).reshape(nl, nw)
np.save(f"{od}/unw.npy", unw.astype(np.float32)); np.save(f"{od}/coh.npy", coh); np.save(f"{od}/cc.npy", cc)
for f in ("crop.int", "crop.unw", "crop.coh", "crop.conncomp", "snaphu.conf"):
    os.remove(f"{od}/{f}")
open(f"{od}/DONE", "w").write("")
print(f"UNW_OK {cal} comp0 {float((cc == 0).mean()):.2f} largest {int(np.bincount(cc.ravel())[1:].max()) if cc.max() > 0 else 0}")
