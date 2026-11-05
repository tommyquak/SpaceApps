"""Stage 3 (variant B): radar-coordinate ifgramStack from the unwrapped crops, then MintPy
without deramp, geocoding after the inversion. Dates, bperp and radar attributes come
from the existing 17x56 stack (same 1650 pairs)."""
import os
import shutil
import sys

import h5py
import numpy as np

from strip_common import OUT, REFSTACK, ROOT, pair_list, window

w = window()
W = f"{OUT}/mintpyB"; os.makedirs(f"{W}/inputs", exist_ok=True)
with h5py.File(REFSTACK) as f:
    dates = f["date"][:]; bperp = f["bperp"][:]; drop = f["dropIfgram"][:]; attrs = {k: v for k, v in f.attrs.items()}
    ref_cals = [f"20{a.decode()[2:]}"[:4] for a in dates[:, 0]]
    cal_of = ["%s-%s-%s_%s-%s-%s" % (a.decode()[:4], a.decode()[4:6], a.decode()[6:], b.decode()[:4], b.decode()[4:6], b.decode()[6:]) for a, b in dates]
have = set(pair_list())
keep = [i for i, c in enumerate(cal_of) if c in have and os.path.exists(f"{ROOT}/pairs/{c}/DONE")]
print("pairs in reference stack", len(cal_of), "with unwrapped crops", len(keep))
n, nl, nw = len(keep), w["length"], w["width"]
for k in ("X_FIRST", "Y_FIRST", "X_STEP", "Y_STEP", "X_UNIT", "Y_UNIT", "REF_Y", "REF_X", "REF_LAT", "REF_LON", "EPSG"):
    attrs.pop(k, None)
attrs.update(LENGTH=str(nl), WIDTH=str(nw), FILE_TYPE="ifgramStack", UNIT="radian",
             processing_note="strip experiment B: same filtered 17x56 interferograms, re-unwrapped in radar coords, inverted before geocoding, no deramp")
with h5py.File(f"{W}/inputs/ifgramStack.h5", "w") as g:
    du = g.create_dataset("unwrapPhase", (n, nl, nw), dtype=np.float32, chunks=(1, nl, nw), compression="lzf")
    dc = g.create_dataset("coherence", (n, nl, nw), dtype=np.float32, chunks=(1, nl, nw), compression="lzf")
    dk = g.create_dataset("connectComponent", (n, nl, nw), dtype=np.int16, chunks=(1, nl, nw), compression="lzf")
    for j, i in enumerate(keep):
        c = cal_of[i]
        du[j] = np.load(f"{ROOT}/pairs/{c}/unw.npy"); dc[j] = np.load(f"{ROOT}/pairs/{c}/coh.npy"); dk[j] = np.load(f"{ROOT}/pairs/{c}/cc.npy").astype(np.int16)
        if j % 200 == 0: print("  loaded", j, c, flush=True)
    g.create_dataset("date", data=dates[keep]); g.create_dataset("bperp", data=bperp[keep].astype(np.float32)); g.create_dataset("dropIfgram", data=np.ones(n, bool))
    for k, v in attrs.items(): g.attrs[k] = v
shutil.copy(f"{OUT}/geometryRadar.h5", f"{W}/inputs/geometryRadar.h5")
# reference pixel: interior, in a component in >= 95 percent of pairs, highest mean proxy coherence
with h5py.File(f"{W}/inputs/ifgramStack.h5") as g:
    frac = np.zeros((nl, nw)); coh = np.zeros((nl, nw))
    for j in range(n):
        frac += (g["connectComponent"][j] > 0); coh += g["coherence"][j]
    frac /= n; coh /= n
m = 30; cand = np.zeros((nl, nw), bool); cand[m:-m, m:-m] = True
print("in-component fraction over interior px: p50 %.2f p90 %.2f p99 %.2f max %.2f" % tuple(np.percentile(frac[cand], [50, 90, 99, 100])))
# the strip crop has no pixel inside a component in every pair (fragmented components):
# take the top 0.1 percent by component presence, then the highest proxy coherence,
# and switch MintPy's reference mask check off (same rule as pick_ref.py)
thr = float(np.percentile(frac[cand], 99.9)); c2 = cand & (frac >= thr)
score = np.where(c2, coh, -1); ry, rx = np.unravel_index(int(np.argmax(score)), score.shape)
print(f"reference y,x = {ry},{rx} in-component {frac[ry, rx]:.2f} proxy coh {coh[ry, rx]:.3f} (top 0.1 percent threshold {thr:.2f}, {int(c2.sum())} candidates)")
thr = frac[ry, rx]
cfg = f"""mintpy.load.processor = isce
mintpy.reference.yx = {ry}, {rx}
mintpy.reference.maskFile = {'auto' if thr >= 1.0 else 'no'}
mintpy.network.coherenceBased = no
mintpy.unwrapError.method = no
mintpy.networkInversion.weightFunc = var
mintpy.networkInversion.minTempCoh = 0.4
mintpy.networkInversion.maskDataset = no
mintpy.troposphericDelay.method = no
mintpy.deramp = no
mintpy.solidEarthTides = yes
mintpy.topographicResidual = no
mintpy.geocode = yes
mintpy.geocode.laloStep = 0.000833333, 0.000833333
mintpy.plot = no
"""
open(f"{W}/smallbaselineApp.cfg", "w").write(cfg)
print("STACKB_OK", n, "pairs")
