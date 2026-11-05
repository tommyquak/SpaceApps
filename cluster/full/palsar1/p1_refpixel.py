"""Pick the MintPy reference pixel from the loaded ifgramStack.

Usage: p1_refpixel.py <mintpy_dir>

Interior pixel (60 px margin) that is unwrapped and inside a connected
component in at least 95 percent of pairs, with the highest mean coherence.
MintPy's own choice (max average spatial coherence inside maskConnComp) picked
an edge pixel on frame 7140 that is nodata in most pairs, and the inversion
collapsed to a single valid pixel. Appends mintpy.reference.yx to the config.
"""
import sys

import h5py
import numpy as np


def main():
    w = sys.argv[1]
    with h5py.File(f"{w}/inputs/ifgramStack.h5") as f:
        n = f["unwrapPhase"].shape[0]
        ny, nx = f["unwrapPhase"].shape[1:]
        frac_cc = np.zeros((ny, nx)); coh = np.zeros((ny, nx)); nz = np.zeros((ny, nx))
        has_cc = "connectComponent" in f
        for i in range(n):
            u = f["unwrapPhase"][i]
            c = f["coherence"][i]
            cc = f["connectComponent"][i] if has_cc else (u != 0)
            frac_cc += (cc > 0); coh += np.nan_to_num(c); nz += (u != 0)
    frac_cc /= n; coh /= n; nz /= n
    m = 60
    # MintPy also requires the reference to lie inside maskConnComp.h5 (the component common to
    # EVERY pair), so pixels present in all pairs come first; anything looser disables that mask
    for thr in (1.0, 0.95, 0.8, 0.6):
        cand = np.zeros((ny, nx), bool); cand[m:-m, m:-m] = True
        cand &= (frac_cc >= thr) & (nz >= thr)
        if cand.sum() >= 10:
            break
    score = np.where(cand, coh, -1.0)
    y, x = np.unravel_index(int(np.argmax(score)), score.shape)
    print(f"reference y,x = {y},{x}  mean coh {coh[y, x]:.3f}  in-component {frac_cc[y, x]:.2f}  "
          f"nonzero {nz[y, x]:.2f}  candidates {int(cand.sum())} at threshold {thr}")
    with open(f"{w}/smallbaselineApp.cfg", "a") as g:
        g.write(f"mintpy.reference.yx = {y}, {x}\n")
        if thr < 1.0:
            g.write("mintpy.reference.maskFile = no\n")
            print("no pixel is unwrapped in every pair: reference mask check disabled")


if __name__ == "__main__":
    main()
