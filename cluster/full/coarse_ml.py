"""Coarse-posting interferogram from the surviving 17x56 filtered wrapped phase.

Usage: coarse_ml.py <filt_int_17x56> <out_filt_prefix> <out_cor_prefix> <rl> <al> <ncorrlooks>

The ALOS-2 coregistered SLC stack was lost on 2026-09-06; what survives per
pair is the ion-corrected, Goldstein-filtered wrapped interferogram at 17x56
looks (770 x 960, radar coordinates). This script complex-averages it by
rl x al (10 x 8 -> 170 x 448 looks, 77 x 120 px), writes ISCE products
(<filt>.int, <cor>.cor as amplitude+coherence BIL, <filt>.unw, <filt>.unw.conncomp)
and unwraps with the standalone snaphu (statistical SMOOTH cost, MCF init).
The "coherence" is the phase consistency of the FILTERED phase inside the block
(|mean z| / mean |z|), an upward-biased proxy that is fine for snaphu costs and
MintPy weights but must not be read as true coherence.
"""
import os
import subprocess
import sys

import numpy as np
import isce  # noqa: F401  (sets up the ISCE component path)
import isceobj


def main():
    src, fpre, cpre, rl, al, ncorr = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
    img = isceobj.createIntImage(); img.load(src + ".xml")
    w, l = img.width, img.length
    a = np.fromfile(src, dtype=np.complex64).reshape(l, w)
    L2, W2 = (l // al) * al, (w // rl) * rl
    b = a[:L2, :W2].reshape(L2 // al, al, W2 // rl, rl)
    z = b.mean(axis=(1, 3)).astype(np.complex64)
    m = np.abs(b).mean(axis=(1, 3))
    coh = np.where(m > 0, np.abs(z) / np.maximum(m, 1e-12), 0).astype(np.float32)
    amp = np.abs(z).astype(np.float32)
    nl, nw = z.shape
    print(f"coarse {nw} x {nl} from {w} x {l} (rl {rl}, al {al}); valid {int((m > 0).sum())} px, median coh {np.median(coh[m > 0]):.3f}")

    z.tofile(fpre + ".int")
    o = isceobj.createIntImage(); o.setFilename(fpre + ".int"); o.setWidth(nw); o.setLength(nl); o.setAccessMode("READ"); o.renderHdr(); o.renderVRT()
    np.stack([amp, coh], axis=1).astype(np.float32).tofile(cpre + ".cor")   # BIL: amp line, coh line
    c = isceobj.createImage(); c.setFilename(cpre + ".cor"); c.setWidth(nw); c.setLength(nl); c.setBands(2)
    c.setDataType("FLOAT"); c.setInterleavedScheme("BIL"); c.setAccessMode("READ"); c.renderHdr(); c.renderVRT()
    coh.tofile(cpre + ".coh1")

    conf = f"""INFILE {fpre}.int
LINELENGTH {nw}
INFILEFORMAT COMPLEX_DATA
OUTFILE {fpre}.unw
OUTFILEFORMAT ALT_LINE_DATA
CORRFILE {cpre}.coh1
CORRFILEFORMAT FLOAT_DATA
STATCOSTMODE SMOOTH
INITMETHOD MCF
NCORRLOOKS {ncorr}
CONNCOMPFILE {fpre}.unw.conncomp
CONNCOMPOUTTYPE UCHAR
"""
    open(fpre + ".snaphu.conf", "w").write(conf)
    r = subprocess.run(["snaphu", "-f", fpre + ".snaphu.conf"], capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(fpre + ".unw"):
        print("SNAPHU_FAIL", r.stdout[-800:], r.stderr[-800:]); sys.exit(1)
    u = isceobj.createImage(); u.setFilename(fpre + ".unw"); u.setWidth(nw); u.setLength(nl); u.setBands(2)
    u.setDataType("FLOAT"); u.setInterleavedScheme("BIL"); u.setImageType("unw"); u.setAccessMode("READ"); u.renderHdr(); u.renderVRT()
    k = isceobj.createImage(); k.setFilename(fpre + ".unw.conncomp"); k.setWidth(nw); k.setLength(nl); k.setBands(1)
    k.setDataType("BYTE"); k.setInterleavedScheme("BIL"); k.setAccessMode("READ"); k.renderHdr(); k.renderVRT()
    unw = np.fromfile(fpre + ".unw", dtype=np.float32).reshape(nl, 2, nw)[:, 1, :]
    cc = np.fromfile(fpre + ".unw.conncomp", dtype=np.uint8).reshape(nl, nw)
    if cc.max() > 0:
        sel = unw[cc > 0]
        print(f"unwrapped: phase range {sel.min():+.1f}..{sel.max():+.1f} rad, largest component {int(np.bincount(cc.ravel())[1:].max())} px")
    else:
        print("unwrapped: NO connected component (fully decorrelated pair); products written, MintPy will drop it")
    print("COARSE_OK")


if __name__ == "__main__":
    main()
