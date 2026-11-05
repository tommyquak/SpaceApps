"""Regenerate a cropped SLC's lower/upper subbands for the ion pass (spec sec 3).

ion_subband.py reads dates_resampled/<date>/f1_*/s*/<date>_{lower,upper}.slc but
does not create them: its rg_filter call is commented out upstream, because in
the stock flow resample_common_grid -subband wrote them. The fleet resamples one
band and crops it, so the ion block job regenerates the subbands from the crop
and restores the demodulation ramp.

PHYSICS (Amendment A.7). rg_filter demodulates sample j of the file it is handed
by exp(-i 2 pi bc j). Upstream filtered each date on its ORIGINAL full swath, so
the authentic subband carries the ramp of j = rg - 1 + off_d(az, rg). Filtering
the crop instead gives the ramp of j = k, the cropped 0-based sample. The ratio
is a per-pixel factor exp(-i 2 pi bc (off_d + sampOff_s)), so each regenerated
subband is multiplied by it. sampOff_s is a per-swath constant that cancels in
every interferogram; off_d does not.

  ion_regen.py <cropped.slc> <sidecar.json> <resampled_ref_dir> <ref_date>

Round-14 finding 20 stands and is bounded, not fixed: rg_filter zero-pads
outside the file, so the outermost hnfilter = 128 samples at a crop edge that is
not also an original swath edge carry a filter transient no ramp can restore.
For the derived window the only such edges are s2's left and s4's right, and
both coincide with the frame crop boundary.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import crop_pack as cp

NFILTER, NFFT, BETA = 257, 2048, 0.1
BLOCK = 4096


def apply_ramp(path, sc, bc):
    """Multiply in place by exp(-i 2 pi bc (off_rg(az, rg) + sampOff)) at the
    sidecar's evaluation coordinates: az = lineOff + i + 1, rg = sampOff + k + 1."""
    w, l = int(sc["width"]), int(sc["length"])
    if os.path.getsize(path) != w * l * 8:
        raise SystemExit(f"ion_regen: {path} is not {w}x{l} complex64")
    poly = sc["offset_polys"]["rg"]
    rg = (int(sc["sampOff"]) + np.arange(w) + 1)[None, :]
    peak = 0.0
    with open(path, "r+b") as f:
        for i0 in range(0, l, BLOCK):
            n = min(BLOCK, l - i0)
            f.seek(i0 * w * 8)
            a = np.fromfile(f, dtype=np.complex64, count=n * w).reshape(n, w)
            az = (int(sc["lineOff"]) + np.arange(i0, i0 + n) + 1)[:, None]
            off = cp.eval_poly(poly, az, rg) + int(sc["sampOff"])
            peak = max(peak, float(np.abs(off).max()))
            a *= np.exp(-2j * np.pi * bc * off).astype(np.complex64)
            f.seek(i0 * w * 8)
            a.tofile(f)
    return peak


def main(slc, sidecar, ref_dir, ref_date):
    sc = json.load(open(sidecar))
    for k in ("width", "length", "sampOff", "lineOff", "offset_polys"):
        if k not in sc:
            raise SystemExit(f"ion_regen: {sidecar} has no {k}")
    loadTrack, subbandParameters = cp._stackpulic("loadTrack", "subbandParameters")
    track = loadTrack(os.path.abspath(ref_dir), ref_date)
    wl, bw, fc, prefix = subbandParameters(track)
    if len(prefix) != 2:
        raise SystemExit(f"ion_regen: expected two subbands, got {prefix}")
    d, base = os.path.split(os.path.abspath(slc))
    if not base.endswith(".slc"):
        raise SystemExit(f"ion_regen: {slc} is not a .slc")
    outs = [base[:-4] + "_" + p + ".slc" for p in prefix]
    from contrib.alos2proc.alos2proc import rg_filter
    cwd = os.getcwd()
    try:
        os.chdir(d)                      # rg_filter writes beside its input
        rg_filter(base, 2, outs, bw, fc, NFILTER, NFFT, BETA, 0, 0.0)
    finally:
        os.chdir(cwd)
    for k, o in enumerate(outs):
        p = os.path.join(d, o)
        if os.path.getsize(p) != int(sc["width"]) * int(sc["length"]) * 8:
            raise SystemExit(f"ion_regen: {o} is not the crop size")
        peak = apply_ramp(p, sc, float(fc[k]))
        print(f"  {o}: bc {fc[k]:+.6f} cycles/sample, max |off + sampOff| {peak:.1f} samples")
    print(f"ion_regen ok for {base} ({prefix[0]}, {prefix[1]})")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    main(*sys.argv[1:5])
