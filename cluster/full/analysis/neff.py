"""Independent measurement of the effective number of looks from speckle.

In a homogeneous, speckle-dominated patch the equivalent number of looks is
ENL = mean(I)^2 / var(I), I = amplitude^2. [E] Oliver & Quegan 2004. Texture
biases ENL DOWN, so upper quantiles across many patches are the fairer read.
The prior 41-look figure came from closure statistics, which count any real
closure signal as noise and can only bias LOW; this is the cross-check.

Amplitude sources, in priority order: a .amp product; band 1 of a .unw (ISCE
writes multilooked amplitude there); |z| of a .int (product of two amplitudes,
so its ENL is NOT the single-image ENL and is reported only as a fallback).
"""
import glob
import os
import numpy as np
from osgeo import gdal

ROOT = os.path.expanduser(os.environ.get("PAIR_ROOT", "~/full/pairs_tay"))
DATES = os.path.expanduser(os.environ.get("DATES_DIR", "/mnt/scratch/t/tommy/full/dates_resampled"))
PATCH = int(os.environ.get("PATCH", "24"))
gdal.UseExceptions()


def open_any(path):
    for p in (path + ".vrt", path):
        try:
            ds = gdal.Open(p)
            if ds is not None:
                return ds
        except Exception:
            pass
    return None


def candidates(tag):
    pats = [
        (os.path.join(ROOT, "*", "insar", f"*{tag}*.amp"), "amp"),
        (os.path.join(ROOT, "*", f"*{tag}*.amp"), "amp"),
        (os.path.join(DATES, "*", f"*{tag}*.amp"), "amp"),
        (os.path.join(ROOT, "*", "insar", f"*{tag}*.unw"), "unw"),
        (os.path.join(ROOT, "*", "insar", f"filt_*{tag}*.unw"), "unw"),
        (os.path.join(ROOT, "*", "insar", f"*{tag}*.int"), "int"),
    ]
    for pat, kind in pats:
        hits = sorted(glob.glob(pat))
        for h in hits[:3]:
            if open_any(h) is not None:
                return h, kind
    return None, None


def amplitude_bands(path, kind):
    ds = open_any(path)
    n = ds.RasterCount
    if kind == "amp":
        return [ds.GetRasterBand(i + 1).ReadAsArray().astype(np.float64) for i in range(n)], "amp bands"
    if kind == "unw":
        return [ds.GetRasterBand(1).ReadAsArray().astype(np.float64)], "unw band 1 (amplitude)"
    z = ds.GetRasterBand(1).ReadAsArray()
    return [np.abs(z).astype(np.float64)], "|int| (product of two amplitudes, fallback)"


def enl_quantiles(amp, patch):
    I = amp * amp
    ny, nx = I.shape
    vals = []
    for r in range(0, ny - patch, patch):
        for c in range(0, nx - patch, patch):
            p = I[r:r + patch, c:c + patch]
            if (p <= 0).mean() > 0.02 or not np.all(np.isfinite(p)):
                continue
            m = p.mean(); v = p.var()
            if m <= 0 or v <= 0:
                continue
            vals.append(m * m / v)
    vals = np.array(vals)
    if vals.size == 0:
        return None
    return {q: float(np.percentile(vals, q)) for q in (10, 25, 50, 75, 90, 95)}, int(vals.size)


for tag, nominal in (("1rlks_14alks", 14), ("17rlks_56alks", 952)):
    path, kind = candidates(tag)
    print(f"\n[{tag}] nominal looks {nominal}")
    if path is None:
        print("   no readable amplitude source found")
        continue
    bands, desc = amplitude_bands(path, kind)
    print(f"   source: {path}\n   using: {desc}; shape {bands[0].shape}")
    for bi, amp in enumerate(bands):
        r = enl_quantiles(amp, PATCH)
        if r is None:
            print(f"   band {bi}: no valid patches"); continue
        q, n = r
        print(f"   band {bi}: {n} patches of {PATCH}x{PATCH} px")
        print("      ENL p10 %.1f  p25 %.1f  p50 %.1f  p75 %.1f  p90 %.1f  p95 %.1f"
              % (q[10], q[25], q[50], q[75], q[90], q[95]))
        print(f"      p75/nominal {q[75]/nominal:.3f}   p90/nominal {q[90]/nominal:.3f}")
        if kind == "int":
            print("      NOTE: |int| ENL is for a product of two speckle fields; single-image ENL is roughly 2x this at low coherence")
print("\nNEFF_DONE")
