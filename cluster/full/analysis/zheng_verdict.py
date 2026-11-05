"""Decide whether the Zheng closure-bias correction changed the spatial pattern
or only removed a near-uniform offset.

Correlation with an external benchmark is invariant under v -> a*v + b, so a
spatially uniform correction cannot move the Tay correlation at all. This
separates "Zheng fixed the bias number" from "Zheng improved the science".
"""
import os
import numpy as np
import h5py

W = os.environ.get("W", os.path.expanduser("~/full/mintpy_tay"))


def load(name, ds="velocity"):
    p = os.path.join(W, name)
    with h5py.File(p, "r") as f:
        if ds not in f:
            ds = list(f.keys())[0]
        return f[ds][:].astype(np.float64)


def mask_tempcoh(thresh):
    with h5py.File(os.path.join(W, "temporalCoherence.h5"), "r") as f:
        tc = f["temporalCoherence"][:].astype(np.float64)
    return tc >= thresh, tc


pre = load("velocity.h5") * 1000.0            # m/yr -> mm/yr
post = load("velocityClosureCor.h5") * 1000.0
corr_field = post - pre

for thr in (0.0, 0.5, 0.7):
    if thr == 0.0:
        m = np.isfinite(pre) & np.isfinite(post)
        label = "all finite"
    else:
        mk, tc = mask_tempcoh(thr)
        m = mk & np.isfinite(pre) & np.isfinite(post)
        label = f"tempCoh>={thr}"
    n = int(m.sum())
    if n < 10:
        print(f"[{label}] n={n} too few")
        continue
    a, b, c = pre[m], post[m], corr_field[m]
    r = np.corrcoef(a, b)[0, 1]
    print(f"[{label}] n={n}")
    print(f"   pre   mean {a.mean():+8.3f}  std {a.std():7.3f} mm/yr")
    print(f"   post  mean {b.mean():+8.3f}  std {b.std():7.3f} mm/yr")
    print(f"   delta mean {c.mean():+8.3f}  std {c.std():7.3f} mm/yr")
    print(f"   corr(pre,post) = {r:.6f}")
    # How much of the correction is a constant?
    frac_uniform = 1.0 - (c.std() ** 2) / max(c.mean() ** 2 + c.std() ** 2, 1e-30)
    print(f"   correction energy that is a pure offset: {100*frac_uniform:.1f}%")
    print()
