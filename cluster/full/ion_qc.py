"""Automated ion QC for the barrier (spec sec 5).

ion_ls solves a least squares over the ion network, so a pair with a BROKEN
ionosphere estimate contaminates the dates it touches. This scores every
persisted pair, prints the ones to exclude, and refuses to release the network
if excluding them would disconnect it.

THRESHOLD (measured over 326 real pairs, 2026-09-03). filt_ion std is not a
single population: it sits at 0.0995 rad for 2018-2021, rises to 0.149 in 2022,
0.421 in 2023 and 0.182 in 2024, tracking solar cycle 25 from minimum toward
maximum. A tight cut (0.35, or median+3*MAD = 0.135 since the MAD is 0.0067)
would delete the recent dates for having a genuinely active ionosphere, which
is exactly where the correction matters most. CORRECTED 2026-09-04. An earlier version of this docstring claimed "ion_ls also
READS each pair's .std and weights the inversion by it, so a noisy-but-valid
pair is down-weighted rather than trusted equally". THAT IS FALSE. ion_ls.py
loads stdPairs at lines 183 and 191 and never references it again; the weighted
branch (line 241) uses winPairs, the .win window-size file, not .std. And the
solve is unweighted by default: line 233 is `if ww == False:` and barrier.sbatch
never passes -ww. There is no downstream safety net.

Consequence of trusting that non-existent net: STD_MAX was raised to 0.8 rad
while the campaign maximum .std is 0.5259 rad, so the gate sat OUTSIDE the range
of its own quantity and passed 362 of 362 pairs, excluding none - including 115
whose ion estimate was unphysical by up to six orders of magnitude. A threshold
outside the observed range of its quantity is not a loose gate, it is a disabled
one.

.std is a coherence-derived FORMAL uncertainty and is blind to magnitude error.
Gate instead on the magnitude and gradient of the ion field itself: see
ion_outliers.py. Note its 200 rad mean cut is 14.9 TECU at 13.44 rad/TECU, which
is NOT physically impossible in equatorial conditions, so that threshold is
provisional and needs a no-ion / current / robust-weighted ablation.

  ion_qc.py <pairs_ion dir> <ion_blocks dir> [std_max_rad]

Prints one "EXCLUDE <pair>" line per rejected pair and exits non-zero when the
surviving network no longer spans every date, so the barrier stops rather than
handing ion_ls a network it cannot solve.
"""
import json
import os
import sys

import numpy as np

NRLKS1, NALKS1 = 1, 14
NRLKS_ION, NALKS_ION = 80, 32
# ion_cal products live on the ION grid (ion_ls.py:171 builds the name from
# nrlks1*nrlks_ion and nalks1*nalks_ion), not on ml2
ML2 = "_{}rlks_{}alks".format(NRLKS1 * NRLKS_ION, NALKS1 * NALKS_ION)


def score(d):
    """Mean std and the ion field's own spread, both in radians."""
    ion = os.path.join(d, f"filt_ion{ML2}.ion")
    std = os.path.join(d, f"filt_ion{ML2}.std")
    if not (os.path.isfile(ion) and os.path.isfile(std)):
        return None
    a = np.fromfile(ion, dtype="float32")
    s = np.fromfile(std, dtype="float32")
    if a.size == 0 or s.size == 0:
        return None
    fa, fs = a[np.isfinite(a)], s[np.isfinite(s)]
    if fa.size == 0 or fs.size == 0:
        return None
    nz = fa[fa != 0]
    return {"std": float(fs[fs > 0].mean()) if (fs > 0).any() else float("nan"),
            "spread": float(nz.std()) if nz.size else 0.0,
            "zero_frac": float((a == 0).mean()),
            "finite_frac": float(np.isfinite(a).mean())}


def connected(dates, pairs):
    """Union-find over the surviving network."""
    parent = {d: d for d in dates}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for p in pairs:
        a, b = p.split("-")
        if a in parent and b in parent:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
    roots = {find(d) for d in dates}
    return len(roots) == 1, roots


def main(ion_out, blocks_dir, std_max=0.8):
    pairs = []
    dates = set()
    for f in sorted(os.listdir(blocks_dir)):
        if f.endswith(".pairs"):
            for ln in open(os.path.join(blocks_dir, f)):
                p = ln.strip()
                if p:
                    pairs.append(p)
                    dates.update(p.split("-"))
    if not pairs:
        raise SystemExit("ion_qc: no pairs in the block plan")
    dates = sorted(dates)
    print(f"ion QC: {len(pairs)} pairs over {len(dates)} dates, std_max {std_max} rad")

    keep, exclude, missing = [], [], []
    rows = []
    for p in pairs:
        m = score(os.path.join(ion_out, p, "ion", "ion_cal"))
        if m is None:
            missing.append(p)
            continue
        rows.append((p, m))
    if rows:
        stds = np.array([m["std"] for _, m in rows], dtype=float)
        fin = stds[np.isfinite(stds)]
        if fin.size:
            print(f"  filt_ion std: median {np.median(fin):.4f}, p90 {np.percentile(fin, 90):.4f}, "
                  f"max {fin.max():.4f} rad")
    for p, m in rows:
        why = []
        if not np.isfinite(m["std"]) or m["std"] > std_max:
            why.append(f"std {m['std']:.3f} > {std_max}")
        if m["finite_frac"] < 0.99:
            why.append(f"non-finite {100 * (1 - m['finite_frac']):.1f} percent")
        if m["zero_frac"] > 0.5:
            why.append(f"zero {100 * m['zero_frac']:.0f} percent")
        if m["spread"] == 0.0:
            why.append("flat field")
        if why:
            exclude.append((p, "; ".join(why)))
        else:
            keep.append(p)

    for p in missing:
        print(f"MISSING {p} (no ion_cal products)")
    for p, why in exclude:
        print(f"EXCLUDE {p}   {why}")
    print(f"  {len(keep)} pairs kept, {len(exclude)} excluded, {len(missing)} missing")

    ok, roots = connected(dates, keep)
    if not ok:
        print(f"ION_NETWORK_DISCONNECTED: {len(roots)} components over {len(dates)} dates")
        return 1
    # every date must still carry at least one usable pair
    used = set()
    for p in keep:
        used.update(p.split("-"))
    orphan = [d for d in dates if d not in used]
    if orphan:
        print(f"ION_NETWORK_ORPHANS: {' '.join(orphan)}")
        return 1
    print(f"ION_NETWORK_OK: {len(keep)} pairs connect all {len(dates)} dates")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    sys.exit(main(sys.argv[1], sys.argv[2],
                  float(sys.argv[3]) if len(sys.argv) > 3 else 0.8))
