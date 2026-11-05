"""Build the ION and MAIN date-block plans (spec sec 5).

Both passes are date-block jobs: a block stages a window of dates node-local
once, then runs every pair whose two dates are both in that window. A pair is
assigned to the FIRST block that contains both of its dates, so no pair is
processed twice.

  ION  (pre-barrier):  nsub_ion=4 network, blocks of 24 dates at stride 20
  MAIN (post-barrier): bandwidth-20 network, blocks of 32 dates at stride 12

Writes ~/full/ion_blocks/<i>.{dates,pairs} and ~/full/main_blocks/<i>.{dates,pairs}
plus a plan summary, and refuses to write a plan that would drop a pair.

  make_blocks.py <dates_file> [outdir]

dates_file: one "YYMMDD YYYY-MM-DD" per line (the wave-1 date list), or any
file whose first column is a YYMMDD date.
"""
import os
import sys

NSUB_ION = 4
ION_BLOCK, ION_STRIDE = 24, 20
MAIN_BANDWIDTH = 20
MAIN_BLOCK, MAIN_STRIDE = 32, 12


def blocks_of(n, size, stride):
    starts = list(range(0, max(n - size, 0) + 1, stride))
    if not starts or starts[-1] + size < n:
        starts.append(max(n - size, 0))
    return [(s, min(s + size, n)) for s in sorted(set(starts))]


def plan(dates, bandwidth, size, stride, label):
    n = len(dates)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, min(i + 1 + bandwidth, n))]
    bl = blocks_of(n, size, stride)
    assign = {}
    for k, (a, b) in enumerate(bl):
        for (i, j) in pairs:
            if (i, j) in assign:
                continue
            if a <= i < b and a <= j < b:
                assign[(i, j)] = k
    missing = [p for p in pairs if p not in assign]
    if missing:
        raise SystemExit(f"{label}: {len(missing)} pairs fit no block "
                         f"(e.g. {dates[missing[0][0]]}-{dates[missing[0][1]]}); "
                         f"block size {size} is too small for bandwidth {bandwidth}")
    per = {}
    for p, k in assign.items():
        per.setdefault(k, []).append(p)
    print(f"{label}: {n} dates, {len(pairs)} pairs, {len(bl)} blocks of {size} at stride {stride}")
    for k, (a, b) in enumerate(bl):
        print(f"  block {k}: dates {dates[a]}..{dates[b - 1]} ({b - a}), "
              f"{len(per.get(k, []))} pairs")
    return bl, per, pairs


def write(outdir, dates, bl, per, label):
    os.makedirs(outdir, exist_ok=True)
    for k, (a, b) in enumerate(bl):
        with open(os.path.join(outdir, f"{k}.dates"), "w") as f:
            f.write("\n".join(dates[a:b]) + "\n")
        with open(os.path.join(outdir, f"{k}.pairs"), "w") as f:
            f.write("\n".join(f"{dates[i]}-{dates[j]}" for i, j in sorted(per.get(k, []))) + "\n")
    print(f"  wrote {len(bl)} {label} block files to {outdir}")


def main(dates_file, outroot):
    dates = []
    for ln in open(dates_file):
        t = ln.split()
        if t and len(t[0]) == 6 and t[0].isdigit():
            dates.append(t[0])
    dates = sorted(set(dates))
    if len(dates) < 2:
        raise SystemExit(f"make_blocks: only {len(dates)} dates in {dates_file}")
    ib, ip, ipairs = plan(dates, NSUB_ION, ION_BLOCK, ION_STRIDE, "ion")
    write(os.path.join(outroot, "ion_blocks"), dates, ib, ip, "ion")
    mb, mp, mpairs = plan(dates, MAIN_BANDWIDTH, MAIN_BLOCK, MAIN_STRIDE, "main")
    write(os.path.join(outroot, "main_blocks"), dates, mb, mp, "main")
    # every pair must appear exactly once across the block files
    for name, root, pairs in (("ion", "ion_blocks", ipairs), ("main", "main_blocks", mpairs)):
        seen = []
        d = os.path.join(outroot, root)
        for f in sorted(os.listdir(d)):
            if f.endswith(".pairs"):
                seen += [x.strip() for x in open(os.path.join(d, f)) if x.strip()]
        if len(seen) != len(set(seen)):
            raise SystemExit(f"{name}: a pair is assigned to more than one block")
        if len(seen) != len(pairs):
            raise SystemExit(f"{name}: {len(seen)} pairs written, {len(pairs)} expected")
        print(f"  {name}: {len(seen)} pairs, each in exactly one block")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser("~/full"))
