#!/bin/bash
#SBATCH --job-name=bw60
#SBATCH --output=/home/t/tommy/full/logs/bw60_%j.out
#SBATCH --error=/home/t/tommy/full/logs/bw60_%j.out
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
# Bandwidth-restricted SBAS test: same 1650-pair stack, network cut to
# temporal baselines <= TBMAX days. closure_strat.py measured 7 percent bad
# triplets under 30 d against 75 percent beyond 400 d, so this is the one
# network that unwraps cleanly. Faithful to Tay otherwise: the cfg is a copy.
set -u
export PROJ_LIB=$HOME/micromamba/envs/insar/share/proj TMPDIR=$HOME/full/tmp PYTHONUNBUFFERED=1
export PATH=$HOME/micromamba/envs/insar/bin:$PATH
TBMAX=${TBMAX:-60}
SRC=${SRC:-$HOME/full/mintpy_tay}                       # any MintPy dir with inputs/ifgramStack.h5 (e.g. the 170x448 block product)
W=${W:-/mnt/scratch/t/tommy/full/mintpy_bw${TBMAX}}
mkdir -p "$W/inputs" "$HOME/full/tmp"
need_kb=$(( $(du -k "$SRC/inputs/ifgramStack.h5" | cut -f1) * 3 / 2 ))
free_kb=$(df -k /mnt/scratch | awk 'NR==2{print $4}')
echo "need ${need_kb} kB, free ${free_kb} kB on /mnt/scratch"
[ "$free_kb" -gt "$need_kb" ] || { echo "BW_ABORT insufficient scratch"; exit 1; }
for f in ifgramStack.h5 geometryGeo.h5 SET.h5; do
  [ -f "$W/inputs/$f" ] || cp "$SRC/inputs/$f" "$W/inputs/$f" || { echo "BW_ABORT copy $f"; exit 1; }
done
python - "$SRC/smallbaselineApp.cfg" "$W/smallbaselineApp.cfg" "$TBMAX" <<'PY'
import re, sys
src, dst, tb = sys.argv[1], sys.argv[2], sys.argv[3]
lines = open(src).read().splitlines()
out, seen = [], False
for L in lines:
    if re.match(r"^\s*mintpy\.network\.tempBaseMax\s*=", L):
        out.append(f"mintpy.network.tempBaseMax = {tb}"); seen = True
    else:
        out.append(L)
if not seen:
    out.append(f"mintpy.network.tempBaseMax = {tb}")
open(dst, "w").write("\n".join(out) + "\n")
print("cfg written; tempBaseMax =", tb)
PY
cd "$W" || exit 1
echo "BW_START $(date)"
smallbaselineApp.py smallbaselineApp.cfg --start modify_network --end velocity
rc=$?
echo "BW_MINTPY_RC=$rc"
[ $rc -eq 0 ] || { echo "BW_FAILED"; exit 1; }
W="$W" python $HOME/full/summarize_bw.py
echo "BW_DONE $(date)"
