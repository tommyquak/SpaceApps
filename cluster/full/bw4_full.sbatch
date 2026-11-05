#!/bin/bash
#SBATCH --job-name=bw4full
#SBATCH --output=/home/t/tommy/full/logs/bw4full_%j.out
#SBATCH --error=/home/t/tommy/full/logs/bw4full_%j.out
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
# Tay's remedy in MintPy's required form: a bandwidth-4 network (connNumMax=4,
# the 362 pairs closure_phase_bias expects), its own SBAS inversion, the FULL
# Zheng et al. 2022 bias estimate on that network, the bias applied, and the
# corrected rate tested against Tay's rate raster. Own directory; the
# bandwidth-20 products are untouched.
set -u
export PATH=$HOME/micromamba/envs/insar/bin:$PATH PROJ_LIB=$HOME/micromamba/envs/insar/share/proj TMPDIR=$HOME/full/tmp PYTHONUNBUFFERED=1
SRC=$HOME/full/mintpy_tay; W=/mnt/scratch/t/tommy/full/mintpy_bw4; mkdir -p $W/inputs
free_kb=$(df -k /mnt/scratch | awk 'NR==2{print $4}'); [ "$free_kb" -gt 30000000 ] || { echo "BW4_ABORT scratch"; exit 1; }
for f in ifgramStack.h5 geometryGeo.h5 SET.h5; do [ -f $W/inputs/$f ] || cp $SRC/inputs/$f $W/inputs/$f || { echo "BW4_ABORT copy $f"; exit 1; }; done
python - "$SRC/smallbaselineApp.cfg" "$W/smallbaselineApp.cfg" <<'PY'
import re, sys
src, dst = sys.argv[1], sys.argv[2]; out = []; seen = False
for L in open(src).read().splitlines():
    if re.match(r"^\s*mintpy\.network\.connNumMax\s*=", L): out.append("mintpy.network.connNumMax = 4"); seen = True
    elif re.match(r"^\s*mintpy\.network\.tempBaseMax\s*=", L): out.append("mintpy.network.tempBaseMax = auto")
    else: out.append(L)
if not seen: out.append("mintpy.network.connNumMax = 4")
open(dst, "w").write("\n".join(out) + "\n"); print("cfg: connNumMax=4")
PY
cd $W || exit 1
echo "BW4_START $(date)"
smallbaselineApp.py smallbaselineApp.cfg --start modify_network --end velocity --dostep modify_network || { echo "BW4_FAILED modify_network"; exit 1; }
n=$(python -c "import h5py,numpy as np;f=h5py.File('inputs/ifgramStack.h5');print(int(np.sum(f['dropIfgram'][:])))"); echo "BW4 kept pairs: $n"
smallbaselineApp.py smallbaselineApp.cfg --start reference_point --end velocity || { echo "BW4_FAILED inversion"; exit 1; }
echo "BW4_INVERTED $(date)"
O=$W/zheng; mkdir -p $O
closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 20 --bw 4 -a mask -o $O; echo "BW4_MASK_RC=$?"
closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 20 --bw 4 -a quick_estimate -c local --num-worker 16 -o $O; echo "BW4_QUICK_RC=$?"
closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 20 --bw 4 -a estimate -c local --num-worker 16 -o $O; rc=$?; echo "BW4_EST_RC=$rc"
ls -la $O/*.h5 2>/dev/null | awk '{print $5, $9}'
for B in timeseriesBias.h5 timeseriesBiasApprox.h5; do
  [ -f $O/$B ] && { echo "=== apply $B to the bandwidth-4 series ==="; BW=$W BIAS=$O/$B python $HOME/full/zheng_apply_bw60.py; }
done
echo "=== also: approximate bias applied to bw60 (for reference) ==="
BW=/mnt/scratch/t/tommy/full/mintpy_bw60 BIAS=$O/timeseriesBiasApprox.h5 python $HOME/full/zheng_apply_bw60.py 2>/dev/null | tail -6
echo "BW4_DONE $(date)"
