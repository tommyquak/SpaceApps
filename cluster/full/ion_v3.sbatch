#!/bin/bash
#SBATCH --job-name=ionv3
#SBATCH --partition=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/home/t/tommy/full/logs/ionv3-%j.log
# Re-solve the per-date ionosphere with the pairs whose own estimate is
# physically impossible removed.
#
# Cunren Liang's insarzd documentation (author of this ionosphere method) states
# the correction "works for approximately 70% of interferograms" and fails where
# "the majority of the imaged area is low coherence area like lake, ocean".
# Measured here: 247 of 362 pairs good = 68.2 percent, matching that figure. The
# method is behaving exactly as designed. The defect was ours: ion_ls.py:233
# solves unweighted normal equations with NO outlier rejection, so the ~30
# percent the author says will fail contaminated all 93 dates. Median per-pair
# input is -2.21 rad; median per-date output was -96,611 rad.
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj
export ISCE_STACK=$HOME/micromamba/envs/insar/share/isce2
export PATH=$ISCE_STACK/alosStack:$PATH
export OMP_NUM_THREADS=4

REF=181123
NRLKS1=1; NALKS1=14; NRLKS2=17; NALKS2=4   # Tay 2025: 17 rg x 56 az, ~180 m
NRLKS_ION=80; NALKS_ION=32
ION_OUT=/mnt/scratch/t/tommy/full/pairs_ion
ION_DIR=/mnt/scratch/t/tommy/full/dates_ion_v3

python3 $HOME/ion_outliers.py > $HOME/full/logs/ion_outliers.txt 2>&1 || { echo "IONFIX_ABORT outlier scan"; exit 1; }
EXC=$(python3 -c "import json;print(' '.join(json.load(open('$HOME/ion_bad_pairs.json'))))")
NEXC=$(echo $EXC | wc -w)
echo "IONFIX excluding $NEXC pairs of $(ls $ION_OUT | wc -l)"
[ "$NEXC" -gt 0 ] || { echo "IONFIX_ABORT empty exclusion list"; exit 1; }

rm -rf $ION_DIR; mkdir -p $ION_DIR
ion_ls.py -idir $ION_OUT -odir $ION_DIR -ref_date_stack $REF \
  -nrlks1 $NRLKS1 -nalks1 $NALKS1 -nrlks2 $NRLKS2 -nalks2 $NALKS2 \
  -nrlks_ion $NRLKS_ION -nalks_ion $NALKS_ION -interp \
  -exc_pair $EXC || { echo "IONFIX_ABORT ion_ls"; exit 1; }

n=$(ls $ION_DIR/*.ion 2>/dev/null | wc -l)
echo "IONFIX produced $n per-date grids, $(du -sm $ION_DIR | cut -f1) MB"
[ "$n" -eq 93 ] || { echo "IONFIX_ABORT expected 93 grids, got $n"; exit 1; }

# Gate on physics: the whole point is that the output must now be plausible.
python3 - <<'PY'
import glob, re, numpy as np, sys, os
d = "/mnt/scratch/t/tommy/full/dates_ion_v3"
fs = sorted(glob.glob(d + "/filt_ion_*.ion"))
t = open(fs[0] + ".vrt").read()
m = re.search(r'rasterXSize="(\d+)"\s+rasterYSize="(\d+)"', t)
W, H = (int(m.group(1)), int(m.group(2))) if m else (0, 0)
mean, gaz, ptp = [], [], []
for f in fs:
    a = np.fromfile(f, dtype=np.float32)
    if a.size != W * H: continue
    a = a.reshape(H, W)
    mean.append(float(np.nanmean(a))); ptp.append(float(np.nanmax(a) - np.nanmin(a)))
    gaz.append(float(np.nanmedian(np.abs(np.diff(a, axis=0)))))
mean, gaz, ptp = map(np.array, (mean, gaz, ptp))
print("IONFIX_STATS dates=%d" % len(mean))
print("  per-date mean:   median %10.2f  min %10.2f  max %10.2f rad" % (np.median(mean), mean.min(), mean.max()))
print("  per-date ptp:    median %10.2f  max %10.2f rad" % (np.median(ptp), ptp.max()))
print("  |grad| azimuth:  median %10.4f  max %10.4f rad/px  (pi=3.1416)" % (np.median(gaz), gaz.max()))
print("  implied TEC:     %.1f TECU (real ionosphere peaks ~200)" % (abs(np.median(mean)) / 13.4))
over = int((gaz > np.pi).sum())
print("  dates over pi/px: %d of %d" % (over, len(gaz)))
if over > 0 or abs(np.median(mean)) > 2000:
    print("IONFIX_FAIL grids still unphysical"); sys.exit(1)
print("IONFIX_OK grids are physically plausible")
PY
[ $? -eq 0 ] && echo "IONFIX_DONE" || echo "IONFIX_ABORT physics gate"
