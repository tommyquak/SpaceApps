#!/bin/bash
#SBATCH --job-name=barrier
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:30:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --no-requeue
#SBATCH --output=/home/t/tommy/full/logs/barrier-%j.log
# alosStack pivot, BARRIER (spec sec 5). Serial, minutes, between the ion pass
# and the main pass:
#
#   1. automated ion QC over every persisted ion_cal product; pairs outside the
#      thresholds go into -exc_pair and are logged ION_EXCLUDED for review
#   2. the surviving ion network must still connect every date, asserted before
#      release, since ion_ls solves a least squares over that network
#   3. compute_baseline over the complete persisted dates/ tree (manifest D)
#   4. ion_ls -> dates_ion/, the per-date ion grids the MAIN pass corrects with
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj
export ISCE_STACK=$HOME/micromamba/envs/insar/share/isce2
export PATH=$ISCE_STACK/alosStack:$PATH
export OMP_NUM_THREADS=8

REF=181123
NRLKS1=1; NALKS1=14; NRLKS2=5; NALKS2=2
NRLKS_ION=80; NALKS_ION=32
STD_MAX=${STD_MAX:-0.8}           # rad; cuts degenerate estimates only, see ion_qc.py:
                                  # filt_ion std tracks solar cycle 25 (0.0995 rad in
                                  # 2018-2021, 0.42 in 2023), and ion_ls weights by it
PERSIST=${PERSIST:-/mnt/scratch/t/tommy/full/stack}
ION_OUT=${ION_OUT:-/mnt/scratch/t/tommy/full/pairs_ion}
ION_DIR=${ION_DIR:-/mnt/scratch/t/tommy/full/dates_ion}
echo "barrier on $(hostname) at $(date -Is)"

WORK=/tmp/barrier_${SLURM_JOB_ID}; mkdir -p $WORK
trap "rm -rf $WORK" EXIT
cd $WORK

# ---------------------------------------------------------------- 1 + 2: QC
python3 $HOME/full/ion_qc.py $ION_OUT $HOME/full/ion_blocks $STD_MAX > $WORK/qc.txt 2>&1
rc=$?
cat $WORK/qc.txt
[ $rc -eq 0 ] || { echo "BARRIER_ABORT ion QC refused to release the network"; exit 1; }
EXC=$(grep "^EXCLUDE " $WORK/qc.txt | cut -d' ' -f2- | tr '\n' ' ')
[ -n "$EXC" ] && echo "ION_EXCLUDED $EXC"
cp $WORK/qc.txt $ION_OUT/ion_qc.txt 2>/dev/null

# --------------------------------------------- 3: the authoritative baseline
# manifest D: written ONCE over the complete persisted dates/ tree
ln -sfn $PERSIST/dates $WORK/dates
DATES=$(ls -d $PERSIST/dates/*/ 2>/dev/null | xargs -n1 basename | grep -v "^$REF$" | tr '\n' ' ')
[ -n "$DATES" ] || { echo "BARRIER_ABORT no persisted dates"; exit 1; }
echo "=== compute_baseline over $(echo $DATES | wc -w) secondary dates at $(date -Is) ==="
compute_baseline.py -idir dates -odir $PERSIST/baseline -ref_date $REF -sec_date $DATES \
  -baseline_center baseline_center.txt || { echo "BARRIER_ABORT compute_baseline"; exit 1; }
[ -s "$PERSIST/baseline/baseline_center.txt" ] || { echo "BARRIER_ABORT no baseline_center.txt"; exit 1; }
# build_stack2._load_baselines reads ~/full/baseline/baseline_center.txt, and the
# MAIN output carries no isce.log, so without a copy here bperp_of returns None
# for every pair and MintPy silently gets bperp=0 (spec sec 5, BPERP BRIDGE)
mkdir -p $HOME/full/baseline
cp $PERSIST/baseline/baseline_center.txt $HOME/full/baseline/ \
  || { echo "BARRIER_ABORT could not publish baseline_center.txt for build_stack2"; exit 1; }
n=$(awk 'NF==4 && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/' $HOME/full/baseline/baseline_center.txt | wc -l)
[ "$n" -ge 2 ] || { echo "BARRIER_ABORT baseline_center.txt has $n parseable rows"; exit 1; }
echo "baseline_center.txt: $(wc -l < $PERSIST/baseline/baseline_center.txt) lines, $n parseable rows, published to ~/full/baseline"

# ------------------------------------------------------------------ 4: ion_ls
mkdir -p $ION_DIR
echo "=== ion_ls at $(date -Is) ==="
ion_ls.py -idir $ION_OUT -odir $ION_DIR -ref_date_stack $REF \
  -nrlks1 $NRLKS1 -nalks1 $NALKS1 -nrlks2 $NRLKS2 -nalks2 $NALKS2 \
  -nrlks_ion $NRLKS_ION -nalks_ion $NALKS_ION -interp \
  ${EXC:+-exc_pair $EXC} || { echo "BARRIER_ABORT ion_ls"; exit 1; }

n=$(ls $ION_DIR/*.ion 2>/dev/null | wc -l)
echo "dates_ion: $n per-date ion grids, $(du -sm $ION_DIR | cut -f1) MB"
[ "$n" -gt 0 ] || { echo "BARRIER_ABORT ion_ls produced no per-date grids"; exit 1; }
python3 $HOME/full/check_run.py dates_ion || { echo "BARRIER_ABORT dates_ion failed the validity check"; exit 1; }
touch $ION_DIR/BARRIER_DONE
echo "BARRIER_OK at $(date -Is)"
