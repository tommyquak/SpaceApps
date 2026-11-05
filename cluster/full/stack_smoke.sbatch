#!/bin/bash
#SBATCH --job-name=stack-smoke
#SBATCH --partition=long
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --tmp=1500000
#SBATCH --time=12:00:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --no-requeue
#SBATCH --output=/home/t/tommy/full/logs/stack-smoke.log
# alosStack pivot, wave 0: 3-date full-frame smoke on one node.
# Validates create_cmds config, per-date read/offset/resample, and interferogram
# formation on the authentic alosStack code path before the cropped int16 waves.
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj
export ISCE_STACK=$HOME/micromamba/envs/insar/share/isce2
export PATH=$ISCE_STACK/alosStack:$PATH
export OMP_NUM_THREADS=48

WORK=/tmp/smoke_${SLURM_JOB_ID}; mkdir -p $WORK/data
avail=$(df -m $WORK | awk 'NR==2{print $4}')
[ "$avail" -lt 1000000 ] && { echo "scratch too small: ${avail}M"; exit 1; }
echo "workdir $WORK (${avail}M free) on $(hostname)"
trap "rm -rf $WORK" EXIT

# three consecutive staged dates; reference in the middle
DATES="2018-10-12 2018-11-23 2018-12-21"
for d in $DATES; do
  yymmdd=$(date -d "$d" +%y%m%d)
  mkdir -p $WORK/data/$yymmdd
  cp ~/full/ceos/$d/* $WORK/data/$yymmdd/ || { echo "stage failed for $d"; exit 1; }
done
mkdir -p $WORK/dem
cp ~/pilot/dem/elevation.dem.wgs84* $WORK/dem/ 2>/dev/null
cp ~/pilot/dem/waterBody.wbd $WORK/dem/waterBody.wbd || { echo WBD_STAGE_FAILED; exit 1; }
cp ~/pilot/dem/waterBody.wbd.xml $WORK/dem/waterBody.wbd.xml || { echo WBD_XML_STAGE_FAILED; exit 1; }
# ISCE requires absolute file_name paths in the image xmls (tutorial's fixImageXml step)
fixImageXml.py -i $WORK/dem/elevation.dem.wgs84 -f && fixImageXml.py -i $WORK/dem/waterBody.wbd -f || { echo DEM_XML_FIX_FAILED; exit 1; }

cd $WORK
cat > alosStack.xml <<XMLEOF
<?xml version="1.0" encoding="UTF-8"?>
<stack>
<component name="stackinsar">
    <property name="data directory">$WORK/data</property>
    <property name="dem for coregistration">$WORK/dem/elevation.dem.wgs84</property>
    <property name="dem for geocoding">$WORK/dem/elevation.dem.wgs84</property>
    <property name="water body">$WORK/dem/waterBody.wbd</property>
    <property name="reference date of the stack">181123</property>
    <property name="polarization">HH</property>
    <property name="starting swath">2</property>
    <property name="ending swath">4</property>
    <property name="number of subsequent dates">1</property>
</component>
</stack>
XMLEOF
echo "=== create_cmds ==="
create_cmds.py -stack_par alosStack.xml || { echo "CREATE_CMDS_FAILED"; exit 1; }
ls cmd_*.sh
echo "=== cmd_1 (per-date) ==="
t0=$(date +%s)
bash ./cmd_1.sh > cmd1.out 2>&1 || { echo "CMD1_FAILED"; tail -40 cmd1.out; exit 1; }
echo "cmd_1 ok in $(( ($(date +%s)-t0)/60 )) min"
echo "=== cmd_2 (per-pair, pre-ion) ==="
t0=$(date +%s)
bash ./cmd_2.sh > cmd2.out 2>&1 || { echo "CMD2_FAILED"; tail -40 cmd2.out; exit 1; }
echo "cmd_2 ok in $(( ($(date +%s)-t0)/60 )) min"
echo "=== inventory ==="
du -sm dates dates_resampled pairs baseline 2>/dev/null
find dates_resampled -maxdepth 3 | head -25
find pairs -maxdepth 4 \( -name '*.int' -o -name 'diff*' \) 2>/dev/null | head -10
# persist the sec-3 gate inputs (one secondary date, all swaths): resampled
# full-band + authentic subbands + xml/vrt + cull.off, df-guarded (~17 GB)
free_home=$(df -m ~ | awk 'NR==2{print $4}')
if [ "$free_home" -gt 40000 ]; then
  mkdir -p ~/full/stack_smoke_artifacts/gate
  rsync -aR --include='*/' --include='181012*.slc' --include='181012*.slc.xml' --include='181012*.slc.vrt' --include='cull.off' --include='*.off.xml' --include='*.track.xml' --include='*.frame.xml' --exclude='*' dates_resampled/181012 dates/181012 ~/full/stack_smoke_artifacts/gate/ 2>/dev/null || echo GATE_PERSIST_PARTIAL
  echo "gate inputs persisted: $(du -sm ~/full/stack_smoke_artifacts/gate 2>/dev/null | cut -f1) MB"
else
  echo "GATE_PERSIST_SKIPPED home headroom ${free_home}M"
fi
# geometric-branch validation: re-run the offset estimation for one secondary
# date with the override armed and keep the variant cull.off for comparison
# NOTE: overwrites the authentic ampcor cull.off in place; runs after the
# gate rsync deliberately. -use_wbd_offset is required: the patched force-gate
# sits inside the useWbdForNumberOffsets branch. No -odir exists in this CLI.
mkdir -p ~/full/stack_smoke_artifacts
ALOS2_FORCE_GEOM_OFFSET=1 python3 $ISCE_STACK/alosStack/estimate_slc_offset.py -idir $WORK/dates -ref_date 181123 -sec_date 181012 -wbd $WORK/dem/waterBody.wbd -dem $WORK/dem/elevation.dem.wgs84 -use_wbd_offset > geomtest.out 2>&1 && echo GEOM_BRANCH_OK || { echo GEOM_BRANCH_FAILED; tail -5 geomtest.out; }
for sw in 2 3 4; do cp $WORK/dates/181012/f1_3700/s$sw/cull.off ~/full/stack_smoke_artifacts/geom_cull_s$sw.off 2>/dev/null || true; done
cp geomtest.out ~/full/stack_smoke_artifacts/ 2>/dev/null || true

# nonzero-data gate on one resampled SLC and one interferogram
python3 - <<'PYX'
import glob, numpy as np, os, sys
ok = True
slcs = glob.glob("dates_resampled/*/f1_*/s*/*.slc") or glob.glob("dates_resampled/*/*/*.slc") or glob.glob("dates_resampled/**/*.slc", recursive=True)
if slcs:
    a = np.fromfile(slcs[0], dtype=np.complex64, count=2_000_000)
    nz = np.count_nonzero(a)
    print(f"slc sample {os.path.basename(slcs[0])}: nonzero {nz}/2000000")
    ok &= nz > 100000
else:
    print("NO RESAMPLED SLC FOUND"); ok = False
ints = glob.glob("pairs/*/**/*.int", recursive=True)
if ints:
    a = np.fromfile(ints[0], dtype=np.complex64, count=2_000_000)
    nz = np.count_nonzero(a)
    print(f"int sample {os.path.basename(ints[0])}: nonzero {nz}/2000000")
    ok &= nz > 100000
else:
    print("NO INTERFEROGRAM FOUND"); ok = False
sys.exit(0 if ok else 1)
PYX
rc=$?
# persist the tiny artefacts we need to design the crop layer
mkdir -p ~/full/stack_smoke_artifacts
cp cmd_*.sh alosStack.xml cmd1.out cmd2.out ~/full/stack_smoke_artifacts/ 2>/dev/null
find dates_resampled -name '*.xml' | head -3 | xargs -I{} cp {} ~/full/stack_smoke_artifacts/ 2>/dev/null
find dates_resampled -name '*.pck' | head -3 | xargs -I{} cp {} ~/full/stack_smoke_artifacts/ 2>/dev/null
[ $rc -eq 0 ] && echo SMOKE_OK || echo SMOKE_DATAGATE_FAILED
