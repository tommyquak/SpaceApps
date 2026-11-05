#!/bin/bash
#SBATCH --job-name=alos2-pairs
#SBATCH --partition=long
#SBATCH --array=1-70%8
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=10:00:00
#SBATCH --exclude=xcng[0-1],xgpj0
#SBATCH --output=/home/t/tommy/pilot/logs/pair-%a.log
set -uo pipefail
grep -q PAIRTEST_OK ~/pilot/logs/pairtest.log 2>/dev/null || { echo "pair test not passed, aborting"; exit 1; }
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/pilot/pairs.txt)
REF=${line%% *}; SEC=${line##* }
OUT=~/pilot/pairs/${REF}_${SEC}
[ -f $OUT/DONE ] && { echo "pair ${REF}_${SEC} already done"; exit 0; }
mkdir -p $OUT
WORK=/temp/${SLURM_JOB_ID}; mkdir -p $WORK 2>/dev/null || { WORK=/tmp/pr_$$; mkdir -p $WORK; }
trap "cp $WORK/isce.log $OUT/ 2>/dev/null; rm -rf $WORK" EXIT
echo "pair ${REF}_${SEC} on $(hostname), workdir $WORK"
cp -r ~/pilot/ceos/$REF ~/pilot/ceos/$SEC $WORK/
cp ~/pilot/dem/elevation.dem.wgs84* ~/pilot/dem/elevation.dem.hdr $WORK/ 2>/dev/null
cp ~/pilot/dem/waterBody.wbd ~/pilot/dem/waterBody.wbd.xml $WORK/ 2>/dev/null
cd $WORK
sed -e "s|REFDIR|$WORK/$REF|; s|SECDIR|$WORK/$SEC|; s|DEMPATH|$WORK/elevation.dem.wgs84|; s|WBDPATH|$WORK/waterBody.wbd|" ~/pilot/alos2App.template.xml > alos2App.xml
export PATH=$HOME/micromamba/envs/insar/bin:$HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj
export OMP_NUM_THREADS=16
set +e
python3 $HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications/alos2App.py alos2App.xml
rc=$?
set -e
echo "alos2App exit: $rc"
# Keep MintPy-ready AOI crops of every geocoded layer (identical grid across
# pairs); full-frame .geo files are ~1 GB each and would blow the home quota.
for G in $WORK/insar/*.geo; do
  [ -e "$G" ] || continue
  B=$(basename $G)
  ~/micromamba/envs/insar/bin/python -c "
from osgeo import gdal
gdal.UseExceptions()
gdal.Translate(\"$OUT/$B.aoi.tif\", \"$G\", projWin=[113.8, -2.0, 115.2, -3.4], creationOptions=[\"COMPRESS=DEFLATE\", \"TILED=YES\"])
" || echo "crop failed for $B"
done
cp $WORK/insar/*.geo.xml $OUT/ 2>/dev/null || true
[ $rc -eq 0 ] && { touch $OUT/DONE; echo "PAIR_DONE ${REF}_${SEC}"; } || { echo "PAIR_FAILED ${REF}_${SEC}"; exit 1; }
