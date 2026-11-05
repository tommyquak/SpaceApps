#!/bin/bash
#SBATCH --job-name=coarse
#SBATCH --partition=normal
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --output=/home/t/tommy/full/logs/coarse-%A_%a.log
# Block-scale ALOS-2 products from the surviving 17x56 filtered interferograms
# (the SLC stack is gone). Per pair: complex multilook (RL x AL) -> snaphu ->
# geocode.py (track xml + 3 arcsec DEM, same bbox as main.sbatch) -> aoi90
# GeoTIFFs under the finale names build_stack2 globs for. Array over 8 blocks.
#   sbatch --array=0-7 --export=ALL,RL=10,AL=8,OUT=$HOME/full/pairs_tay170c coarse_from_filt.sbatch
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj ISCE_STACK=$HOME/micromamba/envs/insar/share/isce2
export PATH=$ISCE_STACK/alosStack:$PATH PYTHONPATH=$ISCE_STACK/alosStack:${PYTHONPATH:-} OMP_NUM_THREADS=1
RL=${RL:-10}; AL=${AL:-8}; NCORR=${NCORR:-800}
NR=$((17*RL)); NA=$((56*AL)); ML2=_${NR}rlks_${NA}alks; ML1=_17rlks_56alks
SRC=${SRC:-$HOME/full/pairs_tay}
OUT=${OUT:-$HOME/full/pairs_tay170c}
TRACK=${TRACK:-/mnt/scratch/t/tommy/full/pairs_ion/190621-190830/181123.track.xml}
DEM3=$HOME/pilot/dem/elevation_3s.dem.wgs84
BBOX=${BBOX:-"-3.4/-2.34/113.8/115.2"}
B=${SLURM_ARRAY_TASK_ID:-0}; NB=${NBLOCKS:-8}; LANES=${LANES:-12}
ONLY=${ONLY:-}                                       # optional single pair "YYYY-MM-DD_YYYY-MM-DD" for a test
WORK=/tmp/coarse_${B}_${SLURM_JOB_ID:-x}; mkdir -p $WORK/dem $OUT; trap "rm -rf $WORK" EXIT
python -c "import shutil,sys; sys.exit(shutil.which('snaphu') is None)" || { echo "COARSE_ABORT no snaphu binary on PATH"; exit 1; }
[ -s "$TRACK" ] || { echo "COARSE_ABORT no track xml $TRACK"; exit 1; }
for f in $DEM3; do cp $f $f.xml $WORK/dem/ && cp $f.vrt $WORK/dem/ 2>/dev/null; done
DEM3W=$WORK/dem/$(basename $DEM3); fixImageXml.py -i $DEM3W -f > /dev/null || { echo "COARSE_ABORT fixImageXml"; exit 1; }
echo "COARSE_LOOKS rl=$RL al=$AL -> $ML2 ; block $B/$NB ; OUT=$OUT"

run_pair() {
  local cal=$1; local r=${cal:2:2}${cal:5:2}${cal:8:2}; local s=${cal:13:2}${cal:16:2}${cal:19:2}; local p=${r}-${s}
  local od=$OUT/$cal; local L=$od/coarse.log
  [ -f $od/DONE ] && { echo "PAIR_SKIP $cal"; return 0; }
  local src=$SRC/$cal/filt_${p}${ML1}.int
  [ -s "$src" ] && [ -s "$src.xml" ] || { echo "PAIR_FAILED $cal no 17x56 filtered int"; return 1; }
  mkdir -p $od $WORK/$p/insar; cp $TRACK $WORK/$p/181123.track.xml; cd $WORK/$p/insar || return 1
  python $HOME/full/coarse_ml.py $src filt_${p}${ML2} ${p}${ML2} $RL $AL $NCORR > $L 2>&1 || { echo "PAIR_FAILED $cal coarse_ml: $(tail -2 $L | tr '\n' ' ')"; return 1; }
  for f in filt_${p}${ML2}.unw ${p}${ML2}.cor filt_${p}${ML2}.unw.conncomp; do
    geocode.py -ref_date_stack_track ../181123.track.xml -dem $DEM3W -input $f -nrlks $NR -nalks $NA -bbox " $BBOX" >> $L 2>&1 \
      || { echo "PAIR_FAILED $cal geocode $f: $(grep -m1 -E 'Error|Traceback' $L)"; return 1; }
    [ -s "$f.geo.vrt" ] || { echo "PAIR_FAILED $cal no $f.geo.vrt"; return 1; }
  done
  conv() { gdal_translate -q -of GTiff -co COMPRESS=DEFLATE -co TILED=YES "$1.vrt" "$2.tmp" >> $L 2>&1 && mv "$2.tmp" "$2"; }
  conv filt_${p}${ML2}.unw.geo $od/filt_${cal}${ML2}.unw.geo.aoi90.tif || { echo "PAIR_FAILED $cal convert unw"; return 1; }
  conv ${p}${ML2}.cor.geo $od/${cal}${ML2}.cor.geo.aoi90.tif || { echo "PAIR_FAILED $cal convert cor"; return 1; }
  conv filt_${p}${ML2}.unw.conncomp.geo $od/filt_${cal}${ML2}.unw.conncomp.geo.aoi90.tif || { echo "PAIR_FAILED $cal convert conncomp"; return 1; }
  cp filt_${p}${ML2}.int filt_${p}${ML2}.int.xml ${p}${ML2}.cor ${p}${ML2}.cor.xml $od/ 2>/dev/null   # radar-coordinate wrapped + coherence kept (small)
  grep -E "^coarse|^unwrapped" $L | tr '\n' ' ' >> $od/coarse.stats; echo >> $od/coarse.stats
  cd $WORK; rm -rf $WORK/$p; touch $od/DONE; echo "PAIR_OK $cal"
}
export -f run_pair; export SRC OUT ML1 ML2 NR NA RL AL NCORR TRACK DEM3W BBOX WORK
if [ -n "$ONLY" ]; then PAIRS=$ONLY; else PAIRS=$(ls -d $SRC/*/ | xargs -n1 basename | awk -v b=$B -v nb=$NB 'NR % nb == b'); fi
echo "=== $(echo $PAIRS | wc -w) pairs, $LANES lanes ==="
printf '%s\n' $PAIRS | xargs -P $LANES -I{} bash -c 'run_pair "$1"' _ {} | tee $WORK/pairs.out
ok=$(grep -c "^PAIR_OK" $WORK/pairs.out || true); skip=$(grep -c "^PAIR_SKIP" $WORK/pairs.out || true); bad=$(grep -c "^PAIR_FAILED" $WORK/pairs.out || true)
echo "block $B: ok=$ok skip=$skip failed=$bad"
[ "$bad" -eq 0 ] && echo "COARSE_BLOCK_OK $B" || { grep "^PAIR_FAILED" $WORK/pairs.out | head -5; echo "COARSE_BLOCK_FAIL $B"; exit 1; }
