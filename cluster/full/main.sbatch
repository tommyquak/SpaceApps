#!/bin/bash
#SBATCH --job-name=main
#SBATCH --partition=long
#SBATCH --cpus-per-task=24
#SBATCH --mem=60G
#SBATCH --tmp=400000
# Sized to the 72-node pool of 24-CPU / 63 GB machines. A 48-CPU 180 GB request
# only fits the 14 large nodes, so blocks queue on "(Resources)" while most of
# the cluster sits idle. Staging is ~3 GB of I/O per swath here (no subband
# regeneration), so the ion pass's NVMe constraint is not needed either.
#SBATCH --time=16:00:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --no-requeue
#SBATCH --output=/home/t/tommy/full/logs/main-%A_%a.log
# alosStack pivot, MAIN PASS (spec sec 5, 5b, 6). One job per date block, every
# pair of the block in one GNU parallel lane:
#
#   form + mosaic + diff (ml1) -> look_coherence (makes the ml2 diff)
#   -> ion_correct (-ion_dir ../../dates_ion) -> filt -> unwrap_snaphu
#   -> geocode.py DIRECTLY with -bbox (create_cmds emits a broken '-bbox /')
#   -> warp onto the aoi90 grid under the finale names, DONE last
#
# NO looks1 pair product is persisted. Blocks come from make_blocks.py.
#   sbatch --array=0-6 main.sbatch
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj
export ISCE_STACK=$HOME/micromamba/envs/insar/share/isce2
export PATH=$ISCE_STACK/alosStack:$PATH
export OMP_NUM_THREADS=2

REF=181123
FRAME=3700
NRLKS1=1; NALKS1=14
# ml2 looks are parameterised. At 5x28 (140 looks) and this AOI's median
# coherence of 0.22, phase noise is ~0.27 rad and snaphu makes frequent 2-pi
# errors: 93 percent of closure triplets inconsistent, temporal coherence 0.03
# median (max 0.215), and MintPy's phase_closure repair cannot run because the
# common connected-component region across 1631 interferograms is empty. The
# WRAPPED closure is clean (+0.2 mm/triplet), so the phase is sound and only the
# unwrapping failed. Phase noise falls as 1/sqrt(looks): 15x84 is 1260 looks,
# 9x the current, taking noise to ~0.09 rad at ~200x286 m ground resolution,
# which still resolves kilometre-scale peat domes.
# DEFECT 25: the defaults were 5 and 2, an abandoned experiment, so a launch that
# forgot to pass the looks silently produced 5x28 products (the era test did).
# Default to the certified Tay looks and say which looks are in force.
NRLKS2=${NRLKS2:-17}; NALKS2=${NALKS2:-4}
echo "MAIN_LOOKS nrlks2=$NRLKS2 nalks2=$NALKS2 (ml2 $((NRLKS1*NRLKS2))x$((NALKS1*NALKS2)))"
FILT_ALPHA=${FILT_ALPHA:-0.3}
ML1=_${NRLKS1}rlks_${NALKS1}alks
ML2=_$((NRLKS1*NRLKS2))rlks_$((NALKS1*NALKS2))alks
BBOX=${BBOX:-"-3.4/-2.34/113.8/115.2"}        # s/n/w/e, spec sec 6. Passed with a
                                              # LEADING SPACE below: the value starts
                                              # with '-', which argparse would take as
                                              # another option flag ("expected one
                                              # argument"). geocode.py:72 does
                                              # float(x) per '/' field, so the space
                                              # is harmless.
LANES=${LANES:-12}
CP=$HOME/full/crop_pack.py
PERSIST=${PERSIST:-/mnt/scratch/t/tommy/full/stack}
ION_DIR=${ION_DIR:-/mnt/scratch/t/tommy/full/dates_ion}
USE_ION=${USE_ION:-1}
OUT=${OUT:-$HOME/full/pairs}                  # finale interface (build_stack2)
DEM=$HOME/pilot/dem/elevation.dem.wgs84
DEM3=$HOME/pilot/dem/elevation_3s.dem.wgs84   # geocoding DEM sets the posting
WBD=$HOME/pilot/dem/waterBody.wbd
B=${SLURM_ARRAY_TASK_ID:-0}
# BLOCKS lets a retry run point at an alternate block dir. The Zheng
# closure-bias gate needs the sequential network complete through level 20,
# and a complete bandwidth-20 network over 93 dates is exactly 1650 pairs,
# so a single missing pair costs the whole correction.
BLOCKS=${BLOCKS:-$HOME/full/main_blocks}
DATES_F=$BLOCKS/$B.dates
PAIRS_F=$BLOCKS/$B.pairs
for f in $DATES_F $PAIRS_F; do [ -s "$f" ] || { echo "MAIN_ABORT no $f"; exit 1; }; done
DATES=$(cat $DATES_F); PAIRS=$(cat $PAIRS_F)
echo "main block $B: $(echo $DATES | wc -w) dates, $(echo $PAIRS | wc -w) pairs on $(hostname) at $(date -Is)"
[ "$USE_ION" = "1" ] && { [ -d "$ION_DIR" ] || { echo "MAIN_ABORT no dates_ion at $ION_DIR (barrier not done)"; exit 1; }; }
[ "$USE_ION" = "1" ] || echo "MAIN_NOTE ionosphere correction DISABLED (USE_ION=0)"

WORK=/tmp/main_${B}_${SLURM_JOB_ID}; mkdir -p $WORK/dem
avail=$(df -m $WORK | awk 'NR==2{print $4}')
[ "$avail" -lt 900000 ] && { echo "MAIN_ABORT scratch too small: ${avail}M"; exit 1; }
trap "rm -rf $WORK" EXIT
cd $WORK

for f in $DEM $DEM3 $WBD; do
  [ -f "$f" ] && [ -f "$f.xml" ] || { echo "MAIN_ABORT missing $f or its xml"; exit 1; }
  cp $f $f.xml $WORK/dem/ || { echo "MAIN_ABORT dem stage failed for $f"; exit 1; }
done
cp $DEM.vrt $DEM3.vrt $WORK/dem/ 2>/dev/null
DEMW=$WORK/dem/$(basename $DEM); DEM3W=$WORK/dem/$(basename $DEM3); WBDW=$WORK/dem/$(basename $WBD)
for f in $DEMW $DEM3W $WBDW; do
  fixImageXml.py -i $f -f || { echo "MAIN_ABORT fixImageXml failed for $f"; exit 1; }
done

# ------------------------------------------------- bundle + block unpack
mkdir -p dates_resampled/$REF/f1_$FRAME dates_resampled/$REF/insar dates/$REF/f1_$FRAME/mosaic
cp $PERSIST/dates_resampled/$REF/$REF.track.xml dates_resampled/$REF/ || exit 1
cp $PERSIST/dates_resampled/$REF/f1_$FRAME/$REF.frame.xml dates_resampled/$REF/f1_$FRAME/ || exit 1
cp $PERSIST/dates_resampled/$REF/insar/* dates_resampled/$REF/insar/ || exit 1
cp $PERSIST/dates/$REF/$REF.track.xml dates/$REF/ || exit 1
cp $PERSIST/dates/$REF/f1_$FRAME/$REF.frame.xml dates/$REF/f1_$FRAME/ || exit 1
cp $PERSIST/dates/$REF/f1_$FRAME/mosaic/swath_offset.txt dates/$REF/f1_$FRAME/mosaic/ || exit 1
for sn in 2 3 4; do mkdir -p dates/$REF/f1_$FRAME/s$sn dates_resampled/$REF/f1_$FRAME/s$sn; done
# ion_correct resolves ../../../dates_ion from inside insar/, i.e. the stack root
ln -sfn $ION_DIR $WORK/dates_ion

STAGE=$(printf '%s\n' $DATES $REF | sort -u)
echo "=== unpack $(echo $STAGE | wc -w) dates at $(date -Is) ==="
unpack_one() {
  d=$1
  src=$PERSIST/dates_resampled/$d/f1_$FRAME
  [ -d "$src" ] || { echo "MAIN_MISSING $d has no crops"; return 1; }
  mkdir -p dates_resampled/$d/f1_$FRAME dates_resampled/$d/insar dates/$d/f1_$FRAME
  cp $PERSIST/dates/$d/$d.track.xml dates/$d/ || return 1
  cp $PERSIST/dates/$d/f1_$FRAME/$d.frame.xml dates/$d/f1_$FRAME/ || return 1
  for sn in 2 3 4; do
    od=dates_resampled/$d/f1_$FRAME/s$sn; mkdir -p $od
    pk=$src/s$sn/$d.slc                       # MAIN needs the full band only
    [ -s "$pk" ] || { echo "MAIN_MISSING $d s$sn"; return 1; }
    python3 $CP unpack $pk $src/s$sn/$d.slc.json $od/$d.slc || return 1
    w=$(python3 -c "import json;s=json.load(open('$src/s$sn/$d.slc.json'));print(s['width'],s['length'])")
    python3 $CP render $od/$d.slc ${w% *} ${w#* } CFLOAT 1 BIL || return 1
  done
  if [ -s "$PERSIST/dates_resampled/$d/insar/${d}${ML1}_rg_rect.off" ]; then
    cp $PERSIST/dates_resampled/$d/insar/${d}${ML1}_rg_rect.off* dates_resampled/$d/insar/ || return 1
  fi
  echo "unpacked $d"
}
export -f unpack_one
export PERSIST FRAME CP ML1
# GNU parallel is not installed on this cluster; xargs -P is
printf '%s\n' $STAGE | xargs -P ${UNPACK_LANES:-16} -I{} bash -c 'unpack_one "$1"' _ {} \
  || { echo "MAIN_ABORT unpack failed"; exit 1; }
mosaic_parameter.py -idir dates_resampled -ref_date $REF -sec_date $REF \
  -nrlks1 $NRLKS1 -nalks1 $NALKS1 >/dev/null || { echo "MAIN_ABORT mosaic_parameter"; exit 1; }

# ml2 geometry. Manifest F persists the ml1 rasters and allows the ml2 set to be
# "regenerated by look_geom at stage time"; this is that step, and it is not
# optional: unwrap_snaphu -wbd_msk opens <ref>_5rlks_28alks.wbd.xml
# (runUnwrapSnaphu.py:81) and the reference los geocode needs <ref>_5rlks_28alks.los.
echo "=== look_geom (ml2 geometry) at $(date -Is) ==="
( cd dates_resampled/$REF && look_geom.py -date $REF -wbd $WBDW \
    -nrlks1 $NRLKS1 -nalks1 $NALKS1 -nrlks2 $NRLKS2 -nalks2 $NALKS2 ) \
  || { echo "MAIN_ABORT look_geom"; exit 1; }
for e in lat lon hgt los wbd; do
  [ -s "dates_resampled/$REF/insar/${REF}${ML2}.$e" ] \
    || { echo "MAIN_ABORT look_geom did not produce ${REF}${ML2}.$e"; exit 1; }
done

echo "=== pair_up at $(date -Is) ==="
pair_up.py -idir1 dates -idir2 dates_resampled -odir pairs -ref_date $REF -pairs $PAIRS \
  || { echo "MAIN_ABORT pair_up"; exit 1; }

pfail() {   # $1=pair $2=step $3=logfile
  echo "PAIR_FAILED $1 $2"
  echo "--- $1 $2 ---"
  grep -viE '^\[|^Using default ISCE|^$' "$3" 2>/dev/null | tail -6
  # product sizes: a "wrong size" failure (snaphu on a non-divisible looks
  # choice, 821292) is undiagnosable from the traceback alone, and the pair
  # directory on node /tmp is gone when the block exits
  echo "--- products (bytes name) ---"
  ls -l "$(dirname "$3")/insar" 2>/dev/null | awk 'NR>1{print $5, $9}' | grep -vE "\.(xml|vrt)$" | head -40
  echo "--- end $1 ---"
}
export -f pfail

run_pair() {
  p=$1
  r=${p%-*}; s=${p#*-}
  cal="20${r:0:2}-${r:2:2}-${r:4:2}_20${s:0:2}-${s:2:2}-${s:4:2}"
  od=$OUT/$cal
  # a pair dir with an isce.log was produced by the banked pairwise campaign
  [ -f "$od/isce.log" ] && { echo "PAIR_SKIP $p (banked pairwise output present)"; return 0; }
  [ -f "$od/DONE" ] && { echo "PAIR_SKIP $p (already done)"; return 0; }
  cd $WORK/pairs/$p || return 1
  L=$WORK/pairs/$p/$p.log
  form_interferogram.py -ref_date $r -sec_date $s -nrlks1 $NRLKS1 -nalks1 $NALKS1 >> $L 2>&1 || { pfail $p form "$L"; return 1; }
  mosaic_interferogram.py -ref_date_stack $REF -ref_date $r -sec_date $s -nrlks1 $NRLKS1 -nalks1 $NALKS1 >> $L 2>&1 || { pfail $p mosaic "$L"; return 1; }
  diff_interferogram.py -idir ../../dates_resampled -ref_date_stack $REF -ref_date $r -sec_date $s -nrlks1 $NRLKS1 -nalks1 $NALKS1 >> $L 2>&1 || { pfail $p diff "$L"; return 1; }
  look_coherence.py -ref_date $r -sec_date $s -nrlks1 $NRLKS1 -nalks1 $NALKS1 -nrlks2 $NRLKS2 -nalks2 $NALKS2 >> $L 2>&1 || { pfail $p look_coherence "$L"; return 1; }
  if [ "${STOP_AFTER:-}" = "look_coherence" ]; then
    # coherence-only mode (era test): geocode the ml2 coherence and stop. No
    # ion, filter, unwrap or conncomp. Writes DONE_COH, never DONE, so a later
    # full run does not mistake this for a finished pair. conv() is defined
    # further down and does not exist yet at this point, so translate inline.
    cd insar || return 1
    f=${p}${ML2}.cor
    [ -s "$f" ] || { echo "PAIR_FAILED $p no $f"; return 1; }
    geocode.py -ref_date_stack_track ../$REF.track.xml -dem $DEM3W -input $f \
      -nrlks $((NRLKS1*NRLKS2)) -nalks $((NALKS1*NALKS2)) -bbox " $BBOX" >> $L 2>&1 \
      || { pfail $p "geocode $f" "$L"; return 1; }
    [ -s "$f.geo.vrt" ] || { echo "PAIR_FAILED $p no $f.geo.vrt"; return 1; }
    mkdir -p $od
    gdal_translate -q -of GTiff -co COMPRESS=DEFLATE -co TILED=YES "$f.geo.vrt" "$od/${cal}${ML2}.cor.geo.aoi90.tif.tmp" >> $L 2>&1 \
      && mv "$od/${cal}${ML2}.cor.geo.aoi90.tif.tmp" "$od/${cal}${ML2}.cor.geo.aoi90.tif" \
      || { pfail $p "convert cor" "$L"; return 1; }
    touch $od/DONE_COH
    echo "PAIR_OK $p (coherence only)"
    return 0
  fi
  if [ "$USE_ION" = "1" ]; then
  ion_correct.py -ion_dir ../../dates_ion -ref_date $r -sec_date $s -nrlks1 $NRLKS1 -nalks1 $NALKS1 -nrlks2 $NRLKS2 -nalks2 $NALKS2 >> $L 2>&1 || { pfail $p ion_correct "$L"; return 1; }
  fi
  filt.py -idir ../../dates_resampled -ref_date_stack $REF -ref_date $r -sec_date $s -nrlks1 $NRLKS1 -nalks1 $NALKS1 -nrlks2 $NRLKS2 -nalks2 $NALKS2 -alpha $FILT_ALPHA -win 32 -step 4 >> $L 2>&1 || { pfail $p filt "$L"; return 1; }
  if [ ! -s "insar/${p}${ML2}.phsig" ]; then
    # ISCE's ICU writes the phase sigma in azimuth patches with a 264-line
    # overlap; on a short image (240 rows at 224 azimuth looks, jobs 821292 and
    # 821322) it writes nothing and snaphu dies with "phsig wrong size".
    # Estimate the phase sigma from the filtered phase instead (Alos2ProcPublic
    # cal_coherence_1, the same phase-only coherence ISCE uses elsewhere) [S].
    # Inert whenever ICU produced the file, so the certified 17x56 path is unchanged.
    echo "PHSIG_FALLBACK $p (ICU wrote an empty phase sigma; estimating from the filtered phase)"
    python3 - "insar/filt_${p}${ML2}.int" "insar/${p}${ML2}.phsig" >> $L 2>&1 <<'PYPHSIG' || { pfail $p phsig "$L"; return 1; }
import sys, numpy as np, xml.etree.ElementTree as ET
from isceobj.Alos2Proc.Alos2ProcPublic import cal_coherence_1
f, o = sys.argv[1], sys.argv[2]
x = ET.parse(f + ".xml").getroot()
w = int(x.find(".//property[@name='width']/value").text); l = int(x.find(".//property[@name='length']/value").text)
a = np.fromfile(f, dtype=np.complex64).reshape(l, w)
cal_coherence_1(a, win=5).astype(np.float32).tofile(o)
print("phsig fallback wrote", l, "x", w)
PYPHSIG
    [ -s "insar/${p}${ML2}.phsig" ] || { echo "PAIR_FAILED $p phsig fallback produced nothing"; return 1; }
  fi
  unwrap_snaphu.py -idir ../../dates_resampled -ref_date_stack $REF -ref_date $r -sec_date $s -nrlks1 $NRLKS1 -nalks1 $NALKS1 -nrlks2 $NRLKS2 -nalks2 $NALKS2 -wbd_msk >> $L 2>&1 || { pfail $p unwrap "$L"; return 1; }
  cd insar || return 1
  for f in filt_${p}${ML2}.unw ${p}${ML2}.cor filt_${p}${ML2}.unw.conncomp; do
    [ -s "$f" ] || { echo "PAIR_FAILED $p no $f"; return 1; }
    geocode.py -ref_date_stack_track ../$REF.track.xml -dem $DEM3W -input $f \
      -nrlks $((NRLKS1*NRLKS2)) -nalks $((NALKS1*NALKS2)) -bbox " $BBOX" >> $L 2>&1 \
      || { pfail $p "geocode $f" "$L"; return 1; }
  done
  mkdir -p $od
  # build_stack2.read_on_grid warps every raster onto its own common grid, so
  # this only has to hand it GeoTIFFs under the names find_pairs globs for.
  # geocode.py writes <file>.geo with a .vrt beside it; translate the .vrt.
  conv() {
    src=$1; dst=$2
    [ -s "$src.vrt" ] || { echo "PAIR_FAILED $p no $src.vrt"; return 1; }
    gdal_translate -q -of GTiff -co COMPRESS=DEFLATE -co TILED=YES "$src.vrt" "$dst.tmp" >> $L 2>&1 || return 1
    mv "$dst.tmp" "$dst"
  }
  conv filt_${p}${ML2}.unw.geo $od/filt_${cal}${ML2}.unw.geo.aoi90.tif || { pfail $p convert unw "$L"; return 1; }
  if [ "${KEEP_WRAPPED:-1}" = "1" ] && [ -s "filt_${p}${ML2}.int" ]; then
    cp filt_${p}${ML2}.int filt_${p}${ML2}.int.xml filt_${p}${ML2}.int.vrt $od/ 2>/dev/null
  fi
  conv ${p}${ML2}.cor.geo $od/${cal}${ML2}.cor.geo.aoi90.tif || { pfail $p convert cor "$L"; return 1; }
  conv filt_${p}${ML2}.unw.conncomp.geo $od/filt_${cal}${ML2}.unw.conncomp.geo.aoi90.tif || { pfail $p convert conncomp "$L"; return 1; }
  touch $od/DONE
  echo "PAIR_OK $p"
}
export -f run_pair
export WORK REF OUT ML1 ML2 DEM3W BBOX NRLKS1 NALKS1 NRLKS2 NALKS2 FILT_ALPHA ION_DIR USE_ION STOP_AFTER

echo "=== $(echo $PAIRS | wc -w) pairs on $LANES lanes at $(date -Is) ==="
printf '%s\n' $PAIRS | xargs -P $LANES -I{} bash -c 'run_pair "$1"' _ {} | tee $WORK/pairs.out
ok=$(grep -c "^PAIR_OK" $WORK/pairs.out || true)
skip=$(grep -c "^PAIR_SKIP" $WORK/pairs.out || true)
bad=$(grep -c "^PAIR_FAILED" $WORK/pairs.out || true)
echo "main block $B: $ok ok, $skip skipped, $bad failed at $(date -Is)"

# the reference LOS geometry is produced by ref_los.sbatch, which writes a
# two-band GeoTIFF into a pair directory under the name build_stack2 globs
# (~/full/pairs/*/ *los.geo.aoi90.tif); copying the raw ISCE .geo here would
# give GDAL a file it cannot open and no .vrt to read it through.
[ "$ok" -gt 0 ] || [ "$skip" -gt 0 ] || { echo "MAIN_ABORT no pair completed"; exit 1; }
echo "MAIN_OK block $B ($ok ok, $skip skipped, $bad failed) at $(date -Is)"
