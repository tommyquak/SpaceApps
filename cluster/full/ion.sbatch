#!/bin/bash
#SBATCH --job-name=ion
#SBATCH --partition=long
#SBATCH --cpus-per-task=24
#SBATCH --mem=60G
#SBATCH --tmp=500000
# Sized to the 72-node pool of 24-CPU / 63 GB machines; a 48-CPU 180 GB request
# only fits the 14 large nodes and queues on "(Resources)".
#SBATCH --time=12:00:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --no-requeue
#SBATCH --output=/home/t/tommy/full/logs/ion-%A_%a.log
# alosStack pivot, ION PASS (spec sec 5, staging protocol 5b). One job per date
# block: stage the block's crops node-local once, then run every ion pair whose
# two dates are both in the block.
#
#   stage bundle E+F+G + manifest-C xmls + the reference skeleton
#   -> unpack crops (full band AND the authentic lower/upper subbands)
#   -> render every .slc.xml/.vrt and the rg_rect.off xmls BEFORE pair_up
#   -> pair_up -odir pairs_ion
#   -> per pair: form + mosaic + diff (cmd_2 scope), ion_subband, ion_unwrap,
#      ion_filt, then the 80x448 _ori.int check artefact
#   -> persist ONLY ion/ion_cal per pair plus one <ref>.track.xml
#
# Blocks come from make_blocks.py: ~/full/ion_blocks/<i>.{dates,pairs}
#   sbatch --array=0-4 ion.sbatch
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj
export ISCE_STACK=$HOME/micromamba/envs/insar/share/isce2
export PATH=$ISCE_STACK/alosStack:$PATH
export OMP_NUM_THREADS=2                     # lanes carry the parallelism

REF=181123
FRAME=3700
NRLKS1=1; NALKS1=14; NRLKS2=5; NALKS2=2
NRLKS_ION=80; NALKS_ION=32
LANES=${LANES:-12}
CP=$HOME/full/crop_pack.py
PERSIST=${PERSIST:-/mnt/scratch/t/tommy/full/stack}
ION_OUT=${ION_OUT:-/mnt/scratch/t/tommy/full/pairs_ion}
DEM=$HOME/pilot/dem/elevation.dem.wgs84
WBD=$HOME/pilot/dem/waterBody.wbd
B=${SLURM_ARRAY_TASK_ID:-0}
DATES_F=$HOME/full/ion_blocks/$B.dates
PAIRS_F=$HOME/full/ion_blocks/$B.pairs
for f in $DATES_F $PAIRS_F; do [ -s "$f" ] || { echo "ION_ABORT no $f"; exit 1; }; done
DATES=$(cat $DATES_F); PAIRS=$(cat $PAIRS_F)
echo "ion block $B: $(echo $DATES | wc -w) dates, $(echo $PAIRS | wc -w) pairs on $(hostname) at $(date -Is)"

WORK=/tmp/ion_${B}_${SLURM_JOB_ID}; mkdir -p $WORK/dem
avail=$(df -m $WORK | awk 'NR==2{print $4}')
[ "$avail" -lt 900000 ] && { echo "ION_ABORT scratch too small: ${avail}M"; exit 1; }
trap "rm -rf $WORK" EXIT
cd $WORK

# ------------------------------------------------------------- DEM / water
for f in $DEM $WBD; do
  [ -f "$f" ] && [ -f "$f.xml" ] || { echo "ION_ABORT missing $f or its xml"; exit 1; }
  cp $f $f.xml $WORK/dem/ || { echo "ION_ABORT dem stage failed for $f"; exit 1; }
done
cp $DEM.vrt $WORK/dem/ 2>/dev/null
DEMW=$WORK/dem/$(basename $DEM); WBDW=$WORK/dem/$(basename $WBD)
for f in $DEMW $WBDW; do
  fixImageXml.py -i $f -f || { echo "ION_ABORT fixImageXml failed for $f"; exit 1; }
done

# ------------------------------------------------- once-bundle E + F + G
mkdir -p dates_resampled/$REF/f1_$FRAME dates_resampled/$REF/insar dates/$REF/f1_$FRAME/mosaic
cp $PERSIST/dates_resampled/$REF/$REF.track.xml dates_resampled/$REF/ || exit 1
cp $PERSIST/dates_resampled/$REF/f1_$FRAME/$REF.frame.xml dates_resampled/$REF/f1_$FRAME/ || exit 1
cp $PERSIST/dates_resampled/$REF/insar/* dates_resampled/$REF/insar/ || exit 1
cp $PERSIST/dates/$REF/$REF.track.xml dates/$REF/ || exit 1
cp $PERSIST/dates/$REF/f1_$FRAME/$REF.frame.xml dates/$REF/f1_$FRAME/ || exit 1
cp $PERSIST/dates/$REF/f1_$FRAME/mosaic/swath_offset.txt dates/$REF/f1_$FRAME/mosaic/ || exit 1
# the reference date's directory skeleton must exist even when the block does
# not contain it: stackDateStatistics enumerates frames and swaths from it
for sn in 2 3 4; do mkdir -p dates/$REF/f1_$FRAME/s$sn dates_resampled/$REF/f1_$FRAME/s$sn; done

# ------------------------------------------------------- unpack the block
# every date the block touches, plus the reference itself
STAGE=$(printf '%s\n' $DATES $REF | sort -u)
echo "=== unpack $(echo $STAGE | wc -w) dates at $(date -Is) ==="
unpack_one() {
  d=$1
  src=$PERSIST/dates_resampled/$d/f1_$FRAME
  [ -d "$src" ] || { echo "ION_MISSING $d has no crops"; return 1; }
  mkdir -p dates_resampled/$d/f1_$FRAME dates_resampled/$d/insar dates/$d/f1_$FRAME
  cp $PERSIST/dates/$d/$d.track.xml dates/$d/ || return 1
  cp $PERSIST/dates/$d/f1_$FRAME/$d.frame.xml dates/$d/f1_$FRAME/ || return 1
  for sn in 2 3 4; do
    od=dates_resampled/$d/f1_$FRAME/s$sn; mkdir -p $od
    pk=$src/s$sn/$d.slc
    [ -s "$pk" ] || { echo "ION_MISSING $d s$sn"; return 1; }
    python3 $CP unpack $pk $src/s$sn/$d.slc.json $od/$d.slc || return 1
    w=$(python3 -c "import json;s=json.load(open('$src/s$sn/$d.slc.json'));print(s['width'],s['length'])")
    python3 $CP render $od/$d.slc ${w% *} ${w#* } CFLOAT 1 BIL || return 1
    if [ -s "$src/s$sn/${d}_lower.slc" ]; then
      # wave 1 persisted the authentic subbands (PERSIST_SUBBANDS=1)
      for band in _lower _upper; do
        python3 $CP unpack $src/s$sn/${d}${band}.slc $src/s$sn/${d}${band}.slc.json $od/${d}${band}.slc || return 1
      done
    else
      # regenerate them from the crop and restore the demodulation ramp (sec 3)
      python3 $HOME/full/ion_regen.py $od/$d.slc $src/s$sn/$d.slc.json \
        $WORK/dates_resampled/$REF $REF || return 1
    fi
    for band in _lower _upper; do
      [ -s "$od/${d}${band}.slc" ] || { echo "ION_MISSING $d s$sn $band"; return 1; }
      python3 $CP render $od/${d}${band}.slc ${w% *} ${w#* } CFLOAT 1 BIL || return 1
    done
  done
  # manifest B: the cropped rg_rect.off and its xml (the reference has none)
  if [ -s "$PERSIST/dates_resampled/$d/insar/${d}_${NRLKS1}rlks_${NALKS1}alks_rg_rect.off" ]; then
    cp $PERSIST/dates_resampled/$d/insar/${d}_${NRLKS1}rlks_${NALKS1}alks_rg_rect.off* \
       dates_resampled/$d/insar/ || return 1
  fi
  echo "unpacked $d"
}
export -f unpack_one
export PERSIST FRAME CP NRLKS1 NALKS1 WORK REF
# GNU parallel is not installed on this cluster; xargs -P is
printf '%s\n' $STAGE | xargs -P ${UNPACK_LANES:-16} -I{} bash -c 'unpack_one "$1"' _ {} \
  || { echo "ION_ABORT unpack failed"; exit 1; }

# the frame/track metadata of the cropped reference drives every mosaic step
mosaic_parameter.py -idir dates_resampled -ref_date $REF -sec_date $REF \
  -nrlks1 $NRLKS1 -nalks1 $NALKS1 >/dev/null || { echo "ION_ABORT mosaic_parameter"; exit 1; }

# ------------------------------------------------------------------ pairs
echo "=== pair_up at $(date -Is) ==="
pair_up.py -idir1 dates -idir2 dates_resampled -odir pairs_ion -ref_date $REF -pairs $PAIRS \
  || { echo "ION_ABORT pair_up"; exit 1; }

run_pair() {
  p=$1
  r=${p%-*}; s=${p#*-}
  cd $WORK/pairs_ion/$p || return 1
  form_interferogram.py -ref_date $r -sec_date $s -nrlks1 $NRLKS1 -nalks1 $NALKS1 >> $p.log 2>&1 || { echo "PAIR_FAILED $p form"; return 1; }
  mosaic_interferogram.py -ref_date_stack $REF -ref_date $r -sec_date $s -nrlks1 $NRLKS1 -nalks1 $NALKS1 >> $p.log 2>&1 || { echo "PAIR_FAILED $p mosaic"; return 1; }
  diff_interferogram.py -idir ../../dates_resampled -ref_date_stack $REF -ref_date $r -sec_date $s -nrlks1 $NRLKS1 -nalks1 $NALKS1 >> $p.log 2>&1 || { echo "PAIR_FAILED $p diff"; return 1; }
  ion_subband.py -idir ../../dates_resampled -ref_date_stack $REF -ref_date $r -sec_date $s -nrlks1 $NRLKS1 -nalks1 $NALKS1 >> $p.log 2>&1 || { echo "PAIR_FAILED $p ion_subband"; return 1; }
  ion_unwrap.py -idir ../../dates_resampled -ref_date_stack $REF -ref_date $r -sec_date $s -wbd $WBDW -nrlks1 $NRLKS1 -nalks1 $NALKS1 -nrlks_ion $NRLKS_ION -nalks_ion $NALKS_ION >> $p.log 2>&1 || { echo "PAIR_FAILED $p ion_unwrap"; return 1; }
  ion_filt.py -idir ../../dates_resampled -idir2 ../../dates -ref_date_stack $REF -ref_date $r -sec_date $s -nrlks1 $NRLKS1 -nalks1 $NALKS1 -nrlks2 $NRLKS2 -nalks2 $NALKS2 -nrlks_ion $NRLKS_ION -nalks_ion $NALKS_ION -win_min 11 -win_max 301 -fit -filt -fit_adaptive -filt_secondary -win_secondary 5 >> $p.log 2>&1 || { echo "PAIR_FAILED $p ion_filt"; return 1; }
  # check artefact: the pair's own diff looked down to the ion grid
  looks.py -i insar/diff_${p}_${NRLKS1}rlks_${NALKS1}alks.int \
    -o ion/ion_cal/diff_${p}_$((NRLKS1*NRLKS_ION))rlks_$((NALKS1*NALKS_ION))alks_ori.int \
    -r $NRLKS_ION -a $NALKS_ION >> $p.log 2>&1 || echo "PAIR_WARN $p ori.int"
  # persist only the ion_cal products (KBs) and, once, the reference track
  out=$ION_OUT/$p; mkdir -p $out/ion
  cp -r ion/ion_cal $out/ion/ || { echo "PAIR_FAILED $p persist"; return 1; }
  cp $REF.track.xml $out/ 2>/dev/null
  touch $out/DONE
  echo "PAIR_OK $p"
}
export -f run_pair
export WORK REF WBDW ION_OUT NRLKS1 NALKS1 NRLKS2 NALKS2 NRLKS_ION NALKS_ION

echo "=== $(echo $PAIRS | wc -w) pairs on $LANES lanes at $(date -Is) ==="
printf '%s\n' $PAIRS | xargs -P $LANES -I{} bash -c 'run_pair "$1"' _ {} | tee $WORK/pairs.out
ok=$(grep -c "^PAIR_OK" $WORK/pairs.out || true)
bad=$(grep -c "^PAIR_FAILED" $WORK/pairs.out || true)
echo "ion block $B: $ok ok, $bad failed at $(date -Is)"
[ "$ok" -gt 0 ] || { echo "ION_ABORT no pair produced ion_cal"; exit 1; }
python3 $HOME/full/check_run.py ion $B || { echo "ION_ABORT block $B failed the validity check"; exit 1; }
touch $ION_OUT/BLOCK_${B}_DONE
echo "ION_OK block $B ($ok/$((ok+bad)) pairs) at $(date -Is)"
