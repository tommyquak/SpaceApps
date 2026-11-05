#!/bin/bash
#SBATCH --job-name=seed
#SBATCH --partition=long
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --tmp=1500000
#SBATCH --time=10:00:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --no-requeue
#SBATCH --output=/home/t/tommy/full/logs/seed-%j.log
# alosStack pivot, SEED (spec sec 5): the one full-frame job that builds the
# reference geometry every later wave depends on, then publishes the cropped
# reference bundle (manifest E, F, G) that gates wave 1.
#
#   read_data(ref + neighbour) -> estimate_slc_offset -> swath/frame offsets
#   -> resample_common_grid(both) -> mosaic_parameter -> rdr2geo
#   -> PUBLISH gate (swath_offset.txt + ref track)
#   -> pair form/mosaic -> radar_dem_offset -> rect_range_offset -> look_geom
#   -> crop layer: geometry -> derive -> pack --reference -> rewrite
#      -> mosaic_parameter re-run -> verify -> bundle F crops -> publish
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj
export ISCE_STACK=$HOME/micromamba/envs/insar/share/isce2   # Amendment A.10
export PATH=$ISCE_STACK/alosStack:$PATH
export OMP_NUM_THREADS=48

REF=181123
NEI=181012                      # neighbour date, for the offsets and the affine
FRAME=3700
NRLKS1=1; NALKS1=14; NRLKS2=5; NALKS2=2
CP=$HOME/full/crop_pack.py
PERSIST=${PERSIST:-/mnt/scratch/t/tommy/full/stack}   # crop store lives on scratch: home has a ~425 GB quota        # persist root; wave-1 and block jobs read this
DEM=$HOME/pilot/dem/elevation.dem.wgs84
WBD=$HOME/pilot/dem/waterBody.wbd

WORK=/tmp/seed_${SLURM_JOB_ID}; mkdir -p $WORK/data $WORK/dem
avail=$(df -m $WORK | awk 'NR==2{print $4}')
[ "$avail" -lt 900000 ] && { echo "SEED_ABORT scratch too small: ${avail}M"; exit 1; }
echo "seed workdir $WORK (${avail}M free) on $(hostname) at $(date -Is)"
trap "rm -rf $WORK" EXIT

# ---------------------------------------------------------------- staging
for d in $REF $NEI; do
  case $d in
    181123) src=$HOME/full/ceos/2018-11-23 ;;
    181012) src=$HOME/full/ceos/2018-10-12 ;;
    *) echo "SEED_ABORT no CEOS mapping for $d"; exit 1 ;;
  esac
  [ -f "$src/.complete" ] || { echo "SEED_ABORT $src not extracted"; exit 1; }
  mkdir -p $WORK/data/$d
  cp $src/* $WORK/data/$d/ || { echo "SEED_ABORT stage failed for $d"; exit 1; }
done
# DEM/water staging rule (spec 5b, r1 finding 1): copy, never '|| true', and
# fixImageXml every image so ISCE resolves the absolute path
for f in $DEM $WBD; do
  [ -f "$f" ] && [ -f "$f.xml" ] || { echo "SEED_ABORT missing $f or its xml"; exit 1; }
  cp $f $f.xml $WORK/dem/ || { echo "SEED_ABORT dem stage failed for $f"; exit 1; }
done
cp $DEM.vrt $WORK/dem/ 2>/dev/null
for f in $WORK/dem/$(basename $DEM) $WORK/dem/$(basename $WBD); do
  fixImageXml.py -i $f -f || { echo "SEED_ABORT fixImageXml failed for $f"; exit 1; }
done
DEMW=$WORK/dem/$(basename $DEM); WBDW=$WORK/dem/$(basename $WBD)

cd $WORK
cat > alosStack.xml <<XMLEOF
<?xml version="1.0" encoding="UTF-8"?>
<stack>
<component name="stackinsar">
    <property name="data directory">$WORK/data</property>
    <property name="dem for coregistration">$DEMW</property>
    <property name="dem for geocoding">$DEMW</property>
    <property name="water body">$WBDW</property>
    <property name="reference date of the stack">$REF</property>
    <property name="polarization">HH</property>
    <property name="starting swath">2</property>
    <property name="ending swath">4</property>
    <property name="number of subsequent dates">1</property>
</component>
</stack>
XMLEOF
create_cmds.py -stack_par alosStack.xml || { echo "SEED_ABORT create_cmds failed"; exit 1; }

# --------------------------------------------------- per-date (cmd_1 scope)
echo "=== cmd_1 (read_data .. rdr2geo) at $(date -Is) ==="
t0=$(date +%s)
bash ./cmd_1.sh > cmd1.out 2>&1 || { echo "SEED_ABORT cmd_1 failed"; tail -40 cmd1.out; exit 1; }
echo "cmd_1 ok in $(( ($(date +%s)-t0)/60 )) min"
REFI=dates_resampled/$REF/insar/${REF}_${NRLKS1}rlks_${NALKS1}alks
for f in $REFI.lat $REFI.lon $REFI.hgt $REFI.los $REFI.wbd; do
  [ -s "$f" ] || { echo "SEED_ABORT rdr2geo did not produce $f"; exit 1; }
done

# ------------------------------------------------------- publish the gate
# wave-1 jobs need only swath_offset.txt and a loadable full-frame ref track
mkdir -p $PERSIST/dates/$REF/f1_$FRAME/mosaic
cp dates/$REF/$REF.track.xml $PERSIST/dates/$REF/ || { echo "SEED_ABORT publish track"; exit 1; }
cp dates/$REF/f1_$FRAME/$REF.frame.xml $PERSIST/dates/$REF/f1_$FRAME/ || { echo "SEED_ABORT publish frame"; exit 1; }
cp dates/$REF/f1_$FRAME/mosaic/swath_offset.txt $PERSIST/dates/$REF/f1_$FRAME/mosaic/ \
  || { echo "SEED_ABORT publish swath_offset"; exit 1; }
python3 -c "
import os, sys
sys.path.insert(0, os.path.expanduser('~/full'))
import crop_pack as cp
loadTrack, = cp._stackpulic('loadTrack')
t = loadTrack('$PERSIST/dates/$REF', '$REF')
print('published ref track loads:', len(t.frames[0].swaths), 'swaths')
" || { echo "SEED_ABORT published ref track does not load"; exit 1; }
touch $PERSIST/GATE_REF_PUBLISHED
echo "SEED_GATE_OPEN at $(date -Is) (wave 1 may start)"

# ------------------------------- pair, affine, rect and ml2 reference geometry
echo "=== pair $NEI-$REF, affine, rect, look_geom at $(date -Is) ==="
PAIR=$NEI-$REF
pair_up.py -idir1 dates -idir2 dates_resampled -odir pairs -ref_date $REF -pairs $PAIR \
  || { echo "SEED_ABORT pair_up"; exit 1; }
( cd pairs/$PAIR && form_interferogram.py -ref_date $NEI -sec_date $REF -nrlks1 $NRLKS1 -nalks1 $NALKS1 ) \
  || { echo "SEED_ABORT form_interferogram"; exit 1; }
( cd pairs/$PAIR && mosaic_interferogram.py -ref_date_stack $REF -ref_date $NEI -sec_date $REF \
    -nrlks1 $NRLKS1 -nalks1 $NALKS1 ) || { echo "SEED_ABORT mosaic_interferogram"; exit 1; }
( cd dates_resampled/$REF && radar_dem_offset.py -track $REF.track.xml -dem $DEMW \
    -wbd insar/${REF}_${NRLKS1}rlks_${NALKS1}alks.wbd -hgt insar/${REF}_${NRLKS1}rlks_${NALKS1}alks.hgt \
    -amp ../../pairs/$PAIR/insar/${PAIR}_${NRLKS1}rlks_${NALKS1}alks.amp \
    -output affine_transform.txt -nrlks1 $NRLKS1 -nalks1 $NALKS1 ) \
  || { echo "SEED_ABORT radar_dem_offset"; exit 1; }
AFF=dates_resampled/$REF/insar/affine_transform.txt
[ -s "$AFF" ] || AFF=dates_resampled/$REF/affine_transform.txt
[ -s "$AFF" ] || { echo "SEED_ABORT no affine_transform.txt"; exit 1; }
( cd dates_resampled/$NEI/insar && rect_range_offset.py -aff $(cd $WORK && readlink -f $AFF) \
    -input ${NEI}_${NRLKS1}rlks_${NALKS1}alks_rg.off \
    -output ${NEI}_${NRLKS1}rlks_${NALKS1}alks_rg_rect.off -nrlks1 $NRLKS1 -nalks1 $NALKS1 ) \
  || { echo "SEED_ABORT rect_range_offset"; exit 1; }
( cd dates_resampled/$REF && look_geom.py -date $REF -wbd $WBDW -nrlks1 $NRLKS1 -nalks1 $NALKS1 \
    -nrlks2 $NRLKS2 -nalks2 $NALKS2 ) || { echo "SEED_ABORT look_geom"; exit 1; }

# ------------------------------------------------------------- crop layer
# Amendment A.4 order: geometry -> derive -> pack --reference -> rewrite ->
# mosaic_parameter re-run -> verify. Nothing may run out of order: the packing
# and the polynomials need the PRE-crop resampled xmls.
echo "=== crop layer at $(date -Is) ==="
python3 $CP geometry dates dates_resampled $REF || { echo "SEED_ABORT crop_pack geometry"; exit 1; }
read ML1W ML1L < <(python3 -c "
import json; g = json.load(open('$HOME/full/ref_geometry.json'))
print(g['frame_ml1']['width'], g['frame_ml1']['length'])")
echo "frame ml1 ${ML1W}x${ML1L}"
python3 $CP derive $REFI.lat $REFI.lon $ML1W $ML1L || { echo "SEED_ABORT crop_pack derive"; exit 1; }

mkdir -p $PERSIST/dates_resampled/$REF/f1_$FRAME
for sn in 2 3 4; do
  sd=$PERSIST/dates_resampled/$REF/f1_$FRAME/s$sn; mkdir -p $sd
  python3 $CP pack dates_resampled/$REF/f1_$FRAME/s$sn/$REF.slc \
    $sd/$REF.slc $sd/$REF.slc.json $sn --reference || { echo "SEED_ABORT pack s$sn"; exit 1; }
done
python3 $CP rewrite dates_resampled/$REF $REF || { echo "SEED_ABORT sec-2 rewrite"; exit 1; }
mosaic_parameter.py -idir dates_resampled -ref_date $REF -sec_date $REF \
  -nrlks1 $NRLKS1 -nalks1 $NALKS1 || { echo "SEED_ABORT mosaic_parameter re-run"; exit 1; }
python3 $CP verify dates_resampled/$REF $REF || { echo "SEED_ABORT crop verify"; exit 1; }

# manifest E: the cropped-rewritten resampled reference xmls
cp dates_resampled/$REF/$REF.track.xml $PERSIST/dates_resampled/$REF/ || exit 1
cp dates_resampled/$REF/f1_$FRAME/$REF.frame.xml $PERSIST/dates_resampled/$REF/f1_$FRAME/ || exit 1
# manifest F: the ml1 geometry cropped at the ml1 window, with xmls
read ROW ROWS COL COLS < <(python3 -c "
import json; w = json.load(open('$HOME/full/crop_window.json'))['ml1']
print(w['rowOff'], w['rows'], w['colOff'], w['cols'])")
echo "ml1 crop window rows ${ROW}+${ROWS} cols ${COL}+${COLS}"
GD=$PERSIST/dates_resampled/$REF/insar; mkdir -p $GD
for spec in "lat double 1" "lon double 1" "hgt double 1" "wbd byte 1" "los float32 2"; do
  set -- $spec; ext=$1; dt=$2; nb=$3
  python3 $CP cropoff $REFI.$ext $GD/${REF}_${NRLKS1}rlks_${NALKS1}alks.$ext \
    $ROW $ROWS $COL $COLS $ML1W $ML1L $nb $dt || { echo "SEED_ABORT crop $ext"; exit 1; }
  python3 $CP render $GD/${REF}_${NRLKS1}rlks_${NALKS1}alks.$ext $COLS $ROWS $dt $nb BIL \
    || { echo "SEED_ABORT render $ext"; exit 1; }
done
cp $(cd $WORK && readlink -f $AFF) $GD/affine_transform.txt || exit 1
cp $HOME/full/ref_geometry.json $HOME/full/crop_window.json $PERSIST/ || exit 1

# validity check on the published bundle and the reference crops
python3 $HOME/full/check_run.py seed || { echo "SEED_ABORT bundle failed the validity check"; exit 1; }
du -sm $PERSIST 2>/dev/null
echo "SEED_OK at $(date -Is)"
