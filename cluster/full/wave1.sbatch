#!/bin/bash
#SBATCH --job-name=wave1
#SBATCH --partition=long
#SBATCH --cpus-per-task=24
#SBATCH --mem=60G
#SBATCH --tmp=250000
# Sized to the cluster's real shape (measured 2026-09-03): 72 idle nodes carry
# 24 CPUs / 63 GB and only 14 carry 128 CPUs / 256 GB, so a 48-CPU 180 GB lane
# could only ever run on a handful of nodes and sat in "(Resources)". A lane's
# measured peak is MaxRSS 51.7 GB and, without the three-band resample, well
# under 100 GB of node-local scratch.
#SBATCH --time=12:00:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --no-requeue
#SBATCH --output=/home/t/tommy/full/logs/wave1-%A_%a.log
# alosStack pivot, WAVE 1 (spec sec 5): one job per date lane set. Full frame
# on node-local scratch, then crop and persist per manifest A/B/C.
#
#   wait for the SEED gate -> stage CEOS + ref bundle -> read_data -virtual
#   -> mosaic_parameter(date) -> estimate_slc_offset -> resample_common_grid
#   -> geo2rdr vs the full-frame ref geometry -> wait for affine
#   -> rect_range_offset -> crop + pack (sec 1-2) -> persist (sec 4)
#
# Submit as an array over lines of ~/full/wave1_dates.txt (one date per line,
# "YYMMDD YYYY-MM-DD"), e.g. sbatch --array=1-20%6 wave1.sbatch
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj
export ISCE_STACK=$HOME/micromamba/envs/insar/share/isce2
export PATH=$ISCE_STACK/alosStack:$PATH
export OMP_NUM_THREADS=24

REF=181123
FRAME=3700
NRLKS1=1; NALKS1=14
ML1=_${NRLKS1}rlks_${NALKS1}alks
CP=$HOME/full/crop_pack.py
PERSIST=${PERSIST:-/mnt/scratch/t/tommy/full/stack}   # crop store lives on scratch: home has a ~425 GB quota
DEM=$HOME/pilot/dem/elevation.dem.wgs84
WBD=$HOME/pilot/dem/waterBody.wbd
GATE_WAIT=7200      # s to wait for the SEED gate before giving up
AFF_WAIT=10800      # s to wait for affine_transform.txt (SEED publishes it late)

WAVE1_LIST=${WAVE1_LIST:-$HOME/full/wave1_dates.txt}   # override to feed a different date list (era test)
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $WAVE1_LIST)
[ -n "$line" ] || { echo "WAVE1_ABORT no line ${SLURM_ARRAY_TASK_ID}"; exit 1; }
D=${line%% *}; CAL=${line##* }
[ "$D" = "$REF" ] && { echo "WAVE1_SKIP reference date is seeded"; exit 0; }
echo "wave1 date $D ($CAL) on $(hostname) at $(date -Is)"

# ------------------------------------------------------------- gate wait
t0=$(date +%s)
until [ -f $PERSIST/GATE_REF_PUBLISHED ]; do
  [ $(( $(date +%s) - t0 )) -gt $GATE_WAIT ] && { echo "WAVE1_ABORT gate never opened"; exit 1; }
  sleep 60
done
for f in $PERSIST/dates/$REF/$REF.track.xml $PERSIST/dates/$REF/f1_$FRAME/$REF.frame.xml \
         $PERSIST/dates/$REF/f1_$FRAME/mosaic/swath_offset.txt; do
  [ -s "$f" ] || { echo "WAVE1_ABORT gate is open but $f is missing"; exit 1; }
done
echo "gate ok at $(date -Is)"

WORK=/tmp/wave1_${D}_${SLURM_JOB_ID}; mkdir -p $WORK/data/$D $WORK/dem
avail=$(df -m $WORK | awk 'NR==2{print $4}')
[ "$avail" -lt 700000 ] && { echo "WAVE1_ABORT scratch too small: ${avail}M"; exit 1; }
trap "rm -rf $WORK" EXIT
cd $WORK

# ---------------------------------------------------------------- staging
src=$HOME/full/ceos/$CAL
[ -f "$src/.complete" ] || { echo "WAVE1_ABORT $src not extracted"; exit 1; }
cp $src/* $WORK/data/$D/ || { echo "WAVE1_ABORT stage failed"; exit 1; }
# read_data.py enumerates the stack's dates from the LED files under -idir and
# refuses to run when the reference is not among them (read_data.py:118-124),
# and it cross-checks frame counts against the reference. Only -sec_date is
# actually processed, so the reference is linked, not copied: no bytes move.
REFCAL=2018-11-23
[ -f "$HOME/full/ceos/$REFCAL/.complete" ] || { echo "WAVE1_ABORT reference CEOS missing"; exit 1; }
mkdir -p $WORK/data/$REF
ln -sf $HOME/full/ceos/$REFCAL/* $WORK/data/$REF/ || { echo "WAVE1_ABORT reference link failed"; exit 1; }
for f in $DEM $WBD; do
  [ -f "$f" ] && [ -f "$f.xml" ] || { echo "WAVE1_ABORT missing $f or its xml"; exit 1; }
  cp $f $f.xml $WORK/dem/ || { echo "WAVE1_ABORT dem stage failed for $f"; exit 1; }
done
cp $DEM.vrt $WORK/dem/ 2>/dev/null
DEMW=$WORK/dem/$(basename $DEM); WBDW=$WORK/dem/$(basename $WBD)
for f in $DEMW $WBDW; do
  fixImageXml.py -i $f -f || { echo "WAVE1_ABORT fixImageXml failed for $f"; exit 1; }
done
# ------------------------------------------------------------ per-date work
echo "=== read_data at $(date -Is) ==="
# BOTH dates: estimate_slc_offset ampcors this date against the reference SLC
# (estimate_slc_offset.py:256 loads dates/<ref>/f1_*/s*/<ref>.slc.xml), so the
# reference needs its raw entry too. -virtual makes it a VRT over the CEOS, so
# no pixels are copied and the reference's own CEOS stays a symlink.
read_data.py -idir $WORK/data -odir dates -ref_date $REF -sec_date $D $REF -pol HH \
  -starting_swath 2 -ending_swath 4 -virtual || { echo "WAVE1_ABORT read_data"; exit 1; }
mosaic_parameter.py -idir dates -ref_date $REF -sec_date $D $REF -nrlks1 $NRLKS1 -nalks1 $NALKS1 \
  || { echo "WAVE1_ABORT mosaic_parameter"; exit 1; }
# the persisted reference metadata is canonical (SEED derived the geometry from
# it) and resample_common_grid reads swath_offset.txt unconditionally
mkdir -p dates/$REF/f1_$FRAME/mosaic dates_resampled/$REF/f1_$FRAME
cp $PERSIST/dates/$REF/$REF.track.xml dates/$REF/ || exit 1
cp $PERSIST/dates/$REF/f1_$FRAME/$REF.frame.xml dates/$REF/f1_$FRAME/ || exit 1
cp $PERSIST/dates/$REF/f1_$FRAME/mosaic/swath_offset.txt dates/$REF/f1_$FRAME/mosaic/ || exit 1
# ampcor fails beyond roughly three months from the reference; the patched
# geometric branch is armed for those dates (spec sec 5, wave 1)
DD=$(python3 -c "
import datetime
a = datetime.datetime.strptime('$D', '%y%m%d'); b = datetime.datetime.strptime('$REF', '%y%m%d')
print(abs((a - b).days))")
FORCE=0; [ "$DD" -gt 90 ] && FORCE=1
echo "=== estimate_slc_offset (|dt| ${DD}d, ALOS2_FORCE_GEOM_OFFSET=$FORCE) at $(date -Is) ==="
ALOS2_FORCE_GEOM_OFFSET=$FORCE estimate_slc_offset.py -idir dates -ref_date $REF -sec_date $D \
  -wbd $WBDW -dem $DEMW -use_wbd_offset || { echo "WAVE1_ABORT estimate_slc_offset"; exit 1; }
for sn in 2 3 4; do
  [ -s dates/$D/f1_$FRAME/s$sn/cull.off ] || { echo "WAVE1_ABORT no cull.off for s$sn"; exit 1; }
done
echo "=== resample_common_grid at $(date -Is) ==="
# NO -subband: with the ion pass regenerating subbands from the crops (sec 3),
# the lower/upper bands resample_common_grid would write here are discarded with
# the node-local tree. Dropping the flag resamples one band instead of three,
# the most expensive step in the job (measured 34 min -> about 12 min per date).
# Set SUBBAND_FLAG=-subband to restore the three-band resample.
resample_common_grid.py -idir dates -odir dates_resampled -ref_date $REF -sec_date $D \
  -nrlks1 $NRLKS1 -nalks1 $NALKS1 ${SUBBAND_FLAG:-} || { echo "WAVE1_ABORT resample_common_grid"; exit 1; }
echo "=== geo2rdr at $(date -Is) ==="
RG=$PERSIST/dates_resampled/$REF/insar/${REF}_${NRLKS1}rlks_${NALKS1}alks
( cd dates_resampled/$D && geo2rdr.py -date $D -date_par_dir ../../dates/$D \
    -lat $RG.lat -lon $RG.lon -hgt $RG.hgt -nrlks1 $NRLKS1 -nalks1 $NALKS1 ) \
  || { echo "WAVE1_ABORT geo2rdr"; exit 1; }

# ------------------------------------------------------- affine wait + rect
t0=$(date +%s)
AFF=$PERSIST/dates_resampled/$REF/insar/affine_transform.txt
until [ -s "$AFF" ]; do
  [ $(( $(date +%s) - t0 )) -gt $AFF_WAIT ] && { echo "WAVE1_ABORT affine never published"; exit 1; }
  sleep 60
done
( cd dates_resampled/$D/insar && rect_range_offset.py -aff $AFF \
    -input ${D}_${NRLKS1}rlks_${NALKS1}alks_rg.off \
    -output ${D}_${NRLKS1}rlks_${NALKS1}alks_rg_rect.off -nrlks1 $NRLKS1 -nalks1 $NALKS1 ) \
  || { echo "WAVE1_ABORT rect_range_offset"; exit 1; }

# ------------------------------------------------------------ crop + persist
# the crop window and the reference geometry are SEED products; they arrive
# with the crop layer, after the gate
t0=$(date +%s)
until [ -s $PERSIST/crop_window.json ] && [ -s $PERSIST/ref_geometry.json ]; do
  [ $(( $(date +%s) - t0 )) -gt $AFF_WAIT ] && { echo "WAVE1_ABORT crop window never published"; exit 1; }
  sleep 60
done
cp $PERSIST/crop_window.json $PERSIST/ref_geometry.json $HOME/full/ 2>/dev/null
echo "=== crop and pack at $(date -Is) ==="
OD=$PERSIST/dates_resampled/$D/f1_$FRAME
for sn in 2 3 4; do
  sd=$OD/s$sn; mkdir -p $sd
  python3 $CP pack dates_resampled/$D/f1_$FRAME/s$sn/$D.slc $sd/$D.slc $sd/$D.slc.json $sn \
    --cull dates/$D/f1_$FRAME/s$sn/cull.off || { echo "WAVE1_ABORT pack s$sn"; exit 1; }
done
# AUTHENTIC SUBBANDS: OFF by default. Persisting them (~5.8 GB/date, 540 GB for
# 93 dates) does not fit the measured ~512 GB scratch quota beside the 270 GB
# crop store, so the ion pass regenerates them per spec sec 3. Round-14 finding
# 20 (the 128-sample rg_filter transient at interior crop edges) is bounded for
# THIS window: the derived crop is s2 1630+2810 (right edge = the swath edge
# 4440), s3 0+5670 (both edges are swath edges) and s4 0+5308, so the only
# interior edges are s2's left and s4's right, and both coincide with the frame
# crop boundary (1630 and 1630+13090 = 14720). The transient therefore lands on
# the outermost ~1.1 km of the product, not inside the AOI, and at nrlks_ion=80
# it is under 2 ion pixels before ion_filt smooths over 11-301. Set
# PERSIST_SUBBANDS=1 if the storage ever allows the exact path.
if [ "${PERSIST_SUBBANDS:-0}" = "1" ] && [ -n "${SUBBAND_FLAG:-}" ]; then
  for sn in 2 3 4; do
    for band in lower upper; do
      src=dates_resampled/$D/f1_$FRAME/s$sn/${D}_${band}.slc
      [ -s "$src" ] || { echo "WAVE1_ABORT $band subband missing for s$sn"; exit 1; }
      sd=$OD/s$sn
      python3 $CP pack $src $sd/${D}_${band}.slc $sd/${D}_${band}.slc.json $sn \
        --cull dates/$D/f1_$FRAME/s$sn/cull.off --date $D \
        || { echo "WAVE1_ABORT pack $band s$sn"; exit 1; }
    done
  done
  echo "authentic subbands persisted ($(du -sm $OD | cut -f1) MB for the date)"
fi

# manifest B: rg_rect.off. geo2rdr ran against the CROPPED reference geometry
# (the only geometry SEED persists), so the offsets are already on the cropped
# ml1 grid and must NOT be cropped again. They are still the right numbers: both
# dates sit on the common resampled grid and take the identical window, so the
# relative offset is unchanged by the crop (sec 4B, "window-relative; valid
# under pure crop"). Assert the dims instead of transforming them.
read ROWS COLS < <(python3 -c "
import json
m = json.load(open('$PERSIST/crop_window.json'))['ml1']; print(m['rows'], m['cols'])")
ID=$PERSIST/dates_resampled/$D/insar; mkdir -p $ID
SRC=dates_resampled/$D/insar/${D}${ML1}_rg_rect.off
[ -s "$SRC" ] || { echo "WAVE1_ABORT no rg_rect.off"; exit 1; }
want=$((ROWS * COLS * 4))
got=$(stat -c %s "$(readlink -f $SRC)")
[ "$got" = "$want" ] || { echo "WAVE1_ABORT rg_rect.off is $got B, expected $want ($COLS x $ROWS float32 on the cropped ml1 grid)"; exit 1; }
cp "$(readlink -f $SRC)" $ID/${D}${ML1}_rg_rect.off || { echo "WAVE1_ABORT persist rg_rect.off"; exit 1; }
python3 $CP render $ID/${D}${ML1}_rg_rect.off $COLS $ROWS float32 1 BIL \
  || { echo "WAVE1_ABORT render rg_rect.off"; exit 1; }

# manifest C: the full-frame date xmls, verbatim
mkdir -p $PERSIST/dates/$D/f1_$FRAME
cp dates/$D/$D.track.xml $PERSIST/dates/$D/ || exit 1
cp dates/$D/f1_$FRAME/$D.frame.xml $PERSIST/dates/$D/f1_$FRAME/ || exit 1

# validity check on the finished date: content, not presence. No DONE without it.
python3 $HOME/full/check_run.py date $D || { echo "WAVE1_ABORT $D failed the validity check"; exit 1; }
du -sm $PERSIST/dates_resampled/$D
touch $PERSIST/dates_resampled/$D/DONE
# the crops are validated and persisted, so this date's CEOS can go: home holds
# a ~23-date window and every later date needs a slot (spec sec 5, wave 1
# logistics, net -13.5 GB/date). Only ever after DONE, and never the reference.
if [ "$D" != "$REF" ] && [ -f $PERSIST/dates_resampled/$D/DONE ]; then
  # ~/full/ceos holds real directories for home-resident dates and symlinks for
  # the ones on scratch, so release the target as well as the link
  tgt=$(readlink -f $HOME/full/ceos/$CAL 2>/dev/null)
  [ -n "$tgt" ] && [ "$tgt" != "/" ] && rm -rf "$tgt"
  rm -rf $HOME/full/ceos/$CAL
  echo "released CEOS $CAL ($(ls $HOME/full/ceos 2>/dev/null | wc -l) dates left)"
fi
echo "WAVE1_OK $D at $(date -Is)"
