#!/bin/bash
#SBATCH --job-name=full-pairs
#SBATCH --partition=long
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --tmp=1800000
#SBATCH --time=20:00:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --no-requeue
#SBATCH --output=/home/t/tommy/full/logs/pack-%j.log
# Packed pair job (v3, audit rounds 1-3). Up to 4 alos2App pairs concurrently.
# Accepted residual: a network-partitioned zombie pack can still replace a
# product with a stale-but-complete file under DONE (scancel cannot reach a
# partitioned node); every publish is whole-file atomic so no torn files.
# Terminal contract: pairs/<R>_<S>/DONE only after data-verified products.
# Infra aborts release QUEUED markers so they never charge a science attempt.
set -uo pipefail
cd ~/full
PACK_IDX=("$@")
infra_abort() {
  echo "INFRA_ABORT job=${SLURM_JOB_ID} host=$(hostname) reason: $1"
  local idx line P n
  for idx in "${PACK_IDX[@]}"; do
    line=$(sed -n "${idx}p" pairs_full.txt); P=pairs/${line%% *}_${line##* }
    [ "$(cat $P/QUEUED 2>/dev/null)" = "$SLURM_JOB_ID" ] && rm -f $P/QUEUED
    n=$(( $(cat $P/INFRA 2>/dev/null || echo 0) + 1 )); mkdir -p $P; echo $n > $P/INFRA
    [ $n -ge 6 ] && { touch $P/FAILED; echo "PERMFAIL_INFRA ${line} after $n infra aborts"; }
  done
  exit 1
}
# job-private /tmp is what --tmp actually guarantees (audit r2 finding 7)
WORK=/tmp/pack_${SLURM_JOB_ID}; mkdir -p $WORK
avail=$(df -m $WORK | awk 'NR==2{print $4}')
[ "$avail" -lt 1500000 ] && infra_abort "scratch too small: ${avail}M at $WORK"
echo "workdir $WORK (${avail}M free)"
trap "rm -rf $WORK" EXIT
export PATH=$HOME/micromamba/envs/insar/bin:$HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj

declare -A NEED
PAIRS=()
for idx in "${PACK_IDX[@]}"; do
  line=$(sed -n "${idx}p" pairs_full.txt); REF=${line%% *}; SEC=${line##* }
  [ -f pairs/${REF}_${SEC}/DONE ] && { echo "pair ${REF}_${SEC} already done"; continue; }
  PAIRS+=("$REF $SEC"); NEED[$REF]=1; NEED[$SEC]=1
  # self-register as the authoritative worker (audit r2 findings 4/6)
  mkdir -p pairs/${REF}_${SEC}; echo ${SLURM_JOB_ID} > pairs/${REF}_${SEC}/QUEUED
done
[ ${#PAIRS[@]} -eq 0 ] && { echo "nothing to do"; exit 0; }
mkdir -p $WORK/ceos
CPFAIL=0
for d in "${!NEED[@]}"; do
  [ -f ceos/$d/.complete ] || infra_abort "scene $d missing or incomplete"
  cp -r ceos/$d $WORK/ceos/ &
done
for j in $(jobs -p); do wait $j || CPFAIL=1; done
[ $CPFAIL -eq 1 ] && infra_abort "scene staging copy failed" 
cp ~/pilot/dem/elevation.dem.wgs84* ~/pilot/dem/elevation.dem.hdr $WORK/ 2>/dev/null || true
cp ~/pilot/dem/waterBody.wbd ~/pilot/dem/waterBody.wbd.xml $WORK/

render_xml() {
  sed -e "s|REFDIR|$WORK/ceos/$1|; s|SECDIR|$WORK/ceos/$2|; s|DEMPATH|$WORK/elevation.dem.wgs84|; s|WBDPATH|$WORK/waterBody.wbd|" ~/full/alos2App.template.xml > alos2App.xml
}

run_pair() {
  local REF=$1 SEC=$2
  local PD=$WORK/p_${REF}_${SEC}
  local OUT=~/full/pairs/${REF}_${SEC}
  [ -f $OUT/DONE ] && { echo "[${REF}_${SEC}] became DONE, skipping"; return 0; }
  mkdir -p $PD $OUT
  cd $PD
  render_xml $REF $SEC
  # Coregistration mode by temporal baseline: amplitude matching fails on swamp
  # beyond ~3 months; long pairs go straight to geometric offsets (precise orbits).
  local dt=$(( ( $(date -d "$SEC" +%s) - $(date -d "$REF" +%s) ) / 86400 ))
  local geom=""
  [ $dt -gt 98 ] && geom=1 && echo "[${REF}_${SEC}] dt=${dt}d: geometric offsets"
  ALOS2_FORCE_GEOM_OFFSET=$geom OMP_NUM_THREADS=12 python3 $HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications/alos2App.py alos2App.xml > alos2app.stdout 2>&1
  local rc=$?
  # Retry geometrically on (a) the deterministic matching crash, (b) the silent
  # cull-fallback that would otherwise ship garbage offsets with rc=0.
  local retry=0
  if [ -z "$geom" ]; then
    if [ $rc -ne 0 ] && grep -qE "zero.size array|not enough offsets left|Bad match at level" alos2app.stdout; then retry=1; fi
    if [ $rc -eq 0 ] && grep -q "use offset without culling" alos2app.stdout; then
      echo "[${REF}_${SEC}] ampcor culling failed (raw offsets) at dt=${dt}d, retrying geometrically"
      retry=1
    fi
  fi
  if [ $retry -eq 1 ]; then
    cp $PD/isce.log $OUT/isce.log.try1 2>/dev/null; cp $PD/alos2app.stdout $OUT/alos2app.stdout.try1 2>/dev/null
    cd $WORK; rm -rf $PD; mkdir -p $PD; cd $PD
    render_xml $REF $SEC
    ALOS2_FORCE_GEOM_OFFSET=1 OMP_NUM_THREADS=12 python3 $HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications/alos2App.py alos2App.xml > alos2app.stdout 2>&1
    rc=$?
  fi
  for F in isce.log alos2app.stdout; do
    [ -f $PD/$F ] && cp $PD/$F $OUT/.tmp.$F.$SLURM_JOB_ID 2>/dev/null && mv -f $OUT/.tmp.$F.$SLURM_JOB_ID $OUT/$F || true
  done
  echo "[${REF}_${SEC}] alos2App exit: $rc (dt=${dt}d geom=${geom:-0} retried=$retry)"
  # AOI products at 3 arcsec: averaged for continuous fields, NEAREST for the
  # integer conncomp labels (never blend labels; 0 is a meaningful class).
  local nprod=0
  for G in $PD/insar/*.geo; do
    [ -e "$G" ] || continue
    local B=$(basename $G)
    local ALG=average NODATA="srcNodata=0, dstNodata=0,"
    case "$B" in *conncomp*) ALG=near; NODATA="";; esac
    python3 -c "
from osgeo import gdal
gdal.UseExceptions()
gdal.Warp(\"$OUT/.tmp.$B.$SLURM_JOB_ID\", \"$G\", format=\"GTiff\",
  outputBounds=(113.8, -3.4, 115.2, -2.34),
  xRes=0.000833333333333333, yRes=0.000833333333333333,
  targetAlignedPixels=True, resampleAlg=\"$ALG\", $NODATA
  creationOptions=[\"COMPRESS=DEFLATE\", \"TILED=YES\"])
" && mv -f "$OUT/.tmp.$B.$SLURM_JOB_ID" "$OUT/$B.aoi90.tif" && nprod=$((nprod+1)) || echo "[${REF}_${SEC}] crop failed for $B"
  done
  for X in $PD/insar/*.geo.xml; do
    [ -e "$X" ] || continue
    local XB=$(basename $X)
    cp $X $OUT/.tmp.$XB.$SLURM_JOB_ID 2>/dev/null && mv -f $OUT/.tmp.$XB.$SLURM_JOB_ID $OUT/$XB || true
  done
  # Data gate, not exit-code gate: unwrapped phase must contain signal.
  local zerocheck=1
  local UNW=$(ls $OUT/filt_*.unw.geo.aoi90.tif 2>/dev/null | grep -v _msk | grep -v conncomp | head -1)
  if [ -n "$UNW" ]; then
    zerocheck=$(python3 -c "
from osgeo import gdal; import numpy as np
gdal.UseExceptions()
ds = gdal.Open(\"$UNW\")
b = ds.GetRasterBand(2) if ds.RasterCount > 1 else ds.GetRasterBand(1)
a = b.ReadAsArray()
print(0 if np.count_nonzero(a) > 1000 else 1)
ds = None" 2>>$OUT/zerocheck.err) || zerocheck=1
    [ -n "$zerocheck" ] || zerocheck=1
  fi
  if [ $rc -eq 0 ] && [ $nprod -ge 4 ] && [ "$zerocheck" = "0" ] && ! grep -q "use offset without culling" alos2app.stdout; then
    touch $OUT/DONE; echo "PAIR_DONE ${REF}_${SEC}"
  else
    echo "PAIR_FAILED ${REF}_${SEC} rc=$rc nprod=$nprod zero=$zerocheck"
  fi
}
for p in "${PAIRS[@]}"; do run_pair $p & done
wait
echo PACK_COMPLETE
