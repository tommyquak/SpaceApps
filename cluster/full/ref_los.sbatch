#!/bin/bash
#SBATCH --job-name=ref-los
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --tmp=100000
#SBATCH --time=01:00:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --output=/home/t/tommy/full/logs/ref-los-%j.log
# The reference line-of-sight geometry MintPy needs (spec sec 5, FINALE
# INTERFACE). build_stack2 globs ~/full/pairs/*/ *los.geo.aoi90.tif and reads
# band 1 as incidence and the last band as azimuth, so this must be a two-band
# GeoTIFF INSIDE a pair directory, not the raw ISCE .geo.
#
#   stage the reference bundle -> look_geom (ml2) -> geocode with -bbox
#   -> two-band GeoTIFF into the first completed pair directory
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj
export ISCE_STACK=$HOME/micromamba/envs/insar/share/isce2
export PATH=$ISCE_STACK/alosStack:$PATH
export OMP_NUM_THREADS=8

REF=181123
FRAME=3700
NRLKS1=1; NALKS1=14; NRLKS2=5; NALKS2=2
ML2=_$((NRLKS1*NRLKS2))rlks_$((NALKS1*NALKS2))alks
BBOX=${BBOX:-"-3.4/-2.34/113.8/115.2"}
PERSIST=${PERSIST:-/mnt/scratch/t/tommy/full/stack}
OUT=${OUT:-$HOME/full/pairs}
DEM=$HOME/pilot/dem/elevation.dem.wgs84
DEM3=$HOME/pilot/dem/elevation_3s.dem.wgs84
WBD=$HOME/pilot/dem/waterBody.wbd

WORK=/tmp/reflos_${SLURM_JOB_ID}; mkdir -p $WORK/dem
trap "rm -rf $WORK" EXIT
cd $WORK
for f in $DEM $DEM3 $WBD; do
  cp $f $f.xml $WORK/dem/ || { echo "REFLOS_ABORT dem stage $f"; exit 1; }
done
cp $DEM.vrt $DEM3.vrt $WORK/dem/ 2>/dev/null
DEM3W=$WORK/dem/$(basename $DEM3); WBDW=$WORK/dem/$(basename $WBD)
for f in $WORK/dem/$(basename $DEM) $DEM3W $WBDW; do
  fixImageXml.py -i $f -f || { echo "REFLOS_ABORT fixImageXml $f"; exit 1; }
done

mkdir -p dates_resampled/$REF/f1_$FRAME dates_resampled/$REF/insar
cp $PERSIST/dates_resampled/$REF/$REF.track.xml dates_resampled/$REF/ || exit 1
cp $PERSIST/dates_resampled/$REF/f1_$FRAME/$REF.frame.xml dates_resampled/$REF/f1_$FRAME/ || exit 1
cp $PERSIST/dates_resampled/$REF/insar/* dates_resampled/$REF/insar/ || exit 1
for sn in 2 3 4; do mkdir -p dates_resampled/$REF/f1_$FRAME/s$sn; done

echo "=== look_geom (ml2) at $(date -Is) ==="
( cd dates_resampled/$REF && look_geom.py -date $REF -wbd $WBDW \
    -nrlks1 $NRLKS1 -nalks1 $NALKS1 -nrlks2 $NRLKS2 -nalks2 $NALKS2 ) \
  || { echo "REFLOS_ABORT look_geom"; exit 1; }
[ -s "dates_resampled/$REF/insar/${REF}${ML2}.los" ] || { echo "REFLOS_ABORT no ml2 los"; exit 1; }

echo "=== geocode at $(date -Is) ==="
( cd dates_resampled/$REF/insar && geocode.py \
    -ref_date_stack_track ../$REF.track.xml -dem $DEM3W -input ${REF}${ML2}.los \
    -nrlks $((NRLKS1*NRLKS2)) -nalks $((NALKS1*NALKS2)) -bbox " $BBOX" ) \
  || { echo "REFLOS_ABORT geocode"; exit 1; }
G=dates_resampled/$REF/insar/${REF}${ML2}.los.geo
[ -s "$G.vrt" ] || { echo "REFLOS_ABORT no $G.vrt"; exit 1; }

# into the first completed pair directory, under the name find_pairs globs
target=$(ls -d $OUT/*/ 2>/dev/null | head -1)
[ -n "$target" ] || { echo "REFLOS_ABORT no pair directory yet"; exit 1; }
dst="${target%/}/${REF}${ML2}.los.geo.aoi90.tif"
gdal_translate -q -of GTiff -co COMPRESS=DEFLATE -co TILED=YES "$G.vrt" "$dst.tmp" \
  || { echo "REFLOS_ABORT gdal_translate"; exit 1; }
mv "$dst.tmp" "$dst"
python3 - <<PYX || { echo "REFLOS_ABORT the written los is not usable"; exit 1; }
from osgeo import gdal
import numpy as np
d = gdal.Open("$dst")
print(f"    {d.RasterXSize}x{d.RasterYSize}, {d.RasterCount} bands")
assert d.RasterCount >= 2, "los must carry incidence and azimuth bands"
inc = d.GetRasterBand(1).ReadAsArray()
azi = d.GetRasterBand(d.RasterCount).ReadAsArray()
v = inc[np.isfinite(inc) & (inc != 0)]
a = azi[np.isfinite(azi) & (azi != 0)]
print(f"    incidence {v.min():.2f}..{v.max():.2f} deg (n={v.size})")
print(f"    azimuth   {a.min():.2f}..{a.max():.2f} deg")
assert 20 < v.min() and v.max() < 60, "incidence outside the plausible ALOS-2 range"
PYX
echo "REFLOS_OK $dst at $(date -Is)"
