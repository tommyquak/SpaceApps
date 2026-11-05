#!/bin/bash
#SBATCH --job-name=strip-geoA
#SBATCH --partition=normal
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --output=/home/t/tommy/full/logs/strip_geoA-%A_%a.log
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj ISCE_STACK=$HOME/micromamba/envs/insar/share/isce2
export PATH=$ISCE_STACK/alosStack:$PATH PYTHONPATH=$HOME/full/strip_test:$ISCE_STACK/alosStack:${PYTHONPATH:-}
cd $HOME/full/strip_test
export OMP_NUM_THREADS=1
B=${SLURM_ARRAY_TASK_ID:-0}; NB=${NBLOCKS:-8}
WORK=/tmp/strip_geoA_$B; mkdir -p $WORK/dem; DEM3=$HOME/pilot/dem/elevation_3s.dem.wgs84
cp $DEM3 $DEM3.xml $WORK/dem/ && cp $DEM3.vrt $WORK/dem/ 2>/dev/null; DEM3W=$WORK/dem/$(basename $DEM3); fixImageXml.py -i $DEM3W -f > /dev/null
TRACK=/mnt/scratch/t/tommy/full/pairs_ion/190621-190830/181123.track.xml; BBOX="-3.4/-2.34/113.8/115.2"
PAIRS=$(python -c "from strip_common import pair_list; print(' '.join(pair_list()))" | tr ' ' '\n' | awk -v b=$B -v nb=$NB 'NR % nb == b')
echo "block $B: $(echo $PAIRS | wc -w) pairs"
printf '%s\n' $PAIRS | xargs -P 12 -I{} python strip_geoA.py {} $TRACK $DEM3W "$BBOX" | tee /tmp/strip_geoA_$B.out
ok=$(grep -c "^GEOA_OK\|^SKIP" /tmp/strip_geoA_$B.out); bad=$(grep -c "GEOCODE_FAIL\|Traceback" /tmp/strip_geoA_$B.out)
rm -rf $WORK; echo "block $B: ok=$ok failed=$bad"; [ "$bad" -eq 0 ] && echo "STRIP_GEOA_BLOCK_OK $B" || { echo "STRIP_GEOA_BLOCK_FAIL $B"; exit 1; }
