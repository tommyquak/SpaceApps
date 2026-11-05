#!/bin/bash
#SBATCH --job-name=alos2-mintpy
#SBATCH --partition=long
#SBATCH --exclude=xcng[0-1],xgpj0
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/home/t/tommy/pilot/logs/mintpy.log
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj
export OMP_NUM_THREADS=16
# Wait for the fanout to finish: all 70 pairs, or 9 hours, whichever first.
for i in $(seq 1 108); do
  done_count=$(ls ~/pilot/pairs/*/DONE 2>/dev/null | wc -l)
  [ "$done_count" -ge 67 ] && break
  echo "waiting: $done_count/70 pairs ($i/108)"
  sleep 300
done
done_count=$(ls ~/pilot/pairs/*/DONE 2>/dev/null | wc -l)
echo "proceeding with $done_count pairs"
[ "$done_count" -lt 55 ] && { echo "too few pairs, aborting for human review"; exit 1; }
W=~/pilot/mintpy
rm -rf $W; mkdir -p $W
python ~/pilot/build_stack.py $W
cd $W
cat > pilotApp.cfg <<CFGEOF
mintpy.load.processor = isce
mintpy.reference.method = maxCoherence
mintpy.reference.coherenceFile = avgSpatialCoh.h5
mintpy.network.coherenceBased = no
mintpy.networkInversion.weightFunc = var
mintpy.networkInversion.maskDataset = coherence
mintpy.networkInversion.maskThreshold = 0.25
mintpy.troposphericDelay.method = no
mintpy.deramp = linear
mintpy.topographicResidual = no
mintpy.velocity.startDate = auto
CFGEOF
smallbaselineApp.py pilotApp.cfg --dir $W --start modify_network --end velocity
echo "--- outputs:"
ls -la $W/*.h5 2>/dev/null | head
for f in velocity temporalCoherence avgSpatialCoh; do
  [ -f $W/$f.h5 ] && save_gdal.py $W/$f.h5 -o $W/$f.tif 2>/dev/null && echo "exported $f.tif"
done
[ -f $W/velocity.h5 ] && echo MINTPY_OK || echo MINTPY_FAILED
