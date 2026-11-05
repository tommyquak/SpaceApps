#!/bin/bash
#SBATCH --job-name=strip-A
#SBATCH --partition=long
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --output=/home/t/tommy/full/logs/strip_A-%j.log
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj ISCE_STACK=$HOME/micromamba/envs/insar/share/isce2
export PATH=$ISCE_STACK/alosStack:$PATH PYTHONPATH=$HOME/full/strip_test:$ISCE_STACK/alosStack:${PYTHONPATH:-}
cd $HOME/full/strip_test
export OMP_NUM_THREADS=16
W=$HOME/full/strip_test/mintpyA; mkdir -p $W; export W PAIR_ROOT=/mnt/scratch/t/tommy/strip_test/pairsA
cd $HOME/full && python build_stack2.py $W || { echo STRIP_A_FAIL_STACK; exit 1; }
cd $W && cat > smallbaselineApp.cfg <<CFG
mintpy.load.processor = isce
mintpy.reference.method = maxCoherence
mintpy.reference.coherenceFile = avgSpatialCoh.h5
mintpy.reference.minCoherence = 0.5
mintpy.reference.maskFile = no
mintpy.network.coherenceBased = no
mintpy.unwrapError.method = no
mintpy.networkInversion.weightFunc = var
mintpy.networkInversion.minTempCoh = 0.4
mintpy.networkInversion.maskDataset = no
mintpy.troposphericDelay.method = no
mintpy.deramp = linear
mintpy.deramp.maskFile = maskTempCoh.h5
mintpy.solidEarthTides = yes
mintpy.topographicResidual = no
mintpy.geocode = no
mintpy.plot = no
CFG
smallbaselineApp.py smallbaselineApp.cfg --start reference_point > smallbaselineApp.log 2>&1; rc=$?
tail -3 smallbaselineApp.log
[ $rc -eq 0 ] && grep -q "Normal end" smallbaselineApp.log || { echo "STRIP_A_FAIL_MINTPY rc=$rc"; grep -nE "Traceback|Error" smallbaselineApp.log | tail -4; exit 1; }
echo "=== compare_tay: variant A' (geocode-first, averaged, deramped; our current order on the same re-unwrapped crops) ==="; cd $HOME/full && SKIP_RUN_GUARD=1 W=$HOME/full/strip_test/mintpyA python compare_tay.py; echo "COMPARE_A_EXIT=$?"
echo STRIP_A_DONE
