#!/bin/bash
#SBATCH --job-name=strip-B
#SBATCH --partition=long
#SBATCH --cpus-per-task=16
#SBATCH --mem=90G
#SBATCH --time=12:00:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --output=/home/t/tommy/full/logs/strip_B-%j.log
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj ISCE_STACK=$HOME/micromamba/envs/insar/share/isce2
export PATH=$ISCE_STACK/alosStack:$PATH PYTHONPATH=$HOME/full/strip_test:$ISCE_STACK/alosStack:${PYTHONPATH:-}
cd $HOME/full/strip_test
export OMP_NUM_THREADS=16
python strip_stackB.py || { echo STRIP_B_FAIL_STACK; exit 1; }
cd mintpyB && smallbaselineApp.py smallbaselineApp.cfg --start reference_point > smallbaselineApp.log 2>&1; rc=$?
tail -3 smallbaselineApp.log
[ $rc -eq 0 ] && grep -q "Normal end" smallbaselineApp.log || { echo "STRIP_B_FAIL_MINTPY rc=$rc"; grep -nE "Traceback|Error" smallbaselineApp.log | tail -4; exit 1; }
mkdir -p taycmp/inputs; ln -sfn $PWD/geo/geo_velocity.h5 taycmp/velocity.h5; ln -sfn $PWD/geo/geo_temporalCoherence.h5 taycmp/temporalCoherence.h5; ln -sfn $PWD/geo/geo_geometryRadar.h5 taycmp/inputs/geometryGeo.h5
echo "=== compare_tay: variant B (radar-coordinate inversion, no deramp) ==="; cd $HOME/full && SKIP_RUN_GUARD=1 W=$HOME/full/strip_test/mintpyB/taycmp python compare_tay.py; echo "COMPARE_B_EXIT=$?"
echo STRIP_B_DONE
