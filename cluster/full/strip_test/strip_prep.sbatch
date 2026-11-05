#!/bin/bash
#SBATCH --job-name=strip-prep
#SBATCH --partition=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --output=/home/t/tommy/full/logs/strip_prep-%j.log
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj ISCE_STACK=$HOME/micromamba/envs/insar/share/isce2
export PATH=$ISCE_STACK/alosStack:$PATH PYTHONPATH=$HOME/full/strip_test:$ISCE_STACK/alosStack:${PYTHONPATH:-}
cd $HOME/full/strip_test
export OMP_NUM_THREADS=4
python strip_prep.py && echo STRIP_PREP_DONE || { echo STRIP_PREP_FAIL; exit 1; }
