#!/bin/bash
#SBATCH --job-name=nisar-refens
#SBATCH --partition=long
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --output=/home/t/tommy/full/logs/nisar_refens-%j.log
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$PATH PROJ_DATA=$HOME/micromamba/envs/insar/share/proj
export OMP_NUM_THREADS=8 GDAL_CACHEMAX=4096
cd $HOME/full/analysis
: "${TRACK:?set TRACK=A or D}"
python nisar_refens.py && echo "REFENS_DONE $TRACK" || { echo "REFENS_FAIL $TRACK"; exit 1; }
