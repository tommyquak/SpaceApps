#!/bin/bash
#SBATCH --job-name=nisar-refsel
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --output=/home/t/tommy/full/logs/nisar_refsel-%j.log
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$PATH PROJ_DATA=$HOME/micromamba/envs/insar/share/proj
export OMP_NUM_THREADS=4 GDAL_CACHEMAX=2048
cd $HOME/full/analysis
python nisar_refsel.py && echo REFSEL_DONE || { echo REFSEL_FAIL; exit 1; }
