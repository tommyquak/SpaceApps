#!/bin/bash
#SBATCH --job-name=elninobox
#SBATCH --output=/home/t/tommy/full/logs/elninobox_%j.out
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
export PROJ_LIB=$HOME/micromamba/envs/insar/share/proj TMPDIR=$HOME/full/tmp PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1
cd $HOME/full && [ -s tay/InSAR_MaxSubsidence-2019ElNino.tif ] || { echo "ELNINO_ABORT no raster"; exit 1; }
$HOME/micromamba/envs/insar/bin/python elnino_box.py; echo ELNINO_EXIT=$?
