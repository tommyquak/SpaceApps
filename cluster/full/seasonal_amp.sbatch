#!/bin/bash
#SBATCH --job-name=seasamp
#SBATCH --output=/home/t/tommy/full/logs/seasamp_%j.out
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
export PROJ_LIB=$HOME/micromamba/envs/insar/share/proj TMPDIR=$HOME/full/tmp PYTHONUNBUFFERED=1 OMP_NUM_THREADS=4
cd $HOME/full && $HOME/micromamba/envs/insar/bin/python seasonal_amp.py; echo SEASAMP_EXIT=$?
