#!/bin/bash
#SBATCH --job-name=season
#SBATCH --output=/home/t/tommy/full/logs/season_%j.out
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
export PROJ_LIB=$HOME/micromamba/envs/insar/share/proj TMPDIR=$HOME/full/tmp PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1
cd $HOME/full && $HOME/micromamba/envs/insar/bin/python season_coh.py; echo SEASON_EXIT=$?
