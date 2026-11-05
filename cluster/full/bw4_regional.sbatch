#!/bin/bash
#SBATCH --job-name=bw4reg
#SBATCH --output=/home/t/tommy/full/logs/bw4reg_%j.out
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
export PROJ_LIB=$HOME/micromamba/envs/insar/share/proj TMPDIR=$HOME/full/tmp PYTHONUNBUFFERED=1 OMP_NUM_THREADS=4
cd $HOME/full && $HOME/micromamba/envs/insar/bin/python bw4_regional.py; echo BW4REG_EXIT=$?
