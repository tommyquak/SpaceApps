#!/bin/bash
#SBATCH --job-name=zapply
#SBATCH --output=/home/t/tommy/full/logs/zapply_%j.out
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
export PROJ_LIB=$HOME/micromamba/envs/insar/share/proj TMPDIR=$HOME/full/tmp PYTHONUNBUFFERED=1 OMP_NUM_THREADS=4
cd $HOME/full && $HOME/micromamba/envs/insar/bin/python zheng_apply_bw60.py; echo ZAPPLY_EXIT=$?
