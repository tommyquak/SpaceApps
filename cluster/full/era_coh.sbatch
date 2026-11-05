#!/bin/bash
#SBATCH --job-name=eracoh
#SBATCH --output=/home/t/tommy/full/logs/eracoh_%j.out
#SBATCH --error=/home/t/tommy/full/logs/eracoh_%j.out
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
set -u
export PROJ_LIB=$HOME/micromamba/envs/insar/share/proj TMPDIR=$HOME/full/tmp PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$HOME/full" || exit 1
$HOME/micromamba/envs/insar/bin/python era_coh.py
echo "ERACOH_EXIT=$?"
