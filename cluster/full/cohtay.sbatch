#!/bin/bash
#SBATCH --job-name=cohtay
#SBATCH --output=/home/t/tommy/full/logs/cohtay_%j.out
#SBATCH --error=/home/t/tommy/full/logs/cohtay_%j.out
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
set -u
export PROJ_LIB=$HOME/micromamba/envs/insar/share/proj
export TMPDIR=$HOME/full/tmp
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$TMPDIR" "$HOME/full/logs"
cd "$HOME/full" || exit 1
$HOME/micromamba/envs/insar/bin/python coh_tbase_tay.py
echo "COHTAY_EXIT=$?"
