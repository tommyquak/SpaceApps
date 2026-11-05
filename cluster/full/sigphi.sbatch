#!/bin/bash
#SBATCH --job-name=sigphi
#SBATCH --output=/home/t/tommy/full/logs/sigphi_%j.out
#SBATCH --error=/home/t/tommy/full/logs/sigphi_%j.out
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
set -u
export PROJ_LIB=$HOME/micromamba/envs/insar/share/proj
export TMPDIR=$HOME/full/tmp
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONUNBUFFERED=1
cd "$HOME/full" || exit 1
$HOME/micromamba/envs/insar/bin/python sigphi_tay.py
echo "SIGPHI_EXIT=$?"
