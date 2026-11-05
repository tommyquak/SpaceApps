#!/bin/bash
#SBATCH --job-name=tropo
#SBATCH --output=/home/t/tommy/full/logs/tropo_%j.out
#SBATCH --error=/home/t/tommy/full/logs/tropo_%j.out
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
set -u
export PROJ_LIB=$HOME/micromamba/envs/insar/share/proj TMPDIR=$HOME/full/tmp PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
cd "$HOME/full" || exit 1
W=${W:-$HOME/full/mintpy_tay} $HOME/micromamba/envs/insar/bin/python tay_tropo.py
echo "TROPO_EXIT=$?"
