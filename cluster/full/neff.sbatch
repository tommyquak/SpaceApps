#!/bin/bash
#SBATCH --job-name=neff
#SBATCH --output=/home/t/tommy/full/logs/neff_%j.out
#SBATCH --error=/home/t/tommy/full/logs/neff_%j.out
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
set -u
export PROJ_LIB=$HOME/micromamba/envs/insar/share/proj TMPDIR=$HOME/full/tmp PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$HOME/full" || exit 1
echo "pair dirs: $(ls -d $HOME/full/pairs_tay/*/ 2>/dev/null | wc -l)"
ls $HOME/full/pairs_tay/$(ls $HOME/full/pairs_tay | head -1)/insar/ 2>/dev/null | head -30
$HOME/micromamba/envs/insar/bin/python neff.py
echo "NEFF_EXIT=$?"
