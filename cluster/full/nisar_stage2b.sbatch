#!/bin/bash
#SBATCH --job-name=nisar2b
#SBATCH --output=/home/t/tommy/full/logs/nisar2b_%j.out
#SBATCH --error=/home/t/tommy/full/logs/nisar2b_%j.out
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
set -u
export PROJ_LIB=$HOME/micromamba/envs/insar/share/proj TMPDIR=$HOME/full/tmp PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$HOME/full" || exit 1
[ -f dam_sites.geojson ] || { echo "NISAR2B_ABORT no dam_sites.geojson"; exit 1; }
$HOME/micromamba/envs/insar/bin/python nisar_stage2b.py
echo "NISAR2B_EXIT=$?"
