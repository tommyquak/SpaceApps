#!/bin/bash
#SBATCH --job-name=p1-tay
#SBATCH --partition=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --output=/home/t/tommy/full/logs/p1_tay-%j.log
# Pattern comparison of the 2007-2009 PALSAR-1 rate against Tay et al. 2025
# (2018-2024) on the overlap. Different eras, so magnitude agreement is NOT
# expected (Zahro 2025: -1.74 cm/yr before restoration vs -0.29 after); the
# question is whether the spatial pattern of subsidence persists (sign,
# correlation with block bootstrap). compare_tay.py expects a MintPy dir
# with geocoded velocity.h5, so the geo/ products are linked under taycmp/.
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$PATH
export PROJ_LIB=$HOME/micromamba/envs/insar/share/proj TMPDIR=$HOME/full/tmp PYTHONUNBUFFERED=1 OMP_NUM_THREADS=4
TAG=${TAG:-p422_f7130}
M=$HOME/full/palsar1_mintpy_$TAG
[ -f $M/geo/geo_velocity.h5 ] || { echo "P1_TAY_FAIL no geocoded velocity in $M"; exit 1; }
T=$M/taycmp; mkdir -p $T/inputs
ln -sfn $M/geo/geo_velocity.h5 $T/velocity.h5
ln -sfn $M/geo/geo_temporalCoherence.h5 $T/temporalCoherence.h5
[ -f $M/geo/geo_geometryRadar.h5 ] && ln -sfn $M/geo/geo_geometryRadar.h5 $T/inputs/geometryGeo.h5
cd $HOME/full
echo "=== compare_tay: PALSAR-1 $TAG (2007-2009) vs Tay 2025 (2018-2024), pattern only ==="
SKIP_RUN_GUARD=1 W=$T python compare_tay.py; echo "COMPARE_EXIT=$?"
echo "P1_TAY_DONE $TAG"
