#!/bin/bash
#SBATCH --job-name=alos2-dem
#SBATCH --partition=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --exclude=xcng[0-1],xgpj0
#SBATCH --output=/home/t/tommy/pilot/logs/dem.log
set -euo pipefail
export MAMBA_ROOT_PREFIX=$HOME/micromamba
[ -f ~/pilot/dem/.dem_complete ] && { echo "DEM marker present, skipping"; echo DEM_DONE; exit 0; }
mkdir -p ~/pilot/dem && cd ~/pilot/dem
export TMPDIR=$HOME/pilot/dem/tmp; mkdir -p $TMPDIR
~/bin/micromamba run -n insar sardem --bbox 112 -5 117 0 --data-source COP -isce --output elevation.dem.wgs84 --output-format ENVI --output-type float32
echo "--- water mask (best effort):"
~/bin/micromamba run -n insar sardem --bbox 112 -5 117 0 --data-source NASA_WATER -o watermask.wbd || echo "water mask failed, continuing without"
touch ~/pilot/dem/.dem_complete
ls -la ~/pilot/dem | grep -v tmp
echo DEM_DONE
