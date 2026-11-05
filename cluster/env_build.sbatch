#!/bin/bash
#SBATCH --job-name=insar-env
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:30:00
#SBATCH --output=%x-%j.log
echo "node: $(hostname)"
echo "outbound test:"; curl -sI --max-time 10 https://conda.anaconda.org 2>&1 | head -1
export MAMBA_ROOT_PREFIX=$HOME/micromamba
~/bin/micromamba create -y -n insar -c conda-forge python=3.11 isce2 mintpy asf_search gdal "numpy<2"
~/bin/micromamba run -n insar python -c "import isce; import mintpy; print(\"ISCE2+MintPy OK\")"
