#!/bin/bash
#SBATCH --job-name=nisar1
#SBATCH --output=/home/t/tommy/full/logs/nisar1_%j.out
#SBATCH --error=/home/t/tommy/full/logs/nisar1_%j.out
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
set -u
export PROJ_LIB=$HOME/micromamba/envs/insar/share/proj TMPDIR=$HOME/full/tmp PYTHONUNBUFFERED=1
D=/mnt/scratch/t/tommy/nisar/gunw; mkdir -p "$D"; cd "$D" || exit 1
[ -f "$HOME/.netrc" ] || { echo "NISAR_ABORT no ~/.netrc"; exit 1; }
n=0
while read -r u; do
  [ -n "$u" ] || continue
  f=$(basename "$u")
  if [ -s "$f" ]; then echo "have $f"; n=$((n+1)); continue; fi
  echo "GET $f"
  wget -q --netrc --auth-no-challenge --tries=3 --timeout=120 -O "$f.part" "$u" && mv "$f.part" "$f" && n=$((n+1)) || { echo "DOWNLOAD_FAILED $f"; rm -f "$f.part"; }
done < "$HOME/full/gunw_urls.txt"
echo "downloaded $n files"; du -sh "$D"
cd "$HOME/full" && $HOME/micromamba/envs/insar/bin/python nisar_stage1.py
echo "NISAR1_EXIT=$?"
