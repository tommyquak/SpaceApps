#!/bin/bash
#SBATCH --job-name=alos2-dl
#SBATCH --partition=normal
#SBATCH --array=1-20%6
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:45:00
#SBATCH --exclude=xcng[0-1],xgpj0
#SBATCH --output=/home/t/tommy/pilot/logs/dl-%a.log
set -uo pipefail
cd ~/pilot
. ~/pilot/lib_checks.sh
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" urls.txt)
date_tag=${line%% *}; url=${line##* }; name=$(basename "$url" .zip)
dest=ceos/$date_tag
if scene_complete "$dest"; then echo "[$date_tag] already extracted"; exit 0; fi
WORK=/temp/${SLURM_JOB_ID:-$$}; mkdir -p $WORK 2>/dev/null || WORK=/tmp/dl_$$; mkdir -p $WORK
zip_src=""
[ -s "zips/$name.zip" ] && zip_src="zips/$name.zip" && echo "[$date_tag] using home zip"
if [ -z "$zip_src" ]; then
  echo "[$date_tag] downloading to node-local $WORK"
  ~/micromamba/envs/insar/bin/python ~/pilot/dl_ranged.py "$url" "$WORK/$name.zip" 12 || { echo "[$date_tag] DOWNLOAD_FAILED"; rm -rf $WORK; exit 1; }
  zip_src="$WORK/$name.zip"
fi
# Staged extraction: $dest appears only once the full member set is verified,
# so a wall-time kill mid-unzip cannot leave a partial scene that passes the
# scene_complete skip check above.
tmp=ceos/.tmp_$date_tag
rm -rf "$dest" "$tmp"; mkdir -p "$tmp"
echo "[$date_tag] selective extract HH F2-F4"
if unzip -q -o "$zip_src" "VOL-*" "LED-*" "TRL-*" "IMG-HH-*-F2" "IMG-HH-*-F3" "IMG-HH-*-F4" -d "$tmp/" && scene_complete "$tmp"; then
  mv "$tmp" "$dest"
  n=$(ls $dest | wc -l)
  echo "[$date_tag] extracted $n files, $(du -sm $dest | cut -f1) MB"
  rm -f "zips/$name.zip"; rm -rf $WORK
  echo "[$date_tag] DONE"
else
  echo "[$date_tag] EXTRACT_FAILED"; rm -rf "$tmp" $WORK; [ "$zip_src" = "zips/$name.zip" ] && rm -f "$zip_src"; exit 1
fi
