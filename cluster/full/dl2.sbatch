#!/bin/bash
#SBATCH --job-name=full-dl
#SBATCH --partition=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:45:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --output=/home/t/tommy/full/logs/dl-%A_%a.log
set -uo pipefail
cd ~/full
# Measured quotas: home ~425 GB, scratch ~512 GB (the RPC quota service is down,
# so only a failed write reveals them). CEOS is split across both and the crop
# store needs 270 GB of scratch, so cap the staged window at 30 dates.
WINDOW=${WINDOW:-18}
SCENES=${SCENES:-scenes_full.txt}   # override to feed a different date list (era test)
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SCENES)
date_tag=${line%% *}; url=${line##* }; name=$(basename "$url" .zip)
dest=ceos/$date_tag
PART="$dest.part"
[ -f "$dest/.complete" ] && { echo "[$date_tag] already extracted"; exit 0; }
# window guard applies only when this scene adds a NEW dir. Wave-1 releases a
# date's CEOS once its crops are validated, so WAIT for a slot rather than
# refusing: refusing here strands the remaining dates once the window fills.
if [ ! -d "$dest" ] && [ ! -d "$PART" ]; then
  t0=$(date +%s)
  while [ "$(ls -d ceos/*/ 2>/dev/null | wc -l)" -ge $WINDOW ]; do
    if [ $(( $(date +%s) - t0 )) -gt 21600 ]; then
      echo "[$date_tag] window still full after 6 h, giving up"; exit 1
    fi
    sleep 300
  done
fi
WORK=/tmp/dl_${SLURM_JOB_ID:-$$}; mkdir -p $WORK
trap 'rm -rf "$PART" "$WORK"' EXIT
# Selective fetch first: every member of these granules is stored uncompressed,
# so the CEOS files we keep (VOL/LED/TRL + IMG-HH F2-F4) are 16.3 GiB of the
# 56.6 GiB granule and can be pulled by byte range. Each member is checked
# against its CRC-32 from the zip directory, so a bad fetch fails loudly.
# Falls back to the whole-zip downloader on any failure.
SELECTIVE=${SELECTIVE:-1}
# fetch to NODE-LOCAL disk, then move: writing 16 GB straight to the NFS home
# from several lanes at once starves everything else on that filesystem
if [ "$SELECTIVE" = "1" ] && ~/micromamba/envs/insar/bin/python ~/full/dl_selective.py get "$url" "$WORK/sel" 12; then
  echo "[$date_tag] selective fetch ok"
  rm -rf "$dest" "$PART"; mkdir -p "$PART"
  # NOTE: the shell glob skips dotfiles, so the marker is written after the move,
  # not carried with it; a date without .complete is invisible to the feeder
  mv $WORK/sel/* "$PART/" && touch "$PART/.complete" && mv "$PART" "$dest"
  echo "[$date_tag] extracted $(ls $dest | wc -l) files, $(du -sm $dest | cut -f1) MB"
  echo "[$date_tag] DONE"
  exit 0
fi
echo "[$date_tag] selective fetch unavailable, falling back to the whole granule"
rm -rf "$PART"
echo "[$date_tag] downloading to node-local $WORK"
~/micromamba/envs/insar/bin/python ~/pilot/dl_ranged.py "$url" "$WORK/$name.zip" 12 || { echo "[$date_tag] DOWNLOAD_FAILED"; exit 1; }
# re-check the window at the moment home space is committed
if [ ! -d "$dest" ] && [ ! -d "$PART" ]; then
  t0=$(date +%s)
  while [ "$(ls -d ceos/*/ 2>/dev/null | wc -l)" -ge $WINDOW ]; do
    if [ $(( $(date +%s) - t0 )) -gt 7200 ]; then
      echo "[$date_tag] window filled during download and stayed full, dropping"; exit 1
    fi
    sleep 300
  done
fi
rm -rf "$dest" "$PART"; mkdir -p "$PART"
echo "[$date_tag] selective extract HH F2-F4 (atomic via $PART)"
if unzip -q -o "$WORK/$name.zip" "VOL-*" "LED-*" "TRL-*" "IMG-HH-*-F2" "IMG-HH-*-F3" "IMG-HH-*-F4" -d "$PART/"; then
  touch "$PART/.complete"
  mv "$PART" "$dest"
  echo "[$date_tag] extracted $(ls $dest | wc -l) files, $(du -sm $dest | cut -f1) MB"
  echo "[$date_tag] DONE"
else
  echo "[$date_tag] EXTRACT_FAILED"; exit 1
fi
