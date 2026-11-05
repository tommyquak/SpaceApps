#!/bin/bash
#SBATCH --job-name=p1-dl
#SBATCH --partition=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:30:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --output=/home/t/tommy/full/logs/p1_dl-%j.log
# ALOS-1 PALSAR L1.0 raw download for the 2007-2009 historical baseline
# (path 422, frames 7120 + 7130, 14 dates each; ASF datapool, ~12 GB).
#
# L1.0 raw, not L1.1 SLC: ISCE2 focuses the raw data itself and sets the
# Doppler centroid; the L1.1 SLC route left it unset and broke geo2rdr
# (insar_salvage/scripts/alos1_process.sh, 2026-05). Earthdata login comes
# from ~/.netrc (same as nisar_pull.sbatch). Every zip is size-checked
# against the CMR manifest and test-read with zipfile before DL_OK is written.
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$PATH
MAN=${MAN:-$HOME/full/palsar1/manifest_p422_l10.csv}
FRAMES=${FRAMES:-7120,7130}
# sbatch --export splits on commas, so multi-frame lists arrive as "7120:7130" or "7120+7130"
FRAMES=$(echo "$FRAMES" | tr ":+" ",,")
ROOT=${ROOT:-/mnt/scratch/t/tommy/palsar1}
Z=$ROOT/zips
mkdir -p $Z
echo "MAN=$MAN FRAMES=$FRAMES Z=$Z"
[ -s "$MAN" ] || { echo "P1_DL_FAIL manifest missing"; exit 1; }
[ -s "$HOME/.netrc" ] || { echo "P1_DL_FAIL no ~/.netrc"; exit 1; }

# rows: date,path,frame,mode,scene,size_mb,url,...
python3 - "$MAN" "$FRAMES" > $Z/want.txt <<'PY'
import csv, sys
man, frames = sys.argv[1], set(sys.argv[2].split(","))
for r in csv.DictReader(open(man)):
    if r["frame"] in frames:
        print(r["url"], r["size_mb"], r["scene"])
PY
n=$(wc -l < $Z/want.txt); echo "granules wanted: $n"
[ "$n" -gt 0 ] || { echo "P1_DL_FAIL nothing selected"; exit 1; }

get() {  # url size_mb scene
  local u=$1 mb=$2 sc=$3 f=$Z/$(basename $1)
  if [ -s "$f" ]; then echo "have $f"; return 0; fi
  for t in 1 2 3; do
    wget -q --netrc --auth-no-challenge --tries=2 --timeout=180 -O "$f.part" "$u" && mv "$f.part" "$f" && { echo "GET ok $sc"; return 0; }
    echo "GET retry $t $sc"; sleep 20
  done
  rm -f "$f.part"; echo "GET FAIL $sc"; return 1
}
export -f get; export Z
cat $Z/want.txt | xargs -P 4 -L 1 bash -c 'get "$0" "$1" "$2"'

# verify: size within 2 percent of the CMR size and the zip is readable
python3 - "$Z" <<'PY' && touch $Z/DL_OK && echo "P1_DL_OK" || { echo "P1_DL_FAIL verification"; exit 1; }
import os, sys, zipfile
Z = sys.argv[1]; bad = 0; n = 0
for line in open(f"{Z}/want.txt"):
    url, mb, scene = line.split()
    f = os.path.join(Z, os.path.basename(url)); n += 1
    if not os.path.exists(f):
        print("MISSING", scene); bad += 1; continue
    got = os.path.getsize(f) / 1048576.0
    if abs(got - float(mb)) > 0.02 * float(mb) + 1:
        print("SIZE", scene, round(got, 1), "vs", mb); bad += 1; continue
    try:
        with zipfile.ZipFile(f) as zf:
            if zf.testzip() is not None:
                print("CORRUPT", scene); bad += 1; continue
            names = zf.namelist()
            if not any("IMG-HH" in x for x in names) or not any("LED-" in x for x in names):
                print("NO IMG/LED", scene); bad += 1; continue
    except zipfile.BadZipFile:
        print("BADZIP", scene); bad += 1; continue
print(f"verified {n - bad}/{n} granules")
sys.exit(1 if bad else 0)
PY
