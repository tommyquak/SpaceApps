#!/bin/bash
#SBATCH --job-name=nisarpull
#SBATCH --output=/home/t/tommy/full/logs/nisarpull_%j.out
#SBATCH --error=/home/t/tommy/full/logs/nisarpull_%j.out
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
# NISAR forward processing on provisional GUNW, 12-day cadence.
#   1. poll CMR for provisional GUNW granules on the two tracks over the AOI
#   2. download what is missing (Earthdata credentials in ~/.netrc)
#   3. rerun the coherence and local-contrast stages (2c, 2d) into a dated dir
#   4. resubmit itself 12 days out, unless the kill switch exists
# Everything written is PROVISIONAL: JPL will reprocess these products.
# Kill switch: touch ~/full/NISAR_PULL_STOP
set -u
export PROJ_LIB=$HOME/micromamba/envs/insar/share/proj TMPDIR=$HOME/full/tmp PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1
cd $HOME/full || exit 1
[ -f NISAR_PULL_STOP ] && { echo "NISAR_PULL_STOP present; not resubmitting"; exit 0; }
D=/mnt/scratch/t/tommy/nisar/gunw; mkdir -p "$D"
STAMP=$(date +%Y%m%d)
echo "=== NISAR pull $STAMP ==="
# 1. poll CMR (provisional collection), both tracks, point inside the AOI
python3 - <<'PY' > gunw_urls_new.txt || { echo "PULL_ABORT cmr"; exit 1; }
import json, ssl, urllib.request, urllib.parse
ctx = ssl.create_default_context()
url = "https://cmr.earthdata.nasa.gov/search/granules.umm_json?" + urllib.parse.urlencode({"collection_concept_id": "C2854335566-ASF", "point": "113.95,-2.85", "page_size": 500})
req = urllib.request.Request(url, headers={"User-Agent": "peatguard-nisar-pull/1.0"})
d = json.load(urllib.request.urlopen(req, timeout=120, context=ctx))
for it in d.get("items", []):
    u = it["umm"]
    h5 = [r["URL"] for r in u.get("RelatedUrls", []) if r.get("Type", "").startswith("GET DATA") and r["URL"].endswith(".h5")]
    if h5: print(h5[0])
PY
echo "CMR granules: $(wc -l < gunw_urls_new.txt)"
# 2. download missing
n_new=0
while read -r u; do
  [ -n "$u" ] || continue; f=$(basename "$u")
  [ -s "$D/$f" ] && continue
  echo "GET $f"; wget -q --netrc --auth-no-challenge --tries=3 --timeout=120 -O "$D/$f.part" "$u" && mv "$D/$f.part" "$D/$f" && n_new=$((n_new+1)) || { echo "DOWNLOAD_FAILED $f"; rm -f "$D/$f.part"; }
done < gunw_urls_new.txt
echo "new granules: $n_new; total on disk: $(ls $D/*.h5 2>/dev/null | wc -l)"
# 3. isolate every attempt, including repeated runs on the same day.
mkdir -p "$HOME/full/nisar_out" || exit 1
OUTD=$(mktemp -d "$HOME/full/nisar_out/${STAMP}_${SLURM_JOB_ID:-manual}_XXXXXX") || exit 1
export GUNW_DIR="$D" NISAR_OUT="$OUTD" HP_KM="${HP_KM:-10}"
INSAR_PYTHON=$HOME/micromamba/envs/insar/bin/python

run_analysis() {
  local stage rc ad suffix
  for stage in stage1 stage2c; do
    "$INSAR_PYTHON" "nisar_${stage}.py" > "$OUTD/${stage}.txt" 2>&1
    rc=$?
    echo "$stage rc=$rc"
    [ "$rc" -eq 0 ] || return "$rc"
  done
  # A fresh directory and direct stage output prevent stale global products
  # being relabeled as results from this attempt. HP_KM is normalized as in Python.
  local hp_label
  hp_label=$("$INSAR_PYTHON" -c 'import os; print(format(float(os.environ["HP_KM"]), "g"))') || return 1
  for ad in D A; do
    for suffix in "disp60_hp${hp_label}km_mm.tif" chain_coh.tif sources.json; do
      [ -s "$OUTD/nisar_${ad}_${suffix}" ] || {
        echo "PULL_ABORT missing current $ad $suffix"
        return 1
      }
    done
  done
  "$INSAR_PYTHON" nisar_stage2d.py > "$OUTD/stage2d.txt" 2>&1
  rc=$?
  echo "stage2d rc=$rc"
  [ "$rc" -eq 0 ] || return "$rc"
  echo "PROVISIONAL: NISAR L2 GUNW provisional products (C2854335566-ASF), pulled $STAMP; subject to JPL reprocessing. Local-contrast diagnostics only; no frame-scale displacement, rate, or validated restoration effect." > "$OUTD/PROVISIONAL.txt" || return 1
  grep -E "site minus control|two-track agreement|r=" "$OUTD/stage2d.txt" | head -6
  return 0
}

analysis_status=0
run_analysis || analysis_status=$?
# 4. resubmit 12 days out
[ -f NISAR_PULL_STOP ] || { J=$(sbatch --parsable --begin=now+12days nisar_pull.sbatch) && echo "resubmitted as $J for $(date -d '+12 days' +%F)"; }
[ "$analysis_status" -eq 0 ] || { echo "PULL_FAILED analysis rc=$analysis_status output=$OUTD"; exit "$analysis_status"; }
echo "PULL_DONE $STAMP"
