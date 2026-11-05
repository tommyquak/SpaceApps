#!/bin/bash
#SBATCH --job-name=p1-stack
#SBATCH --partition=long
#SBATCH --cpus-per-task=24
#SBATCH --mem=60G
#SBATCH --tmp=200000
#SBATCH --time=20:00:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --no-requeue
#SBATCH --output=/home/t/tommy/full/logs/p1_stack-%j.log
# ALOS-1 PALSAR stripmap SBAS stack with ISCE2 stripmapStack (Zahro et al.
# 2025, Sci Rep 15:23749, used the same processor on this ground for
# 2007-2009 and got 76 percent coherent pixels).
#
#   unzip -> unpackFrame_ALOS_raw.py (FBS resampled to FBD bandwidth)
#   -> stackStripMap.py (raw focused by ISCE, interferogram workflow)
#   -> run_files executed in order, each through run.py -p NPROC
#
# FRAMES: along-track frames unzipped under one date directory are
# concatenated by the ALOS raw reader (unpackFrame globs <date>/*/IMG-HH*).
# Frame 7130 alone is the recipe validated on the Vultr VM in May; both
# frames extend the coverage south to the last dam sites (lat -2.985).
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj
export ISCE_STACK=$HOME/micromamba/envs/insar/share/isce2
export PATH=$ISCE_STACK/stripmapStack:$PATH
export PYTHONPATH=$ISCE_STACK:${PYTHONPATH:-}
export OMP_NUM_THREADS=2
NPROC=${NPROC:-12}

MAN=${MAN:-$HOME/full/palsar1/manifest_p422_l10.csv}
FRAMES=${FRAMES:-7120,7130}
# sbatch --export splits on commas, so multi-frame lists arrive as "7120:7130" or "7120+7130"
FRAMES=$(echo "$FRAMES" | tr ":+" ",,")
TAG=${TAG:-p422_f$(echo $FRAMES | tr "," "_")}
ROOT=${ROOT:-/mnt/scratch/t/tommy/palsar1}
Z=$ROOT/zips
W=$ROOT/$TAG
RAW=$W/raw; SLC=$W/SLC; STACK=$W/stack
DEM=${DEM:-$HOME/pilot/dem/elevation.dem.wgs84}   # 1 arcsec, 111-117 E, 0 to -7 S: covers every frame
LEVEL=${LEVEL:-L1.0}                              # L1.0 raw (ISCE focuses) or L1.1 SLC (p1_unpack_slc.py, --nofocus)
RLKS=${RLKS:-6}; ALKS=${ALKS:-28}                 # ~90 m ground posting (slant 9.4 m, az 3.2 m)
TBASE=${TBASE:-800}; BBASE=${BBASE:-3000}         # days, metres; L-band keeps 2007-2009 connected
FILT=${FILT:-0.5}
BBOX=${BBOX:-}                                    # optional "S N W E" crop (stackStripMap -x); halves the footprint of a 3-frame run
mkdir -p $RAW $SLC $STACK $HOME/full/logs
# The scratch quota (~470 GB measured 2026-09-06) is mostly the ALOS-2 stack, so the
# coregistered SLCs and geometry (merged/, ~25 GB per frame) can live in home instead.
if [ "${MERGED_ON_HOME:-0}" = 1 ] && [ ! -L $STACK/merged ]; then
  MH=$HOME/full/palsar1_merged_$TAG; mkdir -p $MH
  [ -d $STACK/merged ] && cp -a $STACK/merged/. $MH/ && rm -rf $STACK/merged
  ln -s $MH $STACK/merged && echo "merged/ -> $MH"
fi
echo "TAG=$TAG FRAMES=$FRAMES LEVEL=$LEVEL W=$W DEM=$DEM looks ${RLKS}x${ALKS} tbase $TBASE bbase $BBASE bbox [${BBOX}]"
[ -f $Z/DL_OK ] || { echo "P1_STACK_FAIL zips not verified (no $Z/DL_OK)"; exit 1; }
[ -f "$DEM" ] && [ -f "$DEM.xml" ] || { echo "P1_STACK_FAIL DEM missing"; exit 1; }
python -c "import shutil,sys; print('stackStripMap:', shutil.which('stackStripMap.py')); sys.exit(shutil.which('stackStripMap.py') is None)" || exit 1

# --- 1. unzip the selected frames under one directory per date -------------
python - "$MAN" "$FRAMES" "$Z" "$RAW" <<'PY' || { echo "P1_STACK_FAIL unzip"; exit 1; }
import csv, os, sys, zipfile
man, frames, Z, RAW = sys.argv[1], set(sys.argv[2].split(",")), sys.argv[3], sys.argv[4]
rows = [r for r in csv.DictReader(open(man)) if r["frame"] in frames]
by_date = {}
for r in rows:
    by_date.setdefault(r["date"].replace("-", ""), []).append(r)
modes = {}
for d, rs in sorted(by_date.items()):
    if len(rs) != len(frames):
        print("SKIP", d, "has", len(rs), "of", len(frames), "frames"); continue
    ms = {r["mode"] for r in rs}
    if len(ms) != 1:
        print("SKIP", d, "mixed modes", ms); continue
    modes[d] = ms.pop()
    dd = os.path.join(RAW, d); os.makedirs(dd, exist_ok=True)
    for r in rs:
        z = os.path.join(Z, os.path.basename(r["url"]))
        gdir = os.path.join(dd, os.path.splitext(os.path.basename(z))[0])
        if (os.path.isdir(gdir) and any("IMG-HH" in x for x in os.listdir(gdir))) or any(x.startswith("IMG-HH") for x in os.listdir(dd)):
            continue
        with zipfile.ZipFile(z) as zf:
            zf.extractall(dd)
        # L1.0 zips carry one top-level directory named like the zip; L1.1 zips
        # extract IMG-HH/LED files straight into the date directory
        if not os.path.isdir(gdir):
            cands = [x for x in os.listdir(dd) if x.startswith(r["scene"][:12])]
            flat = any(x.startswith("IMG-HH") for x in os.listdir(dd))
            if not cands and not flat:
                print("no granule dir for", r["scene"]); sys.exit(1)
    print("raw", d, modes[d], len(rs), "frame(s)")
with open(os.path.join(RAW, "modes.txt"), "w") as f:
    for d in sorted(modes): f.write(f"{d} {modes[d]}\n")
print("dates:", len(modes))
PY
ndates=$(wc -l < $RAW/modes.txt); echo "dates with all frames: $ndates"
[ "$ndates" -ge 8 ] || { echo "P1_STACK_FAIL too few dates ($ndates)"; exit 1; }

# --- 2. unpack: raw ALOS -> ISCE raw + metadata, FBS resampled to FBD --------
unpack() {  # date mode
  local d=$1 m=$2
  if [ "$LEVEL" = L1.1 ]; then
    [ -s "$SLC/$d/$d.slc.xml" ] && ls $SLC/$d/data.* >/dev/null 2>&1 && { echo "unpacked $d (cached)"; return 0; }
    python $HOME/full/palsar1/p1_unpack_slc.py -i $RAW/$d -o $SLC/$d > $SLC/unpack_$d.log 2>&1 \
      && ls $SLC/$d/$d.slc.xml >/dev/null 2>&1 && ls $SLC/$d/data.* >/dev/null 2>&1 && echo "unpacked $d $m (SLC)" || { echo "UNPACK FAIL $d (see $SLC/unpack_$d.log)"; return 1; }
    return 0
  fi
  [ -s "$SLC/$d/$d.raw.xml" ] && { echo "unpacked $d (cached)"; return 0; }
  local extra=""; [ "$m" = FBS ] && extra="-f fbs2fbd"
  [ "$NFRAMES" -gt 1 ] && extra="$extra -m"      # concatenate the along-track frames
  unpackFrame_ALOS_raw.py -i $RAW/$d -o $SLC/$d $extra > $SLC/unpack_$d.log 2>&1 \
    && ls $SLC/$d/*.raw.xml >/dev/null 2>&1 && echo "unpacked $d $m" || { echo "UNPACK FAIL $d (see $SLC/unpack_$d.log)"; return 1; }
}
NFRAMES=$(echo $FRAMES | tr "," "\n" | grep -c .)
export -f unpack; export RAW SLC NFRAMES LEVEL
cat $RAW/modes.txt | xargs -P 6 -L 1 bash -c 'unpack "$0" "$1"'
nslc=$(ls -d $SLC/2* 2>/dev/null | wc -l); echo "unpacked dates: $nslc / $ndates"
[ "$nslc" -eq "$ndates" ] || { echo "P1_STACK_FAIL unpack incomplete"; exit 1; }
rm -rf $RAW/2*/ && echo "unzipped CEOS removed (zips kept)"

# --- 3. stack definition -------------------------------------------------------
# reference = the date closest to the middle of the span (baseline spread)
REF=$(python - $RAW/modes.txt <<'PY'
import sys, datetime as dt
ds = sorted(l.split()[0] for l in open(sys.argv[1]))
t = [dt.datetime.strptime(d, "%Y%m%d") for d in ds]
mid = t[0] + (t[-1] - t[0]) / 2
print(min(ds, key=lambda d: abs(dt.datetime.strptime(d, "%Y%m%d") - mid)))
PY
)
echo "reference date: $REF"
cd $W
if [ ! -d $STACK/run_files ]; then
  XB=""; [ -n "$BBOX" ] && XB="-x $BBOX"
  NF=""; [ "$LEVEL" = L1.1 ] && NF="--nofocus"
  stackStripMap.py -s $SLC -d $DEM -w $STACK -m $REF -W interferogram \
    -t $TBASE -b $BBASE -r $RLKS -a $ALKS -u snaphu -f $FILT $XB $NF > $W/stackStripMap.log 2>&1 \
    || { echo "P1_STACK_FAIL stackStripMap.py"; tail -20 $W/stackStripMap.log; exit 1; }
fi
ls $STACK/run_files/ | grep -vE '\.(job|log|DONE)$' | tr '\n' ' '; echo
npairs=$(cat $STACK/run_files/run_*_igram 2>/dev/null | grep -c .)
echo "pairs in network: $npairs"

# --- 4. execute the run files in order -----------------------------------------
# run.py fans the lines of one run file over NPROC workers (the stripmapStack
# equivalent of the topsStack cmd & ... wait batching; GNU parallel is absent).
for f in $(ls $STACK/run_files/run_* | grep -vE '\.(job|log|DONE)$' | sort -V); do
  b=$(basename $f)
  [ -f $STACK/run_files/$b.DONE ] && { echo "skip $b (done)"; continue; }
  echo "=== $(date '+%H:%M:%S') $b ($(grep -c . $f) lines) ==="
  touch $STACK/run_files/.probe 2>/dev/null || { echo "P1_STACK_FAIL scratch quota exhausted before $b"; exit 1; }
  run.py -i $f -p $NPROC > $STACK/run_files/$b.log 2>&1
  rc=$?
  touch $STACK/run_files/.probe 2>/dev/null || { echo "P1_STACK_FAIL scratch quota exhausted during $b (products may be truncated)"; exit 1; }
  # curl 404s are the water-mask stitcher probing SWBD tiles that do not
  # exist (open water, benign); only Python failures and a non-zero rc count.
  if [ $rc -ne 0 ] || grep -vE "^curl: " $STACK/run_files/$b.log | grep -qE "Traceback|^[A-Za-z_.]+Error:|error while|Killed"; then
    echo "P1_STACK_FAIL $b (rc=$rc); first errors:"; grep -vE "^curl: " $STACK/run_files/$b.log | grep -nE "Traceback|Error" | head -5; exit 1
  fi
  touch $STACK/run_files/$b.DONE
  # Scratch quota is ~512 GB shared with the ALOS-2 stack: a single PALSAR
  # frame run peaks at ~110 GB, so intermediates go as soon as they are consumed.
  # fineResamp reads offsets/ and the focused SLC/<date>/*.raw.slc; after it,
  # merged/SLC holds every coregistered SLC and only the reference shelve is needed.
  case $b in
    *fineResamp)
      n_ok=$(for d in $STACK/merged/SLC/*/; do stat -c %s $d/$(basename $d).slc 2>/dev/null; done | sort | uniq -c | sort -rn | head -1 | awk '{print $1}')
      if [ "${n_ok:-0}" -ge "$ndates" ]; then
        # the igram configs read the REFERENCE from SLC/$REF/$REF.raw.slc, not from
        # merged/SLC, so that one focused SLC must survive (2026-09-06: deleting it
        # cost the 13 reference pairs of the first frame-7130 run)
        rm -rf $STACK/offsets $STACK/coregSLC
        for d in $SLC/2*/; do [ "$(basename $d)" = "$REF" ] && continue; rm -f $d/*.raw $d/*.raw.slc $d/*.raw.aux; [ "$LEVEL" = L1.1 ] && rm -f $d/$(basename $d).slc; done
        echo "intermediates removed after $b (offsets, coregSLC, focused SLCs except the reference $REF)"
      else
        echo "P1_STACK_FAIL merged SLCs incomplete or unequal in size after $b ($n_ok of $ndates)"; exit 1
      fi;;
  esac
done

# --- 5. product check ------------------------------------------------------------
nunw=$(ls $STACK/Igrams/*/filt_*.unw 2>/dev/null | grep -v conncomp | wc -l)
ncor=$(ls $STACK/Igrams/*/filt_*.cor 2>/dev/null | wc -l)
echo "interferograms unwrapped: $nunw  coherence files: $ncor  geometry: $(ls $STACK/geom_reference/ 2>/dev/null | tr '\n' ' ')"
python - $STACK <<'PY'
import glob, os, sys, numpy as np, xml.etree.ElementTree as ET
st = sys.argv[1]
def rd(fn):
    x = ET.parse(fn + ".xml").getroot()
    p = {e.get("name"): e.find("value").text for e in x.findall("property")}
    w, l, nb = int(p["width"]), int(p["length"]), int(p.get("number_bands", 1))
    a = np.fromfile(fn, dtype=np.float32)
    return a.reshape(l, nb, w) if nb > 1 else a.reshape(l, 1, w)
cors = sorted(glob.glob(f"{st}/Igrams/*/filt_*.cor"))
vals = []
for c in cors:
    a = rd(c); coh = a[:, -1, :]
    vals.append((os.path.basename(os.path.dirname(c)), float(np.nanmean(coh[coh > 0])), float((coh > 0.3).mean())))
vals.sort(key=lambda v: v[1])
print("pair, mean coherence, frac>0.3 (worst 5 / best 5):")
for v in vals[:5] + vals[-5:]: print("  %s %.3f %.2f" % v)
print("network mean coherence %.3f over %d pairs" % (np.mean([v[1] for v in vals]), len(vals)))
PY
[ "$nunw" -ge 20 ] && touch $STACK/STACK_OK && echo "P1_STACK_OK $TAG pairs=$nunw" || { echo "P1_STACK_FAIL only $nunw unwrapped pairs"; exit 1; }
