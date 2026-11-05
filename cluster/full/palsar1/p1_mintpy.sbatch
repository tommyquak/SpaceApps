#!/bin/bash
#SBATCH --job-name=p1-mintpy
#SBATCH --partition=long
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=12:00:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --output=/home/t/tommy/full/logs/p1_mintpy-%j.log
# MintPy SBAS on the ALOS-1 stripmapStack stack, then the validation script.
# NaNs in the ISCE products are zeroed before loading: the May attempt on the
# Vultr VM stalled on MintPy NaN handling and worked around it by excluding
# the seven FBS dates, which halves the network; zeroing keeps all 14 dates.
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj
export OMP_NUM_THREADS=16
TAG=${TAG:-p422_f7130}
ROOT=${ROOT:-/mnt/scratch/t/tommy/palsar1}
STACK=$ROOT/$TAG/stack
# NET=short keeps pairs <= TBMAX days / PBMAX m (Zahro et al. 2025: 350 d, 1000 m);
# the 2-year pairs of the full network carry the unwrapping errors that pulled
# temporal coherence down to a 0.49 median on the first frame-7130 run.
NET=${NET:-short}; TBMAX=${TBMAX:-400}; PBMAX=${PBMAX:-2000}
W=$HOME/full/palsar1_mintpy_$TAG; [ "$NET" = full ] && W=${W}_full
CFG_T=${CFG_T:-$HOME/full/palsar1/smallbaselineApp_palsar1.cfg}
mkdir -p $W; cd $W
echo "TAG=$TAG STACK=$STACK W=$W"
[ -f $STACK/STACK_OK ] || { echo "P1_MINTPY_FAIL stack not certified"; exit 1; }

# --- resolve the product names the stack actually wrote ------------------------
META=$(ls $STACK/referenceShelve/data.dat 2>/dev/null || ls $STACK/Igrams/*/referenceShelve/data.dat 2>/dev/null | head -1)
if ls $STACK/Igrams/*/filt_*_snaphu.unw >/dev/null 2>&1; then UNW="$STACK/Igrams/*/filt_*_snaphu.unw"; CONN="$STACK/Igrams/*/filt_*_snaphu.unw.conncomp";
else UNW="$STACK/Igrams/*/filt_*.unw"; CONN="$STACK/Igrams/*/filt_*.unw.conncomp"; fi
[ -n "$META" ] || { echo "P1_MINTPY_FAIL no referenceShelve/data.dat"; exit 1; }
sed "s#__STACK__#$STACK#g; s#__META__#$META#g; s#__UNW__#$UNW#g; s#__CONN__#$CONN#g" $CFG_T > $W/smallbaselineApp.cfg
[ -f $STACK/geom_reference/shadowMask.rdr ] || sed -i '/shadowMaskFile/d' $W/smallbaselineApp.cfg
if [ "$NET" = short ]; then printf "mintpy.network.tempBaseMax = %s\nmintpy.network.perpBaseMax = %s\n" $TBMAX $PBMAX >> $W/smallbaselineApp.cfg; fi
echo "NET=$NET (tempBaseMax $TBMAX d, perpBaseMax $PBMAX m) W=$W"
echo "meta=$META unw=$UNW"

# --- zero NaNs in place (unw, cor, geometry) ------------------------------------
python - $STACK <<'PY'
import glob, sys, numpy as np, xml.etree.ElementTree as ET
st = sys.argv[1]; fixed = 0
files = glob.glob(f"{st}/Igrams/*/filt_*.unw") + glob.glob(f"{st}/Igrams/*/filt_*.cor") + glob.glob(f"{st}/geom_reference/*.rdr")
for fn in files:
    try:
        x = ET.parse(fn + ".xml").getroot()
    except Exception:
        continue
    p = {e.get("name"): e.find("value").text for e in x.findall("property")}
    if p.get("data_type", "FLOAT").upper() not in ("FLOAT", "DOUBLE"):
        continue
    dt = np.float32 if p.get("data_type", "FLOAT").upper() == "FLOAT" else np.float64
    a = np.memmap(fn, dtype=dt, mode="r+")
    bad = ~np.isfinite(a)
    if bad.any():
        a[bad] = 0; a.flush(); fixed += 1
print("files with NaN/Inf zeroed:", fixed, "of", len(files))
PY

# --- reference pixel: interior, coherent, and unwrapped in nearly every pair ------------
# MintPy's automatic choice (max average spatial coherence inside maskConnComp) landed on
# row 67 of frame 7140, an edge pixel whose unwrapped phase is nodata in most pairs; the
# whole inversion then collapsed to one valid pixel. Choose it ourselves from the loaded stack.
smallbaselineApp.py $W/smallbaselineApp.cfg --dostep load_data > $W/load_data.log 2>&1 || { echo "P1_MINTPY_FAIL load_data"; tail -5 $W/load_data.log; exit 1; }
python $HOME/full/palsar1/p1_refpixel.py $W || { echo "P1_MINTPY_FAIL reference pixel"; exit 1; }
grep "mintpy.reference.yx" $W/smallbaselineApp.cfg

# --- run ----------------------------------------------------------------------------
smallbaselineApp.py $W/smallbaselineApp.cfg > $W/smallbaselineApp.log 2>&1
rc=$?
# Frame 7140 (rivers, Palangkaraya): no connected component is common to all 91
# pairs, so MintPy's phase-closure correction cannot place the reference inside
# one and aborts. Fall back to no unwrap-error correction and say so.
if [ $rc -ne 0 ] && grep -qE "reference point is NOT included in the connectComponent|unwrapPhase_phaseClosure.* not found" $W/smallbaselineApp.log; then
  echo "phase_closure correction impossible (no common connected component); rerunning with mintpy.unwrapError.method = no"
  sed -i 's/^mintpy.unwrapError.method .*/mintpy.unwrapError.method = no/' $W/smallbaselineApp.cfg
  echo "UNWRAP_ERROR_CORRECTION=none" > $W/UNWRAP_FALLBACK
  smallbaselineApp.py $W/smallbaselineApp.cfg > $W/smallbaselineApp.log 2>&1
  rc=$?
fi
tail -5 $W/smallbaselineApp.log
if [ $rc -ne 0 ] || ! grep -q "Normal end of smallbaselineApp" $W/smallbaselineApp.log; then
  echo "P1_MINTPY_FAIL rc=$rc"; grep -nE "Traceback|Error" $W/smallbaselineApp.log | tail -5; exit 1
fi
for f in velocity.h5 temporalCoherence.h5 maskTempCoh.h5 geo/geo_velocity.h5; do [ -f $W/$f ] || { echo "P1_MINTPY_FAIL missing $f"; exit 1; }; done
save_gdal.py $W/geo/geo_velocity.h5 -d velocity -o $W/geo_velocity.tif > /dev/null 2>&1
save_gdal.py $W/geo/geo_temporalCoherence.h5 -o $W/geo_temporalCoherence.tif > /dev/null 2>&1
save_gdal.py $W/geo/geo_velocity.h5 -d velocityStd -o $W/geo_velocityStd.tif > /dev/null 2>&1
save_gdal.py $W/geo/geo_geometryRadar.h5 -d incidenceAngle -o $W/geo_incidenceAngle.tif > /dev/null 2>&1
ls -la $W/*.tif | awk '{print $5, $NF}'
echo "P1_MINTPY_OK $TAG"

# --- validation against the literature and the stable-ground null --------------------
python $HOME/full/palsar1/p1_validate.py $W $TAG | tee $W/validate.txt
