#!/bin/bash
#SBATCH --job-name=full-rehearsal
#SBATCH --dependency=singleton
#SBATCH --partition=long
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=08:00:00
#SBATCH --output=/home/t/tommy/full/logs/rehearsal.log

# ============================== SUPERSEDED ==================================
# This rehearsal workflow duplicates the MintPy finale logic and has NOT been
# kept in step with mintpy2.sbatch. As of 2026-09-04 it still contains defects
# that mintpy2.sbatch has fixed:
#   - Zheng apply path looks for timeseries_ramp_demErr.h5 / _demErr_ramp.h5,
#     which this configuration never produces (topographic residual is OFF), so
#     the correction is estimated and then silently not applied (defect 12);
#   - it subtracts the un-deramped bias from an already-deramped series, and
#     deramp(D)-B != deramp(D-B) unless B has zero planar component;
#   - no current-run provenance stamp, so a stale corrected product looks fresh;
#   - it exits zero when the correction was attempted but not produced.
#
# Duplicated correction logic in two scripts is exactly how those defects
# survived, so this file is NOT being re-fixed in parallel. Use mintpy2.sbatch.
# If you genuinely need the rehearsal path, port it to call mintpy2.sbatch
# rather than reimplementing it, then delete this guard.
if [ "${ALLOW_SUPERSEDED_REHEARSAL:-0}" != "1" ]; then
  echo "REHEARSAL_SUPERSEDED: this script contains known-unfixed defects (see header)."
  echo "  Use mintpy2.sbatch. To override deliberately: ALLOW_SUPERSEDED_REHEARSAL=1"
  exit 2
fi
# ============================================================================
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj
export OMP_NUM_THREADS=16
done_count=$(ls -d ~/full/pairs/*/DONE ~/pilot/pairs/*/DONE 2>/dev/null | awk -F/ '{print $(NF-1)}' | sort -u | wc -l)
echo "terminal DONE pairs: $done_count"
[ "$done_count" -lt 200 ] && { echo "too few pairs, aborting for human review"; exit 1; }
W=~/full/mintpy_rehearsal
rm -rf $W; mkdir -p $W
python ~/full/build_stack2.py $W || { echo "build_stack2 FAILED"; echo MINTPY_FAILED; exit 1; }
cd $W
cat > fullApp.cfg <<CFGEOF
mintpy.load.processor = isce
mintpy.reference.method = maxCoherence
mintpy.reference.coherenceFile = avgSpatialCoh.h5
mintpy.network.coherenceBased = no
mintpy.unwrapError.method = phase_closure
mintpy.networkInversion.weightFunc = var
mintpy.networkInversion.maskDataset = coherence
mintpy.networkInversion.maskThreshold = 0.25
mintpy.troposphericDelay.method = no
mintpy.deramp = linear
mintpy.topographicResidual = yes
mintpy.velocity.startDate = auto
CFGEOF
# Two-phase run (rehearsal finding): the phase_closure correction may
# legitimately skip (no integer ambiguities found) without writing the
# unwrapPhase_phaseClosure dataset, and the inversion then crashes asking
# for it. Run through the correction, check the dataset, and fall back to
# plain unwrapPhase if the correction produced nothing.
smallbaselineApp.py fullApp.cfg --dir $W --start modify_network --end correct_unwrap_error || { echo "smallbaselineApp(correct) FAILED"; echo MINTPY_FAILED; exit 1; }
HASPC=$(python - <<'PYX'
import h5py
f = h5py.File("inputs/ifgramStack.h5", "r")
print(1 if "unwrapPhase_phaseClosure" in f else 0)
PYX
)
if [ "$HASPC" != "1" ]; then
  echo "phase_closure produced no corrected dataset (nothing to correct); inverting plain unwrapPhase"
  sed -i 's/mintpy.unwrapError.method = phase_closure/mintpy.unwrapError.method = no/' fullApp.cfg
fi
smallbaselineApp.py fullApp.cfg --dir $W --start invert_network --end velocity || { echo "smallbaselineApp(invert) FAILED"; echo MINTPY_FAILED; exit 1; }
echo "--- closure bias correction (Zheng et al. 2022, MintPy implementation) ---"
# Gate: the estimator raises mid-run if any sequential closure loop is missing
# a pair; --nl 20 is the science assumption and must never be silently lowered.
python - <<'PYEOF2'
import h5py
f = h5py.File("inputs/ifgramStack.h5", "r")
d12 = [(a.decode(), b.decode()) for a, b in f["date"][:]]
dates = sorted(set(x for p in d12 for x in p))
idx = {d: i for i, d in enumerate(dates)}
pairs = set((idx[a], idx[b]) for a, b in d12)
maxL, first_missing = 0, []
for conn in range(1, 21):
    miss = [(dates[i], dates[i + conn]) for i in range(len(dates) - conn)
            if (i, i + conn) not in pairs]
    if miss:
        first_missing = miss[:20]
        break
    maxL = conn
print(f"sequential network complete through connection level {maxL}")
if maxL >= 20:
    print("CLOSURE_NETWORK_COMPLETE")
else:
    print(f"CLOSURE_NETWORK_INCOMPLETE at level {maxL + 1}, missing: {first_missing}")
PYEOF2
python - <<'PYEOF2' || true > closure_gate.txt
import h5py
f = h5py.File("inputs/ifgramStack.h5", "r")
d12 = [(a.decode(), b.decode()) for a, b in f["date"][:]]
dates = sorted(set(x for p in d12 for x in p))
idx = {d: i for i, d in enumerate(dates)}
pairs = set((idx[a], idx[b]) for a, b in d12)
ok = all((i, i + c) in pairs for c in range(1, 21) for i in range(len(dates) - c))
print("COMPLETE" if ok else "INCOMPLETE")
PYEOF2
if grep -q COMPLETE closure_gate.txt && ! grep -q INCOMPLETE closure_gate.txt; then
  set +e
  closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 20 --bw 20 -a mask >> closure_bias.log 2>&1; rc1=$?
  closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 20 --bw 20 -a quick_estimate -c local --num-worker 16 >> closure_bias.log 2>&1; rc2=$?
  closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 20 --bw 20 -a estimate -c local --num-worker 16 >> closure_bias.log 2>&1; rc3=$?
  set -e
  rc_cb=$((rc1 | rc2 | rc3))
  [ $rc_cb -eq 0 ] && echo CLOSURE_BIAS_STEP_OK || echo "CLOSURE_BIAS_STEP_INCOMPLETE rc=$rc1/$rc2/$rc3 (see closure_bias.log)"
  # APPLY the estimated bias (audit r2 finding 1: estimation alone never touches
  # the velocity). Keyed on the estimate step alone, not the aggregate rc.
  if [ $rc3 -eq 0 ] && [ -f $W/timeseriesBias.h5 ]; then
    TS=$(ls $W/timeseries_ramp_demErr.h5 $W/timeseries_demErr_ramp.h5 2>/dev/null | head -1)
    if [ -n "$TS" ] && diff.py $TS $W/timeseriesBias.h5 -o $W/timeseries_closureCor.h5 \
        && timeseries2velocity.py $W/timeseries_closureCor.h5 -t $W/smallbaselineApp.cfg -o $W/velocityClosureCor.h5; then
      echo CLOSURE_BIAS_APPLIED
    else
      echo CLOSURE_BIAS_APPLY_FAILED
    fi
  else
    echo CLOSURE_BIAS_NOT_APPLIED
  fi
else
  echo "CLOSURE_SKIPPED_INCOMPLETE_NETWORK: velocity stands without closure correction; run closure_phase_bias manually on a reduced network (see log above for missing pairs)"
fi
echo "--- outputs:"
ls -la $W/*.h5 2>/dev/null | head -20
for f in velocity velocityClosureCor temporalCoherence avgSpatialCoh; do
  [ -f $W/$f.h5 ] && save_gdal.py $W/$f.h5 -o $W/$f.tif 2>/dev/null && echo "exported $f.tif"
done
if [ -f $W/velocityClosureCor.h5 ]; then echo MINTPY_OK; FINAL_RC=0
elif [ -f $W/velocity.h5 ]; then echo MINTPY_OK_UNCORRECTED; FINAL_RC=0
else echo MINTPY_FAILED; FINAL_RC=1; fi
echo "--- rehearsal diagnostics: is the mini-campaign velocity physical? ---"
python - <<'PYEOF2'
import numpy as np
from osgeo import gdal
gdal.UseExceptions()
import os
W = os.path.expanduser("~/full/mintpy_rehearsal")
for name in ("velocityClosureCor", "velocity"):
    p = f"{W}/{name}.tif"
    if not os.path.exists(p): continue
    v = gdal.Open(p).GetRasterBand(1).ReadAsArray().astype("float64") * 1000
    t = gdal.Open(f"{W}/temporalCoherence.tif").GetRasterBand(1).ReadAsArray().astype("float64")
    ok = np.isfinite(v) & (v != 0)
    for thr in (0.25, 0.7, 0.9):
        m = ok & (t >= thr)
        if m.sum() < 1000: continue
        print(f"{name} tcoh>={thr}: median {np.median(v[m]):7.1f} mm/yr  p5 {np.percentile(v[m],5):7.1f}  p95 {np.percentile(v[m],95):7.1f}  n={int(m.sum()):,}")
print("(pilot failure benchmark: median -91 mm/yr LOS at tcoh>=0.25; physical envelope 0 to -60 mm/yr)")
PYEOF2

exit $FINAL_RC
