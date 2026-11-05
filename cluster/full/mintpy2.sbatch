#!/bin/bash
#SBATCH --job-name=full-mintpy
#SBATCH --dependency=singleton
#SBATCH --partition=long
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --time=48:00:00
#SBATCH --output=/home/t/tommy/full/logs/mintpy.log
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$PATH
PROVENANCE=${PROVENANCE:-$HOME/full/workflow_provenance.py}
[ -s "$PROVENANCE" ] || { echo "missing workflow_provenance.py; deploy the helper with this script"; exit 1; }
ALLOW_UNREPAIRED=${ALLOW_UNREPAIRED:-0}
APPLY_ZHENG=${APPLY_ZHENG:-1}
ZHENG_BW=${ZHENG_BW:-20}
ZHENG_NL=${ZHENG_NL:-20}
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj
export OMP_NUM_THREADS=32
# PAIR_ROOT selects the MAIN pass to invert. It also keys the MintPy work dir,
# so two passes can never share a stack. Before this, done_count counted the OLD
# pair trees and the reuse guard below found the previous run's 23.8 GB
# ifgramStack.h5, so a Tay-config run would have silently inverted the stale
# 5x28 stack and still printed MINTPY_OK.
PAIR_ROOT=${PAIR_ROOT:-$HOME/full/pairs_tay}
TAG=$(basename "$PAIR_ROOT" | sed 's/^pairs_*//')
[ -n "$TAG" ] || TAG=main
export PAIR_ROOT
echo "PAIR_ROOT=$PAIR_ROOT  TAG=$TAG"
done_count=$(ls -d $PAIR_ROOT/*/DONE 2>/dev/null | awk -F/ '{print $(NF-1)}' | sort -u | wc -l)
echo "terminal DONE pairs: $done_count"
[ "$done_count" -lt 1500 ] && { echo "too few pairs, aborting for human review"; exit 1; }
W=${W:-$HOME/full/mintpy_$TAG}
export W
# Rebuilds preserve the previous working directory. Reject broad paths before
# any rename or explicitly requested product reset.
case "$W" in
  "$HOME/full/"mintpy_*) ;;
  *) echo "W must be a named $HOME/full/mintpy_* directory"; exit 1 ;;
esac
[ "$(dirname "$W")" = "$HOME/full" ] || { echo "W must be directly under $HOME/full"; exit 1; }
# reuse an already-built stack: build_stack2 warps 4893 rasters and takes ~30 min,
# and nothing about the pairs has changed. REBUILD=1 forces a fresh build.
STACK_OK=0
if [ -s "$W/inputs/ifgramStack.h5" ]; then
  STACK_OK=$($HOME/micromamba/envs/insar/bin/python - <<'PYCHK'
import h5py, glob, os, re, sys
w = os.environ["W"]; root = os.environ["PAIR_ROOT"]
try:
    f = h5py.File(w + "/inputs/ifgramStack.h5", "r")
    al, rl = int(f.attrs["ALOOKS"]), int(f.attrs["RLOOKS"])
    n = f["unwrapPhase"].shape[0]
    stack_pairs = {tuple(x.decode() if isinstance(x, bytes) else str(x) for x in pair)
                   for pair in f["date"][:]}
except Exception:
    print(0); sys.exit()
looks = set()
done_pairs = set()
for d in glob.glob(root + "/*/"):
    if not os.path.exists(d + "DONE"):
        continue
    parts = os.path.basename(d.rstrip("/")).split("_")
    if len(parts) != 2:
        print(0); sys.exit()
    done_pairs.add(tuple(x.replace("-", "") for x in parts))
    for q in glob.glob(d + "filt_*.unw.geo.aoi90.tif"):
        m = re.search(r"_(\d+)rlks_(\d+)alks", os.path.basename(q))
        if m: looks.add((int(m.group(1)), int(m.group(2))))
# Equal counts alone miss a removed/replaced pair. Compare actual date pairs.
print(1 if (len(looks) == 1 and (rl, al) == looks.pop()
            and n == len(done_pairs) and stack_pairs == done_pairs) else 0)
PYCHK
  )
  [ "$STACK_OK" = "1" ] || echo "existing stack does not match $PAIR_ROOT (looks or pair count); rebuilding"
fi
if [ "${REBUILD:-0}" = "1" ] || [ "$STACK_OK" != "1" ] || [ ! -s "$W/inputs/geometryGeo.h5" ]; then
  if [ -e "$W" ]; then
    BACKUP="${W}.previous.$(date -u +%Y%m%dT%H%M%S).${SLURM_JOB_ID:-$$}"
    [ ! -e "$BACKUP" ] && mv "$W" "$BACKUP" || { echo "cannot preserve old work directory"; exit 1; }
    echo "previous work directory preserved at $BACKUP"
  fi
  mkdir -p "$W" || exit 1
  python ~/full/build_stack2.py "$W" || { echo "build_stack2 FAILED"; echo MINTPY_FAILED; exit 1; }
else
  echo "reusing existing stack: $(du -sm $W/inputs | cut -f1) MB in $W/inputs"
  # MintPy's update_template_file merges the requested custom settings and keeps
  # the effective configuration untouched when values are unchanged. Preserve
  # that file so a normal resume does not invalidate downstream provenance.
  # FRESH=1 explicitly restarts downstream stages while retaining phase caches.
  if [ "${FRESH:-0}" = "1" ]; then
    PRODUCT_BACKUP="$W/products.previous.$(date -u +%Y%m%dT%H%M%S).${SLURM_JOB_ID:-$$}"
    mkdir "$PRODUCT_BACKUP" || exit 1
    for product in "$W"/timeseries*.h5 "$W"/velocity*.h5 "$W"/temporalCoherence.h5 \
      "$W"/maskTempCoh.h5 "$W"/numInvIfgram.h5 "$W"/.inversion_provenance.json \
      "$W"/.bias_provenance.json "$W"/.closureCor_provenance.json "$W"/.closureCor_stamp; do
      [ ! -f "$product" ] || mv "$product" "$PRODUCT_BACKUP/" || exit 1
    done
    echo "FRESH=1: downstream products preserved at $PRODUCT_BACKUP; inversion will re-run"
  else
    n=$(ls $W/*.h5 2>/dev/null | wc -l)
    [ "$n" -gt 0 ] && echo "preserving $n existing product(s); MintPy update mode will skip completed steps"
  fi
fi
cd "$W" || exit 1
# Reference point: chosen on connectComponent coverage, not coherence.
if [ -s "$HOME/full/ref_yx.txt" ]; then
  REF_YX="$(awk '{print $1", "$2}' $HOME/full/ref_yx.txt)"
else
  REF_YX="auto"
fi
export REF_YX
echo "reference pixel: $REF_YX"
# FINDING 3 (audit r5): select_phase_closure_network.py writes a pair list with
# --write-pairs, but nothing consumed it, so the short-baseline network it
# selects was never applied and a rerun still used all 1650 pairs and recreated
# the empty common-region intersection. Translate the list into MintPy's own
# network-exclusion option so the selection actually takes effect.
PAIR_ARGS=()
[ -n "${PHASE_CLOSURE_PAIRS:-}" ] && PAIR_ARGS=(--pairs "$PHASE_CLOSURE_PAIRS")
EXCLUDE_IFG=$(python "$PROVENANCE" selection "$W" ${PAIR_ARGS[@]+"${PAIR_ARGS[@]}"}) || {
  echo "PHASE_CLOSURE_SELECT_FAILED: pair list is missing, invalid or contains unknown pairs"; exit 1;
}
echo "PHASE_CLOSURE_NETWORK_REQUESTED: excludeIfgIndex=$EXCLUDE_IFG"

cat > fullApp.cfg <<CFGEOF
mintpy.load.processor = isce
mintpy.networkInversion.minTempCoh = ${MIN_TCOH:-0.4}
mintpy.reference.yx = ${REF_YX}
# reference.method is retained only as a fallback if REF_YX is unset.
mintpy.reference.method = maxCoherence
mintpy.reference.coherenceFile = avgSpatialCoh.h5
# MintPy defaults to requiring a pixel above 0.85 average spatial coherence to
# anchor the network. Measured on this stack (2026-09-03): median 0.218, p99
# 0.337, single best pixel 0.849, and NOTHING above 0.85, so the default aborts
# at reference_point. 0.5 leaves 1234 candidate pixels and maxCoherence still
# selects the best of them.
mintpy.reference.minCoherence = 0.5
mintpy.network.coherenceBased = no
mintpy.network.excludeIfgIndex = ${EXCLUDE_IFG:-no}
mintpy.unwrapError.method = phase_closure
# [E] Yunjun, Fattahi and Amelung (2019) use phase closure to solve integer
# ambiguities where common connected components support the correction.
# A small wrapped closure residual does not establish deformation accuracy.
# The full stack's empty common-region intersection remains a measured limit.
mintpy.unwrapError.ramp = linear
mintpy.unwrapError.connCompMinArea = 2500
mintpy.reference.maskFile = no
mintpy.networkInversion.weightFunc = var
mintpy.networkInversion.maskDataset = no
mintpy.networkInversion.maskThreshold = 0.25
mintpy.troposphericDelay.method = no
mintpy.deramp = linear
mintpy.deramp.maskFile = maskTempCoh.h5
mintpy.solidEarthTides = yes
mintpy.ionosphericDelay.method = no
mintpy.topographicResidual = no
mintpy.velocity.startDate = auto
CFGEOF
# Two-phase run (rehearsal finding): the phase_closure correction may
# legitimately skip (no integer ambiguities found) without writing the
# unwrapPhase_phaseClosure dataset, and the inversion then crashes asking
# for it. Run through the correction, check the dataset, and fall back to
# plain unwrapPhase if the correction produced nothing.
# BEGIN PHASE_PREREQUISITES
PHASE_STATUS=$(python "$PROVENANCE" prepare-phase "$W" ${PAIR_ARGS[@]+"${PAIR_ARGS[@]}"}) || {
  echo "MINTPY_INCOMPLETE: cached phase provenance is unknown or changed."
  echo "Use an isolated fresh stack; REBUILD=1 preserves the old work directory before rebuilding."
  exit 3
}
# A failed new attempt must not leave the previous corrected result apparently
# current. This marker also protects comparison after interruption.
printf '%s\n' "${SLURM_JOB_ID:-manual} running" > .workflow_incomplete
if [ "$PHASE_STATUS" = "fresh" ]; then
  smallbaselineApp.py fullApp.cfg --dir "$W" --start modify_network --end correct_unwrap_error || { echo "smallbaselineApp(correct) FAILED"; echo MINTPY_FAILED; exit 1; }
  PHASE_STATUS=$(python "$PROVENANCE" record-phase "$W" ${PAIR_ARGS[@]+"${PAIR_ARGS[@]}"}) || { echo MINTPY_FAILED; exit 1; }
else
  # prepare-phase checked the exact stack, retained network, reference and
  # requested phase settings. Re-entering these completed steps would toggle
  # the effective method no -> phase_closure -> no on clean fallback resumes,
  # unnecessarily invalidating MintPy's inversion update check.
  echo "PHASE_PREREQUISITES_CURRENT: $PHASE_STATUS; skipping completed phase steps"
fi
PC_UNREPAIRED=0
if [ "$PHASE_STATUS" != "correction_available" ]; then
  # A missing unwrapPhase_phaseClosure is NOT automatically "nothing to correct".
  # MintPy also produces nothing when ambiguities are widespread but it cannot
  # find a connected region common to every interferogram - which is this stack's
  # situation (maskConnComp is empty by construction; max per-pixel coverage
  # 75.3%). Inverting plain unwrapPhase then silently accepts known-bad phase.
  # Zheng et al. 2022 corrects non-closure BIAS and does not repair gross 2-pi
  # errors, so it is not a substitute. Decide from the measured ambiguity count.
  # Read from $W, not a hardcoded work dir, and FAIL CLOSED: only a nonempty
  # finite dataset that positively shows zero ambiguities may be called clean.
  # Missing, corrupt, empty or unreadable evidence must count as unrepaired.
  if [ "$PHASE_STATUS" = "unrepaired" ]; then
    PC_UNREPAIRED=1
    echo "PC_UNREPAIRED: integer ambiguities remain or could not be verified."
  else
    echo "phase_closure found no ambiguities to correct; inverting plain unwrapPhase"
  fi
  sed -i 's/mintpy.unwrapError.method = phase_closure/mintpy.unwrapError.method = no/' fullApp.cfg
fi
export PC_UNREPAIRED
# END PHASE_PREREQUISITES
# BEGIN UNREPAIRED_STOP
if [ "${PC_UNREPAIRED:-0}" = "1" ] && [ "${ALLOW_UNREPAIRED:-0}" != "1" ]; then
  echo "MINTPY_INCOMPLETE: stopping before inversion because phase is unrepaired."
  echo "ALLOW_UNREPAIRED=1 explicitly permits a diagnostic run; it cannot report success."
  exit 3
fi
# END UNREPAIRED_STOP
# Validate downstream assumptions before paying for the inversion. Exploratory
# phase-only networks may set APPLY_ZHENG=0; nl stays an explicit assumption.
if [ "$APPLY_ZHENG" = "1" ]; then
  python "$PROVENANCE" zheng-network "$W" --bw "$ZHENG_BW" --nl "$ZHENG_NL" || {
    echo "MINTPY_INCOMPLETE: requested Zheng network is incompatible; stopping before inversion"; exit 3;
  }
fi
smallbaselineApp.py fullApp.cfg --dir "$W" --start invert_network --end velocity || { echo "smallbaselineApp(invert) FAILED"; echo MINTPY_FAILED; exit 1; }
INVERSION_ARTIFACTS=(--artifact velocity.h5 --artifact timeseries.h5 --artifact temporalCoherence.h5
  --artifact smallbaselineApp.cfg --artifact inputs/geometryGeo.h5)
[ -f "$W/timeseries_SET.h5" ] && INVERSION_ARTIFACTS+=(--artifact timeseries_SET.h5)
[ -f "$W/maskTempCoh.h5" ] && INVERSION_ARTIFACTS+=(--artifact maskTempCoh.h5)
python "$PROVENANCE" stamp "$W" --stage inversion "${INVERSION_ARTIFACTS[@]}" || exit 3
echo "--- closure bias correction (Zheng et al. 2022, MintPy implementation) ---"
CC_CURRENT=0
if [ "$APPLY_ZHENG" = "1" ]; then
  # DEFECT 13 GUARD. closure_phase_bias caches intermediates under
  # closurePhase/conn*/ and unwrap_closure_phase SKIPS regeneration when the
  # .cor/.unw already exist. Those intermediates carry their OWN .rsc metadata,
  # which isce_utils.unwrap_snaphu reads via readfile.read_attribute - NOT the
  # HDF5 products. So a stack whose attributes were fixed AFTER the intermediates
  # were written silently reuses the stale metadata and dies with KeyError.
  # That is exactly how HEIGHT failed twice. Invalidate the cache whenever the
  # stack is newer than the intermediates.
  # BEGIN CLOSUREPHASE_GUARD
  if [ -d "$W/closurePhase" ]; then
    # DEFECT 17 (found by executing the guard, not by reading it): the original
    # used `stat -c %Y`, which is GNU-only, and `find -printf`, likewise. When
    # either is unavailable or the stack is missing, the variables came back
    # empty, the -n tests failed, and the guard SILENTLY DID NOTHING - reusing
    # exactly the stale cache it exists to discard. Portable, and fails CLOSED.
    _mtime() {
      stat -c %Y "$1" 2>/dev/null || stat -f %m "$1" 2>/dev/null
    }
    OLDEST_RSC=""
    for _r in $(find "$W/closurePhase" -name "*.rsc" 2>/dev/null); do
      _t=$(_mtime "$_r"); [ -z "$_t" ] && continue
      [ -z "$OLDEST_RSC" ] && OLDEST_RSC="$_t"
      [ "${_t%.*}" -lt "${OLDEST_RSC%.*}" ] && OLDEST_RSC="$_t"
    done
    STACK_T=$(_mtime "$W/inputs/ifgramStack.h5")
    if [ -z "$OLDEST_RSC" ] || [ -z "$STACK_T" ]; then
      echo "CLOSUREPHASE_STALE: cannot compare cache and stack timestamps"
      echo "  (rsc='$OLDEST_RSC' stack='$STACK_T'); failing CLOSED and discarding"
      echo "  $W/closurePhase rather than risking stale metadata."
      rm -rf "$W/closurePhase"
    elif [ "${STACK_T%.*}" -gt "${OLDEST_RSC%.*}" ]; then
      echo "CLOSUREPHASE_STALE: ifgramStack.h5 is newer than the cached closure-phase"
      echo "  intermediates, whose .rsc metadata would be reused as-is. Discarding"
      echo "  $W/closurePhase so it is regenerated from current stack metadata."
      rm -rf "$W/closurePhase"
    else
      echo "closure-phase cache is current (stack not newer than intermediates)"
    fi
  fi
  # END CLOSUREPHASE_GUARD
  for ACTION in mask quick_estimate estimate; do
    if ! closure_phase_bias.py -i inputs/ifgramStack.h5 --nl "$ZHENG_NL" --bw "$ZHENG_BW" \
         -a "$ACTION" -c local --num-worker 16 >> closure_bias.log 2>&1; then
      echo "CLOSURE_BIAS_ESTIMATE_FAILED: action=$ACTION (see closure_bias.log)"; exit 3
    fi
  done
  python "$PROVENANCE" stamp "$W" --stage bias --artifact timeseriesBias.h5 || exit 3
  # The shared apply path subtracts timeseriesBias.h5 from the PRE-DERAMP
  # series and uses MintPy's same remove_ramp.py type AND spatial mask.
  if W="$W" PROVENANCE="$PROVENANCE" ALLOW_UNREPAIRED="$ALLOW_UNREPAIRED" RUN_FROM_WORKFLOW=1 \
       bash "$HOME/full/apply_zheng.sh"; then
    CC_CURRENT=1
  else
    echo CLOSURE_BIAS_APPLY_FAILED
  fi
else
  echo "ZHENG_DISABLED: exploratory uncorrected inversion, not the full baseline workflow"
fi
echo "--- outputs:"
ls -la $W/*.h5 2>/dev/null | head -20
for f in velocity velocityClosureCor temporalCoherence avgSpatialCoh; do
  [ "$f" = "velocityClosureCor" ] && [ "$CC_CURRENT" != "1" ] && continue
  [ -f $W/$f.h5 ] && save_gdal.py $W/$f.h5 -o $W/$f.tif 2>/dev/null && echo "exported $f.tif"
done
# CC_CURRENT is set only when THIS invocation of apply_zheng.sh succeeded.
# .closureCor_stamp is retained for human-readable history, not used as proof.
# PC_UNREPAIRED takes precedence over EVERY success branch (audit r6 finding 1).
# Zheng et al. 2022 corrects non-closure BIAS; it does not solve the integer
# ambiguities that Yunjun, Fattahi & Amelung 2019 (phase_closure) targets, so a
# Zheng product sitting on unrepaired 2-pi errors must not print MINTPY_OK.
if [ "${PC_UNREPAIRED:-0}" = "1" ]; then
  echo "MINTPY_INCOMPLETE: gross integer ambiguities were NOT repaired"
  [ "$CC_CURRENT" = "1" ] && echo "  (velocityClosureCor.h5 was produced, but does not license success)"
  FINAL_RC=3
elif [ "$CC_CURRENT" = "1" ]; then
  echo MINTPY_OK; FINAL_RC=0
elif [ "$APPLY_ZHENG" = "1" ]; then
  echo "MINTPY_INCOMPLETE: requested closure bias correction was not applied"
  FINAL_RC=3
elif [ -f $W/velocity.h5 ]; then echo MINTPY_EXPLORATORY_UNCORRECTED; FINAL_RC=0
else echo MINTPY_FAILED; FINAL_RC=1
fi
[ "$FINAL_RC" = "0" ] && rm -f "$W/.workflow_incomplete"
exit $FINAL_RC
