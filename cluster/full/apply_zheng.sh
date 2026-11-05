#!/bin/bash
#SBATCH --job-name=zhapply
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/home/t/tommy/full/logs/zheng_apply.log
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj
export ISCE_STACK=$HOME/micromamba/envs/insar/share/isce2
export PATH=$ISCE_STACK/alosStack:$PATH
W=${W:-$HOME/full/mintpy_tay}
PROVENANCE=${PROVENANCE:-$HOME/full/workflow_provenance.py}
cd "$W" || exit 1
DIAGNOSTIC_ARGS=()
[ "${ALLOW_UNREPAIRED:-0}" = "1" ] && DIAGNOSTIC_ARGS=(--allow-unrepaired)
# Check phase, retained network, base/bias acquisition dates and configuration
# before any existing product can be overwritten or restamped.
SETTINGS=$(python "$PROVENANCE" apply-inputs "$W" ${DIAGNOSTIC_ARGS[@]+"${DIAGNOSTIC_ARGS[@]}"}) || {
  echo APPLY_FAILED_unverified_inputs; exit 3;
}
BASE_NAME=$(printf '%s\n' "$SETTINGS" | sed -n '1p')
RAMP=$(printf '%s\n' "$SETTINGS" | sed -n '2p')
RAMP_MASK=$(printf '%s\n' "$SETTINGS" | sed -n '3p')
BASE="$W/$BASE_NAME"
printf '%s\n' "${SLURM_JOB_ID:-manual} applying closure bias" > .workflow_incomplete

# CORRECTED ORDER (audit r4 finding 2): deramp and bias subtraction do not
# commute, so apply the bias to the PRE-deramp series and re-run the same deramp.
echo "CLOSURE_BIAS_APPLY_BASE $BASE_NAME (PRE-DERAMP)"

diff.py "$BASE" "$W/timeseriesBias.h5" -o "$W/timeseries_closureCor_noramp.h5" || { echo APPLY_FAILED_diff; exit 1; }
# DEFECT 23 (round 5): this read fullApp.cfg, which does not exist at $W, so RAMP
# was always empty and fell back to a hard-coded linear. Read the cfg MintPy
# actually ran with, and refuse to guess if it is missing: a wrong re-ramp makes a
# product that silently disagrees with velocity.h5.
CFG=$W/smallbaselineApp.cfg
[ -f "$CFG" ] || { echo APPLY_FAILED_no_cfg; exit 1; }
if [ "$RAMP" = "no" ]; then
  cp "$W/timeseries_closureCor_noramp.h5" "$W/timeseries_closureCor.h5" || { echo APPLY_FAILED_copy; exit 1; }
else
  remove_ramp.py "$W/timeseries_closureCor_noramp.h5" -s "$RAMP" -m "$RAMP_MASK" \
    -o "$W/timeseries_closureCor.h5" || { echo APPLY_FAILED_ramp; exit 1; }
  echo "CLOSURE_BIAS_RERAMPED with $RAMP using mask $RAMP_MASK"
fi
timeseries2velocity.py "$W/timeseries_closureCor.h5" -t "$CFG" -o "$W/velocityClosureCor.h5" || { echo APPLY_FAILED_vel; exit 1; }
ARTIFACTS=(--artifact velocityClosureCor.h5 --artifact timeseries_closureCor.h5 --artifact "$BASE_NAME"
  --artifact timeseriesBias.h5 --artifact smallbaselineApp.cfg --artifact .inversion_provenance.json
  --artifact .bias_provenance.json --artifact temporalCoherence.h5 --artifact inputs/geometryGeo.h5)
[ "$RAMP" != "no" ] && [ "$RAMP_MASK" != "no" ] && ARTIFACTS+=(--artifact "$RAMP_MASK")
python "$PROVENANCE" stamp "$W" --stage closureCor "${ARTIFACTS[@]}" || { echo APPLY_FAILED_provenance; exit 3; }
date -u +%Y-%m-%dT%H:%M:%SZ > "$W/.closureCor_stamp"
echo "${SLURM_JOB_ID:-manual}" >> "$W/.closureCor_stamp"
echo CLOSURE_BIAS_APPLIED
if [ "${RUN_FROM_WORKFLOW:-0}" != "1" ] && [ "${ALLOW_UNREPAIRED:-0}" != "1" ]; then
  rm -f "$W/.workflow_incomplete"
fi
ls -la "$W/velocityClosureCor.h5"
