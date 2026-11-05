#!/bin/bash
# Reproduce the alosStack processing environment exactly.
#
# Pinned because version drift here is not cosmetic: MintPy 1.6.4 and ISCE 2.6.5
# are what every result in this repo was produced with, and Tay et al. 2025 used
# MintPy 1.5.3 / ISCE 2.6.3, so the difference must be recorded rather than
# discovered later.
set -euo pipefail
ENV_NAME="${ENV_NAME:-insar}"
PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
MM="${MM:-$HOME/bin/micromamba}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

[ -x "$MM" ] || { echo "micromamba not found at $MM (set MM=)"; exit 1; }
export MAMBA_ROOT_PREFIX="$PREFIX"

echo "== creating $ENV_NAME from the lockfile =="
# Run this on a COMPUTE NODE: the login node has ulimit -v 1 GB and the solver
# will die with 'Out of memory allocating ...' partway through.
"$MM" create -y -n "$ENV_NAME" -c conda-forge --file "$HERE/insar-env.lock"

E="$PREFIX/envs/$ENV_NAME"

echo "== linking the snaphu BINARY onto PATH =="
# The conda 'snaphu' package (build py311...) installs a PYTHON WRAPPER and puts
# NO snaphu command on PATH. ISCE decides whether to run the full statistical
# solver with shutil.which('snaphu'); without it, Alos2Proc/runUnwrapSnaphu.py
# falls back to snaphuUnwrap(..., initOnly=True), which runs minimum-cost-flow
# INITIALISATION ONLY and never invokes the statistical cost solver. Every
# interferogram this project produced before 2026-09-03 was unwrapped that way.
# The real Stanford binary (v2.0.7) ships INSIDE the python package.
BIN="$E/lib/python3.11/site-packages/snaphu/snaphu"
[ -x "$BIN" ] || { echo "bundled snaphu binary not found at $BIN"; exit 1; }
ln -sfn "$BIN" "$E/bin/snaphu"

echo "== verification =="
# Check with shutil.which, NOT `which`: shutil.which is the test ISCE performs.
"$E/bin/python" - <<'PY'
import shutil, sys
w = shutil.which("snaphu")
print("shutil.which('snaphu') ->", w)
if not w:
    print("FAIL: ISCE would silently take the initOnly branch"); sys.exit(1)
import mintpy, isce
print("mintpy", mintpy.__version__)
try:
    import isce as _i; print("isce ok")
except Exception as e:
    print("isce import:", e)
print("ENV_OK")
PY
echo
echo "Recorded provenance: this environment is MintPy 1.6.4 / ISCE 2.6.5 /"
echo "snaphu 2.0.7. Tay et al. 2025 used MintPy 1.5.3 / ISCE 2.6.3."
