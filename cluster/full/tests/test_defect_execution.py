"""EXECUTABLE tests for the fixes that static checks cannot verify.

tests/test_defect_regressions.py greps the scripts: it proves a fix is PRESENT,
not that it WORKS. That distinction matters here, because three of the sixteen
defects fixed on 2026-09-04/05 were found only by running the pipeline and
reading a traceback, and no static check would have surfaced them.

These tests build synthetic fixtures and actually execute the logic.
"""
import os
import re
import subprocess
import sys
import time
from pathlib import Path
import importlib.util

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
_SPEC = importlib.util.spec_from_file_location("workflow_provenance", Path(ROOT) / "workflow_provenance.py")
workflow = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(workflow)


# --- defect 13: the stale closure-phase cache guard -------------------------
def _guard_snippet():
    """Extract the CLOSUREPHASE_STALE guard so it can be run standalone."""
    s = open(os.path.join(ROOT, "mintpy2.sbatch")).read()
    i = s.index("# BEGIN CLOSUREPHASE_GUARD")
    j = s.index("# END CLOSUREPHASE_GUARD")
    return s[i:j]


def _run_guard(tmp_path, stack_newer):
    w = tmp_path / "work"
    (w / "closurePhase" / "conn2").mkdir(parents=True)
    (w / "inputs").mkdir(parents=True)
    rsc = w / "closurePhase" / "conn2" / "filt_01.int.rsc"
    rsc.write_text("WIDTH 10\n")
    stack = w / "inputs" / "ifgramStack.h5"
    stack.write_bytes(b"x")
    old = time.time() - 3600
    os.utime(rsc, (old, old))
    if stack_newer:
        os.utime(stack, (time.time(), time.time()))
    else:
        os.utime(stack, (old - 3600, old - 3600))
    script = tmp_path / "g.sh"
    script.write_text('W="%s"\n' % w + _guard_snippet())
    out = subprocess.run(["bash", str(script)], capture_output=True, text=True)
    return out.stdout, (w / "closurePhase").exists()


def test_stale_closure_cache_is_discarded(tmp_path):
    out, still_there = _run_guard(tmp_path, stack_newer=True)
    assert "CLOSUREPHASE_STALE" in out, out
    assert not still_there, "stale cache was detected but not removed"


def test_current_closure_cache_is_kept(tmp_path):
    out, still_there = _run_guard(tmp_path, stack_newer=False)
    assert "CLOSUREPHASE_STALE" not in out, out
    assert still_there, "a current cache must not be discarded"


# --- defect 14: the pair-list consumer --------------------------------------
def test_pair_list_translates_to_exclusion_indices(tmp_path):
    h5py = pytest.importorskip("h5py")
    import numpy as np
    w = tmp_path / "work"
    (w / "inputs").mkdir(parents=True)
    dates = [("20200101", "20200113"), ("20200113", "20200125"),
             ("20200125", "20200206"), ("20200101", "20200125")]
    with h5py.File(w / "inputs" / "ifgramStack.h5", "w") as f:
        f.create_dataset("date", data=np.array(
            [[a.encode(), b.encode()] for a, b in dates]))
    keep = tmp_path / "pairs.txt"
    keep.write_text("20200101_20200113\n20200113_20200125\n")

    out = subprocess.run([sys.executable, str(Path(ROOT) / "workflow_provenance.py"),
                          "selection", str(w), "--pairs", str(keep)],
                         capture_output=True, text=True)
    got = out.stdout.strip()
    # indices 2 and 3 are not in the keep list
    assert got == "2,3", "expected '2,3', got %r (stderr %s)" % (got, out.stderr[:200])


def test_pair_list_matching_nothing_is_reported(tmp_path):
    h5py = pytest.importorskip("h5py")
    import numpy as np
    w = tmp_path / "work"
    (w / "inputs").mkdir(parents=True)
    with h5py.File(w / "inputs" / "ifgramStack.h5", "w") as f:
        f.create_dataset("date", data=np.array([[b"20200101", b"20200113"]]))
    keep = tmp_path / "pairs.txt"
    keep.write_text("19990101_19990113\n")
    out = subprocess.run([sys.executable, str(Path(ROOT) / "workflow_provenance.py"),
                          "selection", str(w), "--pairs", str(keep)],
                         capture_output=True, text=True)
    assert out.stdout.strip() == "", "must emit no exclusion list when nothing matches"
    assert out.returncode == 3
    assert "absent from stack" in out.stderr


# --- defect 16: closure statistics must never crash on None -----------------
def test_fmt_handles_unavailable_statistics():
    sys.path.insert(0, ROOT)
    src = open(os.path.join(ROOT, "select_phase_closure_network.py")).read()
    ns = {}
    m = re.search(r"def _fmt\(.*?\n(?:.*\n)*?    return format\(value, spec\)\n", src)
    assert m, "_fmt helper not found"
    exec(m.group(0), ns)
    # the exact crash that occurred: TypeError on NoneType.__format__
    assert "n/a" in ns["_fmt"](None)
    assert ns["_fmt"](0.5) == "0.500000"
    assert ns["_fmt"](0.5, ".2f") == "0.50"


# Defect 18 now uses network/reference/configuration manifests. See the real
# HDF5 acceptance, missing-input, raw-timeseries-change and legacy-stamp cases
# in test_workflow_provenance.py; timestamps alone no longer establish validity.


def test_defect_23_deramp_read_executes_correctly(tmp_path):
    """Round 6: execute the shared config parser used by the apply path against
    three cfgs. A trailing comment must be stripped, 'no' must survive, and a
    missing key must yield an empty string so the guard can refuse to guess."""
    import os, re, subprocess
    here = os.path.dirname(os.path.abspath(__file__)); root = os.path.dirname(here)
    cases = {
        "commented": ("mintpy.deramp = linear  #[no / linear / quadratic], auto for no\n", "linear"),
        "plain_no": ("mintpy.deramp = no\n", "no"),
        "quadratic_tight": ("mintpy.deramp=quadratic\n", "quadratic"),
        "absent": ("mintpy.other = 1\n", ""),
    }
    for case, (text, want) in cases.items():
        cfg = tmp_path / f"{case}.cfg"
        cfg.write_text(text)
        got = workflow.read_config(cfg).get("mintpy.deramp", "")
        assert got == want, f"[{case}]: parser gave {got!r}, want {want!r}"
