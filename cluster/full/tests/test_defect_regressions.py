"""Regression tests for the defects found in the 2026-09-03/04 audit rounds.

These are pure static checks on the shipped scripts: no cluster, no ISCE, no
data. They exist because every one of these defects passed syntax checking, ran
to completion, and reported success while producing a wrong or missing result.
Each test names the defect it guards.
"""
import os
import re

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def _read(name):
    p = os.path.join(ROOT, name)
    if not os.path.exists(p):
        pytest.skip("%s not present" % name)
    return open(p).read()


# --- defect 4: run_pair executes in an `xargs bash -c` subshell -------------
def test_main_exports_every_var_run_pair_reads():
    s = _read("main.sbatch")
    # There are SEVERAL bare `export A B C` lines (one per exported function's
    # environment). Union them all: any of them puts the name in the child env.
    exported = set()
    for ln in s.splitlines():
        if ln.startswith("export ") and "=" not in ln and not ln.startswith("export -f"):
            exported.update(ln.split()[1:])
    assert exported, "no bare export line found in main.sbatch"
    body = s[s.index("run_pair()"):]
    if "export -f run_pair" in body:
        body = body[:body.index("export -f run_pair")]
    used = set(re.findall(r"\$\{?([A-Z][A-Z0-9_]{2,})\b", body))
    assigned = set(re.findall(r"^([A-Z][A-Z0-9_]*)=", s, re.M))
    local = {"SLURM_ARRAY_TASK_ID", "SLURM_JOB_ID", "PATH", "HOME", "TMPDIR",
             "ISCE_STACK", "PROJ_DATA", "OMP_NUM_THREADS"}
    missing = sorted((used & assigned) - exported - local)
    assert not missing, (
        "read inside run_pair but never exported (defect 4 class): %s" % missing)


# --- defect 5: looks metadata must not be hardcoded ------------------------
def test_build_stack2_derives_looks_not_hardcoded():
    s = _read("build_stack2.py")
    assert "_looks_from_products" in s, "looks must be derived from product names"
    assert not re.search(r'"ALOOKS"\s*:\s*"\d+"', s), "ALOOKS is hardcoded again"
    assert not re.search(r'"RLOOKS"\s*:\s*"\d+"', s), "RLOOKS is hardcoded again"


# --- defects 9 and 11: metadata that downstream MintPy steps require -------
@pytest.mark.parametrize("attr", ["CENTER_LINE_UTC", "HEIGHT",
                                  "CENTER_INCIDENCE_ANGLE", "EARTH_RADIUS"])
def test_build_stack2_writes_required_metadata(attr):
    # correct_SET died on CENTER_LINE_UTC and closure_phase_bias on HEIGHT,
    # each AFTER hours of successful work.
    assert attr in _read("build_stack2.py"), "%s missing from stack attrs" % attr


# --- defect 6: the pair root must be selectable ----------------------------
def test_build_stack2_pair_root_is_parameterised():
    s = _read("build_stack2.py")
    assert "PAIR_ROOT" in s
    assert 'os.path.expanduser("~/full/pairs")' not in s.replace("PAIR_ROOT", ""), \
        "hardcoded ~/full/pairs would silently rebuild from the old tree"


# --- defect 12: Zheng apply path ------------------------------------------
def test_zheng_apply_does_not_hardcode_demErr_names():
    s = _read("apply_zheng.sh")
    assert "timeseriesBias.h5" in s, "closure bias apply step is absent"
    bad = re.search(r"TS=\$\(ls \$W/timeseries_ramp_demErr\.h5", s)
    assert not bad, "apply path hardcodes files this config never produces"


def test_zheng_applied_before_deramp():
    # deramp(D) - B != deramp(D - B) unless B has zero planar component.
    s = _read("apply_zheng.sh")
    assert "PRE-DERAMP" in s, "correction order not documented"
    assert "remove_ramp.py" in s, "bias must be applied pre-deramp then re-ramped"
    # check CODE, not prose: comments legitimately name the old broken filename
    code = "\n".join(l for l in s.splitlines() if not l.lstrip().startswith("#"))
    assert code.index('diff.py "$BASE"') < code.index('remove_ramp.py ')
    assert 'apply-inputs "$W"' in code, "base must come from verified correction configuration"
    assert '-m "$RAMP_MASK"' in code, "re-ramp must use the original fitted mask"


# --- defect 10 and stale-product provenance --------------------------------
def test_resume_does_not_delete_products_by_default():
    s = _read("mintpy2.sbatch")
    assert "FRESH" in s, "product deletion must be gated"
    assert not re.search(r"^\s*rm -f \$W/\*\.h5 \$W/smallbaselineApp\.cfg", s, re.M), \
        "unconditional product deletion destroys a completed inversion on resume"


def test_success_requires_current_run_provenance():
    s = _read("mintpy2.sbatch")
    assert ".closureCor_stamp" in s, "no provenance stamp on the corrected product"
    assert "CC_CURRENT" in s, "success gate does not check the stamp"


def test_unapplied_correction_is_not_reported_as_success():
    s = _read("mintpy2.sbatch")
    assert "FINAL_RC=3" in s, "an estimated-but-unapplied bias must exit non-zero"


# --- defect 7: stack reuse must verify it matches the products -------------
def test_stack_reuse_is_guarded():
    s = _read("mintpy2.sbatch")
    assert "STACK_OK" in s, "stale-stack guard missing"
    assert "ALOOKS" in s and "RLOOKS" in s, "guard must compare looks"


# --- Tay comparison: zero-fill and framing ---------------------------------
def test_compare_tay_masks_zero_fill_and_is_labelled():
    s = _read("compare_tay.py")
    assert "srcNodata" in s, "their raster declares no NoData; zeros must be masked"
    assert "valid > 0.999" in s or "valid>0.999" in s, "validity mask not enforced"
    assert "NEVER CALIBRATE" in s.upper()
    assert "not independent validation" in s.lower(), "same-sensor caveat missing"
    assert "demeaned" in s.lower(), "datum-contaminated RMSE reported unqualified"


# --- audit r6 finding 1: unrepaired phase must outrank every success --------
def test_unrepaired_phase_outranks_success():
    s = _read("mintpy2.sbatch")
    i_pc = s.find('if [ "${PC_UNREPAIRED:-0}" = "1" ]; then')
    i_cc = s.find('elif [ "$CC_CURRENT" = "1" ]; then')
    assert i_pc != -1 and i_cc != -1, "final status gate not found"
    assert i_pc < i_cc, (
        "CC_CURRENT is checked before PC_UNREPAIRED: a Zheng product could report "
        "MINTPY_OK while gross 2-pi ambiguities remain unrepaired")


# --- audit r6 finding 2: the ambiguity probe must fail CLOSED --------------
def test_ambiguity_probe_uses_workdir_and_fails_closed():
    s = _read("mintpy2.sbatch")
    assert 'expanduser("~/full/mintpy_tay/numTri' not in s, \
        "ambiguity probe reads a hardcoded work dir instead of $W"
    assert 'record-phase "$W"' in s and '"$PHASE_STATUS" = "unrepaired"' in s, \
        "probe must use the shared fail-closed provenance helper"
    assert "print(-1)" not in s, \
        "a -1 sentinel is silently treated as 'no ambiguities'"


# --- defect 13: stale closure-phase intermediates carry their own metadata ---
def test_closure_phase_cache_is_invalidated_when_stack_changes():
    # closure_phase_bias caches closurePhase/conn*/ and unwrap_closure_phase
    # skips regeneration when .cor/.unw exist. Those intermediates carry their
    # OWN .rsc, which isce_utils.unwrap_snaphu reads instead of the HDF5
    # products, so a stack fixed after they were written is silently ignored.
    s = _read("mintpy2.sbatch")
    assert "CLOSUREPHASE_STALE" in s, "no staleness guard on the closure-phase cache"
    i_guard = s.find("CLOSUREPHASE_STALE")
    i_run = s.find("closure_phase_bias.py -i inputs/ifgramStack.h5 --nl")
    assert i_guard < i_run, "staleness check must run BEFORE closure_phase_bias"
    assert 'rm -rf "$W/closurePhase"' in s, "stale cache is detected but not discarded"


def test_rsc_backfill_uses_a_separator():
    # A key of 20+ chars with "%-20s%s" produces no separator, so
    # readfile.read_roipac_rsc splits one token and raises IndexError on c[1].
    # CENTER_INCIDENCE_ANGLE is 22 chars, which is how this actually broke.
    for name in ("mintpy2.sbatch",):
        s = _read(name)
        assert '"%-20s%s' not in s and "'%-20s%s" not in s, (
            "%-20s padding is too narrow for CENTER_INCIDENCE_ANGLE (22 chars) "
            "and yields an unparseable .rsc line")


# --- audit r5 finding 3: the selected pair list must actually be applied -----
def test_phase_closure_pair_list_has_a_consumer():
    s = _read("mintpy2.sbatch")
    assert "PHASE_CLOSURE_PAIRS" in s, (
        "select_phase_closure_network.py --write-pairs has no consumer, so the "
        "selected short-baseline network is never applied and a rerun uses all "
        "1650 pairs again")
    assert "mintpy.network.excludeIfgIndex" in s, \
        "selection is read but never turned into a MintPy network restriction"
    i_read = s.index("PHASE_CLOSURE_PAIRS")
    i_cfg = s.index("cat > fullApp.cfg")
    assert i_read < i_cfg, "exclusion list must be computed before the config is written"


def test_failed_selection_does_not_silently_use_full_network():
    s = _read("mintpy2.sbatch")
    assert "PHASE_CLOSURE_SELECT_FAILED" in s, (
        "a pair list that yields no exclusions must fail loudly; silently "
        "inverting the full network is the failure this option exists to avoid")


# --- audit r6 finding 4: viability must consider WHERE the region is --------
def test_selector_viability_requires_aoi_overlap():
    s = _read("select_phase_closure_network.py")
    assert "--aoi-bbox" in s, "no operational-AOI option on the selector"
    assert "min_aoi_pixels" in s, "no minimum overlap requirement"
    i = s.index("viable = (")
    seg = s[i:i + 500]
    assert "aoi_ok" in seg, (
        "viability is decided on region size alone; a region sitting entirely "
        "on the mineral hills would be called viable")


# --- audit r6 finding 5: closure sampling must validate and never format None -
def test_selector_validates_reference_and_handles_no_samples():
    s = _read("select_phase_closure_network.py")
    assert "def _fmt(" in s, "None statistics would crash on float formatting"
    assert "interferograms_with_invalid_reference" in s, \
        "a zero or non-finite reference phase makes every sample meaningless"
    assert "finite_coverage_fraction" in s, \
        "fraction silently describes only the surviving samples"


def test_defect_23_deramp_type_read_from_the_cfg_mintpy_ran_with():
    """Defect 23: both Zheng apply paths read mintpy.deramp from fullApp.cfg,
    which does not exist at $W (only inputs/fullApp.cfg does), so the read was
    always empty and a hard-coded 'linear' was used. Correct by coincidence
    today; silently wrong for any other deramp setting."""
    import os, re
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    src = _read("workflow_provenance.py")
    apply = src[src.index("def apply_inputs("):src.index("def select_velocity(")]
    assert 'read_config(work / "smallbaselineApp.cfg")' in apply
    assert 'cfg.get("mintpy.deramp", "")' in apply
    assert "unsupported or unresolved deramp setting" in apply


def test_era_hooks_are_inert_by_default_and_exported():
    """The era test adds three env overrides. Each must default to the old
    literal (so the certified campaign is unchanged) and STOP_AFTER must be in
    run_pair's export list, or it is unset inside the xargs subshell (defect 1)."""
    import os, re
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dl2 = open(os.path.join(root, "dl2.sbatch")).read()
    assert 'SCENES=${SCENES:-scenes_full.txt}' in dl2 and '$SCENES)' in dl2
    w1 = open(os.path.join(root, "wave1.sbatch")).read()
    assert 'WAVE1_LIST=${WAVE1_LIST:-$HOME/full/wave1_dates.txt}' in w1 and '$WAVE1_LIST)' in w1
    m = open(os.path.join(root, "main.sbatch")).read()
    assert '[ "${STOP_AFTER:-}" = "look_coherence" ]' in m
    exp = [l for l in m.splitlines() if l.startswith("export WORK REF OUT")]
    assert len(exp) == 1 and "STOP_AFTER" in exp[0], "STOP_AFTER not exported to the pair subshell"
    # the hook must write DONE_COH and never DONE
    hook = m.split('STOP_AFTER:-}" = "look_coherence"')[1].split("return 0")[0]
    assert "DONE_COH" in hook and "touch $od/DONE\n" not in hook


def test_defect_25_main_defaults_to_the_certified_looks_and_announces_them():
    """Defect 25: main.sbatch defaulted NRLKS2/NALKS2 to 5 and 2 (an abandoned
    experiment), so a launch that omitted the looks silently produced 5x28
    products. The default must be the certified 17x4 (ml2 17x56) and the run
    must print the looks in force."""
    import os, re
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    m = open(os.path.join(root, "main.sbatch")).read()
    assert "NRLKS2=${NRLKS2:-17}; NALKS2=${NALKS2:-4}" in m
    assert not re.search(r"NRLKS2=\$\{NRLKS2:-5\}", m)
    assert "MAIN_LOOKS nrlks2=$NRLKS2 nalks2=$NALKS2" in m
    # the announcement must come after NRLKS1/NALKS1 are set, or it prints garbage
    assert m.index("NRLKS1=1; NALKS1=14") < m.index("MAIN_LOOKS nrlks2")


def test_pfail_lists_product_sizes():
    """After 821292 (snaphu 'wrong size' with the pair dir gone on node /tmp),
    pfail must list the pair's product sizes so a size mismatch is visible in
    the block log."""
    import os
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    m = open(os.path.join(root, "main.sbatch")).read()
    body = m.split("pfail() {")[1].split("\n}")[0]
    assert "products (bytes name)" in body and "/insar" in body and "ls -l" in body


def test_phsig_fallback_is_guarded_and_fails_closed():
    """ICU writes an empty phase sigma on short images (240 rows). The fallback
    must run only when the file is empty, and must fail the pair if it still
    produces nothing rather than let snaphu die later with 'wrong size'."""
    import os
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    m = open(os.path.join(root, "main.sbatch")).read()
    i = m.index('if [ ! -s "insar/${p}${ML2}.phsig" ]')
    block = m[i:m.index("\n  fi\n", i)]
    assert "cal_coherence_1" in block and "PHSIG_FALLBACK" in block
    assert 'phsig fallback produced nothing"; return 1' in block
    assert m.index("filt.py -idir") < i < m.index("unwrap_snaphu.py")
