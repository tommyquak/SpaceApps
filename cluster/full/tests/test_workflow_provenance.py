"""Behavioral tests for phase cache reuse and downstream product acceptance."""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("workflow_provenance", ROOT / "workflow_provenance.py")
wf = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(wf)


def make_work(tmp_path, ambiguous=False):
    work = tmp_path / "work"
    (work / "inputs").mkdir(parents=True)
    dates = ["20200101", "20200113", "20200125", "20200206"]
    pairs = [(a, b) for i, a in enumerate(dates) for b in dates[i + 1:]]
    attrs = dict(REF_Y="1", REF_X="1", X_FIRST="114", X_STEP="0.001",
                 Y_FIRST="-2", Y_STEP="-0.001", UNIT="m", PLATFORM="ALOS-2")
    with h5py.File(work / "inputs/ifgramStack.h5", "w") as f:
        f.create_dataset("date", data=np.array(pairs, dtype="S8"))
        f.create_dataset("dropIfgram", data=np.ones(len(pairs), dtype=bool))
        f.create_dataset("unwrapPhase", data=np.ones((len(pairs), 3, 3), dtype="float32"))
        f.create_dataset("connectComponent", data=np.ones((len(pairs), 3, 3), dtype="uint8"))
        f.attrs.update(attrs)
    cfg = ("mintpy.unwrapError.method = phase_closure\n"
           "mintpy.reference.yx = 1, 1\n"
           "mintpy.solidEarthTides = yes\n"
           "mintpy.deramp = linear  # same fitted domain\n"
           "mintpy.deramp.maskFile = maskTempCoh.h5\n"
           "mintpy.troposphericDelay.method = no\n"
           "mintpy.ionosphericDelay.method = no\n"
           "mintpy.topographicResidual = no\n")
    for name in ("fullApp.cfg", "smallbaselineApp.cfg"):
        (work / name).write_text(cfg)
    for name in ("timeseries.h5", "timeseries_SET.h5", "timeseriesBias.h5"):
        with h5py.File(work / name, "w") as f:
            f.create_dataset("date", data=np.array(dates, dtype="S8"))
            f.create_dataset("timeseries", data=np.zeros((len(dates), 3, 3), dtype="float32"))
            f.attrs.update(attrs)
    for name, dataset in (("velocity.h5", "velocity"), ("temporalCoherence.h5", "temporalCoherence"),
                          ("maskTempCoh.h5", "mask"), ("inputs/geometryGeo.h5", "incidenceAngle"),
                          ("numTriNonzeroIntAmbiguity.h5", "mask")):
        with h5py.File(work / name, "w") as f:
            data = np.zeros((3, 3), dtype="float32")
            if name == "numTriNonzeroIntAmbiguity.h5" and ambiguous:
                data[0, 0] = 4
            f.create_dataset(dataset, data=data)
            f.attrs.update(attrs)
    wf.record_phase(work)
    wf.stamp_stage(work, "inversion", ["velocity.h5", "timeseries.h5", "timeseries_SET.h5",
                                     "smallbaselineApp.cfg", "temporalCoherence.h5",
                                     "maskTempCoh.h5", "inputs/geometryGeo.h5"])
    wf.stamp_stage(work, "bias", ["timeseriesBias.h5"])
    return work, pairs


def corrected_product(work):
    (work / "velocityClosureCor.h5").write_bytes((work / "velocity.h5").read_bytes())
    wf.stamp_stage(work, "closureCor", ["velocityClosureCor.h5", "timeseries_SET.h5",
                                       "timeseriesBias.h5", "smallbaselineApp.cfg",
                                       ".inversion_provenance.json", ".bias_provenance.json"])


def test_current_phase_cache_can_be_reused(tmp_path):
    work, _ = make_work(tmp_path)
    assert wf.prepare_phase(work) == "no_integer_ambiguities"


def test_fresh_phase_prerequisites_are_reported_explicitly(tmp_path):
    work, _ = make_work(tmp_path)
    for name in ("numTriNonzeroIntAmbiguity.h5", "timeseries.h5"):
        (work / name).unlink()
    assert wf.prepare_phase(work) == "fresh"


def test_changed_requested_reference_rejects_cached_phase(tmp_path):
    work, _ = make_work(tmp_path)
    path = work / "fullApp.cfg"
    path.write_text(path.read_text().replace("mintpy.reference.yx = 1, 1", "mintpy.reference.yx = 1, 2"))
    with pytest.raises(wf.ProvenanceError, match="phase/reference/network settings changed"):
        wf.prepare_phase(work)


def test_repeated_clean_resume_keeps_effective_config_and_skips_reinversion(tmp_path):
    """Execute the actual phase-resume block against a MintPy update-check stub.

    MintPy itself is a cluster dependency. The stub models its relevant
    contract: unchanged custom settings do not rewrite the effective config,
    while a newer effective config invalidates the inversion update check.
    """
    work, _ = make_work(tmp_path)
    custom = work / "fullApp.cfg"
    requested_phase_config = custom.read_text()
    effective = work / "smallbaselineApp.cfg"
    effective.write_text(requested_phase_config.replace(
        "mintpy.unwrapError.method = phase_closure", "mintpy.unwrapError.method = no"))
    old_ns = 1_600_000_000_000_000_000
    os.utime(effective, ns=(old_ns, old_ns))
    inversion = work / "inversion_complete"
    inversion.write_text("already current\n")
    before_config = effective.stat().st_mtime_ns
    before_inversion = inversion.stat().st_mtime_ns
    commands = tmp_path / "bin"
    commands.mkdir()
    wrappers = {
        "python": f'#!/bin/sh\nexec "{sys.executable}" "$@"\n',
        # GNU sed is a cluster prerequisite; emulate only this literal
        # substitution to run the production shell block on macOS too.
        "sed": f'''#!{sys.executable}
import pathlib, sys
assert sys.argv[1] == "-i"
assert sys.argv[2] == "s/mintpy.unwrapError.method = phase_closure/mintpy.unwrapError.method = no/"
p = pathlib.Path(sys.argv[3])
p.write_text(p.read_text().replace("mintpy.unwrapError.method = phase_closure", "mintpy.unwrapError.method = no"))
''',
        "smallbaselineApp.py": f'''#!{sys.executable}
import pathlib, sys
args = sys.argv[1:]
work = pathlib.Path(args[args.index("--dir") + 1])
requested = pathlib.Path(args[0]).read_text()
effective = work / "smallbaselineApp.cfg"
with (work / "mintpy_calls.log").open("a") as stream:
    stream.write(" ".join(args) + "\\n")
if effective.read_text() != requested:
    effective.write_text(requested)
if args[args.index("--start") + 1] == "invert_network":
    marker = work / "inversion_complete"
    if effective.stat().st_mtime_ns > marker.stat().st_mtime_ns:
        marker.write_text("reinverted\\n")
''',
    }
    for name, text in wrappers.items():
        executable = commands / name
        executable.write_text(text)
        executable.chmod(0o755)
    source = (ROOT / "mintpy2.sbatch").read_text()
    block = source.split("# BEGIN PHASE_PREREQUISITES", 1)[1].split("# END PHASE_PREREQUISITES", 1)[0]
    script = ('set -uo pipefail\ncd "$W"\nPAIR_ARGS=()\n' + block
              + '\nsmallbaselineApp.py fullApp.cfg --dir "$W" --start invert_network --end velocity\n')
    for _ in range(2):
        # The production script regenerates this requested phase config.
        custom.write_text(requested_phase_config)
        result = subprocess.run(["bash", "-c", script], text=True, capture_output=True,
                                env={**os.environ, "W": str(work),
                                     "PROVENANCE": str(ROOT / "workflow_provenance.py"),
                                     "PATH": str(commands) + os.pathsep + os.environ["PATH"]})
        assert result.returncode == 0, result.stdout + result.stderr
        assert "PHASE_PREREQUISITES_CURRENT: no_integer_ambiguities" in result.stdout
        assert effective.stat().st_mtime_ns == before_config
        assert inversion.stat().st_mtime_ns == before_inversion
    calls = (work / "mintpy_calls.log").read_text().splitlines()
    assert len(calls) == 2
    assert all("--start invert_network" in call for call in calls)
    assert all("modify_network" not in call for call in calls)


def test_changed_requested_network_is_rejected_before_mintpy(tmp_path):
    work, pairs = make_work(tmp_path)
    selection = tmp_path / "pairs.txt"
    selection.write_text("\n".join("_".join(pair) for pair in pairs[:-1]))
    with pytest.raises(wf.ProvenanceError, match="requested network differs"):
        wf.prepare_phase(work, selection)


def test_network_content_change_detected_even_with_restored_timestamp(tmp_path):
    work, _ = make_work(tmp_path)
    path = work / "inputs/ifgramStack.h5"
    st = path.stat()
    with h5py.File(path, "r+") as f:
        f["dropIfgram"][1] = False
    os.utime(path, ns=(st.st_atime_ns, st.st_mtime_ns))
    with pytest.raises(wf.ProvenanceError, match="network/reference"):
        wf.require_phase(work)


@pytest.mark.parametrize("flags", [np.array([1, 1, 2, 1, 1, 1]),
                                    np.array([1, 1, np.nan, 1, 1, 1]),
                                    np.ones((6, 1), dtype=bool)])
def test_malformed_retention_flags_are_rejected(tmp_path, flags):
    work, _ = make_work(tmp_path)
    with h5py.File(work / "inputs/ifgramStack.h5", "r+") as f:
        del f["dropIfgram"]
        f.create_dataset("dropIfgram", data=flags)
    with pytest.raises(wf.ProvenanceError, match="dropIfgram must"):
        wf.stack_state(work)


def test_duplicate_pair_and_wrong_phase_dimensions_are_rejected(tmp_path):
    work, _ = make_work(tmp_path)
    with h5py.File(work / "inputs/ifgramStack.h5", "r+") as f:
        f["date"][1] = f["date"][0]
    with pytest.raises(wf.ProvenanceError, match="duplicate"):
        wf.stack_state(work)


def test_phase_dimensions_must_match_date_pairs(tmp_path):
    work, _ = make_work(tmp_path)
    with h5py.File(work / "inputs/ifgramStack.h5", "r+") as f:
        del f["unwrapPhase"]
        f.create_dataset("unwrapPhase", data=np.ones((5, 3, 3)))
    with pytest.raises(wf.ProvenanceError, match="dimensions"):
        wf.stack_state(work)


def test_changed_reference_rejects_corrected_output(tmp_path):
    work, _ = make_work(tmp_path)
    corrected_product(work)
    with h5py.File(work / "inputs/ifgramStack.h5", "r+") as f:
        f.attrs["REF_X"] = "2"
    with pytest.raises(wf.ProvenanceError):
        wf.select_velocity(work)


def test_unknown_legacy_cache_is_not_silently_adopted(tmp_path):
    work, _ = make_work(tmp_path)
    (work / ".phase_provenance.json").unlink()
    with pytest.raises(wf.ProvenanceError, match="legacy evidence is unverified"):
        wf.prepare_phase(work)


def test_missing_pair_list_and_partial_unknown_selection_fail(tmp_path):
    work, pairs = make_work(tmp_path)
    with pytest.raises(OSError):
        wf.network_selection(work, tmp_path / "missing.txt")
    selection = tmp_path / "pairs.txt"
    selection.write_text("_".join(pairs[0]) + "\n19990101_19990113\n")
    with pytest.raises(wf.ProvenanceError, match="absent from stack"):
        wf.network_selection(work, selection)


def test_keep_all_pair_list_is_valid_and_emits_no_exclusions(tmp_path, capsys):
    work, pairs = make_work(tmp_path)
    selection = tmp_path / "pairs.txt"
    selection.write_text("# full selected network\n" + "\n".join("_".join(p) for p in pairs))
    assert wf.main(["selection", str(work), "--pairs", str(selection)]) == 0
    assert capsys.readouterr().out == "no\n"


@pytest.mark.parametrize("bad_counts", [np.full((3, 3), np.nan), -np.ones((3, 3)), np.full((3, 3), 0.5)])
def test_invalid_ambiguity_evidence_is_unrepaired(tmp_path, bad_counts):
    work, _ = make_work(tmp_path)
    with h5py.File(work / "numTriNonzeroIntAmbiguity.h5", "r+") as f:
        f["mask"][:] = bad_counts
    assert wf.record_phase(work)["status"] == "unrepaired"


def test_missing_ambiguity_evidence_is_unrepaired(tmp_path):
    work, _ = make_work(tmp_path)
    (work / "numTriNonzeroIntAmbiguity.h5").unlink()
    assert wf.record_phase(work)["status"] == "unrepaired"


def test_unrepaired_stage_cannot_be_compared_or_applied_by_default(tmp_path):
    work, _ = make_work(tmp_path, ambiguous=True)
    with pytest.raises(wf.ProvenanceError):
        wf.select_velocity(work)
    with pytest.raises(wf.ProvenanceError):
        wf.apply_inputs(work)
    product, record = wf.select_velocity(work, allow_unverified=True)
    assert product.name == "velocity.h5"
    assert record["status"] == "UNVERIFIED_DIAGNOSTIC"


def test_valid_corrected_product_is_selected(tmp_path):
    work, _ = make_work(tmp_path)
    corrected_product(work)
    assert wf.select_velocity(work)[0].name == "velocityClosureCor.h5"


def test_missing_bias_rejects_corrected_but_verified_inversion_is_available(tmp_path):
    work, _ = make_work(tmp_path)
    corrected_product(work)
    (work / "timeseriesBias.h5").unlink()
    assert wf.select_velocity(work)[0].name == "velocity.h5"
    with pytest.raises(OSError):
        wf.apply_inputs(work)


def test_changed_raw_timeseries_invalidates_corrected_product_recursively(tmp_path):
    work, _ = make_work(tmp_path)
    corrected_product(work)
    with h5py.File(work / "timeseries.h5", "r+") as f:
        f["timeseries"][1, 1, 1] = 0.5
    with pytest.raises(wf.ProvenanceError):
        wf.select_velocity(work)


def test_incomplete_attempt_blocks_old_success(tmp_path):
    work, _ = make_work(tmp_path)
    corrected_product(work)
    (work / ".workflow_incomplete").write_text("later failed job\n")
    with pytest.raises(wf.ProvenanceError, match="latest workflow attempt"):
        wf.select_velocity(work)
    product, record = wf.select_velocity(work, allow_unverified=True)
    assert product.name == "velocity.h5"
    assert record["status"] == "UNVERIFIED_DIAGNOSTIC"


def test_zheng_checks_retained_bandwidth_separately_from_source_loops(tmp_path):
    work, _ = make_work(tmp_path)
    wf.zheng_network(work, 3, 3)
    with h5py.File(work / "inputs/ifgramStack.h5", "r+") as f:
        # Retain bandwidth 2; conn-3 source pair stays available for nl=3.
        f["dropIfgram"][2] = False
    wf.zheng_network(work, 2, 3)
    with pytest.raises(wf.ProvenanceError, match="bandwidth-3"):
        wf.zheng_network(work, 3, 3)


def test_apply_requires_exact_dates_and_same_mask(tmp_path):
    work, _ = make_work(tmp_path)
    assert wf.apply_inputs(work) == ("timeseries_SET.h5", "linear", "maskTempCoh.h5")
    with h5py.File(work / "timeseriesBias.h5", "r+") as f:
        f["date"][1] = b"20200114"
    wf.stamp_stage(work, "bias", ["timeseriesBias.h5"])
    with pytest.raises(wf.ProvenanceError, match="acquisition dates"):
        wf.apply_inputs(work)


def test_set_enabled_does_not_fall_back_to_plain_timeseries(tmp_path):
    work, _ = make_work(tmp_path)
    (work / "timeseries_SET.h5").unlink()
    with pytest.raises(OSError):
        wf.apply_inputs(work)


@pytest.mark.parametrize("unrepaired,override,expected", [(1, 0, 3), (1, 1, 0), (0, 0, 0)])
def test_unrepaired_gate_stops_before_costly_stage(unrepaired, override, expected):
    source = (ROOT / "mintpy2.sbatch").read_text()
    block = source.split("# BEGIN UNREPAIRED_STOP", 1)[1].split("# END UNREPAIRED_STOP", 1)[0]
    result = subprocess.run(["bash", "-c", block + '\necho INVERSION_REACHED\n'],
                            text=True, capture_output=True,
                            env={**os.environ, "PC_UNREPAIRED": str(unrepaired), "ALLOW_UNREPAIRED": str(override)})
    assert result.returncode == expected
    assert ("INVERSION_REACHED" in result.stdout) == (expected == 0)


def test_standalone_apply_passes_the_original_ramp_mask(tmp_path):
    work, _ = make_work(tmp_path)
    commands = tmp_path / "bin"
    commands.mkdir()
    (commands / "python").write_text(f'#!/bin/sh\nexec "{sys.executable}" "$@"\n')
    (commands / "python").chmod(0o755)
    stub = '''#!{python}
import pathlib, shutil, sys
args = sys.argv[1:]
with open({log!r}, "a") as stream:
    stream.write(pathlib.Path(sys.argv[0]).name + " " + " ".join(args) + "\\n")
shutil.copyfile(args[0], args[args.index("-o") + 1])
'''.format(python=sys.executable, log=str(work / "commands.log"))
    for name in ("diff.py", "remove_ramp.py", "timeseries2velocity.py"):
        path = commands / name
        path.write_text(stub)
        path.chmod(0o755)
    result = subprocess.run(["bash", str(ROOT / "apply_zheng.sh")], text=True, capture_output=True,
                            env={**os.environ, "W": str(work), "PROVENANCE": str(ROOT / "workflow_provenance.py"),
                                 "PATH": str(commands) + os.pathsep + os.environ["PATH"]})
    assert result.returncode == 0, result.stdout + result.stderr
    assert "-s linear -m maskTempCoh.h5" in (work / "commands.log").read_text()
    assert not (work / ".workflow_incomplete").exists()
    assert wf.select_velocity(work)[0].name == "velocityClosureCor.h5"
