"""Adversarial checks for NISAR chain selection and dated publication isolation."""

import importlib.util
import os
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np
import pytest

FULL = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("nisar_safety", FULL / "analysis/nisar_safety.py")
safety = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(safety)


def pair(d1, d2, **changes):
    return dict(dict(path=f"{d1}_{d2}.h5", track=125, frame=176, ad="A",
                     crid="P05023", d1=d1, d2=d2), **changes)


def test_identity_uses_hdf5_and_rejects_beta(tmp_path):
    path = tmp_path / "arbitrary_filename.h5"
    with h5py.File(path, "w") as source:
        ident = source.create_group("science/LSAR/identification")
        for key, value in dict(trackNumber=125, frameNumber=176,
                               orbitPassDirection=b"Ascending", compositeReleaseId=b"P05023",
                               referenceZeroDopplerStartTime=b"2026-06-22T22:07:10",
                               secondaryZeroDopplerStartTime=b"2026-07-04T22:07:10").items():
            ident[key] = value
    assert safety.read_pair_identity(str(path))["d1"] == "20260622"
    with h5py.File(path, "a") as source:
        source["science/LSAR/identification/compositeReleaseId"][()] = b"X05010"
    with pytest.raises(ValueError, match="provisional"):
        safety.read_pair_identity(str(path))


def test_chain_rejects_mixed_groups_and_honors_selection():
    one = [pair("20260622", "20260704"), pair("20260704", "20260716")]
    two = [dict(record, frame=177) for record in one]
    with pytest.raises(ValueError, match="one track/frame/CRID"):
        safety.select_chain(one + two, "A")
    assert safety.select_chain(one + two, "A", frame="176") == one
    with pytest.raises(ValueError, match="one track/frame/CRID"):
        safety.select_chain(one + [dict(one[1], crid="P05024")], "A")


def test_chain_rejects_duplicates_and_does_not_silently_drop_gap():
    chain = [pair("20260622", "20260704"), pair("20260704", "20260716")]
    with pytest.raises(ValueError, match="duplicate"):
        safety.select_chain(chain + [chain[0]], "A")
    with pytest.raises(ValueError, match="gap or branching"):
        safety.select_chain(chain + [pair("20260809", "20260821")], "A")


def grid():
    values = np.ones((30, 30))
    return dict(coh=values, unw=values, cc=values, ion=values, ionu=values,
                x=np.arange(30) * 80., y=np.arange(30) * -80., epsg=32750, wl=0.242)


@pytest.mark.parametrize("change", [dict(x=np.arange(30) * 80. + 40), dict(epsg=32749),
                                      dict(wl=0.2384), dict(ionu=None)])
def test_grids_reject_mixed_support_or_missing_uncertainty(change):
    with pytest.raises(ValueError):
        safety.validate_grids([grid(), dict(grid(), **change)])


def test_deramp_never_fits_target_when_stable_mask_is_insufficient():
    phase = np.ones((30, 30))
    stable = np.zeros(phase.shape, bool)
    stable[:10] = True
    with pytest.raises(ValueError, match="Insufficient stable-ground"):
        safety.deramp(phase, stable)


def test_deramp_preserves_target_anomaly():
    yy, xx = np.mgrid[:40, :40]
    phase = 4 + 0.2 * xx + 0.1 * yy + 0.01 * xx * yy
    stable = np.ones(phase.shape, bool)
    stable[15:25, 15:25] = False
    phase[~stable] += 20
    result = safety.deramp(phase, stable)
    np.testing.assert_allclose(result[stable], 0, atol=1e-10)
    np.testing.assert_allclose(result[~stable], 20, atol=1e-10)


def run_analysis(tmp_path, fail_stage=None, omit_direction=None):
    """Execute the actual runner function with tiny local stage stand-ins."""
    text = (FULL / "nisar_pull.sbatch").read_text()
    function = "run_analysis() {" + text.split("run_analysis() {", 1)[1].split("\n}\n", 1)[0] + "\n}\n"
    out = tmp_path / "fresh"
    out.mkdir()
    # Old global artifacts must never satisfy checks for this fresh attempt.
    for direction in ("A", "D"):
        (tmp_path / f"nisar_{direction}_disp60_hp10km_mm.tif").write_bytes(b"old")
    for stage in ("stage1", "stage2c", "stage2d"):
        code = "import os, pathlib, sys\np = pathlib.Path(os.environ['NISAR_OUT'])\n"
        if stage == fail_stage:
            code += "sys.exit(7)\n"
        elif stage == "stage2c":
            for direction in ("A", "D"):
                if direction == omit_direction:
                    continue
                for suffix in ("disp60_hp10km_mm.tif", "chain_coh.tif", "sources.json"):
                    code += f"(p / 'nisar_{direction}_{suffix}').write_bytes(b'new')\n"
        elif stage == "stage2d":
            code += "(p / 'stage2d_ran').touch()\n"
        (tmp_path / f"nisar_{stage}.py").write_text(code)
    env = dict(os.environ, INSAR_PYTHON=sys.executable, OUTD=str(out), NISAR_OUT=str(out),
               STAMP="20260910", HP_KM="10.0")
    result = subprocess.run(["bash", "-c", function + "run_analysis"], cwd=tmp_path,
                            env=env, capture_output=True, text=True)
    return result, out


@pytest.mark.parametrize("failed", ["stage1", "stage2c", "stage2d"])
def test_runner_rejects_failed_stage_without_success_marker(tmp_path, failed):
    result, out = run_analysis(tmp_path, fail_stage=failed)
    assert result.returncode == 7, result.stderr
    assert not (out / "PROVISIONAL.txt").exists()
    if failed != "stage2d":
        assert not (out / "stage2d_ran").exists()


def test_runner_cannot_relabel_old_global_files(tmp_path):
    result, out = run_analysis(tmp_path, omit_direction="D")
    assert result.returncode != 0
    assert "missing current D" in result.stdout
    assert not (out / "PROVISIONAL.txt").exists()
    assert not (out / "stage2d_ran").exists()


def test_runner_marks_only_complete_new_run(tmp_path):
    result, out = run_analysis(tmp_path)
    assert result.returncode == 0, result.stderr
    assert (out / "PROVISIONAL.txt").exists()
    assert (out / "nisar_D_disp60_hp10km_mm.tif").read_bytes() == b"new"
