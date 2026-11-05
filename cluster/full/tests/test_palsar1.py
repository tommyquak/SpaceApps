"""Static checks for the ALOS-1 PALSAR historical-baseline campaign."""
import csv
import pathlib
import re
import subprocess

HERE = pathlib.Path(__file__).resolve().parent.parent
P1 = HERE / "palsar1"


def test_manifest_has_both_frames_and_https_urls():
    rows = list(csv.DictReader(open(P1 / "manifest_p422_l10.csv")))
    frames = {r["frame"] for r in rows}
    assert frames == {"7120", "7130", "7140"}
    for fr in frames:
        dates = {r["date"] for r in rows if r["frame"] == fr}
        assert len(dates) == 14, fr
    assert all(r["url"].startswith("https://") and r["url"].endswith("-L1.0.zip") for r in rows)
    assert all(float(r["size_mb"]) > 100 for r in rows)
    assert {r["mode"] for r in rows} == {"FBD", "FBS"}


def test_sbatch_scripts_parse():
    for name in ("p1_dl.sbatch", "p1_stack.sbatch", "p1_mintpy.sbatch"):
        subprocess.run(["bash", "-n", str(P1 / name)], check=True)


def test_cfg_placeholders_are_rendered_by_the_job():
    cfg = (P1 / "smallbaselineApp_palsar1.cfg").read_text()
    job = (P1 / "p1_mintpy.sbatch").read_text()
    for ph in set(re.findall(r"__[A-Z]+__", cfg)):
        assert ph in job, ph
    assert "excludeDate" not in cfg            # all 14 dates stay in the network


def test_validate_module_imports():
    subprocess.run([__import__("sys").executable, "-m", "py_compile", str(P1 / "p1_validate.py")], check=True)


def test_stack_uses_the_real_stripmapstack_flags():
    job = (P1 / "p1_stack.sbatch").read_text()
    assert "-a $ALKS" in job and "-z $ALKS" not in job   # -z is the zero-doppler switch, -a is azimuth looks
    assert 'extra="$extra -m"' in job                    # unpackFrame_ALOS_raw.py -m for multi-frame dates


def test_stack_loop_skips_markers_and_logs():
    job = (P1 / "p1_stack.sbatch").read_text()
    assert "grep -vE '\\.(job|log|DONE)$' | sort -V" in job   # run.py must never be fed a .log or .DONE file


def test_stack_cleans_intermediates_and_probes_quota():
    job = (P1 / "p1_stack.sbatch").read_text()
    assert "rm -rf $STACK/offsets $STACK/coregSLC" in job
    assert '[ "$(basename $d)" = "$REF" ] && continue' in job   # the reference focused SLC feeds the igram step
    assert "scratch quota exhausted" in job
    assert 'XB="-x $BBOX"' in job


def test_mosaic_module_compiles():
    subprocess.run([__import__("sys").executable, "-m", "py_compile", str(P1 / "p1_mosaic.py")], check=True)


def test_frames_list_accepts_colon_separator():
    for name in ("p1_dl.sbatch", "p1_stack.sbatch"):
        assert 'tr ":+" ",,"' in (P1 / name).read_text()   # --export=ALL,FRAMES=a,b would truncate to a
