"""Smoke installed science interfaces from outside the source checkout."""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    root = args.root.resolve()
    probe = """
import importlib.resources
import pathlib
import sys
import peatguard
from peatguard.timeseries.local_contrast import extract_local_contrast, solve_contrast_network
from peatguard.timeseries.nisar_local_contrast import audit_local_contrasts
from peatguard.insar.nisar_validation import compare_fine_to_mission
import peatguard.dashboard.app
source = pathlib.Path(sys.argv[1]).resolve()
installed = pathlib.Path(peatguard.__file__).resolve()
assert not installed.is_relative_to(source), f'Source checkout imported instead of wheel: {installed}'
resources = importlib.resources.files('peatguard')
assert resources.joinpath('data/default.yaml').is_file(), 'Packaged default configuration missing'
assert resources.joinpath('dashboard/templates/index.html').is_file(), 'Packaged dashboard missing'
print(f'Installed package smoke passed: {installed}')
"""
    with tempfile.TemporaryDirectory(prefix="peatguard-wheel-smoke-") as directory:
        # -I ignores PYTHONPATH, the current directory and the user site. These
        # commands must resolve the installed distribution, not the source tree.
        subprocess.run(
            [sys.executable, "-I", "-c", probe, str(root / "src")],
            cwd=directory,
            check=True,
        )
        subprocess.run(
            [sys.executable, "-I", "-m", "peatguard.cli", "--help"],
            cwd=directory,
            check=True,
        )
        for filename in (
            "export_nisar_gunw.py",
            "audit_nisar_local_contrasts.py",
            "validate_nisar_gslc_gunw.py",
            "nisar_gslc_interferogram.py",
        ):
            subprocess.run(
                [sys.executable, "-I", str(root / "scripts" / filename), "--help"],
                cwd=directory,
                check=True,
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
