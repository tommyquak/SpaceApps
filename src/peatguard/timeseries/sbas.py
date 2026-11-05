"""MintPy SBAS time-series inversion.

Orchestrates the MintPy processing chain for Small BAseline Subset
(SBAS) time-series analysis. This replaces the naive per-pair
displacement summation with a proper inversion that accounts for
atmospheric delays, unwrapping errors, and orbital ramps.

The MintPy chain (smallbaselineApp step order):
    load_data -> modify_network -> reference_point ->
    correct_unwrap_error -> invert_network -> correct_troposphere ->
    deramp -> correct_topography -> velocity
"""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import Optional

from peatguard.config import PeatGuardConfig

logger = logging.getLogger(__name__)

# smallbaselineApp step names in execution order, matching STEP_LIST in
# mintpy.defaults.template (the subset this pipeline runs). Steps are
# dispatched through smallbaselineApp.py --dostep; the standalone per-step
# scripts take the ifgramStack/timeseries file as their positional argument,
# not the template, so they cannot be invoked with the config path.
# Corrections (troposphere, deramp, topography) operate on the inverted
# timeseries.h5 and therefore run after invert_network.
MINTPY_STEPS = [
    ("load_data", "Loading interferogram stack"),
    ("modify_network", "Modifying network (removing low-coherence pairs)"),
    ("reference_point", "Setting reference point"),
    ("quick_overview", "Generating quick overview"),
    ("correct_unwrap_error", "Correcting unwrapping errors"),
    ("invert_network", "Running SBAS network inversion"),
    ("correct_troposphere", "Correcting tropospheric delay (ERA5)"),
    ("deramp", "Removing orbital ramps"),
    ("correct_topography", "Correcting topographic residual"),
    ("residual_RMS", "Computing residual RMS"),
    ("reference_date", "Selecting reference date"),
    ("velocity", "Estimating velocity"),
]

_STEP_NAMES = [name for name, _ in MINTPY_STEPS]


def _run_mintpy_step(
    step: str,
    config_path: Path,
    work_dir: Path,
) -> None:
    """Run a single smallbaselineApp step via --dostep.

    Args:
        step: smallbaselineApp step name (e.g., 'invert_network').
        config_path: Path to the MintPy template file.
        work_dir: MintPy working directory.
    """
    cmd = ["smallbaselineApp.py", str(config_path), "--dostep", step]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(work_dir),
    )

    if result.returncode != 0:
        logger.error(
            "smallbaselineApp step %s failed:\nSTDOUT: %s\nSTDERR: %s",
            step,
            result.stdout[-1000:],
            result.stderr[-1000:],
        )
        raise RuntimeError(f"MintPy step {step} failed with return code {result.returncode}")


def run_sbas_inversion(
    config_path: Path,
    work_dir: Path,
    skip_steps: Optional[list[str]] = None,
    start_step: Optional[str] = None,
) -> Path:
    """Run the MintPy SBAS chain step by step via smallbaselineApp --dostep.

    Args:
        config_path: Path to the MintPy template file.
        work_dir: MintPy working directory.
        skip_steps: List of step names to skip (smallbaselineApp names,
            without a .py extension).
        start_step: Step name to start from (skips earlier steps).

    Returns:
        Path to the velocity file produced by MintPy.

    Raises:
        ValueError: If a step name is not a known smallbaselineApp step.
    """
    if skip_steps is None:
        skip_steps = []

    requested = list(skip_steps) + ([start_step] if start_step else [])
    unknown = [s for s in requested if s not in _STEP_NAMES]
    if unknown:
        raise ValueError(
            f"Unknown MintPy step name(s) {unknown}; valid steps: {_STEP_NAMES}"
        )

    started = start_step is None

    for step, description in MINTPY_STEPS:
        if not started:
            if step == start_step:
                started = True
            else:
                continue

        if step in skip_steps:
            logger.info("Skipping: %s (%s)", step, description)
            continue

        logger.info("Running: %s -- %s", step, description)
        _run_mintpy_step(step, config_path, work_dir)

    velocity_file = work_dir / "velocity.h5"
    if velocity_file.exists():
        logger.info("SBAS inversion complete. Velocity file: %s", velocity_file)
    else:
        logger.warning("Velocity file not found at expected path: %s", velocity_file)

    return velocity_file


def run_smallbaselineApp(
    config_path: Path,
    work_dir: Path,
) -> Path:
    """Run MintPy's smallbaselineApp.py as a single command.

    This is the all-in-one MintPy entry point that runs the entire
    SBAS chain. Use this instead of run_sbas_inversion when you
    don't need per-step control.

    Args:
        config_path: Path to the MintPy template file.
        work_dir: MintPy working directory.

    Returns:
        Path to the velocity file.
    """
    logger.info("Running smallbaselineApp.py with config: %s", config_path)

    cmd = ["smallbaselineApp.py", str(config_path)]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(work_dir),
    )

    if result.returncode != 0:
        logger.error(
            "smallbaselineApp.py failed:\nSTDOUT: %s\nSTDERR: %s",
            result.stdout[-2000:],
            result.stderr[-2000:],
        )
        raise RuntimeError(f"smallbaselineApp.py failed with return code {result.returncode}")

    velocity_file = work_dir / "velocity.h5"
    logger.info("smallbaselineApp complete. Velocity: %s", velocity_file)
    return velocity_file


def run_smallbaselineApp_from_step(
    config_path: Path,
    work_dir: Path,
    start_step: str = "modify_network",
    end_step: Optional[str] = None,
) -> Path:
    """Run smallbaselineApp.py starting from a specific step.

    Use this to skip earlier steps (e.g. load_data) when the HDF5
    input files have already been created manually.

    Args:
        config_path: Path to the MintPy template file.
        work_dir: MintPy working directory.
        start_step: Step name to start from (without .py extension).
        end_step: Optional step name to stop after (inclusive).
            If None, runs through to the end.

    Returns:
        Path to the velocity file.
    """
    if end_step:
        logger.info(
            "Running smallbaselineApp.py --start %s --end %s",
            start_step, end_step,
        )
    else:
        logger.info("Running smallbaselineApp.py --start %s", start_step)

    cmd = ["smallbaselineApp.py", str(config_path), "--start", start_step]
    if end_step:
        cmd.extend(["--end", end_step])

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(work_dir),
    )

    if result.returncode != 0:
        logger.error(
            "smallbaselineApp.py (from %s) failed:\nSTDOUT: %s\nSTDERR: %s",
            start_step,
            result.stdout[-3000:],
            result.stderr[-3000:],
        )
        raise RuntimeError(
            f"smallbaselineApp.py failed (start={start_step}) with rc {result.returncode}"
        )

    velocity_file = work_dir / "velocity.h5"
    logger.info("smallbaselineApp complete (from %s). Velocity: %s", start_step, velocity_file)
    return velocity_file
