#!/usr/bin/env python3
"""Prepare and optionally invert transparent NISAR GUNW phase variants."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

from peatguard.timeseries.nisar_gunw import (
    MINTPY_NISAR_COMMIT,
    assess_temporal_readiness,
    prepare_phase_variant,
    select_reference_patch,
)
from peatguard.timeseries.nisar_commissioning import sha256_file


def _file_record(path: Path) -> dict[str, object]:
    return {
        "sha256": sha256_file(path),
        "bytes": path.stat().st_size,
    }


def _invert(
    executable: Path,
    stack: Path,
    water_mask: Path,
    output_dir: Path,
) -> dict[str, Path]:
    stack = stack.resolve()
    water_mask = water_mask.resolve()
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    matplotlib_cache = output_dir / ".matplotlib"
    general_cache = output_dir / ".cache"
    matplotlib_cache.mkdir(parents=True, exist_ok=True)
    general_cache.mkdir(parents=True, exist_ok=True)
    outputs = {
        "timeseries": (output_dir / "timeseries.h5").resolve(),
        "temporal_coherence": (
            output_dir / "temporalCoherence.h5"
        ).resolve(),
        "num_inverted": (output_dir / "numInvIfgram.h5").resolve(),
    }
    command = [
        str(executable),
        str(stack),
        "--skip-reference",
        "--weight-func",
        "var",
        "--water-mask",
        str(water_mask),
        "--mask-dset",
        "connectComponent",
        "--mask-thres",
        "0.5",
        "--min-redun",
        "1",
        "--output",
        str(outputs["timeseries"]),
        str(outputs["temporal_coherence"]),
        str(outputs["num_inverted"]),
    ]
    environment = os.environ.copy()
    environment["MPLCONFIGDIR"] = str(matplotlib_cache)
    environment["XDG_CACHE_HOME"] = str(general_cache)
    subprocess.run(
        command,
        check=True,
        cwd=output_dir,
        env=environment,
    )
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mintpy-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--peat-geojson", type=Path)
    parser.add_argument("--canal-mask", type=Path)
    parser.add_argument("--invert", action="store_true")
    args = parser.parse_args()

    inputs = args.mintpy_dir / "inputs"
    source_stack = inputs / "ifgramStack.h5"
    geometry = inputs / "geometryGeo.h5"
    water_mask = args.mintpy_dir / "waterMask.h5"
    corrections = {
        "ionosphere": inputs / "ionStack.h5",
        "troposphere": inputs / "tropoStack.h5",
        "solid_earth_tides": inputs / "setStack.h5",
    }
    required = [source_stack, geometry, water_mask, *corrections.values()]
    reference_support = [
        path
        for path in (args.peat_geojson, args.canal_mask)
        if path is not None
    ]
    missing = [str(path) for path in required if not path.exists()]
    missing.extend(
        str(path) for path in reference_support if not path.exists()
    )
    if missing:
        raise FileNotFoundError(f"Missing required NISAR MintPy files: {missing}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    readiness = assess_temporal_readiness(source_stack)
    reference_parameters = {
        "minimum_coherence": 0.70,
        "maximum_local_phase_std_rad": 0.80,
        "minimum_canal_distance_m": 1_000.0,
        "minimum_peat_distance_m": 500.0,
        "minimum_edge_distance_m": 1_000.0,
        "patch_radius_pixels": 2,
        "local_window_pixels": 7,
    }
    reference = select_reference_patch(
        source_stack,
        geometry,
        peat_geojson=args.peat_geojson,
        canal_mask=args.canal_mask,
        **reference_parameters,
    )

    variants_dir = args.output_dir / "variants"
    variants = {
        "raw": prepare_phase_variant(
            source_stack,
            variants_dir / "ifgramStack_raw_referenced.h5",
            reference,
        ),
        "ionosphere": prepare_phase_variant(
            source_stack,
            variants_dir / "ifgramStack_ionosphere_referenced.h5",
            reference,
            correction_stacks={"ionosphere": corrections["ionosphere"]},
        ),
        "mission_full": prepare_phase_variant(
            source_stack,
            variants_dir / "ifgramStack_mission_full_referenced.h5",
            reference,
            correction_stacks=corrections,
        ),
    }

    inversions: dict[str, dict[str, str]] = {}
    if args.invert:
        executable = Path(sys.executable).parent / "ifgram_inversion.py"
        if not executable.exists():
            raise FileNotFoundError(
                f"Run this script in the pinned MintPy environment: {executable}"
            )
        for name, stack in variants.items():
            outputs = _invert(
                executable,
                stack,
                water_mask,
                args.output_dir / "inversions" / name,
            )
            inversions[name] = {
                key: str(path) for key, path in outputs.items()
            }

    manifest = {
        "status": "commissioning_not_velocity",
        "mintpy_commit": MINTPY_NISAR_COMMIT,
        "temporal_readiness": readiness.to_dict(),
        "reference": reference.to_dict(),
        "source_files": {
            str(path.resolve()): _file_record(path)
            for path in [*required, *reference_support]
        },
        "derived_files": {
            str(path.resolve()): _file_record(path)
            for path in [
                *variants.values(),
                *(
                    Path(path)
                    for outputs in inversions.values()
                    for path in outputs.values()
                ),
            ]
        },
        "reference_support": {
            "peat_geojson": (
                str(args.peat_geojson.resolve())
                if args.peat_geojson
                else None
            ),
            "canal_mask": (
                str(args.canal_mask.resolve()) if args.canal_mask else None
            ),
            "selection_parameters": reference_parameters,
        },
        "variants": {name: str(path) for name, path in variants.items()},
        "inversions": inversions,
    }
    manifest_path = args.output_dir / "commissioning_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
