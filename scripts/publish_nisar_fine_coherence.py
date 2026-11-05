#!/usr/bin/env python3
"""Publish validated 5 m-posted NISAR coherence as a compact COG."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from peatguard.insar.nisar_fine_product import build_minimum_coherence_cog


def _load_quality_gate(report_path: Path) -> dict[str, object]:
    report = json.loads(report_path.read_text())
    gate = report.get("release_gates", {}).get("supporting_coherence", {})
    if gate.get("passed") is not True:
        raise RuntimeError(
            f"Supporting-coherence validation gate failed: {report_path}"
        )
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pair-coherence",
        type=Path,
        nargs="+",
        required=True,
    )
    parser.add_argument(
        "--validation-report",
        type=Path,
        nargs="+",
        required=True,
    )
    parser.add_argument("--closure-report", type=Path, required=True)
    parser.add_argument("--working-tif", type=Path, required=True)
    parser.add_argument("--output-cog", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()

    if len(args.pair_coherence) != len(args.validation_report):
        raise ValueError(
            "Each pair-coherence raster needs one validation report"
        )
    validations = [
        _load_quality_gate(path) for path in args.validation_report
    ]
    validated_inputs = {
        str(Path(report["inputs"]["fine_coherence"]).resolve())
        for report in validations
    }
    supplied_inputs = {
        str(path.resolve()) for path in args.pair_coherence
    }
    if validated_inputs != supplied_inputs:
        raise ValueError(
            "Validation reports do not exactly match pair-coherence inputs"
        )

    closure = json.loads(args.closure_report.read_text())
    if closure.get("gate", {}).get("passed") is not True:
        raise RuntimeError("Three-date phase-closure gate did not pass")

    result = build_minimum_coherence_cog(
        args.pair_coherence,
        args.working_tif,
        args.output_cog,
        extra_tags={
            "validation_status": "coherence_pass_phase_displacement_fail",
            "validation_report_count": str(len(validations)),
            "closure_gate": "passed_internal_integrity_only",
            "acquisition_start": "2026-06-22",
            "acquisition_end": "2026-07-16",
            "span_days": "24",
            "track_number": "125",
            "frame_number": "176",
            "polarization": "HH",
        },
    )
    result["validation"] = {
        "reports": [str(path) for path in args.validation_report],
        "coherence_gates_passed": True,
        "phase_displacement_gates_passed": all(
            report["release_gates"]["displacement_phase"]["passed"]
            for report in validations
        ),
        "closure_report": str(args.closure_report),
        "closure_gate_passed": True,
        "interpretation": (
            "Releasable as supporting radar-observability quality only; "
            "not releasable as displacement, subsidence, or velocity."
        ),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
