"""Validate a fine GSLC interferogram against the matching mission GUNW."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from peatguard.insar.nisar_validation import compare_fine_to_mission


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fine-phase", required=True, type=Path)
    parser.add_argument("--fine-coherence", required=True, type=Path)
    parser.add_argument("--mission-phase", required=True, type=Path)
    parser.add_argument("--mission-coherence", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--min-coherence", type=float, default=0.25)
    parser.add_argument("--min-phasor-strength", type=float, default=0.25)
    parser.add_argument(
        "--aggregation", choices=("unit_phasor", "complex_weighted"), default="unit_phasor",
        help="Legacy unit-phasor diagnostic or amplitude-preserving complex average",
    )
    parser.add_argument(
        "--fine-cross-product", type=Path,
        help="Two-band GSLC cross-product export, required only for complex_weighted",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    report = compare_fine_to_mission(
        fine_phase_path=args.fine_phase,
        fine_coherence_path=args.fine_coherence,
        mission_phase_path=args.mission_phase,
        mission_coherence_path=args.mission_coherence,
        output_json=args.output_json,
        min_coherence=args.min_coherence,
        min_phasor_strength=args.min_phasor_strength,
        aggregation=args.aggregation,
        fine_cross_product_path=args.fine_cross_product,
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
