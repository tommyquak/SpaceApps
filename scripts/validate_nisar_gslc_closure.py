"""Validate phase closure for a three-date NISAR GSLC network."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from peatguard.insar.nisar_validation import audit_three_date_phase_closure


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase-ab", required=True, type=Path)
    parser.add_argument("--coherence-ab", required=True, type=Path)
    parser.add_argument("--phase-bc", required=True, type=Path)
    parser.add_argument("--coherence-bc", required=True, type=Path)
    parser.add_argument("--phase-ac", required=True, type=Path)
    parser.add_argument("--coherence-ac", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--min-coherence", type=float, default=0.25)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    report = audit_three_date_phase_closure(
        phase_ab_path=args.phase_ab,
        coherence_ab_path=args.coherence_ab,
        phase_bc_path=args.phase_bc,
        coherence_bc_path=args.coherence_bc,
        phase_ac_path=args.phase_ac,
        coherence_ac_path=args.coherence_ac,
        output_json=args.output_json,
        min_coherence=args.min_coherence,
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
