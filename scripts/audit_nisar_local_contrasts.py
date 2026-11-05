#!/usr/bin/env python3
"""Audit native NISAR target-minus-background changes without a fixed station."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from peatguard.timeseries.nisar_local_contrast import audit_local_contrasts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gunw", nargs="+", type=Path, required=True)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target-radius-m", type=float, default=500)
    parser.add_argument("--background-inner-m", type=float, default=1000)
    parser.add_argument("--background-outer-m", type=float, default=2000)
    parser.add_argument("--min-coherence", type=float, default=0.3)
    parser.add_argument("--min-pixels", type=int, default=10)
    parser.add_argument("--min-fraction", type=float, default=0.5)
    parser.add_argument("--polarization", choices=["HH", "VV"], default="HH")
    args = parser.parse_args(argv)
    if args.output.exists() or args.output.resolve() in {p.resolve() for p in [args.sites, *args.gunw]}:
        parser.error("Output must be a new path distinct from every input")
    report = audit_local_contrasts(
        args.gunw, args.sites, target_radius_m=args.target_radius_m,
        background_inner_m=args.background_inner_m, background_outer_m=args.background_outer_m,
        min_coherence=args.min_coherence, min_pixels=args.min_pixels,
        min_fraction=args.min_fraction, polarization=args.polarization,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(f"Relative diagnostic support: {report['available_sites']}/{report['total_sites']} sites")
    print(f"Observed network cycles: {report['network']['cycle_count']}; release_ready=false")
    print(f"Report: {args.output}")
    return 0 if report["available_sites"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
