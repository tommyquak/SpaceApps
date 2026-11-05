"""Audit independent baseline products used by the NISAR integration."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from peatguard.analysis.baseline_audit import (
    TAY_EXPECTED_MD5,
    audit_tay_long_term_rate,
    write_audit_json,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Verify PeatGuard's Tay long-term-rate derivative against the "
            "official source without calibrating either product."
        )
    )
    parser.add_argument("--tay-source", required=True, type=Path)
    parser.add_argument("--tay-derivative", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--expected-md5", default=TAY_EXPECTED_MD5)
    parser.add_argument("--scale", type=float, default=1_000.0)
    parser.add_argument("--absolute-limit-mm-yr", type=float, default=100.0)
    args = parser.parse_args()

    report = audit_tay_long_term_rate(
        args.tay_source,
        args.tay_derivative,
        derivative_scale=args.scale,
        derivative_absolute_limit=args.absolute_limit_mm_yr,
        expected_source_md5=args.expected_md5,
    )
    output = write_audit_json(report, args.output_json)
    print(json.dumps({"output": str(output), "audit_passed": report["audit_passed"]}))
    if not report["audit_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
