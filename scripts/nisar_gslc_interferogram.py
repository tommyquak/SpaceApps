"""Create a fine-resolution wrapped interferogram from two NISAR GSLCs."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from peatguard.insar.nisar_gslc import form_gslc_interferogram


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Cross-multiply two matching NISAR L2 GSLC products and export "
            "block-wise wrapped phase and local coherence GeoTIFFs."
        )
    )
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--secondary", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--bbox",
        nargs=4,
        type=float,
        metavar=("WEST", "SOUTH", "EAST", "NORTH"),
        default=None,
        help="Optional EPSG:4326 crop. PeatGuard operational AOI: 114.21 -3.0 114.83 -2.42",
    )
    parser.add_argument("--frequency", default="frequencyA")
    parser.add_argument("--polarization", default="HH")
    parser.add_argument("--coherence-window", type=int, default=5)
    parser.add_argument("--block-size", type=int, default=1024)
    parser.add_argument(
        "--export-cross-product",
        action="store_true",
        help="Also preserve native complex amplitudes in a two-band GeoTIFF for aggregation tests",
    )
    parser.add_argument(
        "--allow-missing-quality-mask",
        action="store_true",
        help="Diagnostic only: process a GSLC lacking the mission valid-data mask",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    products = form_gslc_interferogram(
        reference_path=args.reference,
        secondary_path=args.secondary,
        output_dir=args.output_dir,
        bbox_wgs84=tuple(args.bbox) if args.bbox else None,
        frequency=args.frequency,
        polarization=args.polarization,
        coherence_window=args.coherence_window,
        block_size=args.block_size,
        require_quality_mask=not args.allow_missing_quality_mask,
        export_cross_product=args.export_cross_product,
    )
    for label, path in products.items():
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
