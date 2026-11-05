"""Export native-resolution QA layers from NISAR L2 GUNW products."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from peatguard.insar.nisar_gunw import DEFAULT_BBOX, export_gunw_layers


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Crop and export transparent 20 m and 80 m NISAR GUNW QA layers. "
            "No atmospheric correction or velocity conversion is performed."
        )
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--bbox",
        nargs=4,
        type=float,
        default=DEFAULT_BBOX,
        metavar=("WEST", "SOUTH", "EAST", "NORTH"),
    )
    parser.add_argument(
        "--layers",
        nargs="+",
        default=None,
        help="Optional layer keys. By default all supported QA layers are exported.",
    )
    parser.add_argument("--block-size", type=int, default=512)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    outputs = export_gunw_layers(
        path=args.input,
        output_dir=args.output_dir,
        bbox_wgs84=tuple(args.bbox),
        layers=args.layers,
        block_size=args.block_size,
    )
    for key, path in outputs.items():
        print(f"{key}: {path}")


if __name__ == "__main__":
    main()
