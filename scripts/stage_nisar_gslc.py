"""Search and stage cropped 40 MHz NISAR GSLC products for PeatGuard."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from peatguard.ingest.nisar_gslc import search_gslc_products, stage_gslc_products


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", default="2026-06-01")
    parser.add_argument("--end", default=None)
    parser.add_argument("--path", type=int, default=125)
    parser.add_argument("--frame", type=int, default=176)
    parser.add_argument("--mode-token", default="4005")
    parser.add_argument(
        "--bbox",
        nargs=4,
        type=float,
        default=(114.21, -3.0, 114.83, -2.42),
        metavar=("WEST", "SOUTH", "EAST", "NORTH"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("grid_products/_deletable_work/nisar/gslc_subset"),
    )
    parser.add_argument("--block-rows", type=int, default=512)
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Concurrent GSLC products to stage; each retains its own checkpoint",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    logging.getLogger("fsspec").setLevel(logging.WARNING)
    bbox = tuple(args.bbox)
    products = search_gslc_products(
        bbox_wgs84=bbox,
        start_date=args.start,
        end_date=args.end,
        path_number=args.path,
        frame_number=args.frame,
        mode_token=args.mode_token,
    )
    if not products:
        raise SystemExit("No matching GSLC products found")
    logging.info("Found %d matching GSLC products", len(products))
    outputs = stage_gslc_products(
        products,
        output_dir=args.output_dir,
        bbox_wgs84=bbox,
        block_rows=args.block_rows,
        workers=args.workers,
    )
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
