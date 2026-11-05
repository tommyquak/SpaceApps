"""MODIS MCD64A1 burned-area ingestion.

Downloads monthly MCD64A1 v061 burned-area granules from NASA LP DAAC via
earthaccess (Earthdata credentials read from ~/.netrc). MCD64A1 is the
standard global burned-area product: 500 m resolution, one granule per
calendar month, with each pixel tagged by the day-of-year it burned that
month (0 = unburned, -1 = unmapped, -2 = water).

Burned-area delineation is preferred over active-fire point detections for
a fire-history layer because it captures actual burn extent rather than
just overpass-time thermal anomalies.

Downstream: analysis/fire_history.py turns these granules into the
fire_history.tif risk layer.

Usage:
    python -m peatguard.ingest.fire \
        --bbox 114.21 -2.71 114.53 -2.41 \
        --start-year 2015 --end-year 2024 \
        --out Datasets/fire_mcd64a1
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("fire_ingest")

# Sebangau Block C bounding box (W, S, E, N) -- matches the canal_risk COG.
SEBANGAU_BBOX = (114.21, -2.71, 114.53, -2.41)


def download_mcd64a1(
    bbox: tuple[float, float, float, float],
    start_year: int,
    end_year: int,
    out_dir: Path,
) -> list[str]:
    """Download every MCD64A1 v061 granule intersecting bbox over the year
    range. Returns the list of local file paths."""
    import earthaccess

    auth = earthaccess.login(strategy="netrc")
    if not auth.authenticated:
        raise RuntimeError(
            "NASA Earthdata authentication failed -- check ~/.netrc for a "
            "urs.earthdata.nasa.gov entry."
        )

    results = earthaccess.search_data(
        short_name="MCD64A1",
        version="061",
        bounding_box=bbox,
        temporal=(f"{start_year}-01-01", f"{end_year}-12-31"),
        count=4000,
    )
    log.info(
        "MCD64A1 v061: %d granules over bbox %s, %d-%d",
        len(results), bbox, start_year, end_year,
    )
    if not results:
        raise RuntimeError("no MCD64A1 granules found for the requested AOI")

    out_dir.mkdir(parents=True, exist_ok=True)
    files = earthaccess.download(results, str(out_dir))
    log.info("downloaded %d files to %s", len(files), out_dir)
    return files


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bbox", nargs=4, type=float, default=list(SEBANGAU_BBOX),
                   metavar=("W", "S", "E", "N"))
    p.add_argument("--start-year", type=int, default=2015)
    p.add_argument("--end-year", type=int, default=2024)
    p.add_argument("--out", type=Path, required=True,
                   help="output directory for downloaded HDF granules")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    download_mcd64a1(
        tuple(args.bbox), args.start_year, args.end_year, args.out,
    )


if __name__ == "__main__":
    main()
