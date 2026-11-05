"""ECOSTRESS land-surface-temperature ingestion.

Downloads ECOSTRESS ECO_L2T_LSTE granules from NASA LP DAAC via
earthaccess (Earthdata credentials from ~/.netrc). ECO_L2T_LSTE is the
Level-2 *tiled* LST product: 70 m, Cloud-Optimized GeoTIFFs on the
Military Grid Reference System, one set of single-band rasters per
overpass per tile. Bands include LST (Kelvin), cloud, water and QC.

ECOSTRESS rides the ISS, so overpass local-times drift -- which is what
makes it useful for thermal work, but also means coverage is irregular
and cloud-limited. We pull the dry-season window and composite.

Downstream: analysis/thermal_stress.py turns these into the
thermal_stress.tif risk layer (daytime LST anomaly vs intact peat).

Usage:
    python -m peatguard.ingest.thermal \
        --bbox 114.21 -2.71 114.53 -2.41 \
        --start 2024-06-01 --end 2024-10-31 \
        --out Datasets/ecostress_lste
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("thermal_ingest")

SEBANGAU_BBOX = (114.21, -2.71, 114.53, -2.41)
# ECOSTRESS land-surface-temperature, Level-2 tiled.
SHORT_NAME = "ECO_L2T_LSTE"


def download_ecostress_lst(
    bbox: tuple[float, float, float, float],
    start: str,
    end: str,
    out_dir: Path,
) -> list[str]:
    """Download every ECO_L2T_LSTE granule intersecting bbox over the
    date window. Returns the list of local file paths."""
    import earthaccess

    auth = earthaccess.login(strategy="netrc")
    if not auth.authenticated:
        raise RuntimeError(
            "NASA Earthdata authentication failed -- check ~/.netrc."
        )

    results = earthaccess.search_data(
        short_name=SHORT_NAME,
        bounding_box=bbox,
        temporal=(start, end),
        count=4000,
    )
    log.info("%s: %d granules over bbox %s, %s..%s",
             SHORT_NAME, len(results), bbox, start, end)
    if not results:
        raise RuntimeError("no ECO_L2T_LSTE granules found for the AOI")

    out_dir.mkdir(parents=True, exist_ok=True)
    files = earthaccess.download(results, str(out_dir))
    log.info("downloaded %d files to %s", len(files), out_dir)
    return files


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bbox", nargs=4, type=float, default=list(SEBANGAU_BBOX),
                   metavar=("W", "S", "E", "N"))
    p.add_argument("--start", default="2024-06-01")
    p.add_argument("--end", default="2024-10-31")
    p.add_argument("--out", type=Path, required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    download_ecostress_lst(tuple(args.bbox), args.start, args.end, args.out)


if __name__ == "__main__":
    main()
