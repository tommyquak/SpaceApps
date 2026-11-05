"""Audit PeatGuard restoration-priority sensitivity to project assumptions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject

from peatguard.analysis.risk_sensitivity import audit_priority_sensitivity


def _pixel_size_m(transform: rasterio.Affine, height: int) -> tuple[float, float]:
    center_latitude = transform.f + (height / 2.0) * transform.e
    y_m = abs(transform.e) * 110_574.0
    x_m = abs(transform.a) * 111_320.0 * np.cos(np.radians(center_latitude))
    return float(y_m), float(x_m)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--canal", required=True, type=Path)
    parser.add_argument("--subsidence", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    args = parser.parse_args()

    with rasterio.open(args.canal) as canal_source:
        canal = canal_source.read(1) > 0
        shape = canal_source.shape
        transform = canal_source.transform
        crs = canal_source.crs

    with rasterio.open(args.subsidence) as subsidence_source:
        subsidence = np.full(shape, np.nan, dtype=np.float32)
        reproject(
            source=rasterio.band(subsidence_source, 1),
            destination=subsidence,
            src_transform=subsidence_source.transform,
            src_crs=subsidence_source.crs,
            src_nodata=subsidence_source.nodata,
            dst_transform=transform,
            dst_crs=crs,
            dst_nodata=np.nan,
            resampling=Resampling.bilinear,
        )

    report = audit_priority_sensitivity(
        canal,
        subsidence,
        pixel_size_m=_pixel_size_m(transform, shape[0]),
    )
    report["inputs"] = {
        "canal": str(args.canal),
        "subsidence": str(args.subsidence),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
