#!/usr/bin/env python3
"""Audit NISAR coherence around PlanetScope-mapped drainage canals."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import rasterio
from rasterio.features import rasterize
from rasterio.warp import Resampling, reproject, transform_geom

from peatguard.analysis.nisar_planet import (
    block_bootstrap_near_far_difference,
    canal_distance_m,
    distance_band_statistics,
)

AOIS = {
    "SBG-NW": {
        "name": "Sebangau North-West",
        "area_ha": 121_200,
        "latitude": -2.565,
        "longitude": 114.38,
    },
    "SBG-NE": {
        "name": "Sebangau North-East",
        "area_ha": 68_900,
        "latitude": -2.52,
        "longitude": 114.69,
    },
    "SBG-CE": {
        "name": "Sebangau Central-East",
        "area_ha": 31_000,
        "latitude": -2.665,
        "longitude": 114.69,
    },
    "SBG-S": {
        "name": "Sebangau South",
        "area_ha": 221_000,
        "latitude": -2.855,
        "longitude": 114.52,
    },
}


def _reproject_canal_mask(
    path: Path,
    shape: tuple[int, int],
    transform: rasterio.Affine,
    crs: rasterio.crs.CRS,
) -> np.ndarray:
    output = np.zeros(shape, dtype=np.uint8)
    with rasterio.open(path) as source:
        reproject(
            rasterio.band(source, 1),
            output,
            src_transform=source.transform,
            src_crs=source.crs,
            src_nodata=source.nodata,
            dst_transform=transform,
            dst_crs=crs,
            dst_nodata=0,
            resampling=Resampling.max,
        )
    return output > 0


def _rasterize_geojson(
    path: Path,
    shape: tuple[int, int],
    transform: rasterio.Affine,
    crs: rasterio.crs.CRS,
) -> np.ndarray:
    collection = json.loads(path.read_text())
    geometries = [
        transform_geom("EPSG:4326", crs, feature["geometry"])
        for feature in collection.get("features", [])
        if feature.get("geometry")
    ]
    return rasterize(
        [(geometry, 1) for geometry in geometries],
        out_shape=shape,
        transform=transform,
        fill=0,
        dtype=np.uint8,
    ).astype(bool)


def _aoi_mask(
    metadata: dict[str, float | str],
    shape: tuple[int, int],
    transform: rasterio.Affine,
    crs: rasterio.crs.CRS,
) -> np.ndarray:
    side_m = max(500.0, np.sqrt(float(metadata["area_ha"]) * 10_000.0))
    latitude = float(metadata["latitude"])
    longitude = float(metadata["longitude"])
    half_lat = side_m / 2.0 / 111_320.0
    half_lon = side_m / 2.0 / (
        111_320.0 * np.cos(np.deg2rad(latitude))
    )
    west = longitude - half_lon
    east = longitude + half_lon
    south = latitude - half_lat
    north = latitude + half_lat
    polygon = {
        "type": "Polygon",
        "coordinates": [
            [
                [west, south],
                [east, south],
                [east, north],
                [west, north],
                [west, south],
            ]
        ],
    }
    projected = transform_geom("EPSG:4326", crs, polygon)
    return rasterize(
        [(projected, 1)],
        out_shape=shape,
        transform=transform,
        fill=0,
        dtype=np.uint8,
    ).astype(bool)


def _bootstrap_or_status(
    coherence: np.ndarray,
    distance: np.ndarray,
    mask: np.ndarray,
    *,
    block_size_pixels: int,
) -> dict[str, object]:
    try:
        return block_bootstrap_near_far_difference(
            coherence,
            distance,
            mask=mask,
            block_size_pixels=block_size_pixels,
        )
    except ValueError as error:
        return {
            "status": "not_available",
            "reason": str(error),
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coherence", type=Path, required=True)
    parser.add_argument("--canal-mask", type=Path, required=True)
    parser.add_argument("--peat-geojson", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with rasterio.open(args.coherence) as source:
        coherence = source.read(1).astype(np.float64)
        raw_valid = np.isfinite(coherence)
        if source.nodata is not None:
            raw_valid &= coherence != source.nodata
        scale = float(source.scales[0])
        offset = float(source.offsets[0])
        coherence[raw_valid] = coherence[raw_valid] * scale + offset
        coherence[~raw_valid] = np.nan
        coherence[(coherence < 0.0) | (coherence > 1.0)] = np.nan
        shape = source.shape
        transform = source.transform
        crs = source.crs
        spacing = (abs(transform.e), abs(transform.a))
        bounds = tuple(source.bounds)
        coherence_tags = source.tags()
    if crs is None:
        raise ValueError("NISAR coherence raster has no CRS")

    canals = _reproject_canal_mask(
        args.canal_mask, shape, transform, crs
    )
    peat = _rasterize_geojson(
        args.peat_geojson, shape, transform, crs
    )
    distance = canal_distance_m(canals, spacing)
    valid = np.isfinite(coherence)
    # Preserve the earlier 4 km spatial block width when moving from the
    # mission 20 m product to the 5 m-posted fine product.
    block_size_pixels = max(1, round(4_000.0 / max(spacing)))

    distance_bands: dict[str, list[dict[str, object]]] = {}
    spatial_sensitivity: dict[str, dict[str, object]] = {}

    def analyze_stratum(name: str, mask: np.ndarray) -> None:
        if not mask.any():
            return
        distance_bands[name] = [
            statistic.to_dict()
            for statistic in distance_band_statistics(
                coherence,
                distance,
                mask=mask,
            )
        ]
        if mask.sum() >= 1_000:
            spatial_sensitivity[name] = _bootstrap_or_status(
                coherence,
                distance,
                mask,
                block_size_pixels=block_size_pixels,
            )

    analyze_stratum("all_valid", valid)
    analyze_stratum("mapped_peat", valid & peat)
    analyze_stratum("outside_mapped_peat", valid & ~peat)
    for code, metadata in AOIS.items():
        aoi = _aoi_mask(metadata, shape, transform, crs)
        analyze_stratum(f"{code}_all", valid & aoi)
        analyze_stratum(f"{code}_peat", valid & aoi & peat)

    report = {
        "status": "coherence_observability_not_canal_causality",
        "source": {
            "coherence": str(args.coherence),
            "canal_mask": str(args.canal_mask),
            "peat_extent": str(args.peat_geojson),
            "coherence_posting_m": list(spacing),
            "coherence_bounds": bounds,
            "coherence_scale": scale,
            "coherence_offset": offset,
            "effective_support_m": coherence_tags.get(
                "effective_support_m"
            ),
        },
        "coverage": {
            "valid_coherence_pixels": int(valid.sum()),
            "canal_pixels_on_coherence_grid": int(canals.sum()),
            "mapped_peat_pixels": int(peat.sum()),
        },
        "distance_bands": distance_bands,
        "spatial_sensitivity": spatial_sensitivity,
        "spatial_bootstrap": {
            "target_block_width_m": 4_000,
            "block_size_pixels": block_size_pixels,
        },
        "interpretation": {
            "allowed": [
                "Use PlanetScope for canal geometry and continuity.",
                "Use NISAR coherence to map where future radar deformation is observable.",
                "Use low coherence to trigger manual review or alternative validation.",
            ],
            "not_allowed": [
                "Do not treat coherence as proof of subsidence or drainage impact.",
                "Do not use provisional NISAR displacement as Planet training truth.",
                "Do not alter restoration priority from this short commissioning stack.",
            ],
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
