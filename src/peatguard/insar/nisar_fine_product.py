"""Publish a compact, quality-only product from fine NISAR GSLC coherence.

The source rasters are floating-point, 5 m-posted local coherence products.
They are large enough that combining them must be block-wise. The public
product stores coherence as UInt8 with a scale of 1/254 and reserves 255 for
NoData. This retains approximately 0.004 coherence precision while reducing
the final COG substantially.

This module does not process phase or displacement. Its output is explicitly a
radar observability layer.
"""

from __future__ import annotations

import hashlib
import math
import os
import tempfile
from contextlib import ExitStack
from pathlib import Path
from typing import Sequence

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.shutil import copy as raster_copy

COHERENCE_SCALE = 1.0 / 254.0
COHERENCE_NODATA = 255


def quantize_coherence(data: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Encode 0..1 coherence as 0..254 and reserve 255 for NoData."""
    values = np.asarray(data)
    mask = np.asarray(valid, dtype=bool)
    if values.shape != mask.shape:
        raise ValueError("Coherence values and validity mask differ")
    output = np.full(values.shape, COHERENCE_NODATA, dtype=np.uint8)
    usable = mask & np.isfinite(values) & (values >= 0.0) & (values <= 1.0)
    output[usable] = np.rint(values[usable] * 254.0).astype(np.uint8)
    return output


def _same_grid(
    first: rasterio.DatasetReader,
    second: rasterio.DatasetReader,
) -> bool:
    return (
        first.crs == second.crs
        and first.width == second.width
        and first.height == second.height
        and first.transform.almost_equals(second.transform)
    )


def _valid(data: np.ndarray, nodata: float | None) -> np.ndarray:
    output = np.isfinite(data) & (data >= 0.0) & (data <= 1.0)
    if nodata is not None and np.isfinite(nodata):
        output &= data != nodata
    return output


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(16 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def build_minimum_coherence_cog(
    coherence_paths: Sequence[Path],
    working_tif: Path,
    output_cog: Path,
    *,
    extra_tags: dict[str, str] | None = None,
) -> dict[str, object]:
    """Build a scaled UInt8 minimum-pair coherence COG block-wise.

    Args:
        coherence_paths: Two or more same-grid 0..1 coherence rasters.
        working_tif: Large temporary tiled raster. Keep this under the
            deletable-work root.
        output_cog: Compact final COG.
        extra_tags: Additional provenance tags.

    Returns:
        Reproducible output metadata and coverage statistics.
    """
    paths = [Path(path) for path in coherence_paths]
    if len(paths) < 2:
        raise ValueError("At least two coherence pairs are required")
    if len(set(paths)) != len(paths):
        raise ValueError("Coherence input paths must be unique")

    working_tif = Path(working_tif)
    output_cog = Path(output_cog)
    working_tif.parent.mkdir(parents=True, exist_ok=True)
    output_cog.parent.mkdir(parents=True, exist_ok=True)

    thresholds = (0.25, 0.35, 0.50, 0.70)
    valid_pixels = 0
    value_sum = 0.0
    minimum = math.inf
    maximum = -math.inf
    threshold_counts = {threshold: 0 for threshold in thresholds}

    with ExitStack() as stack:
        sources = [
            stack.enter_context(rasterio.open(path)) for path in paths
        ]
        reference = sources[0]
        if reference.crs is None:
            raise ValueError("Fine coherence raster has no CRS")
        grid_report = {
            "crs": str(reference.crs),
            "width": reference.width,
            "height": reference.height,
            "transform": tuple(reference.transform),
            "pixel_posting_m": 5,
            "effective_support_m": 25,
        }
        for path, source in zip(paths[1:], sources[1:]):
            if not _same_grid(reference, source):
                raise ValueError(f"Fine coherence grid mismatch: {path}")

        profile = reference.profile.copy()
        profile.update(
            driver="GTiff",
            count=1,
            dtype="uint8",
            nodata=COHERENCE_NODATA,
            tiled=True,
            blockxsize=512,
            blockysize=512,
            compress="DEFLATE",
            predictor=2,
            BIGTIFF="IF_SAFER",
        )
        with rasterio.open(working_tif, "w", **profile) as destination:
            destination.scales = (COHERENCE_SCALE,)
            destination.offsets = (0.0,)
            destination.units = ("1",)
            destination.set_band_description(
                1, "minimum_pair_local_interferometric_coherence"
            )
            destination.update_tags(
                sensor="NISAR L-SAR",
                source_product="NISAR L2 GSLC",
                product_maturity="PROVISIONAL",
                processing_status="supporting_quality_only_not_displacement",
                measurement="minimum_pair_local_interferometric_coherence",
                pair_count=str(len(paths)),
                pixel_posting_m="5",
                effective_support_m="25",
                coherence_window_pixels="5",
                quantization="uint8_round_coherence_times_254",
                quantization_scale=str(COHERENCE_SCALE),
                quantization_nodata=str(COHERENCE_NODATA),
                not_displacement="true",
                not_subsidence="true",
                not_velocity="true",
                not_ground_truth="true",
                **(extra_tags or {}),
            )

            for _, window in destination.block_windows(1):
                blocks = [
                    source.read(1, window=window).astype(np.float32)
                    for source in sources
                ]
                common = np.ones(blocks[0].shape, dtype=bool)
                for block, source in zip(blocks, sources):
                    common &= _valid(block, source.nodata)
                minimum_block = np.full(
                    blocks[0].shape, np.nan, dtype=np.float32
                )
                if np.any(common):
                    stack_values = np.stack(blocks)
                    values = np.min(stack_values[:, common], axis=0)
                    minimum_block[common] = values
                    count = int(values.size)
                    valid_pixels += count
                    value_sum += float(np.sum(values, dtype=np.float64))
                    minimum = min(minimum, float(np.min(values)))
                    maximum = max(maximum, float(np.max(values)))
                    for threshold in thresholds:
                        threshold_counts[threshold] += int(
                            np.count_nonzero(values >= threshold)
                        )
                destination.write(
                    quantize_coherence(minimum_block, common),
                    1,
                    window=window,
                )

    if valid_pixels == 0:
        raise RuntimeError("No pixels are valid in every coherence pair")

    # Convert the large work raster into a formally ordered COG. Write to a
    # sibling temporary path and replace atomically only after GDAL succeeds.
    with tempfile.NamedTemporaryFile(
        prefix=f".{output_cog.stem}.",
        suffix=".tif",
        dir=output_cog.parent,
        delete=False,
    ) as temporary:
        temporary_cog = Path(temporary.name)
    temporary_cog.unlink()
    try:
        raster_copy(
            working_tif,
            temporary_cog,
            driver="COG",
            compress="DEFLATE",
            blocksize=512,
            overview_resampling=Resampling.average.name,
            overview_count=4,
            bigtiff="IF_SAFER",
        )
        os.replace(temporary_cog, output_cog)
    finally:
        if temporary_cog.exists():
            temporary_cog.unlink()

    with rasterio.open(output_cog) as result:
        image_structure = result.tags(ns="IMAGE_STRUCTURE")
        if result.dtypes != ("uint8",):
            raise RuntimeError("Fine coherence COG has unexpected dtype")
        if result.nodata != COHERENCE_NODATA:
            raise RuntimeError("Fine coherence COG has unexpected NoData")
        if not math.isclose(result.scales[0], COHERENCE_SCALE):
            raise RuntimeError("Fine coherence COG lost its scale metadata")
        if image_structure.get("LAYOUT") != "COG":
            raise RuntimeError("GDAL did not produce a COG-ordered output")
        overview_levels = result.overviews(1)

    return {
        "output": str(output_cog),
        "output_bytes": output_cog.stat().st_size,
        "output_sha256": _sha256(output_cog),
        "working_tif": str(working_tif),
        "source_paths": [str(path) for path in paths],
        "grid": grid_report,
        "encoding": {
            "dtype": "uint8",
            "nodata": COHERENCE_NODATA,
            "scale": COHERENCE_SCALE,
            "decode": "coherence = raw_value * scale",
            "overview_levels": overview_levels,
            "layout": "COG",
        },
        "coverage": {
            "valid_pixels": valid_pixels,
            "mean": value_sum / valid_pixels,
            "minimum": minimum,
            "maximum": maximum,
            "fraction_at_or_above": {
                str(threshold): threshold_counts[threshold] / valid_pixels
                for threshold in thresholds
            },
        },
    }
