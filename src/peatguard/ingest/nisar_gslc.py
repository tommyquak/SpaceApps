"""Search and AOI-stage NISAR GSLC products from ASF Earthdata Cloud.

Full-frame 40 MHz GSLC products over PeatGuard are approximately 20 GiB each.
The complex rasters are HDF5 chunked, so authenticated HTTP range access can
copy only the chunks intersecting the operational AOI. The resulting local
subset preserves the GSLC paths consumed by :mod:`peatguard.insar.nisar_gslc`.

Subset writes are restartable. A ``.part`` HDF5 records the number of rows
completed and is atomically renamed only after all selected datasets validate.
"""

from __future__ import annotations

import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable, Optional

import asf_search as asf
import numpy as np
from pyproj import Transformer

from peatguard.ingest.download import authenticate_earthdata

logger = logging.getLogger(__name__)

DEFAULT_GSLC_BBOX = (114.21, -3.0, 114.83, -2.42)
DEFAULT_PATH = 125
DEFAULT_FRAME = 176
DEFAULT_MODE_TOKEN = "4005"
COMPLEX_CACHE_BLOCK_BYTES = 32 * 1024 * 1024
MASK_CACHE_BLOCK_BYTES = 1 * 1024 * 1024


def search_gslc_products(
    bbox_wgs84: tuple[float, float, float, float] = DEFAULT_GSLC_BBOX,
    start_date: str = "2026-06-01",
    end_date: Optional[str] = None,
    path_number: int = DEFAULT_PATH,
    frame_number: int = DEFAULT_FRAME,
    mode_token: str = DEFAULT_MODE_TOKEN,
    max_results: int = 500,
) -> list[asf.ASFProduct]:
    """Find full-cover 40 MHz NISAR GSLC products for the PeatGuard AOI."""
    west, south, east, north = bbox_wgs84
    wkt = (
        f"POLYGON(({west} {south},{east} {south},{east} {north},"
        f"{west} {north},{west} {south}))"
    )
    products = asf.search(
        platform="NISAR",
        processingLevel="GSLC",
        intersectsWith=wkt,
        start=start_date,
        end=end_date,
        maxResults=max_results,
    )
    selected = []
    for product in products:
        properties = product.properties
        filename = properties.get("fileName", "")
        if properties.get("pathNumber") != path_number:
            continue
        if properties.get("frameNumber") != frame_number:
            continue
        if f"_{mode_token}_" not in filename:
            continue
        selected.append(product)
    return sorted(selected, key=lambda item: item.properties.get("startTime", ""))


def _copy_attrs(source: object, destination: object) -> None:
    """Copy HDF5 attributes except object-reference dimension metadata."""
    for key, value in source.attrs.items():
        if key in {"DIMENSION_LIST", "REFERENCE_LIST"}:
            continue
        destination.attrs[key] = value


def _copy_small_group(source: object, destination: object) -> None:
    """Recursively copy metadata groups without spatial raster datasets."""
    import h5py

    _copy_attrs(source, destination)
    for name, source_item in source.items():
        if isinstance(source_item, h5py.Group):
            child = destination.require_group(name)
            _copy_small_group(source_item, child)
            continue
        data = source_item[()]
        target = destination.create_dataset(name, data=data, dtype=source_item.dtype)
        _copy_attrs(source_item, target)


def _window_from_coordinates(
    x_coordinates: np.ndarray,
    y_coordinates: np.ndarray,
    epsg: int,
    bbox_wgs84: tuple[float, float, float, float],
) -> tuple[int, int, int, int]:
    """Convert a WGS84 bbox to an inclusive GSLC coordinate-vector window."""
    west, south, east, north = bbox_wgs84
    if west >= east or south >= north:
        raise ValueError("bbox_wgs84 must be (west, south, east, north)")

    transformer = Transformer.from_crs(4326, epsg, always_xy=True)
    min_x, min_y, max_x, max_y = transformer.transform_bounds(
        west, south, east, north, densify_pts=21
    )
    col_start = int(np.searchsorted(x_coordinates, min_x, side="left"))
    col_stop = int(np.searchsorted(x_coordinates, max_x, side="right"))
    negative_y = -y_coordinates
    row_start = int(np.searchsorted(negative_y, -max_y, side="left"))
    row_stop = int(np.searchsorted(negative_y, -min_y, side="right"))

    row_start = max(0, min(len(y_coordinates), row_start))
    row_stop = max(0, min(len(y_coordinates), row_stop))
    col_start = max(0, min(len(x_coordinates), col_start))
    col_stop = max(0, min(len(x_coordinates), col_stop))
    if row_stop <= row_start or col_stop <= col_start:
        raise ValueError("Requested AOI does not intersect the GSLC grid")
    return row_start, row_stop, col_start, col_stop


def _initialize_subset(
    source: object,
    part_path: Path,
    source_name: str,
    bbox_wgs84: tuple[float, float, float, float],
    frequency: str,
    polarization: str,
) -> tuple[int, int, int, int]:
    """Create the restartable subset HDF5 and return the source window."""
    import h5py

    grid_path = f"science/LSAR/GSLC/grids/{frequency}"
    source_grid = source[grid_path]
    x_coordinates = source_grid["xCoordinates"][:].astype(np.float64)
    y_coordinates = source_grid["yCoordinates"][:].astype(np.float64)
    epsg = int(source_grid["projection"][()])
    row_start, row_stop, col_start, col_stop = _window_from_coordinates(
        x_coordinates, y_coordinates, epsg, bbox_wgs84
    )
    output_shape = (row_stop - row_start, col_stop - col_start)

    with h5py.File(part_path, "w") as destination:
        destination.attrs["peatguard_subset_version"] = "1"
        destination.attrs["peatguard_source_filename"] = source_name
        destination.attrs["peatguard_bbox_wgs84"] = np.asarray(
            bbox_wgs84, dtype=np.float64
        )
        destination.attrs["peatguard_source_window"] = np.asarray(
            (row_start, row_stop, col_start, col_stop), dtype=np.int64
        )
        destination.attrs["peatguard_rows_completed"] = 0
        destination.attrs["peatguard_complete"] = False

        identification_path = "science/LSAR/identification"
        destination_identification = destination.require_group(identification_path)
        _copy_small_group(source[identification_path], destination_identification)

        destination_grid = destination.require_group(grid_path)
        _copy_attrs(source_grid, destination_grid)
        for name in (
            "azimuthBandwidth",
            "centerFrequency",
            "listOfPolarizations",
            "numberOfSubSwaths",
            "projection",
            "rangeBandwidth",
            "slantRangeSpacing",
            "xCoordinateSpacing",
            "yCoordinateSpacing",
            "zeroDopplerTimeSpacing",
        ):
            if name not in source_grid:
                continue
            source_dataset = source_grid[name]
            target = destination_grid.create_dataset(
                name, data=source_dataset[()], dtype=source_dataset.dtype
            )
            _copy_attrs(source_dataset, target)

        x_target = destination_grid.create_dataset(
            "xCoordinates", data=x_coordinates[col_start:col_stop]
        )
        y_target = destination_grid.create_dataset(
            "yCoordinates", data=y_coordinates[row_start:row_stop]
        )
        _copy_attrs(source_grid["xCoordinates"], x_target)
        _copy_attrs(source_grid["yCoordinates"], y_target)

        source_complex = source_grid[polarization]
        chunks = (
            min(source_complex.chunks[0] if source_complex.chunks else 512, output_shape[0]),
            min(source_complex.chunks[1] if source_complex.chunks else 512, output_shape[1]),
        )
        complex_target = destination_grid.create_dataset(
            polarization,
            shape=output_shape,
            dtype=source_complex.dtype,
            chunks=chunks,
            compression="gzip",
            compression_opts=4,
        )
        _copy_attrs(source_complex, complex_target)

        # Preserve mission quality masks when they share the GSLC raster grid.
        for name in ("mask", "inputDataExceptionMask"):
            if name not in source_grid:
                continue
            source_mask = source_grid[name]
            if source_mask.shape != source_complex.shape:
                continue
            mask_target = destination_grid.create_dataset(
                name,
                shape=output_shape,
                dtype=source_mask.dtype,
                chunks=chunks,
                compression="gzip",
                compression_opts=4,
            )
            _copy_attrs(source_mask, mask_target)

    return row_start, row_stop, col_start, col_stop


def _dataset_resume_row(
    attrs: object,
    dataset_name: str,
    polarization: str,
) -> int:
    """Return a dataset checkpoint without leaking the legacy HH row to masks."""
    checkpoint_key = f"peatguard_rows_completed__{dataset_name}"
    if checkpoint_key in attrs:
        return int(attrs[checkpoint_key])
    if dataset_name == polarization:
        return int(attrs.get("peatguard_rows_completed", 0))
    return 0


def subset_open_gslc(
    source: object,
    output_path: Path,
    source_name: str,
    bbox_wgs84: tuple[float, float, float, float] = DEFAULT_GSLC_BBOX,
    frequency: str = "frequencyA",
    polarization: str = "HH",
    block_rows: int = 512,
) -> Path:
    """Copy an AOI subset from an open GSLC HDF5 object with restart support."""
    import h5py

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    part_path = output_path.with_suffix(output_path.suffix + ".part")

    if output_path.exists():
        with h5py.File(output_path, "r") as existing:
            complete = bool(existing.attrs.get("peatguard_complete", False))
            existing_source = str(existing.attrs.get("peatguard_source_filename", ""))
        if complete and existing_source == source_name:
            logger.info("GSLC AOI subset already complete: %s", output_path)
            return output_path
        raise RuntimeError(f"Existing GSLC subset is not valid for {source_name}: {output_path}")

    if not part_path.exists():
        _initialize_subset(
            source,
            part_path,
            source_name,
            bbox_wgs84,
            frequency,
            polarization,
        )

    grid_path = f"science/LSAR/GSLC/grids/{frequency}"
    with h5py.File(part_path, "r+") as destination:
        if str(destination.attrs["peatguard_source_filename"]) != source_name:
            raise RuntimeError(f"Partial subset belongs to another source: {part_path}")
        expected_bbox = np.asarray(destination.attrs["peatguard_bbox_wgs84"])
        if not np.allclose(expected_bbox, np.asarray(bbox_wgs84), atol=1e-10):
            raise RuntimeError(f"Partial subset uses another AOI: {part_path}")

        row_start, row_stop, col_start, col_stop = (
            int(value)
            for value in destination.attrs["peatguard_source_window"]
        )
        rows_completed = int(destination.attrs.get("peatguard_rows_completed", 0))
        output_rows = row_stop - row_start
        source_grid = source[grid_path]
        destination_grid = destination[grid_path]
        raster_names = [
            name
            for name in (polarization, "mask", "inputDataExceptionMask")
            if name in destination_grid
        ]

        logger.info(
            "Staging %s: %dx%d pixels, legacy checkpoint row %d/%d",
            source_name,
            col_stop - col_start,
            output_rows,
            rows_completed,
            output_rows,
        )
        # Read one source dataset to completion before seeking to the next.
        # NISAR's large chunked datasets occupy distant regions of the HDF5.
        # Alternating HH -> mask -> exception mask for every row block defeats
        # HTTP readahead and repeats expensive range-request seeks.
        for name in raster_names:
            checkpoint_key = f"peatguard_rows_completed__{name}"
            dataset_rows_completed = _dataset_resume_row(
                destination.attrs,
                name,
                polarization,
            )
            logger.info(
                "  dataset %s: resuming row %d/%d",
                name,
                dataset_rows_completed,
                output_rows,
            )
            for output_row in range(
                dataset_rows_completed, output_rows, block_rows
            ):
                output_row_stop = min(output_rows, output_row + block_rows)
                source_slice = np.s_[
                    row_start + output_row : row_start + output_row_stop,
                    col_start:col_stop,
                ]
                output_slice = np.s_[output_row:output_row_stop, :]
                destination_grid[name][output_slice] = source_grid[name][source_slice]
                destination.attrs[checkpoint_key] = output_row_stop
                destination.flush()
                logger.info(
                    "  %s rows %d/%d (%.1f%%)",
                    name,
                    output_row_stop,
                    output_rows,
                    output_row_stop / output_rows * 100.0,
                )

        destination.attrs["peatguard_rows_completed"] = output_rows
        destination.attrs["peatguard_complete"] = True
        destination.flush()

    part_path.replace(output_path)
    logger.info("GSLC AOI subset complete: %s", output_path)
    return output_path


def _temporary_signed_url(
    product: asf.ASFProduct,
    session: asf.ASFSession,
) -> str:
    """Resolve an authenticated ASF URL to a short-lived cloud URL."""
    source_url = product.properties.get("url")
    if not source_url:
        raise ValueError(f"ASF product has no URL: {product.properties.get('fileName')}")
    response = session.get(
        source_url,
        headers={"Range": "bytes=0-0"},
        stream=True,
        timeout=(30, 60),
    )
    response.raise_for_status()
    signed_url = response.url
    response.close()
    return signed_url


def _is_expired_signed_url_error(error: BaseException) -> bool:
    """Return whether a range read failed because its temporary URL expired."""
    status = getattr(error, "status", None)
    if status in {401, 403}:
        return True
    message = str(error).lower()
    return (
        ("401" in message or "403" in message)
        and ("unauthorized" in message or "forbidden" in message)
    )


def _remote_cache_settings(
    output_path: Path,
    polarization: str,
) -> tuple[int, str, dict[str, int]]:
    """Choose a remote cache based on the next unfinished raster.

    Large sequential readahead is efficient for complex GSLC chunks. Once the
    complex raster is checkpointed, the remaining uint8 masks are tiny but
    their compressed HDF5 chunks are interleaved with metadata across a
    roughly 20 GiB object. A smaller multi-block cache avoids repeatedly
    fetching tens of megabytes around each mask chunk.
    """
    import h5py

    part_path = Path(output_path).with_suffix(output_path.suffix + ".part")
    if part_path.exists():
        try:
            with h5py.File(part_path, "r") as partial:
                row_start, row_stop, _, _ = (
                    int(value)
                    for value in partial.attrs["peatguard_source_window"]
                )
                output_rows = row_stop - row_start
                primary_rows = int(
                    partial.attrs.get(
                        f"peatguard_rows_completed__{polarization}", 0
                    )
                )
            if primary_rows >= output_rows:
                return (
                    MASK_CACHE_BLOCK_BYTES,
                    "blockcache",
                    {"maxblocks": 128},
                )
        except (KeyError, OSError, ValueError):
            logger.warning(
                "Could not inspect the partial GSLC cache checkpoint: %s",
                part_path,
            )
    return COMPLEX_CACHE_BLOCK_BYTES, "readahead", {}


def stage_gslc_product(
    product: asf.ASFProduct,
    output_dir: Path,
    bbox_wgs84: tuple[float, float, float, float] = DEFAULT_GSLC_BBOX,
    session: Optional[asf.ASFSession] = None,
    frequency: str = "frequencyA",
    polarization: str = "HH",
    block_rows: int = 512,
    max_url_refreshes: int = 3,
) -> Path:
    """Range-read and stage one ASF NISAR GSLC product to a local AOI HDF5."""
    import certifi
    import fsspec
    import h5py

    if session is None:
        session = authenticate_earthdata()
    filename = product.properties.get("fileName")
    if not filename:
        raise ValueError("ASF product is missing fileName")
    output_path = Path(output_dir) / filename.replace(".h5", "_peatguard_aoi.h5")
    if output_path.exists():
        # Avoid resolving a signed URL when a validated subset already exists.
        with h5py.File(output_path, "r") as existing:
            if bool(existing.attrs.get("peatguard_complete", False)):
                return output_path

    # aiohttp does not always find macOS system roots under standalone Python.
    os.environ.setdefault("SSL_CERT_FILE", certifi.where())
    if max_url_refreshes < 0:
        raise ValueError("max_url_refreshes must be non-negative")
    for refresh_count in range(max_url_refreshes + 1):
        signed_url = _temporary_signed_url(product, session)
        cache_block_size, cache_type, cache_options = (
            _remote_cache_settings(output_path, polarization)
        )
        logger.info(
            "Opening remote GSLC with authenticated range requests: %s "
            "(cache=%s, block=%d bytes)",
            filename,
            cache_type,
            cache_block_size,
        )
        try:
            with fsspec.open(
                signed_url,
                "rb",
                block_size=cache_block_size,
                cache_type=cache_type,
                cache_options=cache_options,
                expand=False,
            ) as remote_file:
                with h5py.File(remote_file, "r") as source:
                    return subset_open_gslc(
                        source=source,
                        output_path=output_path,
                        source_name=filename,
                        bbox_wgs84=bbox_wgs84,
                        frequency=frequency,
                        polarization=polarization,
                        block_rows=block_rows,
                    )
        except Exception as error:
            if (
                not _is_expired_signed_url_error(error)
                or refresh_count >= max_url_refreshes
            ):
                raise
            logger.warning(
                "Temporary GSLC URL expired for %s; refreshing and resuming "
                "from the HDF5 row checkpoint (%d/%d)",
                filename,
                refresh_count + 1,
                max_url_refreshes,
            )
    raise RuntimeError(f"Failed to stage GSLC after URL refreshes: {filename}")


def stage_gslc_products(
    products: Iterable[asf.ASFProduct],
    output_dir: Path,
    bbox_wgs84: tuple[float, float, float, float] = DEFAULT_GSLC_BBOX,
    workers: int = 1,
    **kwargs: object,
) -> list[Path]:
    """Stage multiple GSLC products with restartable optional concurrency."""
    product_list = list(products)
    if workers < 1:
        raise ValueError("workers must be at least 1")
    if workers > 1:
        outputs_by_index: dict[int, Path] = {}
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(
                    stage_gslc_product,
                    product,
                    output_dir=output_dir,
                    bbox_wgs84=bbox_wgs84,
                    session=None,
                    **kwargs,
                ): index
                for index, product in enumerate(product_list)
            }
            for future in as_completed(futures):
                index = futures[future]
                outputs_by_index[index] = future.result()
        return [
            outputs_by_index[index] for index in range(len(product_list))
        ]

    session = authenticate_earthdata()
    outputs = []
    for product in product_list:
        outputs.append(
            stage_gslc_product(
                product,
                output_dir=output_dir,
                bbox_wgs84=bbox_wgs84,
                session=session,
                **kwargs,
            )
        )
    return outputs
