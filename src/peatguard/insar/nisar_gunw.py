"""Inspect and export transparent QA layers from NISAR L2 GUNW products.

Mission GUNW products contain a 20 m wrapped complex interferogram and
coherence layer, plus 80 m unwrapped phase, coherence, connected components,
and ionosphere estimates.  The atmospheric and solid-Earth correction cubes
are stored on a separate three-dimensional radar grid and are intentionally
not interpolated here.  PeatGuard uses the upstream MintPy NISAR preparation
workflow for that scientifically sensitive step.

All phase exported by this module remains in radians.  No two-epoch layer is
labelled as velocity, and no correction is silently applied.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.enums import Resampling
from rasterio.transform import Affine
from rasterio.windows import Window

logger = logging.getLogger(__name__)

DEFAULT_BBOX = (114.21, -3.0, 114.83, -2.42)
DEFAULT_NODATA = -9999.0
SPEED_OF_LIGHT_M_S = 299_792_458.0


@dataclass(frozen=True)
class GUNWProductInfo:
    """Identity and release metadata for one NISAR GUNW pair."""

    path: Path
    granule_id: str
    reference_start: str
    secondary_start: str
    track_number: int
    frame_number: int
    orbit_direction: str
    composite_release_id: str
    product_version: str
    product_specification_version: str
    maturity: str
    center_frequency_hz: float
    wavelength_m: float

    @property
    def reference_date(self) -> str:
        """Return the reference date as YYYYMMDD."""
        return _date_token(self.reference_start)

    @property
    def secondary_date(self) -> str:
        """Return the secondary date as YYYYMMDD."""
        return _date_token(self.secondary_start)


@dataclass(frozen=True)
class GUNWGridInfo:
    """Map-grid metadata for one GUNW resolution group."""

    group_path: str
    epsg: int
    x_coordinates: np.ndarray
    y_coordinates: np.ndarray
    x_spacing: float
    y_spacing: float

    @property
    def shape(self) -> tuple[int, int]:
        """Return grid shape as rows, columns."""
        return self.y_coordinates.size, self.x_coordinates.size

    @property
    def transform(self) -> Affine:
        """Return the pixel-corner affine transform for centre coordinates."""
        return Affine(
            self.x_spacing,
            0.0,
            float(self.x_coordinates[0]) - self.x_spacing / 2.0,
            0.0,
            self.y_spacing,
            float(self.y_coordinates[0]) - self.y_spacing / 2.0,
        )


@dataclass(frozen=True)
class _LayerSpec:
    key: str
    resolution_group: str
    dataset_name: str
    filename_suffix: str
    description: str
    units: str
    complex_phase: bool = False
    categorical: bool = False


_LAYER_SPECS = {
    "wrapped_phase_20m": _LayerSpec(
        key="wrapped_phase_20m",
        resolution_group="wrappedInterferogram",
        dataset_name="wrappedInterferogram",
        filename_suffix="wrapped_phase_20m",
        description="mission_wrapped_interferometric_phase",
        units="radians",
        complex_phase=True,
    ),
    "coherence_20m": _LayerSpec(
        key="coherence_20m",
        resolution_group="wrappedInterferogram",
        dataset_name="coherenceMagnitude",
        filename_suffix="coherence_20m",
        description="mission_normalized_interferometric_coherence",
        units="1",
    ),
    "unwrapped_phase_80m": _LayerSpec(
        key="unwrapped_phase_80m",
        resolution_group="unwrappedInterferogram",
        dataset_name="unwrappedPhase",
        filename_suffix="unwrapped_phase_raw_80m",
        description="mission_unwrapped_phase_raw_uncorrected",
        units="radians",
    ),
    "coherence_80m": _LayerSpec(
        key="coherence_80m",
        resolution_group="unwrappedInterferogram",
        dataset_name="coherenceMagnitude",
        filename_suffix="coherence_80m",
        description="mission_normalized_interferometric_coherence",
        units="1",
    ),
    "connected_components_80m": _LayerSpec(
        key="connected_components_80m",
        resolution_group="unwrappedInterferogram",
        dataset_name="connectedComponents",
        filename_suffix="connected_components_80m",
        description="mission_unwrapping_connected_components",
        units="1",
        categorical=True,
    ),
    "ionosphere_phase_80m": _LayerSpec(
        key="ionosphere_phase_80m",
        resolution_group="unwrappedInterferogram",
        dataset_name="ionospherePhaseScreen",
        filename_suffix="ionosphere_phase_screen_80m",
        description="mission_ionosphere_phase_screen_not_applied",
        units="radians",
    ),
    "ionosphere_uncertainty_80m": _LayerSpec(
        key="ionosphere_uncertainty_80m",
        resolution_group="unwrappedInterferogram",
        dataset_name="ionospherePhaseScreenUncertainty",
        filename_suffix="ionosphere_phase_uncertainty_80m",
        description="mission_ionosphere_phase_screen_uncertainty",
        units="radians",
    ),
}


def _decode_scalar(value: object) -> object:
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.generic):
        return value.item()
    return value


def _date_token(value: str) -> str:
    dashed = re.search(r"(\d{4})-(\d{2})-(\d{2})", value)
    if dashed:
        return "".join(dashed.groups())
    compact = re.search(r"(?<!\d)(\d{8})T", value)
    return compact.group(1) if compact else "unknown"


def _maturity_from_granule(granule_id: str, release_id: str) -> str:
    if "_PR_" in granule_id or release_id.startswith("P"):
        return "PROVISIONAL"
    if "_BETA_" in granule_id:
        return "BETA"
    return "UNKNOWN"


def inspect_gunw(
    path: Path,
    frequency: str = "frequencyA",
) -> GUNWProductInfo:
    """Read GUNW identity metadata without loading raster arrays."""
    import h5py

    path = Path(path)
    identification_path = "science/LSAR/identification"
    with h5py.File(path, "r") as h5:
        if identification_path not in h5:
            raise KeyError(f"Missing GUNW identification group in {path}")
        identification = h5[identification_path]

        def required(name: str) -> object:
            if name not in identification:
                raise KeyError(f"Missing GUNW identification dataset: {name}")
            return _decode_scalar(identification[name][()])

        granule_id = str(required("granuleId"))
        release_id = str(required("compositeReleaseId"))
        frequency_path = f"science/LSAR/GUNW/grids/{frequency}"
        if frequency_path not in h5:
            raise KeyError(f"Missing GUNW frequency group: {frequency_path}")
        frequency_group = h5[frequency_path]
        if "centerFrequency" not in frequency_group:
            raise KeyError(
                f"Missing GUNW frequency dataset: {frequency_path}/centerFrequency"
            )
        center_frequency_hz = float(frequency_group["centerFrequency"][()])
        if not np.isfinite(center_frequency_hz) or center_frequency_hz <= 0:
            raise ValueError(
                f"Invalid GUNW center frequency: {center_frequency_hz}"
            )
        return GUNWProductInfo(
            path=path,
            granule_id=granule_id,
            reference_start=str(required("referenceZeroDopplerStartTime")),
            secondary_start=str(required("secondaryZeroDopplerStartTime")),
            track_number=int(required("trackNumber")),
            frame_number=int(required("frameNumber")),
            orbit_direction=str(required("orbitPassDirection")),
            composite_release_id=release_id,
            product_version=str(required("productVersion")),
            product_specification_version=str(
                required("productSpecificationVersion")
            ),
            maturity=_maturity_from_granule(granule_id, release_id),
            center_frequency_hz=center_frequency_hz,
            wavelength_m=SPEED_OF_LIGHT_M_S / center_frequency_hz,
        )


def _regular_spacing(coordinates: np.ndarray, label: str) -> float:
    if coordinates.ndim != 1 or coordinates.size < 2:
        raise ValueError(f"{label} must contain at least two coordinates")
    differences = np.diff(coordinates.astype(np.float64))
    spacing = float(np.median(differences))
    tolerance = max(abs(spacing) * 1e-6, 1e-8)
    if spacing == 0 or not np.allclose(
        differences, spacing, rtol=0.0, atol=tolerance
    ):
        raise ValueError(f"{label} is not regularly spaced")
    return spacing


def inspect_gunw_grid(
    path: Path,
    resolution_group: str,
    frequency: str = "frequencyA",
    polarization: str = "HH",
) -> GUNWGridInfo:
    """Inspect one 20 m or 80 m GUNW grid."""
    import h5py

    group_path = (
        f"science/LSAR/GUNW/grids/{frequency}/"
        f"{resolution_group}/{polarization}"
    )
    with h5py.File(path, "r") as h5:
        if group_path not in h5:
            raise KeyError(f"Missing GUNW grid group: {group_path}")
        group = h5[group_path]
        for name in ("projection", "xCoordinates", "yCoordinates"):
            if name not in group:
                raise KeyError(f"Missing GUNW grid dataset: {group_path}/{name}")
        x_coordinates = group["xCoordinates"][:].astype(np.float64)
        y_coordinates = group["yCoordinates"][:].astype(np.float64)
        epsg = int(group["projection"][()])

    x_spacing = _regular_spacing(x_coordinates, "xCoordinates")
    y_spacing = _regular_spacing(y_coordinates, "yCoordinates")
    if x_spacing <= 0 or y_spacing >= 0:
        raise ValueError(
            f"Unexpected GUNW coordinate direction: dx={x_spacing}, dy={y_spacing}"
        )
    return GUNWGridInfo(
        group_path=group_path,
        epsg=epsg,
        x_coordinates=x_coordinates,
        y_coordinates=y_coordinates,
        x_spacing=x_spacing,
        y_spacing=y_spacing,
    )


def _window_for_bbox(
    grid: GUNWGridInfo,
    bbox_wgs84: Optional[tuple[float, float, float, float]],
) -> tuple[int, int, int, int]:
    height, width = grid.shape
    if bbox_wgs84 is None:
        return 0, height, 0, width

    west, south, east, north = bbox_wgs84
    if west >= east or south >= north:
        raise ValueError("bbox_wgs84 must be (west, south, east, north)")
    transformer = Transformer.from_crs(4326, grid.epsg, always_xy=True)
    min_x, min_y, max_x, max_y = transformer.transform_bounds(
        west, south, east, north, densify_pts=21
    )

    col_start = int(np.searchsorted(grid.x_coordinates, min_x, side="left"))
    col_stop = int(np.searchsorted(grid.x_coordinates, max_x, side="right"))
    negative_y = -grid.y_coordinates
    row_start = int(np.searchsorted(negative_y, -max_y, side="left"))
    row_stop = int(np.searchsorted(negative_y, -min_y, side="right"))

    row_start = max(0, min(height, row_start))
    row_stop = max(0, min(height, row_stop))
    col_start = max(0, min(width, col_start))
    col_stop = max(0, min(width, col_stop))
    if row_stop <= row_start or col_stop <= col_start:
        raise ValueError("Requested AOI does not intersect the GUNW grid")
    return row_start, row_stop, col_start, col_stop


def _dataset_path(
    spec: _LayerSpec,
    frequency: str,
    polarization: str,
) -> str:
    return (
        f"science/LSAR/GUNW/grids/{frequency}/{spec.resolution_group}/"
        f"{polarization}/{spec.dataset_name}"
    )


def _build_overviews(path: Path, categorical: bool, wrapped_phase: bool = False) -> None:
    # Wrapped phase must never be averaged across the -pi/+pi branch cut:
    # averaging +3.1 and -3.1 rad fabricates ~0 instead of ~pi. Use nearest
    # for categorical AND phase layers, matching insar/nisar_gslc.py.
    resampling = (
        Resampling.nearest if categorical or wrapped_phase else Resampling.average
    )
    with rasterio.open(path, "r+") as dataset:
        minimum_dimension = min(dataset.width, dataset.height)
        levels = [
            level for level in (2, 4, 8, 16) if minimum_dimension // level >= 1
        ]
        if levels:
            dataset.build_overviews(levels, resampling)
            dataset.update_tags(ns="rio_overview", resampling=resampling.name)


def export_gunw_layers(
    path: Path,
    output_dir: Path,
    bbox_wgs84: Optional[tuple[float, float, float, float]] = DEFAULT_BBOX,
    layers: Optional[list[str]] = None,
    frequency: str = "frequencyA",
    polarization: str = "HH",
    block_size: int = 512,
    nodata: float = DEFAULT_NODATA,
) -> dict[str, Path]:
    """Export cropped GUNW QA layers as tiled GeoTIFFs.

    The outputs retain the native 20 m or 80 m projected grids. Atmospheric
    correction cubes are excluded because they require DEM-aware 3-D
    interpolation; use MintPy's current NISAR preparation workflow for those.
    """
    import h5py

    if block_size < 1:
        raise ValueError("block_size must be positive")
    selected = layers or list(_LAYER_SPECS)
    unknown = sorted(set(selected) - set(_LAYER_SPECS))
    if unknown:
        raise ValueError(f"Unknown GUNW layers: {unknown}")

    path = Path(path)
    info = inspect_gunw(path, frequency=frequency)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"nisar_gunw_{info.reference_date}_{info.secondary_date}"
    outputs: dict[str, Path] = {}
    grid_cache: dict[str, GUNWGridInfo] = {}

    with h5py.File(path, "r") as h5:
        for key in selected:
            spec = _LAYER_SPECS[key]
            grid = grid_cache.setdefault(
                spec.resolution_group,
                inspect_gunw_grid(
                    path,
                    spec.resolution_group,
                    frequency=frequency,
                    polarization=polarization,
                ),
            )
            dataset_path = _dataset_path(spec, frequency, polarization)
            if dataset_path not in h5:
                raise KeyError(f"Missing GUNW layer: {dataset_path}")
            source = h5[dataset_path]
            if source.shape != grid.shape:
                raise ValueError(
                    f"GUNW layer shape {source.shape} does not match "
                    f"{spec.resolution_group} grid {grid.shape}"
                )

            row_start, row_stop, col_start, col_stop = _window_for_bbox(
                grid, bbox_wgs84
            )
            height = row_stop - row_start
            width = col_stop - col_start
            transform = grid.transform * Affine.translation(col_start, row_start)
            output_path = output_dir / f"{prefix}_{spec.filename_suffix}.tif"

            source_fill = source.attrs.get("_FillValue")
            if spec.categorical:
                dtype = str(source.dtype)
                output_nodata = int(
                    source_fill
                    if source_fill is not None
                    else np.iinfo(source.dtype).max
                )
                predictor = 2
            else:
                dtype = "float32"
                output_nodata = nodata
                predictor = 3

            profile = {
                "driver": "GTiff",
                "dtype": dtype,
                "width": width,
                "height": height,
                "count": 1,
                "crs": f"EPSG:{grid.epsg}",
                "transform": transform,
                "nodata": output_nodata,
                "tiled": True,
                "blockxsize": min(512, max(16, (width // 16) * 16)),
                "blockysize": min(512, max(16, (height // 16) * 16)),
                "compress": "DEFLATE",
                "predictor": predictor,
                "BIGTIFF": "IF_SAFER",
            }
            tags = {
                "sensor": "NISAR L-SAR",
                "source_product": "L2 GUNW",
                "source_granule": info.granule_id,
                "reference_date": info.reference_start,
                "secondary_date": info.secondary_start,
                "track_number": str(info.track_number),
                "frame_number": str(info.frame_number),
                "orbit_direction": info.orbit_direction,
                "composite_release_id": info.composite_release_id,
                "product_version": info.product_version,
                "product_specification_version": info.product_specification_version,
                "product_maturity": info.maturity,
                "center_frequency_hz": str(info.center_frequency_hz),
                "wavelength_m": str(info.wavelength_m),
                "wavelength_source": "source_product_centerFrequency",
                "pixel_posting_x": str(abs(grid.x_spacing)),
                "pixel_posting_y": str(abs(grid.y_spacing)),
                "processing_status": spec.description,
                "corrections_applied": "none",
                "not_ground_truth": "true",
                "units": spec.units,
            }

            logger.info(
                "Exporting %s: %dx%d pixels at %.1f m",
                key,
                width,
                height,
                abs(grid.x_spacing),
            )
            with rasterio.open(output_path, "w", **profile) as destination:
                destination.update_tags(**tags)
                destination.set_band_description(1, spec.description)
                for output_row in range(0, height, block_size):
                    output_row_stop = min(height, output_row + block_size)
                    source_row = row_start + output_row
                    source_row_stop = row_start + output_row_stop
                    for output_col in range(0, width, block_size):
                        output_col_stop = min(width, output_col + block_size)
                        source_col = col_start + output_col
                        source_col_stop = col_start + output_col_stop
                        data = source[
                            source_row:source_row_stop,
                            source_col:source_col_stop,
                        ]

                        if spec.complex_phase:
                            # np.angle(0+0j) is 0.0, indistinguishable from a
                            # genuine zero-phase sample: exclude zero-magnitude
                            # and declared-fill samples from validity.
                            valid = (
                                np.isfinite(data.real)
                                & np.isfinite(data.imag)
                                & (np.abs(data) > 0)
                            )
                            if source_fill is not None and np.isfinite(
                                complex(source_fill).real
                            ):
                                valid &= data != source_fill
                            block = np.full(data.shape, nodata, dtype=np.float32)
                            block[valid] = np.angle(data[valid]).astype(np.float32)
                        elif spec.categorical:
                            block = np.asarray(data, dtype=source.dtype)
                        else:
                            block = np.asarray(data, dtype=np.float32)
                            valid = np.isfinite(block)
                            if source_fill is not None and np.isfinite(source_fill):
                                valid &= block != source_fill
                            block = np.where(valid, block, nodata).astype(
                                np.float32, copy=False
                            )

                        destination.write(
                            block,
                            1,
                            window=Window(
                                output_col,
                                output_row,
                                output_col_stop - output_col,
                                output_row_stop - output_row,
                            ),
                        )

            _build_overviews(
                output_path,
                categorical=spec.categorical,
                wrapped_phase=spec.complex_phase,
            )
            outputs[key] = output_path

    return outputs
