"""Fine-resolution interferograms from NISAR L2 GSLC products.

NISAR 40 MHz Frequency A GSLC products are delivered on a 5 m map grid.
Their complex phase is geocoded and topographically flattened, so matching
same-track products can be cross-multiplied to form a fine-resolution wrapped
interferogram.

This module deliberately stops at wrapped phase and local coherence. It does
not call a two-epoch phase difference a velocity product. Unwrapping,
ionosphere/troposphere correction, reference-area selection, and time-series
inversion are separate downstream stages.

The implementation is block-wise because a full PeatGuard AOI at 5 m contains
well over 100 million pixels and should not be loaded into memory at once.
"""

from __future__ import annotations

import logging
import re
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.enums import Resampling
from rasterio.transform import Affine
from rasterio.windows import Window
from scipy.ndimage import uniform_filter

logger = logging.getLogger(__name__)

DEFAULT_FREQUENCY = "frequencyA"
DEFAULT_POLARIZATION = "HH"
DEFAULT_NODATA = -9999.0
SPEED_OF_LIGHT_M_S = 299_792_458.0


@dataclass(frozen=True)
class GSLCProductInfo:
    """Metadata needed to align and process one NISAR GSLC product."""

    path: Path
    dataset_path: str
    frequency: str
    polarization: str
    epsg: int
    shape: tuple[int, int]
    x_coordinates: np.ndarray
    y_coordinates: np.ndarray
    x_spacing: float
    y_spacing: float
    center_frequency_hz: float
    wavelength_m: float
    number_of_subswaths: int
    acquisition_start: str
    track_number: Optional[int]
    frame_number: Optional[int]
    granule_id: str

    @property
    def transform(self) -> Affine:
        """Return a pixel-corner affine transform for centre coordinates."""
        return Affine(
            self.x_spacing,
            0.0,
            float(self.x_coordinates[0]) - self.x_spacing / 2.0,
            0.0,
            self.y_spacing,
            float(self.y_coordinates[0]) - self.y_spacing / 2.0,
        )


def _decode_scalar(value: object) -> object:
    """Decode HDF5 scalar byte strings while preserving numeric values."""
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.generic):
        return value.item()
    return value


def _read_optional_scalar(group: object, name: str) -> object | None:
    """Read an optional scalar dataset from an h5py group."""
    if name not in group:
        return None
    dataset = group[name]
    if dataset.shape != ():
        return None
    return _decode_scalar(dataset[()])


def _regular_spacing(coordinates: np.ndarray, label: str) -> float:
    """Validate a coordinate vector and return its regular spacing."""
    if coordinates.ndim != 1 or coordinates.size < 2:
        raise ValueError(f"{label} must be a one-dimensional vector with at least 2 values")

    differences = np.diff(coordinates.astype(np.float64))
    spacing = float(np.median(differences))
    tolerance = max(abs(spacing) * 1e-6, 1e-8)
    if spacing == 0 or not np.allclose(differences, spacing, rtol=0.0, atol=tolerance):
        raise ValueError(f"{label} is not regularly spaced")
    return spacing


def _fallback_date_from_name(name: str) -> str:
    """Extract the first acquisition timestamp from a NISAR granule name."""
    match = re.search(r"_(\d{8}T\d{6})", name)
    return match.group(1) if match else ""


def inspect_gslc(
    path: Path,
    frequency: str = DEFAULT_FREQUENCY,
    polarization: str = DEFAULT_POLARIZATION,
) -> GSLCProductInfo:
    """Inspect a NISAR GSLC without reading its complex raster into memory.

    Args:
        path: NISAR L2 GSLC HDF5 file.
        frequency: Frequency group, normally ``frequencyA``.
        polarization: Complex polarization dataset, normally ``HH``.

    Returns:
        Product and grid metadata required for interferogram formation.
    """
    import h5py

    path = Path(path)
    grid_path = f"science/LSAR/GSLC/grids/{frequency}"
    dataset_path = f"{grid_path}/{polarization}"
    identification_path = "science/LSAR/identification"

    with h5py.File(path, "r") as h5:
        if dataset_path not in h5:
            available = []
            if grid_path in h5:
                available = sorted(h5[grid_path].keys())
            raise KeyError(
                f"GSLC dataset not found: {dataset_path}; available entries: {available}"
            )

        grid = h5[grid_path]
        for coordinate_name in ("xCoordinates", "yCoordinates", "projection"):
            if coordinate_name not in grid:
                raise KeyError(f"Missing GSLC grid dataset: {grid_path}/{coordinate_name}")

        x_coordinates = grid["xCoordinates"][:].astype(np.float64)
        y_coordinates = grid["yCoordinates"][:].astype(np.float64)
        epsg = int(grid["projection"][()])
        if "centerFrequency" not in grid:
            raise KeyError(f"Missing GSLC grid dataset: {grid_path}/centerFrequency")
        center_frequency_hz = float(grid["centerFrequency"][()])
        if not np.isfinite(center_frequency_hz) or center_frequency_hz <= 0:
            raise ValueError(
                f"Invalid GSLC center frequency: {center_frequency_hz}"
            )
        if "numberOfSubSwaths" not in grid:
            raise KeyError(
                f"Missing GSLC grid dataset: {grid_path}/numberOfSubSwaths"
            )
        number_of_subswaths = int(grid["numberOfSubSwaths"][()])
        if not 1 <= number_of_subswaths <= 254:
            raise ValueError(
                f"Invalid GSLC numberOfSubSwaths: {number_of_subswaths}"
            )
        shape = tuple(int(value) for value in h5[dataset_path].shape)
        if shape != (y_coordinates.size, x_coordinates.size):
            raise ValueError(
                f"GSLC raster shape {shape} does not match coordinate vectors "
                f"({y_coordinates.size}, {x_coordinates.size})"
            )

        identification = h5[identification_path]
        acquisition_start = str(
            _read_optional_scalar(identification, "zeroDopplerStartTime")
            or _read_optional_scalar(identification, "referenceZeroDopplerStartTime")
            or _fallback_date_from_name(path.name)
        )
        granule_id = str(
            _read_optional_scalar(identification, "granuleId") or path.stem
        )
        track_number = _read_optional_scalar(identification, "trackNumber")
        frame_number = _read_optional_scalar(identification, "frameNumber")

    x_spacing = _regular_spacing(x_coordinates, "xCoordinates")
    y_spacing = _regular_spacing(y_coordinates, "yCoordinates")
    if x_spacing <= 0:
        raise ValueError(f"xCoordinates must increase from west to east, got {x_spacing}")
    if y_spacing >= 0:
        raise ValueError(f"yCoordinates must decrease from north to south, got {y_spacing}")

    return GSLCProductInfo(
        path=path,
        dataset_path=dataset_path,
        frequency=frequency,
        polarization=polarization,
        epsg=epsg,
        shape=shape,
        x_coordinates=x_coordinates,
        y_coordinates=y_coordinates,
        x_spacing=x_spacing,
        y_spacing=y_spacing,
        center_frequency_hz=center_frequency_hz,
        wavelength_m=SPEED_OF_LIGHT_M_S / center_frequency_hz,
        number_of_subswaths=number_of_subswaths,
        acquisition_start=acquisition_start,
        track_number=int(track_number) if track_number is not None else None,
        frame_number=int(frame_number) if frame_number is not None else None,
        granule_id=granule_id,
    )


def _validate_pair(reference: GSLCProductInfo, secondary: GSLCProductInfo) -> None:
    """Require matching grids and acquisition geometry for cross-multiplication."""
    if reference.epsg != secondary.epsg:
        raise ValueError(f"GSLC EPSG mismatch: {reference.epsg} != {secondary.epsg}")
    if reference.shape != secondary.shape:
        raise ValueError(f"GSLC shape mismatch: {reference.shape} != {secondary.shape}")
    if reference.track_number != secondary.track_number:
        raise ValueError(
            f"GSLC track mismatch: {reference.track_number} != {secondary.track_number}"
        )
    if reference.frame_number != secondary.frame_number:
        raise ValueError(
            f"GSLC frame mismatch: {reference.frame_number} != {secondary.frame_number}"
        )
    if reference.number_of_subswaths != secondary.number_of_subswaths:
        raise ValueError(
            "GSLC subswath-count mismatch: "
            f"{reference.number_of_subswaths} != "
            f"{secondary.number_of_subswaths}"
        )
    if not np.isclose(
        reference.center_frequency_hz,
        secondary.center_frequency_hz,
        rtol=0.0,
        atol=1.0,
    ):
        raise ValueError(
            "GSLC center-frequency mismatch: "
            f"{reference.center_frequency_hz} != {secondary.center_frequency_hz}"
        )
    # Coordinate magnitude must not relax alignment: NumPy's default relative
    # tolerance can admit whole-pixel shifts at UTM eastings and northings.
    if not np.allclose(
        reference.x_coordinates, secondary.x_coordinates, rtol=0.0, atol=1e-6
    ):
        raise ValueError("GSLC x-coordinate grids do not match")
    if not np.allclose(
        reference.y_coordinates, secondary.y_coordinates, rtol=0.0, atol=1e-6
    ):
        raise ValueError("GSLC y-coordinate grids do not match")


def _as_complex(array: np.ndarray) -> np.ndarray:
    """Normalize HDF5 complex encodings to a NumPy complex64 array."""
    if np.iscomplexobj(array):
        return array.astype(np.complex64, copy=False)

    fields = array.dtype.fields or {}
    real_name = next((name for name in ("r", "real") if name in fields), None)
    imag_name = next((name for name in ("i", "imag") if name in fields), None)
    if real_name and imag_name:
        return (
            array[real_name].astype(np.float32)
            + 1j * array[imag_name].astype(np.float32)
        ).astype(np.complex64)

    raise TypeError(f"Unsupported GSLC complex dtype: {array.dtype}")


def _window_for_bbox(
    info: GSLCProductInfo,
    bbox_wgs84: Optional[tuple[float, float, float, float]],
) -> tuple[int, int, int, int]:
    """Return row/column bounds covering a WGS84 bbox."""
    height, width = info.shape
    if bbox_wgs84 is None:
        return 0, height, 0, width

    west, south, east, north = bbox_wgs84
    if west >= east or south >= north:
        raise ValueError("bbox_wgs84 must be (west, south, east, north)")

    transformer = Transformer.from_crs(4326, info.epsg, always_xy=True)
    min_x, min_y, max_x, max_y = transformer.transform_bounds(
        west, south, east, north, densify_pts=21
    )

    col_start = int(np.searchsorted(info.x_coordinates, min_x, side="left"))
    col_stop = int(np.searchsorted(info.x_coordinates, max_x, side="right"))
    negative_y = -info.y_coordinates
    row_start = int(np.searchsorted(negative_y, -max_y, side="left"))
    row_stop = int(np.searchsorted(negative_y, -min_y, side="right"))

    row_start = max(0, min(height, row_start))
    row_stop = max(0, min(height, row_stop))
    col_start = max(0, min(width, col_start))
    col_stop = max(0, min(width, col_stop))
    if row_stop <= row_start or col_stop <= col_start:
        raise ValueError(
            f"Requested AOI does not intersect GSLC grid for {info.path.name}"
        )
    return row_start, row_stop, col_start, col_stop


def _local_coherence(
    reference: np.ndarray,
    secondary: np.ndarray,
    window_size: int,
    reference_mask: Optional[np.ndarray] = None,
    secondary_mask: Optional[np.ndarray] = None,
    reference_subswaths: Optional[int] = None,
    secondary_subswaths: Optional[int] = None,
    reference_exception_mask: Optional[np.ndarray] = None,
    secondary_exception_mask: Optional[np.ndarray] = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return wrapped phase and normalized local interferometric coherence.

    Positive NISAR L2 GSLC mask values identify fully focused samples from the
    corresponding sub-swath. Value 0 identifies partially focused or invalid
    data and 255 is fill. When supplied, ``numberOfSubSwaths`` bounds the
    accepted sub-swath identifiers.
    """
    valid = (
        np.isfinite(reference.real)
        & np.isfinite(reference.imag)
        & np.isfinite(secondary.real)
        & np.isfinite(secondary.imag)
        & (np.abs(reference) > 0)
        & (np.abs(secondary) > 0)
    )
    for mask, number_of_subswaths in (
        (reference_mask, reference_subswaths),
        (secondary_mask, secondary_subswaths),
    ):
        if mask is None:
            continue
        if mask.shape != reference.shape:
            raise ValueError("GSLC quality mask shape does not match complex data")
        maximum_valid = (
            number_of_subswaths
            if number_of_subswaths is not None
            else 254
        )
        if not 1 <= maximum_valid <= 254:
            raise ValueError(
                f"Invalid GSLC numberOfSubSwaths: {maximum_valid}"
            )
        valid &= (mask >= 1) & (mask <= maximum_valid)
    for exception_mask in (
        reference_exception_mask,
        secondary_exception_mask,
    ):
        if exception_mask is None:
            continue
        if exception_mask.shape != reference.shape:
            raise ValueError(
                "GSLC input-data exception mask shape does not match "
                "complex data"
            )
        valid &= exception_mask == 0

    with np.errstate(invalid="ignore", over="ignore", under="ignore"):
        cross = reference * np.conj(secondary)
    valid &= np.isfinite(cross.real) & np.isfinite(cross.imag) & (np.abs(cross) > 0)
    cross = np.where(valid, cross, 0.0j)
    reference_power = np.where(valid, np.abs(reference) ** 2, 0.0).astype(np.float32)
    secondary_power = np.where(valid, np.abs(secondary) ** 2, 0.0).astype(np.float32)

    mean_cross = uniform_filter(
        cross.real.astype(np.float32), size=window_size, mode="constant", cval=0.0
    ) + 1j * uniform_filter(
        cross.imag.astype(np.float32), size=window_size, mode="constant", cval=0.0
    )
    mean_reference_power = uniform_filter(
        reference_power, size=window_size, mode="constant", cval=0.0
    )
    mean_secondary_power = uniform_filter(
        secondary_power, size=window_size, mode="constant", cval=0.0
    )
    valid_fraction = uniform_filter(
        valid.astype(np.float32), size=window_size, mode="constant", cval=0.0
    )

    denominator = np.sqrt(mean_reference_power * mean_secondary_power)
    coherence = np.full(reference.shape, np.nan, dtype=np.float32)
    # Do not spatially infill coherence into a centre pixel that the mission
    # mask itself marks invalid, even if enough valid neighbours exist.
    usable = valid & (denominator > 0) & (valid_fraction >= 0.5)
    coherence[usable] = np.clip(
        np.abs(mean_cross[usable]) / denominator[usable], 0.0, 1.0
    )

    phase = np.full(reference.shape, np.nan, dtype=np.float32)
    phase[valid] = np.angle(cross[valid]).astype(np.float32)
    return phase, coherence


def _date_token(value: str) -> str:
    """Return YYYYMMDD for stable filenames."""
    dashed = re.search(r"(\d{4})-(\d{2})-(\d{2})", value)
    if dashed:
        return "".join(dashed.groups())
    compact = re.search(r"(?<!\d)(\d{8})T", value)
    return compact.group(1) if compact else "unknown"


def _build_overviews(path: Path, resampling: Resampling) -> None:
    """Build internal overviews after block-wise raster writing."""
    with rasterio.open(path, "r+") as dataset:
        minimum_dimension = min(dataset.width, dataset.height)
        levels = [level for level in (2, 4, 8, 16) if minimum_dimension // level >= 1]
        if levels:
            dataset.build_overviews(levels, resampling)
            dataset.update_tags(ns="rio_overview", resampling=resampling.name)


def form_gslc_interferogram(
    reference_path: Path,
    secondary_path: Path,
    output_dir: Path,
    bbox_wgs84: Optional[tuple[float, float, float, float]] = None,
    frequency: str = DEFAULT_FREQUENCY,
    polarization: str = DEFAULT_POLARIZATION,
    coherence_window: int = 5,
    block_size: int = 1024,
    nodata: float = DEFAULT_NODATA,
    require_quality_mask: bool = True,
    export_cross_product: bool = False,
) -> dict[str, Path]:
    """Form wrapped phase and coherence from a matching NISAR GSLC pair.

    The phase convention is ``reference * conjugate(secondary)``. Before
    converting phase to displacement, compare an aggregated version against
    the corresponding mission GUNW to verify the date/sign convention.

    Args:
        reference_path: Earlier NISAR GSLC HDF5 file.
        secondary_path: Later NISAR GSLC HDF5 file.
        output_dir: Directory for tiled GeoTIFF products.
        bbox_wgs84: Optional ``(west, south, east, north)`` crop.
        frequency: Frequency group, normally ``frequencyA``.
        polarization: Polarization dataset, normally ``HH``.
        coherence_window: Odd moving-window size in GSLC pixels.
        block_size: Output processing block size.
        nodata: Output NoData value.
        require_quality_mask: Require and apply the mission GSLC valid-data
            mask from both acquisitions. Disabling this is diagnostic only.
        export_cross_product: Also export the amplitude-preserving complex
            cross-product as real and imaginary float32 bands. This requires
            the original GSLCs, not phase or coherence reconstructions.

    Returns:
        Mapping with ``wrapped_phase`` and ``coherence`` GeoTIFF paths, plus
        ``complex_cross_product`` when requested.
    """
    import h5py

    if coherence_window < 3 or coherence_window % 2 == 0:
        raise ValueError("coherence_window must be an odd integer >= 3")
    if block_size < coherence_window:
        raise ValueError("block_size must be at least coherence_window")

    reference = inspect_gslc(reference_path, frequency, polarization)
    secondary = inspect_gslc(secondary_path, frequency, polarization)
    _validate_pair(reference, secondary)
    if reference.acquisition_start >= secondary.acquisition_start:
        raise ValueError(
            "reference_path must contain the earlier acquisition: "
            f"{reference.acquisition_start} >= {secondary.acquisition_start}"
        )

    row_start, row_stop, col_start, col_stop = _window_for_bbox(
        reference, bbox_wgs84
    )
    output_height = row_stop - row_start
    output_width = col_stop - col_start
    output_transform = reference.transform * Affine.translation(col_start, row_start)

    reference_date = _date_token(reference.acquisition_start)
    secondary_date = _date_token(secondary.acquisition_start)
    prefix = f"nisar_gslc_{reference_date}_{secondary_date}"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    phase_path = output_dir / f"{prefix}_wrapped_phase_5m.tif"
    coherence_path = output_dir / f"{prefix}_coherence_5m.tif"
    cross_path = output_dir / f"{prefix}_complex_cross_product_5m.tif"

    profile = {
        "driver": "GTiff",
        "dtype": "float32",
        "width": output_width,
        "height": output_height,
        "count": 1,
        "crs": f"EPSG:{reference.epsg}",
        "transform": output_transform,
        "nodata": nodata,
        "tiled": True,
        "blockxsize": 512,
        "blockysize": 512,
        "compress": "DEFLATE",
        "predictor": 3,
        "BIGTIFF": "IF_SAFER",
    }
    common_tags = {
        "sensor": "NISAR L-SAR",
        "source_product": "L2 GSLC",
        "reference_granule": reference.granule_id,
        "secondary_granule": secondary.granule_id,
        "reference_date": reference.acquisition_start,
        "secondary_date": secondary.acquisition_start,
        "frequency": frequency,
        "center_frequency_hz": str(reference.center_frequency_hz),
        "wavelength_m": str(reference.wavelength_m),
        "wavelength_source": "source_product_centerFrequency",
        "polarization": polarization,
        "track_number": str(reference.track_number),
        "frame_number": str(reference.frame_number),
        "coherence_window_pixels": str(coherence_window),
        "pixel_spacing_x_crs_units": str(abs(reference.x_spacing)),
        "pixel_spacing_y_crs_units": str(abs(reference.y_spacing)),
        "phase_convention": "reference_times_conjugate_secondary",
        "mission_valid_data_mask": (
            "subswath_ids_1_to_product_number_applied"
            if require_quality_mask
            else "applied_when_available"
        ),
        "reference_number_of_subswaths": str(reference.number_of_subswaths),
        "secondary_number_of_subswaths": str(secondary.number_of_subswaths),
        "mission_valid_data_mask_source": (
            "NISAR_L2_GSLC_Product_Specification_Rev_E_section_4.4"
        ),
        "mission_input_data_exception_mask": (
            "nonzero_values_excluded"
            if require_quality_mask
            else "applied_when_available"
        ),
        "processing_status": "wrapped_phase_not_displacement_or_velocity",
    }
    output_crs = rasterio.crs.CRS.from_epsg(reference.epsg)
    metre_grid = (
        output_crs.is_projected
        and output_crs.linear_units_factor[1] == 1.0
    )

    logger.info(
        "Forming GSLC interferogram %s--%s: %dx%d pixels at %.1fx%.1f m",
        reference_date,
        secondary_date,
        output_width,
        output_height,
        abs(reference.x_spacing),
        abs(reference.y_spacing),
    )

    halo = coherence_window // 2
    quality_mask_path = f"science/LSAR/GSLC/grids/{frequency}/mask"
    exception_mask_path = (
        f"science/LSAR/GSLC/grids/{frequency}/inputDataExceptionMask"
    )
    if require_quality_mask:
        with (
            h5py.File(reference.path, "r") as reference_check,
            h5py.File(secondary.path, "r") as secondary_check,
        ):
            missing = [
                f"{label}:{path}"
                for label, source in (
                    ("reference", reference_check),
                    ("secondary", secondary_check),
                )
                for path in (quality_mask_path, exception_mask_path)
                if path not in source
            ]
            if missing:
                raise KeyError(
                    "Both GSLC products must contain the mission valid-data "
                    "mask and input-data exception mask; missing "
                    f"{missing}. Use require_quality_mask=False only for "
                    "diagnostics"
                )

    with (
        ExitStack() as optional_outputs,
        h5py.File(reference.path, "r") as reference_h5,
        h5py.File(secondary.path, "r") as secondary_h5,
        rasterio.open(phase_path, "w", **profile) as phase_output,
        rasterio.open(coherence_path, "w", **profile) as coherence_output,
    ):
        cross_output = (
            optional_outputs.enter_context(
                rasterio.open(cross_path, "w", **{**profile, "count": 2, "nodata": np.nan})
            )
            if export_cross_product
            else None
        )
        reference_dataset = reference_h5[reference.dataset_path]
        secondary_dataset = secondary_h5[secondary.dataset_path]
        reference_mask_dataset = reference_h5.get(quality_mask_path)
        secondary_mask_dataset = secondary_h5.get(quality_mask_path)
        reference_exception_dataset = reference_h5.get(exception_mask_path)
        secondary_exception_dataset = secondary_h5.get(exception_mask_path)
        phase_tags = {
            **common_tags,
            "units": "radians",
            "effective_support": "native_mode_dependent",
        }
        coherence_tags = {
            **common_tags,
            "units": "1",
            "effective_support": (
                f"{coherence_window}_by_{coherence_window}_source_pixels"
            ),
        }
        if metre_grid:
            pixel_posting_m = max(
                abs(reference.x_spacing),
                abs(reference.y_spacing),
            )
            phase_tags["pixel_posting_m"] = str(pixel_posting_m)
            coherence_tags.update(
                pixel_posting_m=str(pixel_posting_m),
                effective_support_x_m=str(
                    abs(reference.x_spacing) * coherence_window
                ),
                effective_support_y_m=str(
                    abs(reference.y_spacing) * coherence_window
                ),
            )
        phase_output.update_tags(**phase_tags)
        coherence_output.update_tags(**coherence_tags)
        phase_output.set_band_description(1, "wrapped_interferometric_phase")
        coherence_output.set_band_description(1, "local_interferometric_coherence")
        if cross_output is not None:
            cross_output.update_tags(
                **{**phase_tags, "units": "source_complex_amplitude_product"},
                complex_representation="real_imaginary_bands",
                amplitude_source="unscaled_reference_times_conjugate_secondary",
                spatial_aggregation="none_native_geocoded_gslc_grid",
                operator_caveat="geocoded_cross_product_not_radar_multilook_then_geocode",
            )
            cross_output.set_band_description(1, "cross_product_real")
            cross_output.set_band_description(2, "cross_product_imaginary")

        for output_row in range(0, output_height, block_size):
            output_row_stop = min(output_height, output_row + block_size)
            global_row = row_start + output_row
            global_row_stop = row_start + output_row_stop
            read_row = max(0, global_row - halo)
            read_row_stop = min(reference.shape[0], global_row_stop + halo)

            for output_col in range(0, output_width, block_size):
                output_col_stop = min(output_width, output_col + block_size)
                global_col = col_start + output_col
                global_col_stop = col_start + output_col_stop
                read_col = max(0, global_col - halo)
                read_col_stop = min(reference.shape[1], global_col_stop + halo)

                source_slice = np.s_[read_row:read_row_stop, read_col:read_col_stop]
                reference_block = _as_complex(reference_dataset[source_slice])
                secondary_block = _as_complex(secondary_dataset[source_slice])
                reference_mask_block = (
                    reference_mask_dataset[source_slice]
                    if reference_mask_dataset is not None
                    else None
                )
                secondary_mask_block = (
                    secondary_mask_dataset[source_slice]
                    if secondary_mask_dataset is not None
                    else None
                )
                reference_exception_block = (
                    reference_exception_dataset[source_slice]
                    if reference_exception_dataset is not None
                    else None
                )
                secondary_exception_block = (
                    secondary_exception_dataset[source_slice]
                    if secondary_exception_dataset is not None
                    else None
                )
                phase, coherence = _local_coherence(
                    reference_block,
                    secondary_block,
                    coherence_window,
                    reference_mask=reference_mask_block,
                    secondary_mask=secondary_mask_block,
                    reference_subswaths=reference.number_of_subswaths,
                    secondary_subswaths=secondary.number_of_subswaths,
                    reference_exception_mask=reference_exception_block,
                    secondary_exception_mask=secondary_exception_block,
                )

                crop_row = global_row - read_row
                crop_col = global_col - read_col
                block_height = output_row_stop - output_row
                block_width = output_col_stop - output_col
                crop = np.s_[
                    crop_row : crop_row + block_height,
                    crop_col : crop_col + block_width,
                ]
                phase_block = np.where(np.isfinite(phase[crop]), phase[crop], nodata)
                coherence_block = np.where(
                    np.isfinite(coherence[crop]), coherence[crop], nodata
                )
                output_window = Window(
                    output_col, output_row, block_width, block_height
                )
                phase_output.write(
                    phase_block.astype(np.float32, copy=False), 1, window=output_window
                )
                coherence_output.write(
                    coherence_block.astype(np.float32, copy=False),
                    1,
                    window=output_window,
                )
                if cross_output is not None:
                    # Use the original samples and exactly the phase validity
                    # mask. Coherence is not an amplitude proxy.
                    with np.errstate(invalid="ignore", over="ignore", under="ignore"):
                        cross_block = reference_block[crop] * np.conj(secondary_block[crop])
                    cross_valid = np.isfinite(phase[crop])
                    cross_output.write(
                        np.where(cross_valid, cross_block.real, np.nan).astype(np.float32),
                        1,
                        window=output_window,
                    )
                    cross_output.write(
                        np.where(cross_valid, cross_block.imag, np.nan).astype(np.float32),
                        2,
                        window=output_window,
                    )

    # Nearest avoids circular-mean artefacts across +/-pi phase wraps.
    _build_overviews(phase_path, Resampling.nearest)
    _build_overviews(coherence_path, Resampling.average)
    logger.info("Wrote %s and %s", phase_path, coherence_path)
    products = {"wrapped_phase": phase_path, "coherence": coherence_path}
    if export_cross_product:
        products["complex_cross_product"] = cross_path
    return products
