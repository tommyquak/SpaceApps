"""Scientific QA helpers for short NISAR GUNW commissioning stacks.

The functions in this module are deliberately separate from annual velocity
estimation.  They validate MintPy's phase-to-displacement algebra, quantify
mission correction sensitivity, propagate the delivered ionosphere
uncertainty, and test how strongly the relative datum depends on the selected
reference patch.
"""

from __future__ import annotations

import hashlib
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Optional, Sequence

import numpy as np
from pyproj import Transformer
from rasterio.transform import xy
from scipy.ndimage import binary_erosion, distance_transform_edt

from peatguard.timeseries.nisar_gunw import (
    ReferenceSelection,
    _dominant_positive_label,
    _mintpy_transform,
    _nan_local_std,
    _rasterize_geojson,
    _reproject_binary_raster,
)


def sha256_file(path: Path) -> str:
    """Return a streaming SHA-256 digest without loading the file in memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_file_manifest(
    entries: Mapping[str, Mapping[str, object]],
    *,
    replacements: Optional[Mapping[str, Path]] = None,
) -> dict[str, dict[str, object]]:
    """Verify byte size and SHA-256 for every file in a manifest mapping.

    ``replacements`` allows a reference-support file to be relocated while
    still requiring byte-for-byte identity with the file used to select the
    relative reference.
    """
    replacements = replacements or {}
    verified: dict[str, dict[str, object]] = {}
    errors: list[str] = []
    for recorded_path, metadata in entries.items():
        actual_path = Path(replacements.get(recorded_path, Path(recorded_path)))
        expected_size = int(metadata["bytes"])
        expected_hash = str(metadata["sha256"])
        if not actual_path.exists():
            errors.append(f"missing: {actual_path}")
            continue
        actual_size = actual_path.stat().st_size
        if actual_size != expected_size:
            errors.append(
                f"size mismatch: {actual_path} "
                f"({actual_size} != {expected_size})"
            )
            continue
        actual_hash = sha256_file(actual_path)
        if actual_hash != expected_hash:
            errors.append(f"SHA-256 mismatch: {actual_path}")
            continue
        verified[recorded_path] = {
            "actual_path": str(actual_path.resolve()),
            "bytes": actual_size,
            "sha256": actual_hash,
            "verified": True,
        }
    if errors:
        raise ValueError(
            "Commissioning provenance verification failed: "
            + "; ".join(errors)
        )
    return verified


@dataclass(frozen=True)
class InversionValidation:
    """Numerical comparison between MintPy output and input phase."""

    variant: str
    pixel_count: int
    date_count: int
    maximum_absolute_error_m: float
    root_mean_square_error_m: float
    absolute_tolerance_m: float
    passed: bool
    phase_to_displacement_sign: str

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class DatumSensitivity:
    """Spread among plausible, spatially separated relative references."""

    candidate_count: int
    evaluated_count: int
    minimum_separation_m: float
    offsets_los_mm: tuple[float, ...]
    median_los_mm: float
    median_absolute_deviation_los_mm: float
    percentile_05_los_mm: float
    percentile_95_los_mm: float
    percentile_05_95_span_los_mm: float
    reference_status: str

    def to_dict(self) -> dict:
        return asdict(self)


def phase_to_los_displacement(
    phase_rad: np.ndarray,
    wavelength_m: float,
) -> np.ndarray:
    """Convert NISAR/MintPy interferometric phase to LOS displacement.

    NISAR GUNW uses reference times conjugate secondary phase.  MintPy's
    convention therefore multiplies phase by ``-wavelength / (4*pi)``.
    Negative displacement represents motion away from the sensor.
    """
    if not np.isfinite(wavelength_m) or wavelength_m <= 0:
        raise ValueError(f"Invalid radar wavelength: {wavelength_m}")
    return np.asarray(phase_rad) * (-wavelength_m / (4.0 * np.pi))


def los_to_vertical_displacement(
    los_displacement: np.ndarray,
    incidence_angle_degrees: np.ndarray,
) -> np.ndarray:
    """Apply a vertical-only projection using the per-pixel incidence angle."""
    incidence = np.asarray(incidence_angle_degrees, dtype=np.float64)
    if incidence.shape != np.asarray(los_displacement).shape:
        raise ValueError("LOS displacement and incidence-angle shapes differ")
    valid = np.isfinite(incidence) & (incidence > 0.0) & (incidence < 90.0)
    output = np.full(incidence.shape, np.nan, dtype=np.float64)
    output[valid] = np.asarray(los_displacement)[valid] / np.cos(
        np.deg2rad(incidence[valid])
    )
    return output


def validate_tree_inversion(
    variant: str,
    stack_path: Path,
    timeseries_path: Path,
    *,
    absolute_tolerance_m: float = 1e-6,
) -> InversionValidation:
    """Validate an exactly determined, acyclic MintPy inversion.

    This commissioning check intentionally rejects redundant or disconnected
    networks.  Once NISAR has enough epochs for network redundancy, residual
    and closure tests replace this exact tree-algebra check.
    """
    import h5py

    with h5py.File(stack_path, "r") as stack:
        phase = stack["unwrapPhase"][:].astype(np.float64)
        connected = stack["connectComponent"][:]
        pair_dates = [
            tuple(
                value.decode() if isinstance(value, bytes) else str(value)
                for value in pair
            )
            for pair in stack["date"][:]
        ]
        wavelength_m = float(stack.attrs["WAVELENGTH"])
    with h5py.File(timeseries_path, "r") as timeseries:
        observed = timeseries["timeseries"][:].astype(np.float64)
        dates = [
            value.decode() if isinstance(value, bytes) else str(value)
            for value in timeseries["date"][:]
        ]

    acquisitions = sorted({date for pair in pair_dates for date in pair})
    if dates != acquisitions:
        raise ValueError(
            f"MintPy dates differ from sorted stack acquisitions: {dates} != "
            f"{acquisitions}"
        )
    if len(pair_dates) != len(acquisitions) - 1:
        raise ValueError("Exact inversion validation requires an acyclic tree")

    date_index = {date: index for index, date in enumerate(acquisitions)}
    adjacency: dict[str, list[tuple[str, int, float]]] = {
        date: [] for date in acquisitions
    }
    for pair_index, (reference, secondary) in enumerate(pair_dates):
        adjacency[reference].append((secondary, pair_index, 1.0))
        adjacency[secondary].append((reference, pair_index, -1.0))

    expected = np.full(observed.shape, np.nan, dtype=np.float64)
    expected[0] = 0.0
    visited = {acquisitions[0]}
    queue = [acquisitions[0]]
    while queue:
        current = queue.pop(0)
        current_index = date_index[current]
        for neighbour, pair_index, direction in adjacency[current]:
            if neighbour in visited:
                continue
            neighbour_index = date_index[neighbour]
            delta = phase_to_los_displacement(
                phase[pair_index], wavelength_m
            )
            expected[neighbour_index] = (
                expected[current_index] + direction * delta
            )
            visited.add(neighbour)
            queue.append(neighbour)
    if len(visited) != len(acquisitions):
        raise ValueError("Interferogram network is disconnected")

    valid_pixels = np.all(
        np.isfinite(phase) & (phase != 0.0) & (connected > 0), axis=0
    )
    valid = (
        np.broadcast_to(valid_pixels, observed.shape)
        & np.isfinite(expected)
        & np.isfinite(observed)
    )
    differences = observed[valid] - expected[valid]
    if differences.size == 0:
        raise ValueError("No valid pixels available for inversion validation")

    maximum_error = float(np.max(np.abs(differences)))
    rmse = float(np.sqrt(np.mean(differences * differences)))
    return InversionValidation(
        variant=variant,
        pixel_count=int(valid_pixels.sum()),
        date_count=len(acquisitions),
        maximum_absolute_error_m=maximum_error,
        root_mean_square_error_m=rmse,
        absolute_tolerance_m=absolute_tolerance_m,
        passed=maximum_error <= absolute_tolerance_m,
        phase_to_displacement_sign="-wavelength/(4*pi)",
    )


def build_commissioning_quality_mask(
    coherence: np.ndarray,
    connected_components: np.ndarray,
    land_mask: np.ndarray,
    incidence_angle_degrees: np.ndarray,
    *,
    minimum_pair_coherence: float = 0.50,
) -> np.ndarray:
    """Build a conservative common-pair quality mask."""
    if coherence.shape != connected_components.shape:
        raise ValueError("Coherence and connected-component shapes differ")
    if coherence.shape[1:] != land_mask.shape:
        raise ValueError("Interferograms and land mask have different grids")
    if land_mask.shape != incidence_angle_degrees.shape:
        raise ValueError("Land mask and incidence-angle grids differ")
    incidence = np.asarray(incidence_angle_degrees)
    return (
        np.asarray(land_mask, dtype=bool)
        & np.all(np.isfinite(coherence), axis=0)
        & np.all(connected_components > 0, axis=0)
        & (np.min(coherence, axis=0) >= minimum_pair_coherence)
        & np.isfinite(incidence)
        & (incidence > 0.0)
        & (incidence < 90.0)
    )


def propagate_ionosphere_uncertainty(
    pair_uncertainty_rad: np.ndarray,
    wavelength_m: float,
    incidence_angle_degrees: np.ndarray,
    *,
    reference: Optional[ReferenceSelection] = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Propagate delivered pair uncertainties to final-date vertical mm.

    The first output is the pixel-only root-sum-square lower bound.  The
    second conservatively includes each pair's median reference-patch
    uncertainty without assuming that averaging removes spatially correlated
    ionosphere error.  Both exclude tropospheric, unwrapping, orbit, DEM,
    horizontal-motion, and absolute-datum uncertainty.
    """
    uncertainty = np.asarray(pair_uncertainty_rad, dtype=np.float64)
    if uncertainty.ndim != 3:
        raise ValueError("Pair uncertainty must have shape (pairs, rows, cols)")
    if uncertainty.shape[1:] != incidence_angle_degrees.shape:
        raise ValueError("Uncertainty and incidence-angle grids differ")
    if np.any(uncertainty[np.isfinite(uncertainty)] < 0):
        raise ValueError("Ionosphere uncertainty cannot be negative")

    # Strict sum: a pixel missing any pair's uncertainty gets NaN sigma and
    # is excluded downstream, never a silently understated finite sigma.
    lower_phase = np.sqrt(np.sum(uncertainty * uncertainty, axis=0))
    conservative_phase = lower_phase.copy()
    if reference is not None:
        radius = reference.patch_radius_pixels
        rows = slice(reference.y - radius, reference.y + radius + 1)
        columns = slice(reference.x - radius, reference.x + radius + 1)
        reference_sigma = np.array(
            [
                float(np.nanmedian(layer[rows, columns]))
                for layer in uncertainty
            ],
            dtype=np.float64,
        )
        if not np.all(np.isfinite(reference_sigma)):
            raise ValueError("Reference patch has invalid ionosphere uncertainty")
        conservative_phase = np.sqrt(
            np.sum(
                uncertainty * uncertainty
                + reference_sigma[:, np.newaxis, np.newaxis] ** 2,
                axis=0,
            )
        )

    lower_los = np.abs(
        phase_to_los_displacement(lower_phase, wavelength_m)
    )
    conservative_los = np.abs(
        phase_to_los_displacement(conservative_phase, wavelength_m)
    )
    return (
        los_to_vertical_displacement(
            lower_los, incidence_angle_degrees
        )
        * 1_000.0,
        los_to_vertical_displacement(
            conservative_los, incidence_angle_degrees
        )
        * 1_000.0,
    )


def rank_reference_patches(
    stack_path: Path,
    geometry_path: Path,
    *,
    peat_geojson: Optional[Path] = None,
    canal_mask: Optional[Path] = None,
    limit: int = 12,
    minimum_separation_m: float = 4_000.0,
    minimum_coherence: float = 0.70,
    maximum_local_phase_std_rad: float = 0.80,
    minimum_canal_distance_m: float = 1_000.0,
    minimum_peat_distance_m: float = 500.0,
    minimum_edge_distance_m: float = 1_000.0,
    patch_radius_pixels: int = 2,
    local_window_pixels: int = 7,
) -> tuple[ReferenceSelection, ...]:
    """Rank spatially separated alternative relative reference patches.

    The candidate gates and ranking weights intentionally match
    :func:`select_reference_patch`.  Spatial separation prevents a dense
    cluster of adjacent pixels from understating datum sensitivity.
    """
    import h5py

    if limit < 1:
        raise ValueError("Reference candidate limit must be positive")
    if minimum_separation_m < 0:
        raise ValueError("Reference separation cannot be negative")

    with h5py.File(stack_path, "r") as stack:
        attrs = dict(stack.attrs)
        coherence = stack["coherence"][:].astype(np.float32)
        phase = stack["unwrapPhase"][:].astype(np.float32)
        connected = stack["connectComponent"][:]
    with h5py.File(geometry_path, "r") as geometry:
        land = geometry["waterMask"][:].astype(bool)

    shape = land.shape
    transform = _mintpy_transform(attrs)
    target_crs = f"EPSG:{attrs['EPSG']}"
    spacing_y = abs(float(attrs["Y_STEP"]))
    spacing_x = abs(float(attrs["X_STEP"]))
    peat = (
        _rasterize_geojson(peat_geojson, shape, transform, target_crs)
        if peat_geojson
        else np.zeros(shape, dtype=bool)
    )
    canals = (
        _reproject_binary_raster(canal_mask, shape, transform, target_crs)
        if canal_mask
        else np.zeros(shape, dtype=bool)
    )
    peat_distance = (
        distance_transform_edt(~peat, sampling=(spacing_y, spacing_x))
        if peat.any()
        else np.full(shape, np.inf, dtype=np.float32)
    )
    canal_distance = (
        distance_transform_edt(~canals, sampling=(spacing_y, spacing_x))
        if canals.any()
        else np.full(shape, np.inf, dtype=np.float32)
    )
    rows, columns = np.indices(shape)
    edge_distance = np.minimum.reduce(
        [
            rows * spacing_y,
            (shape[0] - 1 - rows) * spacing_y,
            columns * spacing_x,
            (shape[1] - 1 - columns) * spacing_x,
        ]
    )
    dominant_labels = tuple(
        _dominant_positive_label(layer) for layer in connected
    )
    common_component = np.ones(shape, dtype=bool)
    for layer, label in zip(connected, dominant_labels):
        common_component &= layer == label

    local_roughness = np.maximum.reduce(
        [_nan_local_std(layer, local_window_pixels) for layer in phase]
    )
    minimum_pair_coherence = np.min(coherence, axis=0)
    finite = np.all(np.isfinite(phase), axis=0) & np.all(
        np.isfinite(coherence), axis=0
    )
    candidates = (
        land
        & finite
        & common_component
        & ~peat
        & (minimum_pair_coherence >= minimum_coherence)
        & (local_roughness <= maximum_local_phase_std_rad)
        & (canal_distance >= minimum_canal_distance_m)
        & (peat_distance >= minimum_peat_distance_m)
        & (edge_distance >= minimum_edge_distance_m)
    )
    if patch_radius_pixels > 0:
        candidates = binary_erosion(
            candidates,
            structure=np.ones(
                (2 * patch_radius_pixels + 1, 2 * patch_radius_pixels + 1),
                dtype=bool,
            ),
            border_value=0,
        )
    candidate_count = int(candidates.sum())
    if candidate_count == 0:
        raise ValueError("No alternative reference candidate passed QA gates")

    roughness_score = 1.0 - np.clip(
        local_roughness / maximum_local_phase_std_rad, 0.0, 1.0
    )
    peat_score = (
        np.ones(shape, dtype=np.float32)
        if not peat.any()
        else np.clip(peat_distance / 5_000.0, 0.0, 1.0)
    )
    canal_score = (
        np.ones(shape, dtype=np.float32)
        if not canals.any()
        else np.clip(canal_distance / 5_000.0, 0.0, 1.0)
    )
    score = (
        0.55 * minimum_pair_coherence
        + 0.20 * roughness_score
        + 0.15 * peat_score
        + 0.10 * canal_score
    )
    score[~candidates] = -np.inf
    available = candidates.copy()
    selections: list[ReferenceSelection] = []

    while available.any() and len(selections) < limit:
        available_score = np.where(available, score, -np.inf)
        selected_y, selected_x = np.unravel_index(
            np.argmax(available_score), shape
        )
        projected_x, projected_y = xy(
            transform, selected_y, selected_x, offset="center"
        )
        longitude, latitude = Transformer.from_crs(
            target_crs, 4326, always_xy=True
        ).transform(projected_x, projected_y)
        selections.append(
            ReferenceSelection(
                x=int(selected_x),
                y=int(selected_y),
                longitude=float(longitude),
                latitude=float(latitude),
                patch_radius_pixels=patch_radius_pixels,
                pixel_spacing_m=(spacing_x, spacing_y),
                candidate_count=candidate_count,
                minimum_pair_coherence=float(
                    minimum_pair_coherence[selected_y, selected_x]
                ),
                pair_coherence=tuple(
                    float(value)
                    for value in coherence[:, selected_y, selected_x]
                ),
                maximum_local_phase_std_rad=float(
                    local_roughness[selected_y, selected_x]
                ),
                distance_to_mapped_canal_m=float(
                    canal_distance[selected_y, selected_x]
                ),
                distance_to_mapped_peat_m=float(
                    peat_distance[selected_y, selected_x]
                ),
                dominant_connected_components=dominant_labels,
                reference_status="relative_reference_not_ground_control",
            )
        )
        if minimum_separation_m == 0:
            available[selected_y, selected_x] = False
        else:
            distance_squared = (
                (rows - selected_y) * spacing_y
            ) ** 2 + ((columns - selected_x) * spacing_x) ** 2
            available[distance_squared < minimum_separation_m**2] = False

    return tuple(selections)


def assess_datum_sensitivity(
    corrected_stack_path: Path,
    candidates: Sequence[ReferenceSelection],
    *,
    minimum_reference_coherence: float = 0.70,
    minimum_separation_m: float = 4_000.0,
) -> DatumSensitivity:
    """Measure final-date LOS datum offsets among plausible references."""
    import h5py

    if not candidates:
        raise ValueError("At least one reference candidate is required")
    with h5py.File(corrected_stack_path, "r") as stack:
        phase = stack["unwrapPhase"][:].astype(np.float64)
        coherence = stack["coherence"][:]
        connected = stack["connectComponent"][:]
        wavelength_m = float(stack.attrs["WAVELENGTH"])
        pair_dates = [
            tuple(
                value.decode() if isinstance(value, bytes) else str(value)
                for value in row
            )
            for row in stack["date"][:]
        ]

    # The unconditional sum of per-pair phases below telescopes to the
    # final-date phase only for a date-ordered nearest-neighbour chain.
    # Refuse any other network topology loudly rather than summing phases
    # with wrong signs or non-telescoping spans.
    for reference_date, secondary_date in pair_dates:
        if reference_date >= secondary_date:
            raise ValueError(
                "assess_datum_sensitivity requires reference < secondary "
                f"in every pair, got {reference_date} >= {secondary_date}"
            )
    for index in range(len(pair_dates) - 1):
        if pair_dates[index][1] != pair_dates[index + 1][0]:
            raise ValueError(
                "assess_datum_sensitivity requires a sequential "
                "nearest-neighbour chain of pairs; found a gap between "
                f"{pair_dates[index]} and {pair_dates[index + 1]}"
            )

    offsets: list[float] = []
    for candidate in candidates:
        radius = candidate.patch_radius_pixels
        rows = slice(candidate.y - radius, candidate.y + radius + 1)
        columns = slice(candidate.x - radius, candidate.x + radius + 1)
        cumulative_phase = 0.0
        usable = True
        for pair_index in range(phase.shape[0]):
            valid = (
                np.isfinite(phase[pair_index, rows, columns])
                & (
                    coherence[pair_index, rows, columns]
                    >= minimum_reference_coherence
                )
                & (connected[pair_index, rows, columns] > 0)
            )
            values = phase[pair_index, rows, columns][valid]
            minimum_pixels = max(3, (2 * radius + 1) ** 2 // 2)
            if values.size < minimum_pixels:
                usable = False
                break
            cumulative_phase += float(np.median(values))
        if usable:
            offset_m = float(
                phase_to_los_displacement(
                    np.asarray(cumulative_phase), wavelength_m
                )
            )
            offsets.append(offset_m * 1_000.0)
    if not offsets:
        raise ValueError("No reference candidate was valid after corrections")

    values = np.asarray(offsets, dtype=np.float64)
    median = float(np.median(values))
    mad = float(np.median(np.abs(values - median)))
    p05, p95 = np.percentile(values, [5.0, 95.0])
    return DatumSensitivity(
        candidate_count=int(candidates[0].candidate_count),
        evaluated_count=len(offsets),
        minimum_separation_m=minimum_separation_m,
        offsets_los_mm=tuple(float(value) for value in values),
        median_los_mm=median,
        median_absolute_deviation_los_mm=mad,
        percentile_05_los_mm=float(p05),
        percentile_95_los_mm=float(p95),
        percentile_05_95_span_los_mm=float(p95 - p05),
        reference_status="relative_reference_not_ground_control",
    )


def summary_statistics(
    data: np.ndarray,
    mask: np.ndarray,
) -> dict[str, object]:
    """Return stable finite-value statistics for a report."""
    values = np.asarray(data)[np.asarray(mask, dtype=bool)]
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {"count": 0}
    percentiles = np.percentile(
        values, [2.5, 25.0, 50.0, 75.0, 97.5]
    )
    return {
        "count": int(values.size),
        "mean": float(np.mean(values)),
        "standard_deviation": float(np.std(values)),
        "minimum": float(np.min(values)),
        "percentile_02_5": float(percentiles[0]),
        "percentile_25": float(percentiles[1]),
        "median": float(percentiles[2]),
        "percentile_75": float(percentiles[3]),
        "percentile_97_5": float(percentiles[4]),
        "maximum": float(np.max(values)),
    }
