"""Auditable NISAR GUNW commissioning utilities.

NISAR GUNW phase corrections are delivered as phase-screen contributions
that are not applied to the unwrapped interferogram. Corrected phase is:

    raw phase - ionosphere - troposphere - solid Earth tides

This module prepares explicit comparison variants and selects a reproducible
relative reference patch. It does not promote a short commissioning stack to
an annual velocity product.
"""

from __future__ import annotations

import json
import shutil
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Mapping, Optional

import numpy as np
from pyproj import Transformer
from rasterio.features import rasterize
from rasterio.transform import Affine, xy
from rasterio.warp import Resampling, reproject, transform_geom
from scipy.ndimage import binary_erosion, distance_transform_edt, uniform_filter

MINTPY_NISAR_COMMIT = "53030ef4ae390ff17e845a342e086993eca81115"
CORRECTION_SIGN_CONVENTION = "corrected_phase=raw-ionosphere-troposphere-set"


@dataclass(frozen=True)
class TemporalReadiness:
    """Temporal network status and conservative release gates."""

    pair_count: int
    acquisition_count: int
    first_date: str
    last_date: str
    span_days: int
    connected_components: int
    network_cycle_count: int
    minimum_acquisition_degree: int
    bridge_edge_count: int
    maximum_connection_bandwidth: int
    complete_connection_bandwidth: int
    tay_target_bandwidth: int
    tay_required_pair_count: int
    tay_present_pair_count: int
    tay_required_pair_fraction: float
    redundant_network_ready: bool
    tay_network_geometry_ready: bool
    commissioning_displacement_ready: bool
    experimental_trend_ready: bool
    annual_velocity_ready: bool
    release_status: str

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class ReferenceSelection:
    """Provenance and diagnostics for a relative reference patch."""

    x: int
    y: int
    longitude: float
    latitude: float
    patch_radius_pixels: int
    pixel_spacing_m: tuple[float, float]
    candidate_count: int
    minimum_pair_coherence: float
    pair_coherence: tuple[float, ...]
    maximum_local_phase_std_rad: float
    distance_to_mapped_canal_m: float
    distance_to_mapped_peat_m: float
    dominant_connected_components: tuple[int, ...]
    reference_status: str

    def to_dict(self) -> dict:
        return asdict(self)


def _decode_dates(date_dataset: np.ndarray) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    for reference, secondary in date_dataset:
        ref = reference.decode() if isinstance(reference, bytes) else str(reference)
        sec = secondary.decode() if isinstance(secondary, bytes) else str(secondary)
        pairs.append((ref, sec))
    return pairs


def _count_bridge_edges(adjacency: Mapping[str, set[str]]) -> int:
    """Count graph edges whose removal disconnects the acquisition network."""
    discovery: dict[str, int] = {}
    low: dict[str, int] = {}
    parent: dict[str, str] = {}
    timestamp = 0
    bridge_count = 0

    def visit(node: str) -> None:
        nonlocal timestamp, bridge_count
        timestamp += 1
        discovery[node] = timestamp
        low[node] = timestamp
        for neighbour in sorted(adjacency[node]):
            if neighbour not in discovery:
                parent[neighbour] = node
                visit(neighbour)
                low[node] = min(low[node], low[neighbour])
                if low[neighbour] > discovery[node]:
                    bridge_count += 1
            elif parent.get(node) != neighbour:
                low[node] = min(low[node], discovery[neighbour])

    for node in sorted(adjacency):
        if node not in discovery:
            visit(node)
    return bridge_count


def assess_temporal_readiness(stack_path: Path) -> TemporalReadiness:
    """Assess what the current GUNW network can support scientifically.

    The annual gate is deliberately conservative for peatlands because a short
    interval can be dominated by seasonal hydrology and atmospheric residuals.
    """
    import h5py

    with h5py.File(stack_path, "r") as h5:
        pairs = _decode_dates(h5["date"][:])

    acquisitions = sorted({date for pair in pairs for date in pair})
    if not acquisitions:
        raise ValueError(f"No acquisitions found in {stack_path}")

    adjacency = {date: set() for date in acquisitions}
    for reference, secondary in pairs:
        adjacency[reference].add(secondary)
        adjacency[secondary].add(reference)
    unique_edges = {
        tuple(sorted((reference, secondary)))
        for reference, secondary in pairs
        if reference != secondary
    }

    unvisited = set(acquisitions)
    component_count = 0
    while unvisited:
        component_count += 1
        stack = [unvisited.pop()]
        while stack:
            node = stack.pop()
            neighbours = adjacency[node] & unvisited
            unvisited -= neighbours
            stack.extend(neighbours)

    first = datetime.strptime(acquisitions[0], "%Y%m%d")
    last = datetime.strptime(acquisitions[-1], "%Y%m%d")
    span_days = (last - first).days
    cycle_count = max(
        0, len(unique_edges) - len(acquisitions) + component_count
    )
    minimum_degree = min(len(adjacency[date]) for date in acquisitions)
    bridge_count = _count_bridge_edges(adjacency)
    acquisition_rank = {
        acquisition: index for index, acquisition in enumerate(acquisitions)
    }
    maximum_connection_bandwidth = max(
        (
            abs(
                acquisition_rank[secondary]
                - acquisition_rank[reference]
            )
            for reference, secondary in unique_edges
        ),
        default=0,
    )
    complete_connection_bandwidth = 0
    for lag in range(1, len(acquisitions)):
        lag_edges = {
            (acquisitions[index], acquisitions[index + lag])
            for index in range(len(acquisitions) - lag)
        }
        if not lag_edges.issubset(unique_edges):
            break
        complete_connection_bandwidth = lag

    tay_target_bandwidth = min(20, len(acquisitions) - 1)
    tay_required_edges = {
        (acquisitions[index], acquisitions[index + lag])
        for lag in range(1, tay_target_bandwidth + 1)
        for index in range(len(acquisitions) - lag)
    }
    tay_present_pair_count = len(tay_required_edges & unique_edges)
    tay_required_pair_count = len(tay_required_edges)
    tay_required_pair_fraction = (
        tay_present_pair_count / tay_required_pair_count
        if tay_required_pair_count
        else 0.0
    )
    redundant_network_ready = (
        len(acquisitions) >= 3
        and component_count == 1
        and cycle_count >= 1
        and bridge_count == 0
        and minimum_degree >= 2
    )

    commissioning_ready = (
        len(acquisitions) >= 2 and len(pairs) >= 1 and component_count == 1
    )
    experimental_ready = (
        len(acquisitions) >= 12
        and span_days >= 180
        and redundant_network_ready
    )
    annual_ready = (
        len(acquisitions) >= 24
        and span_days >= 365
        and redundant_network_ready
    )
    tay_network_ready = (
        annual_ready
        and complete_connection_bandwidth >= tay_target_bandwidth
        and tay_present_pair_count == tay_required_pair_count
    )

    if tay_network_ready:
        status = "tay_network_candidate_requires_correction_and_quality_gates"
    elif annual_ready:
        status = "annual_velocity_candidate_requires_quality_gates"
    elif experimental_ready:
        status = "experimental_trend_only"
    elif commissioning_ready:
        status = "relative_displacement_commissioning_only"
    else:
        status = "insufficient_temporal_network"

    return TemporalReadiness(
        pair_count=len(pairs),
        acquisition_count=len(acquisitions),
        first_date=acquisitions[0],
        last_date=acquisitions[-1],
        span_days=span_days,
        connected_components=component_count,
        network_cycle_count=cycle_count,
        minimum_acquisition_degree=minimum_degree,
        bridge_edge_count=bridge_count,
        maximum_connection_bandwidth=maximum_connection_bandwidth,
        complete_connection_bandwidth=complete_connection_bandwidth,
        tay_target_bandwidth=tay_target_bandwidth,
        tay_required_pair_count=tay_required_pair_count,
        tay_present_pair_count=tay_present_pair_count,
        tay_required_pair_fraction=tay_required_pair_fraction,
        redundant_network_ready=redundant_network_ready,
        tay_network_geometry_ready=tay_network_ready,
        commissioning_displacement_ready=commissioning_ready,
        experimental_trend_ready=experimental_ready,
        annual_velocity_ready=annual_ready,
        release_status=status,
    )


def _mintpy_transform(attrs: Mapping[str, object]) -> Affine:
    return Affine(
        float(attrs["X_STEP"]),
        0.0,
        float(attrs["X_FIRST"]),
        0.0,
        float(attrs["Y_STEP"]),
        float(attrs["Y_FIRST"]),
    )


def _nan_local_std(data: np.ndarray, size: int) -> np.ndarray:
    valid = np.isfinite(data)
    values = np.where(valid, data, 0.0).astype(np.float32)
    fraction = uniform_filter(valid.astype(np.float32), size=size, mode="nearest")
    mean = np.divide(
        uniform_filter(values, size=size, mode="nearest"),
        fraction,
        out=np.full(data.shape, np.nan, dtype=np.float32),
        where=fraction > 0,
    )
    second = np.divide(
        uniform_filter(values * values, size=size, mode="nearest"),
        fraction,
        out=np.full(data.shape, np.nan, dtype=np.float32),
        where=fraction > 0,
    )
    return np.sqrt(np.maximum(second - mean * mean, 0.0))


def _dominant_positive_label(array: np.ndarray) -> int:
    values = array[np.isfinite(array) & (array > 0)].astype(np.int64)
    if values.size == 0:
        raise ValueError("Connected-component layer has no positive labels")
    return int(np.bincount(values).argmax())


def _rasterize_geojson(
    path: Path,
    shape: tuple[int, int],
    transform: Affine,
    target_crs: str,
) -> np.ndarray:
    with Path(path).open() as stream:
        collection = json.load(stream)
    geometries = [
        transform_geom("EPSG:4326", target_crs, feature["geometry"])
        for feature in collection.get("features", [])
        if feature.get("geometry")
    ]
    if not geometries:
        return np.zeros(shape, dtype=bool)
    return rasterize(
        [(geometry, 1) for geometry in geometries],
        out_shape=shape,
        transform=transform,
        fill=0,
        dtype=np.uint8,
    ).astype(bool)


def _reproject_binary_raster(
    path: Path,
    shape: tuple[int, int],
    transform: Affine,
    target_crs: str,
) -> np.ndarray:
    import rasterio

    output = np.zeros(shape, dtype=np.uint8)
    with rasterio.open(path) as source:
        reproject(
            rasterio.band(source, 1),
            output,
            src_transform=source.transform,
            src_crs=source.crs,
            src_nodata=source.nodata,
            dst_transform=transform,
            dst_crs=target_crs,
            dst_nodata=0,
            resampling=Resampling.max,
        )
    return output.astype(bool)


def select_reference_patch(
    stack_path: Path,
    geometry_path: Path,
    *,
    peat_geojson: Optional[Path] = None,
    canal_mask: Optional[Path] = None,
    minimum_coherence: float = 0.70,
    maximum_local_phase_std_rad: float = 0.80,
    minimum_canal_distance_m: float = 1_000.0,
    minimum_peat_distance_m: float = 500.0,
    minimum_edge_distance_m: float = 1_000.0,
    patch_radius_pixels: int = 2,
    local_window_pixels: int = 7,
) -> ReferenceSelection:
    """Select a high-coherence relative reference away from mapped peat/canals.

    This is a relative datum selection, not proof that the patch is absolutely
    stable. A GNSS or surveyed stable-ground control point is still required
    for an absolute displacement claim.
    """
    import h5py

    with h5py.File(stack_path, "r") as stack:
        attrs = dict(stack.attrs)
        coherence = stack["coherence"][:].astype(np.float32)
        phase = stack["unwrapPhase"][:].astype(np.float32)
        connected = stack["connectComponent"][:]
    with h5py.File(geometry_path, "r") as geometry:
        land = geometry["waterMask"][:].astype(bool)

    if phase.shape != coherence.shape or phase.shape != connected.shape:
        raise ValueError("NISAR phase, coherence, and connected-component shapes differ")
    if phase.shape[1:] != land.shape:
        raise ValueError("NISAR geometry mask does not match interferogram grid")

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
        raise ValueError(
            "No NISAR reference candidate passed the configured coherence, "
            "roughness, land, peat, canal, and edge gates"
        )

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
    selected_y, selected_x = np.unravel_index(np.argmax(score), shape)

    projected_x, projected_y = xy(
        transform, selected_y, selected_x, offset="center"
    )
    longitude, latitude = Transformer.from_crs(
        target_crs, 4326, always_xy=True
    ).transform(projected_x, projected_y)

    return ReferenceSelection(
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
            float(value) for value in coherence[:, selected_y, selected_x]
        ),
        maximum_local_phase_std_rad=float(
            local_roughness[selected_y, selected_x]
        ),
        distance_to_mapped_canal_m=float(
            canal_distance[selected_y, selected_x]
        ),
        distance_to_mapped_peat_m=float(peat_distance[selected_y, selected_x]),
        dominant_connected_components=dominant_labels,
        reference_status="relative_reference_not_ground_control",
    )


def _validate_correction_stack(
    source: object,
    correction: object,
    correction_name: str,
) -> None:
    if source["unwrapPhase"].shape != correction["unwrapPhase"].shape:
        raise ValueError(f"{correction_name} correction shape does not match source")
    if not np.array_equal(source["date"][:], correction["date"][:]):
        raise ValueError(f"{correction_name} correction dates do not match source")
    for key in ("EPSG", "X_FIRST", "Y_FIRST", "X_STEP", "Y_STEP", "WAVELENGTH"):
        if str(source.attrs.get(key)) != str(correction.attrs.get(key)):
            raise ValueError(
                f"{correction_name} correction metadata {key} does not match source"
            )


def prepare_phase_variant(
    source_stack: Path,
    output_stack: Path,
    reference: ReferenceSelection,
    *,
    correction_stacks: Optional[Mapping[str, Path]] = None,
    minimum_reference_coherence: float = 0.70,
) -> Path:
    """Create a patch-referenced raw or corrected MintPy interferogram stack."""
    import h5py

    correction_stacks = correction_stacks or {}
    allowed = {"ionosphere", "troposphere", "solid_earth_tides"}
    unknown = sorted(set(correction_stacks) - allowed)
    if unknown:
        raise ValueError(f"Unknown NISAR phase corrections: {unknown}")

    source_stack = Path(source_stack)
    output_stack = Path(output_stack)
    output_stack.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_stack, output_stack)

    correction_files = {
        name: h5py.File(path, "r") for name, path in correction_stacks.items()
    }
    try:
        with h5py.File(source_stack, "r") as source:
            for name, correction in correction_files.items():
                _validate_correction_stack(source, correction, name)

            with h5py.File(output_stack, "r+") as output:
                phase_dataset = output["unwrapPhase"]
                coherence = source["coherence"]
                connected = source["connectComponent"]
                radius = reference.patch_radius_pixels
                row_slice = slice(reference.y - radius, reference.y + radius + 1)
                col_slice = slice(reference.x - radius, reference.x + radius + 1)
                reference_values: list[float] = []

                for pair_index in range(phase_dataset.shape[0]):
                    phase = source["unwrapPhase"][pair_index].astype(np.float32)
                    valid = np.isfinite(phase)
                    for correction in correction_files.values():
                        correction_phase = correction["unwrapPhase"][
                            pair_index
                        ].astype(np.float32)
                        valid &= np.isfinite(correction_phase)
                        phase -= correction_phase

                    patch_valid = (
                        valid[row_slice, col_slice]
                        & (
                            coherence[pair_index, row_slice, col_slice]
                            >= minimum_reference_coherence
                        )
                        & (connected[pair_index, row_slice, col_slice] > 0)
                    )
                    patch_values = phase[row_slice, col_slice][patch_valid]
                    minimum_pixels = max(3, (2 * radius + 1) ** 2 // 2)
                    if patch_values.size < minimum_pixels:
                        raise ValueError(
                            f"Reference patch has only {patch_values.size} valid "
                            f"pixels for pair {pair_index}; need {minimum_pixels}"
                        )
                    reference_value = float(np.median(patch_values))
                    reference_values.append(reference_value)
                    phase[valid] -= reference_value
                    phase[~valid] = np.nan
                    phase_dataset[pair_index] = phase

                output.attrs["PROCESSOR"] = "nisar"
                output.attrs["REF_X"] = reference.x
                output.attrs["REF_Y"] = reference.y
                output.attrs["REF_LON"] = reference.longitude
                output.attrs["REF_LAT"] = reference.latitude
                output.attrs["PEATGUARD_REFERENCE_METHOD"] = (
                    f"median_{2 * radius + 1}x{2 * radius + 1}_patch"
                )
                output.attrs["PEATGUARD_REFERENCE_STATUS"] = (
                    reference.reference_status
                )
                output.attrs["PEATGUARD_REFERENCE_PHASE_RAD"] = json.dumps(
                    reference_values
                )
                output.attrs["PEATGUARD_CORRECTIONS_APPLIED"] = ",".join(
                    correction_files
                ) or "none"
                output.attrs["PEATGUARD_CORRECTION_SIGN"] = (
                    CORRECTION_SIGN_CONVENTION
                )
                output.attrs["PEATGUARD_MINTPY_COMMIT"] = MINTPY_NISAR_COMMIT
    finally:
        for correction in correction_files.values():
            correction.close()

    return output_stack
