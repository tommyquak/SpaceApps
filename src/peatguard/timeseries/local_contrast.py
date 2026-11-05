"""Conservative, relative phase contrasts with fixed spatial support.

[S] This is an engineered local-observable adaptation, not a reproduction of
Tay et al. (2025), "Satellite radar advances carbon emissions accountability
over tropical peat", or the ambiguity correction of Yunjun, Fattahi and
Amelung (2019), "Small baseline InSAR time series analysis: Unwrapping error
correction and noise reduction". A target-minus-background contrast cancels
terms shared by the two regions, including a common component's arbitrary
phase offset. Background motion and spatially varying atmospheric or
unwrapping errors do not cancel in general.

The caller must supply aligned, unwrapped phase in radians and apply product
validity masks and documented component fill values before calling. Connected
islands are checked within the supplied window, not beyond it. Neither this
conservative support check nor network self-consistency proves correct internal
unwrapping, a stable datum, ground accuracy, or annual-velocity readiness.
"""

from __future__ import annotations

from datetime import datetime
from numbers import Integral
from typing import Sequence

import numpy as np
from scipy.ndimage import label


def _real_array(value: object, name: str) -> np.ndarray:
    array = np.asarray(value)
    if array.dtype.kind not in "buif":
        raise ValueError(f"{name} must contain real numeric values")
    return array


def _mask(value: object, shape: tuple[int, int], name: str) -> np.ndarray:
    array = _real_array(value, name)
    if array.shape != shape or not np.all(np.isin(array, (0, 1))):
        raise ValueError(f"{name} must be a matching boolean or binary mask")
    return array.astype(bool)


def _fraction(value: float, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
        raise ValueError(f"{name} must be a finite number in [0, 1]")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite number in [0, 1]") from exc
    if not np.isfinite(result) or not 0 <= result <= 1:
        raise ValueError(f"{name} must be a finite number in [0, 1]")
    return result


def _spatial_islands(components: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Split each original component into deterministic four-connected islands."""
    islands = np.zeros(components.shape, dtype=np.int64)
    cross = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=bool)
    for pair_index in range(len(components)):
        next_id = 0
        for component in np.unique(components[pair_index][valid[pair_index]]):
            same = valid[pair_index] & (components[pair_index] == component)
            parts, count = label(same, structure=cross)
            islands[pair_index][same] = parts[same] + next_id
            next_id += count
    return islands


def extract_local_contrast(
    phase: np.ndarray,
    coherence: np.ndarray,
    components: np.ndarray,
    target_mask: np.ndarray,
    background_mask: np.ndarray,
    *,
    min_coherence: float = 0.3,
    min_pixels: int = 10,
    min_target_fraction: float = 0.5,
    min_background_fraction: float = 0.5,
) -> dict:
    """Extract pairwise target-minus-background means using frozen support.

    Shape is ``(pairs, rows, columns)``. Finite zero phase is valid. A pixel
    must pass every pair's coherence and component checks. Each original
    component is further split at invalid pixels into four-connected spatial
    islands. Select the complete island signature with the most target pixels;
    ties use the first target pixel in row-major order, independent of numeric
    component labels. Background pixels must
    share that entire signature. No phase magnitude, desired sign, benchmark,
    or background availability is used to choose a different signature.

    Fractions use the full requested masks, not only surviving valid pixels.
    Insufficient support returns an unavailable report without observations.
    Malformed inputs raise ValueError. Spatial standard deviations use ddof=0
    and are descriptive scatter, not standard errors or observation covariance.
    Thresholds are engineering screening choices, not published accuracy bounds.
    """
    phase = _real_array(phase, "phase").astype(np.float64)
    coherence = _real_array(coherence, "coherence").astype(np.float64)
    components = _real_array(components, "components")
    if phase.ndim != 3 or any(size == 0 for size in phase.shape):
        raise ValueError("phase must have nonempty (pairs, rows, columns) shape")
    if coherence.shape != phase.shape or components.shape != phase.shape:
        raise ValueError("phase, coherence, and components must have identical shapes")
    target = _mask(target_mask, phase.shape[1:], "target_mask")
    background = _mask(background_mask, phase.shape[1:], "background_mask")
    if not target.any() or not background.any():
        raise ValueError("target and background masks must both be nonempty")
    if np.any(target & background):
        raise ValueError("target and background masks must be disjoint")
    min_coherence = _fraction(min_coherence, "min_coherence")
    min_target_fraction = _fraction(min_target_fraction, "min_target_fraction")
    min_background_fraction = _fraction(min_background_fraction, "min_background_fraction")
    if isinstance(min_pixels, (bool, np.bool_)) or not isinstance(min_pixels, Integral):
        raise ValueError("min_pixels must be a positive integer")
    if min_pixels < 1:
        raise ValueError("min_pixels must be a positive integer")

    valid = (
        np.isfinite(phase)
        & np.isfinite(coherence)
        & (coherence >= min_coherence)
        & (coherence <= 1)
        & np.isfinite(components)
        & (components > 0)
        & (components == np.floor(components))
    )
    complete = np.all(valid, axis=0)
    islands = _spatial_islands(components, valid)
    target_complete = target & complete
    background_complete = background & complete
    requested_target_count = int(target.sum())
    requested_background_count = int(background.sum())
    report = {
        "available": False,
        "status": "unavailable",
        "reasons": [],
        "contrast_definition": "target_minus_background_unwrapped_phase_rad",
        "contrast_rad": None,
        "per_pair": [],
        "component_labels": None,
        "spatial_island_ids": None,
        "support_shape": list(phase.shape[1:]),
        "target_support_flat_indices": [],
        "background_support_flat_indices": [],
        "pair_count": int(len(phase)),
        "requested_target_count": requested_target_count,
        "requested_background_count": requested_background_count,
        "complete_target_candidate_count": int(target_complete.sum()),
        "complete_background_candidate_count": int(background_complete.sum()),
        "target_count": 0,
        "background_count": 0,
        "target_fraction": 0.0,
        "background_fraction": 0.0,
        "thresholds": {
            "min_coherence": min_coherence,
            "min_pixels": int(min_pixels),
            "min_target_fraction": min_target_fraction,
            "min_background_fraction": min_background_fraction,
        },
        "selection": "largest_complete_target_island_signature_then_first_target_pixel",
        "support_policy": "same_pixels_and_same_four_connected_island_in_every_pair",
        "spatial_scatter_is_uncertainty": False,
        "limitations": [
            "Relative to the chosen background, which may move.",
            "Connectivity is checked only within the supplied window.",
            "Spatial islands do not prove internal unwrapping correctness.",
            "Spatially varying propagation errors may remain in the contrast.",
        ],
        "release_ready": False,
        "annual_velocity_ready": False,
    }
    if not target_complete.any():
        report["reasons"] = ["no_complete_target_support"]
        return report

    signatures, first_indices, counts = np.unique(
        islands[:, target_complete].T, axis=0, return_index=True, return_counts=True
    )
    # Component numbers are arbitrary. Use geometry, never numeric label order,
    # to break a support-count tie between distinct complete signatures.
    tied = np.flatnonzero(counts == counts.max())
    chosen = signatures[tied[int(np.argmin(first_indices[tied]))]]
    same_signature = complete & np.all(islands == chosen[:, None, None], axis=0)
    target_support = target & same_signature
    background_support = background & same_signature
    target_count = int(target_support.sum())
    background_count = int(background_support.sum())
    target_fraction = target_count / requested_target_count
    background_fraction = background_count / requested_background_count
    original_labels = components[:, target_support][:, 0]
    report.update(
        component_labels=[int(value) for value in original_labels],
        spatial_island_ids=[int(value) for value in chosen],
        target_support_flat_indices=np.flatnonzero(target_support).tolist(),
        background_support_flat_indices=np.flatnonzero(background_support).tolist(),
        target_count=target_count,
        background_count=background_count,
        target_fraction=target_fraction,
        background_fraction=background_fraction,
    )
    reasons = []
    if target_count < min_pixels:
        reasons.append("insufficient_target_pixels")
    if background_count == 0:
        reasons.append("no_matching_background_support")
    elif background_count < min_pixels:
        reasons.append("insufficient_background_pixels")
    if target_fraction < min_target_fraction:
        reasons.append("insufficient_target_fraction")
    if background_fraction < min_background_fraction:
        reasons.append("insufficient_background_fraction")
    if reasons:
        report["reasons"] = reasons
        return report

    target_values = phase[:, target_support]
    background_values = phase[:, background_support]
    # Center at one shared anchor before differencing to avoid needless loss of
    # precision from a large arbitrary per-pair unwrapping/reference offset.
    anchor = target_values[:, :1]
    target_centered = target_values - anchor
    background_centered = background_values - anchor
    contrast = target_centered.mean(axis=1) - background_centered.mean(axis=1)
    target_mean = target_centered.mean(axis=1) + anchor[:, 0]
    background_mean = background_centered.mean(axis=1) + anchor[:, 0]
    target_std = target_centered.std(axis=1, ddof=0)
    background_std = background_centered.std(axis=1, ddof=0)
    if not np.isfinite([contrast, target_mean, background_mean, target_std, background_std]).all():
        raise ValueError("phase magnitudes overflow the contrast calculation")
    report.update(
        available=True,
        status="relative_contrast_only",
        contrast_rad=contrast.tolist(),
        per_pair=[
            {
                "pair_index": index,
                "component_label": int(original_labels[index]),
                "spatial_island_id": int(chosen[index]),
                "target_count": target_count,
                "background_count": background_count,
                "target_mean_rad": float(target_mean[index]),
                "background_mean_rad": float(background_mean[index]),
                "contrast_rad": float(contrast[index]),
                "target_spatial_std_rad": float(target_std[index]),
                "background_spatial_std_rad": float(background_std[index]),
            }
            for index in range(len(phase))
        ],
    )
    return report


def _date(value: object) -> str:
    if not isinstance(value, str) or len(value) != 8 or not value.isascii() or not value.isdigit():
        raise ValueError("dates must be YYYYMMDD strings")
    try:
        datetime.strptime(value, "%Y%m%d")
    except ValueError as exc:
        raise ValueError(f"invalid calendar date: {value}") from exc
    return value


def _cycle_matrix(pairs: list[tuple[str, str]], dates: list[str]) -> np.ndarray:
    """A deterministic fundamental cycle basis, oriented by the closing edge."""
    tree: dict[str, list[tuple[str, int, int]]] = {date: [] for date in dates}
    cycles = []
    for edge_index, (first, second) in enumerate(pairs):
        # A path in the current spanning forest makes this edge a cycle.
        paths = {first: []}
        queue = [first]
        for node in queue:
            for neighbour, index, orientation in tree[node]:
                if neighbour not in paths:
                    paths[neighbour] = paths[node] + [(index, orientation)]
                    queue.append(neighbour)
        if second in paths:
            cycle = np.zeros(len(pairs), dtype=int)
            cycle[edge_index] = 1
            for index, orientation in paths[second]:
                cycle[index] = -orientation
            cycles.append(cycle)
        else:
            tree[first].append((second, edge_index, 1))
            tree[second].append((first, edge_index, -1))
    return np.asarray(cycles, dtype=int).reshape(-1, len(pairs))


def solve_contrast_network(
    date_pairs: Sequence[tuple[str, str]],
    observations: Sequence[float],
    *,
    observation_covariance: np.ndarray | None = None,
) -> dict:
    """Solve ``observation(ref, sec) = epoch(sec) - epoch(ref)`` in radians.

    The earliest epoch is fixed at zero, a temporal gauge, not a stable spatial
    reference. Edges must be strictly earlier-to-later and unique. A supplied
    covariance must be finite, symmetric and positive definite, in radians
    squared, for these actual contrast observations. GLS uses its full matrix.
    Covariance propagation is ``L @ C @ L.T`` through the reported operator.
    Without a supplied covariance, unweighted least squares reports no inferred
    uncertainty. A tree fits exactly but supplies no redundant closure test.
    """
    try:
        pair_iterator = iter(date_pairs)
    except TypeError as exc:
        raise ValueError("date_pairs must be a sequence of date pairs") from exc
    pairs = []
    seen = set()
    for pair in pair_iterator:
        if isinstance(pair, (str, bytes)):
            raise ValueError("each date pair must contain exactly two dates")
        try:
            if len(pair) != 2:
                raise ValueError("each date pair must contain exactly two dates")
            first, second = (_date(value) for value in pair)
        except TypeError as exc:
            raise ValueError("each date pair must contain exactly two dates") from exc
        if first >= second:
            raise ValueError("date pairs must be strictly earlier-to-later, never same or reversed")
        if (first, second) in seen:
            raise ValueError("duplicate date pair")
        seen.add((first, second))
        pairs.append((first, second))
    if not pairs:
        raise ValueError("at least one date pair is required")
    y = _real_array(observations, "observations").astype(np.float64)
    if y.shape != (len(pairs),) or not np.isfinite(y).all():
        raise ValueError("observations must be a finite vector matching date_pairs")
    dates = sorted({date for pair in pairs for date in pair})
    indices = {date: index for index, date in enumerate(dates)}
    incidence = np.zeros((len(pairs), len(dates)))
    for row, (first, second) in enumerate(pairs):
        incidence[row, indices[first]] = -1
        incidence[row, indices[second]] = 1
    design = incidence[:, 1:]
    rank = int(np.linalg.matrix_rank(design))
    if rank != len(dates) - 1:
        raise ValueError("date network is disconnected")

    covariance = None
    if observation_covariance is None:
        reduced_operator = np.linalg.lstsq(design, np.eye(len(pairs)), rcond=None)[0]
        solver = "unweighted_least_squares"
    else:
        covariance = _real_array(observation_covariance, "observation_covariance").astype(
            np.float64
        )
        if covariance.shape != (len(pairs), len(pairs)) or not np.isfinite(covariance).all():
            raise ValueError("observation_covariance must be a finite square matrix matching pairs")
        scale = float(np.max(np.abs(covariance)))
        if not np.allclose(covariance, covariance.T, rtol=0, atol=scale * 1e-12):
            raise ValueError("observation_covariance must be symmetric positive definite")
        covariance = (covariance + covariance.T) / 2
        try:
            cholesky = np.linalg.cholesky(covariance)
        except np.linalg.LinAlgError as exc:
            raise ValueError("observation_covariance must be symmetric positive definite") from exc
        whitening = np.linalg.solve(cholesky, np.eye(len(pairs)))
        whitened_design = whitening @ design
        reduced_operator, _, weighted_rank, _ = np.linalg.lstsq(
            whitened_design, whitening, rcond=None
        )
        if weighted_rank != rank:
            raise ValueError("observation_covariance makes the solve numerically rank deficient")
        solver = "generalized_least_squares"

    operator = np.vstack([np.zeros(len(pairs)), reduced_operator])
    epochs = operator @ y
    predicted = incidence @ epochs
    residual = y - predicted
    epoch_covariance = None if covariance is None else operator @ covariance @ operator.T
    cycles = _cycle_matrix(pairs, dates)
    closure = cycles @ y
    values = [operator, epochs, predicted, residual, closure]
    if epoch_covariance is not None:
        values.append(epoch_covariance)
    if any(not np.isfinite(value).all() for value in values):
        raise ValueError("network inputs overflow the numerical calculation")
    redundancy = len(pairs) - rank
    return {
        "status": "relative_contrast_network_only",
        "dates": dates,
        "date_pairs": [list(pair) for pair in pairs],
        "gauge_date": dates[0],
        "gauge_definition": "earliest_epoch_contrast_fixed_to_zero_not_absolute_motion",
        "solver": solver,
        "epoch_contrast_rad": epochs.tolist(),
        "observations_rad": y.tolist(),
        "predicted_observations_rad": predicted.tolist(),
        "residual_rad": residual.tolist(),
        "linear_operator": operator.tolist(),
        "epoch_covariance_rad2": None if epoch_covariance is None else epoch_covariance.tolist(),
        "pair_count": len(pairs),
        "acquisition_count": len(dates),
        "rank": rank,
        "redundancy": redundancy,
        "cycle_count": len(cycles),
        "cycle_matrix": cycles.tolist(),
        "cycle_closure_rad": closure.tolist(),
        "has_redundant_closure_test": redundancy > 0,
        "uncertainty_source": (
            "not_supplied" if covariance is None else "supplied_observation_covariance"
        ),
        "limitations": [
            "All epochs measure target motion relative to potentially moving background.",
            "Network closure is internal consistency, not independent accuracy.",
            "A tree has no redundant observations even when residuals are zero.",
            "No spatial-scatter-derived or residual-derived uncertainty is assumed.",
        ],
        "release_ready": False,
        "annual_velocity_ready": False,
    }
