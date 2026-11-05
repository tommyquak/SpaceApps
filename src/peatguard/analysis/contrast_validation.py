"""Spatially grouped diagnostics for paired satellite contrasts.

[S] Spatial groups reduce pseudoreplication (Roberts et al., 2017,
doi:10.1111/ecog.02881). These diagnostics never certify displacement,
annual velocity, ground accuracy, or restoration attribution.
"""
from __future__ import annotations

import numpy as np


def paired_metrics(first: np.ndarray, second: np.ndarray) -> dict:
    first, second = np.asarray(first, float), np.asarray(second, float)
    if first.ndim != 1 or first.shape != second.shape or first.size < 2:
        raise ValueError("paired observations must be matching vectors with at least two entries")
    if not np.isfinite(first).all() or not np.isfinite(second).all():
        raise ValueError("paired observations must all be finite")
    correlation = None
    if np.std(first) > 0 and np.std(second) > 0:
        correlation = float(np.corrcoef(first, second)[0, 1])
    nonzero = (first != 0) & (second != 0)
    return {
        "n": int(first.size),
        "pearson_r": correlation,
        "median_difference_mm": float(np.median(first - second)),
        "rmse_difference_mm": float(np.sqrt(np.mean((first - second) ** 2))),
        "same_sign_fraction_nonzero": (
            float(np.mean(np.sign(first[nonzero]) == np.sign(second[nonzero])))
            if nonzero.any() else None
        ),
        "nonzero_pairs": int(nonzero.sum()),
    }


def spatial_group_diagnostics(
    first: np.ndarray, second: np.ndarray, xy_metres: np.ndarray,
    *, block_metres: float = 5000, draws: int = 2000, seed: int = 7,
    min_groups: int = 8,
) -> dict:
    """Bootstrap whole spatial groups, retaining pairing and unequal group sizes.

    Grid groups use a fixed projected-coordinate origin, never outcome values.
    min_groups is an engineering screening choice, not a universal theorem.
    Block-size sensitivity is required because no block size proves independence.
    """
    first, second = np.asarray(first, float), np.asarray(second, float)
    xy = np.asarray(xy_metres, float)
    metrics = paired_metrics(first, second)
    if xy.shape != (first.size, 2) or not np.isfinite(xy).all():
        raise ValueError("coordinates must be finite projected metre pairs")
    if not np.isfinite(block_metres) or block_metres <= 0 or draws < 100 or min_groups < 3:
        raise ValueError("invalid spatial bootstrap parameters")
    groups, inverse = np.unique(np.floor(xy / block_metres).astype(np.int64), axis=0,
                                return_inverse=True)
    members = [np.flatnonzero(inverse == i) for i in range(len(groups))]
    report = {
        "block_metres": float(block_metres), "groups": len(groups),
        "group_sizes": [int(len(m)) for m in members],
        "metrics": metrics, "bootstrap_draws_requested": draws,
        "bootstrap_valid_correlation_draws": 0, "pearson_ci95": None,
        "median_difference_ci95_mm": None, "leave_one_group_out_r_range": None,
        "status": "insufficient_spatial_groups", "release_ready": False,
    }
    if len(groups) < min_groups:
        return report
    rng = np.random.default_rng(seed)
    correlations, differences = [], []
    for _ in range(draws):
        idx = np.concatenate([members[i] for i in rng.integers(0, len(groups), len(groups))])
        sample = paired_metrics(first[idx], second[idx])
        if sample["pearson_r"] is not None:
            correlations.append(sample["pearson_r"])
        differences.append(sample["median_difference_mm"])
    report["bootstrap_valid_correlation_draws"] = len(correlations)
    report["median_difference_ci95_mm"] = np.percentile(differences, [2.5, 97.5]).tolist()
    if len(correlations) < 0.95 * draws:
        report["status"] = "degenerate_correlation_resamples"
        return report
    report["pearson_ci95"] = np.percentile(correlations, [2.5, 97.5]).tolist()
    leave_out = []
    for i in range(len(groups)):
        idx = inverse != i
        if idx.sum() >= 2:
            value = paired_metrics(first[idx], second[idx])["pearson_r"]
            if value is not None:
                leave_out.append(value)
    if leave_out:
        report["leave_one_group_out_r_range"] = [min(leave_out), max(leave_out)]
    report["status"] = "descriptive_cross_geometry_evidence_only"
    return report
