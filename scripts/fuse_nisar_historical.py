#!/usr/bin/env python3
"""Build an internal Tay historical and NISAR commissioning fusion audit.

The quantitative output is a Bayesian sensitivity experiment. It must not be
published as current subsidence because the NISAR input spans only a short
commissioning window,
uses a relative non-ground-controlled datum, and has incomplete uncertainty.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject
from scipy.ndimage import label
from scipy.stats import spearmanr

from peatguard.analysis.baseline_audit import hash_file
from peatguard.analysis.historical_nisar import (
    DAYS_PER_YEAR,
    EVIDENCE_LABELS,
    annualize_short_span_displacement,
    bayesian_rate_sensitivity,
    classify_historical_nisar_evidence,
)
from peatguard.export.cog import write_cog

NODATA = -9999.0
TAY_PERIOD = "2018-10-12/2024-08-09"
TAY_PAPER = (
    "Tay et al. (2025), Satellite radar advances carbon emissions "
    "accountability over tropical peat, doi:10.1038/s43247-025-02926-6"
)


def _iso_date_tag(tags: dict[str, Any], key: str) -> str:
    # nisar_period must match the ISO format of historical_period (TAY_PERIOD)
    # on the same tag block, and a missing acquisition date must fail loudly
    # rather than ship a product with an unknown period.
    value = tags.get(key)
    if value is None:
        raise ValueError(f"NISAR displacement has no {key} tag")
    text = str(value)
    for pattern in ("%Y%m%d", "%Y-%m-%d"):
        try:
            return datetime.strptime(text, pattern).date().isoformat()
        except ValueError:
            continue
    raise ValueError(
        f"NISAR displacement {key} tag is not a date: {text!r}"
    )


def _valid(data: np.ndarray, nodata: float | None) -> np.ndarray:
    mask = np.isfinite(data)
    if nodata is not None and np.isfinite(nodata):
        mask &= data != nodata
    return mask


def _read_float(path: Path) -> tuple[np.ndarray, dict[str, Any]]:
    with rasterio.open(path) as source:
        data = source.read(1).astype(np.float64)
        valid = _valid(data, source.nodata)
        data[~valid] = np.nan
        return data, {
            "shape": source.shape,
            "crs": source.crs,
            "transform": source.transform,
            "bounds": tuple(source.bounds),
            "nodata": source.nodata,
            "tags": source.tags(),
            "width": source.width,
            "height": source.height,
        }


def _same_grid(first: dict[str, Any], second: dict[str, Any]) -> bool:
    return (
        first["shape"] == second["shape"]
        and first["crs"] == second["crs"]
        and first["transform"].almost_equals(second["transform"])
    )


def _reproject_path(
    path: Path,
    destination: dict[str, Any],
    *,
    resampling: Resampling,
) -> np.ndarray:
    output = np.full(destination["shape"], np.nan, dtype=np.float64)
    with rasterio.open(path) as source:
        reproject(
            source=rasterio.band(source, 1),
            destination=output,
            src_transform=source.transform,
            src_crs=source.crs,
            src_nodata=source.nodata,
            dst_transform=destination["transform"],
            dst_crs=destination["crs"],
            dst_nodata=np.nan,
            resampling=resampling,
            init_dest_nodata=True,
        )
    return output


def _reproject_array(
    data: np.ndarray,
    source: dict[str, Any],
    destination: dict[str, Any],
    *,
    resampling: Resampling,
) -> np.ndarray:
    output = np.full(destination["shape"], np.nan, dtype=np.float64)
    reproject(
        source=data,
        destination=output,
        src_transform=source["transform"],
        src_crs=source["crs"],
        src_nodata=np.nan,
        dst_transform=destination["transform"],
        dst_crs=destination["crs"],
        dst_nodata=np.nan,
        resampling=resampling,
        init_dest_nodata=True,
    )
    return output


def _summary(data: np.ndarray, mask: np.ndarray | None = None) -> dict[str, Any]:
    valid = np.isfinite(data)
    if mask is not None:
        valid &= mask
    values = data[valid]
    if values.size == 0:
        return {"count": 0}
    return {
        "count": int(values.size),
        "minimum": float(np.min(values)),
        "percentile_02_5": float(np.percentile(values, 2.5)),
        "percentile_25": float(np.percentile(values, 25.0)),
        "median": float(np.median(values)),
        "mean": float(np.mean(values)),
        "standard_deviation": float(np.std(values)),
        "percentile_75": float(np.percentile(values, 75.0)),
        "percentile_97_5": float(np.percentile(values, 97.5)),
        "maximum": float(np.max(values)),
    }


def _correlations(first: np.ndarray, second: np.ndarray) -> dict[str, Any]:
    valid = np.isfinite(first) & np.isfinite(second)
    if valid.sum() < 3:
        return {
            "count": int(valid.sum()),
            "pearson": None,
            "spearman": None,
        }
    x = first[valid]
    y = second[valid]
    return {
        "count": int(valid.sum()),
        "pearson": float(np.corrcoef(x, y)[0, 1]),
        "spearman": float(spearmanr(x, y).statistic),
    }


def _component_summary(mask: np.ndarray) -> dict[str, Any]:
    components, count = label(
        np.asarray(mask, dtype=bool),
        structure=np.ones((3, 3), dtype=np.uint8),
    )
    sizes = np.bincount(components.ravel())[1:]
    return {
        "pixel_count": int(np.sum(mask)),
        "component_count": int(count),
        "maximum_component_pixels": int(sizes.max()) if sizes.size else 0,
        "component_sizes_descending": sorted(
            (int(value) for value in sizes),
            reverse=True,
        ),
        "connectivity": 8,
    }


def _masked_float(data: np.ndarray) -> np.ndarray:
    output = np.full(data.shape, NODATA, dtype=np.float32)
    valid = np.isfinite(data)
    output[valid] = data[valid].astype(np.float32)
    return output


def _write_float(
    data: np.ndarray,
    path: Path,
    grid: dict[str, Any],
    *,
    description: str,
    units: str,
    tags: dict[str, str],
) -> Path:
    output = write_cog(
        _masked_float(data),
        path,
        crs=grid["crs"],
        transform=grid["transform"],
        nodata=NODATA,
        band_names=[description],
        dtype="float32",
    )
    with rasterio.open(output, "r+") as dataset:
        dataset.update_tags(
            processing_status="experimental_historical_nisar_fusion",
            public_layer_recommendation="do_not_publish_as_current_subsidence",
            not_ground_truth="true",
            units=units,
            **tags,
        )
    return output


def _write_evidence(
    evidence: np.ndarray,
    path: Path,
    grid: dict[str, Any],
    *,
    tags: dict[str, str],
) -> Path:
    output = write_cog(
        evidence.astype(np.uint8),
        path,
        crs=grid["crs"],
        transform=grid["transform"],
        nodata=0,
        band_names=["historical_nisar_evidence_class"],
        dtype="uint8",
        overview_resampling=Resampling.nearest,
    )
    with rasterio.open(output, "r+") as dataset:
        dataset.update_tags(
            processing_status="experimental_historical_nisar_fusion",
            public_layer_recommendation=(
                "context_only_not_current_subsidence"
            ),
            not_ground_truth="true",
            units="class_code",
            class_labels=json.dumps(EVIDENCE_LABELS, sort_keys=True),
            **tags,
        )
    return output


def _input_record(path: Path, metadata: dict[str, Any]) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": hash_file(path, "sha256"),
        "crs": metadata["crs"].to_string(),
        "shape": list(metadata["shape"]),
        "bounds": list(metadata["bounds"]),
        "tags": metadata["tags"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--historical-rate", type=Path, required=True)
    parser.add_argument("--nisar-displacement", type=Path, required=True)
    parser.add_argument("--nisar-uncertainty", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--span-days", type=float)
    parser.add_argument(
        "--historical-measurement-sigma",
        type=float,
        default=3.6,
        help="Tay eight-site rate RMSE in mm/year.",
    )
    parser.add_argument(
        "--historical-change-sigma",
        type=float,
        nargs="+",
        default=[20.0, 50.0, 100.0],
        help="Hypothetical between-period change sigma scenarios in mm/year.",
    )
    parser.add_argument(
        "--reference-change-sigma",
        type=float,
        default=50.0,
    )
    parser.add_argument(
        "--historical-subsidence-threshold",
        type=float,
        default=-5.0,
    )
    parser.add_argument("--significance-sigma", type=float, default=2.0)
    args = parser.parse_args()

    historical, historical_grid = _read_float(args.historical_rate)
    displacement, nisar_grid = _read_float(args.nisar_displacement)
    uncertainty, uncertainty_grid = _read_float(args.nisar_uncertainty)
    if not _same_grid(nisar_grid, uncertainty_grid):
        raise ValueError("NISAR displacement and uncertainty grids differ")

    span_days = args.span_days
    if span_days is None:
        raw_span = nisar_grid["tags"].get("span_days")
        if raw_span is None:
            raise ValueError("NISAR displacement has no span_days tag")
        span_days = float(raw_span)
    scenarios = sorted(set(float(value) for value in args.historical_change_sigma))
    if args.reference_change_sigma not in scenarios:
        scenarios.append(float(args.reference_change_sigma))
        scenarios.sort()

    historical_on_nisar = _reproject_path(
        args.historical_rate,
        nisar_grid,
        resampling=Resampling.bilinear,
    )
    nisar_rate, nisar_rate_sigma = annualize_short_span_displacement(
        displacement,
        uncertainty,
        span_days,
    )
    evidence, z_score = classify_historical_nisar_evidence(
        historical_on_nisar,
        displacement,
        uncertainty,
        historical_subsidence_threshold_mm_per_year=(
            args.historical_subsidence_threshold
        ),
        significance_sigma=args.significance_sigma,
    )

    scenario_results = {}
    reference_result = None
    for process_sigma in scenarios:
        result = bayesian_rate_sensitivity(
            historical_on_nisar,
            nisar_rate,
            nisar_rate_sigma,
            historical_measurement_sigma_mm_per_year=(
                args.historical_measurement_sigma
            ),
            historical_change_sigma_mm_per_year=process_sigma,
        )
        if process_sigma == args.reference_change_sigma:
            reference_result = result
        valid = np.isfinite(result.posterior_rate_mm_per_year)
        difference = result.posterior_rate_mm_per_year - historical_on_nisar
        historical_state = (
            historical_on_nisar < args.historical_subsidence_threshold
        )
        posterior_state = (
            result.posterior_rate_mm_per_year
            < args.historical_subsidence_threshold
        )
        scenario_results[str(process_sigma)] = {
            "assumptions": result.assumptions(),
            "valid_fused_pixels_on_80m_posted_grid": int(valid.sum()),
            "nisar_weight_upper_bound": _summary(
                result.nisar_weight_upper_bound
            ),
            "posterior_minus_historical_mm_per_year": _summary(
                difference
            ),
            "posterior_rate_mm_per_year": _summary(
                result.posterior_rate_mm_per_year
            ),
            "historical_subsidence_class_changes": int(
                np.sum(valid & (historical_state != posterior_state))
            ),
            "posterior_vs_historical": _correlations(
                result.posterior_rate_mm_per_year,
                historical_on_nisar,
            ),
        }
    if reference_result is None:
        raise RuntimeError("Reference process-change scenario was not evaluated")

    # Evaluate comparisons at Tay's native support instead of counting
    # upsampled 80 m pixels as independent historical observations.
    displacement_on_historical = _reproject_array(
        displacement,
        nisar_grid,
        historical_grid,
        resampling=Resampling.average,
    )
    uncertainty_on_historical = _reproject_array(
        uncertainty,
        nisar_grid,
        historical_grid,
        resampling=Resampling.max,
    )
    expected_historical_displacement = historical * span_days / DAYS_PER_YEAR
    native_overlap = (
        np.isfinite(historical)
        & np.isfinite(displacement_on_historical)
        & np.isfinite(uncertainty_on_historical)
        & (uncertainty_on_historical > 0)
    )
    native_significant_downward = native_overlap & (
        displacement_on_historical
        <= -args.significance_sigma * uncertainty_on_historical
    )
    native_significant_upward = native_overlap & (
        displacement_on_historical
        >= args.significance_sigma * uncertainty_on_historical
    )
    native_difference = (
        displacement_on_historical - expected_historical_displacement
    )
    overlap_count = int(native_overlap.sum())
    rmse = (
        float(np.sqrt(np.mean(native_difference[native_overlap] ** 2)))
        if overlap_count
        else None
    )
    sign_agreement = (
        float(
            np.mean(
                np.sign(displacement_on_historical[native_overlap])
                == np.sign(expected_historical_displacement[native_overlap])
            )
        )
        if overlap_count
        else None
    )

    class_counts = {
        str(code): {
            "label": label,
            "pixels_on_80m_posted_grid": int(np.sum(evidence == code)),
        }
        for code, label in EVIDENCE_LABELS.items()
    }
    reference_scenario_label = f"{args.reference_change_sigma:g}".replace(
        ".", "p"
    )
    common_tags = {
        "historical_source": TAY_PAPER,
        "historical_period": TAY_PERIOD,
        "nisar_period": (
            f"{_iso_date_tag(nisar_grid['tags'], 'reference_date')}/"
            f"{_iso_date_tag(nisar_grid['tags'], 'secondary_date')}"
        ),
        "nisar_span_days": str(span_days),
        "nisar_product_maturity": nisar_grid["tags"].get(
            "product_maturity", "unknown"
        ),
        "nisar_reference_status": nisar_grid["tags"].get(
            "reference_status", "unknown"
        ),
        "no_fit_to_tay": "true",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    products = {
        "historical_rate_on_nisar_grid": _write_float(
            historical_on_nisar,
            args.output_dir / "tay_historical_rate_on_nisar_80m.tif",
            nisar_grid,
            description="tay_historical_rate_mm_per_year",
            units="mm/year",
            tags={
                **common_tags,
                "effective_support": (
                    "Tay approximately 180m; 80m is output posting only"
                ),
            },
        ),
        "nisar_equivalent_rate_diagnostic": _write_float(
            nisar_rate,
            args.output_dir
            / f"nisar_{int(span_days)}day_equivalent_rate_diagnostic_80m.tif",
            nisar_grid,
            description="nisar_equivalent_rate_not_velocity",
            units="mm/year",
            tags={
                **common_tags,
                "annualization_status": (
                    "unit_conversion_only_not_annual_velocity"
                ),
            },
        ),
        "bayesian_candidate_rate": _write_float(
            reference_result.posterior_rate_mm_per_year,
            args.output_dir
            / (
                "bayesian_candidate_rate_process"
                f"{reference_scenario_label}_diagnostic_80m.tif"
            ),
            nisar_grid,
            description="bayesian_candidate_rate_not_current_velocity",
            units="mm/year",
            tags={
                **common_tags,
                "historical_change_sigma_mm_per_year": str(
                    args.reference_change_sigma
                ),
                "fusion_status": "sensitivity_experiment_only",
            },
        ),
        "nisar_weight_upper_bound": _write_float(
            reference_result.nisar_weight_upper_bound,
            args.output_dir
            / (
                "nisar_weight_upper_bound_process"
                f"{reference_scenario_label}_diagnostic_80m.tif"
            ),
            nisar_grid,
            description="nisar_weight_upper_bound",
            units="fraction",
            tags={
                **common_tags,
                "historical_change_sigma_mm_per_year": str(
                    args.reference_change_sigma
                ),
                "upper_bound_reason": (
                    "NISAR uncertainty excludes several error sources"
                ),
            },
        ),
        "nisar_displacement_z_score": _write_float(
            z_score,
            args.output_dir / "nisar_displacement_z_score_diagnostic_80m.tif",
            nisar_grid,
            description="nisar_displacement_over_partial_uncertainty",
            units="sigma",
            tags={
                **common_tags,
                "uncertainty_scope": (
                    "ionosphere_plus_reference_excludes_other_error_sources"
                ),
            },
        ),
        "evidence_class": _write_evidence(
            evidence,
            args.output_dir / "historical_nisar_evidence_class_80m.tif",
            nisar_grid,
            tags=common_tags,
        ),
    }

    reference_weight = scenario_results[str(args.reference_change_sigma)][
        "nisar_weight_upper_bound"
    ]
    median_weight = float(reference_weight.get("median", np.nan))
    native_significant_count = int(
        native_significant_downward.sum() + native_significant_upward.sum()
    )
    significant_fraction = (
        native_significant_count / overlap_count if overlap_count else 0.0
    )
    dominated = median_weight < 0.05
    report = {
        "schema_version": 1,
        "status": (
            "fusion_dominated_by_historical_prior_not_current_rate"
            if dominated
            else "experimental_fusion_requires_independent_validation"
        ),
        "release_decision": "internal_evidence_fusion_only",
        "current_subsidence_claim": False,
        "method": {
            "historical_source": TAY_PAPER,
            "historical_role": (
                "independent historical benchmark and explicit prior scenario"
            ),
            "nisar_role": (
                f"{int(span_days)}-day commissioning displacement with partial uncertainty"
            ),
            "quantitative_fusion": (
                "independent Gaussian inverse-variance sensitivity test"
            ),
            "categorical_fusion": (
                "historical state crossed with NISAR directional significance"
            ),
            "no_bias_scale_or_datum_fit_to_tay": True,
        },
        "inputs": {
            "historical_rate": _input_record(
                args.historical_rate, historical_grid
            ),
            "nisar_displacement": _input_record(
                args.nisar_displacement, nisar_grid
            ),
            "nisar_uncertainty": _input_record(
                args.nisar_uncertainty, uncertainty_grid
            ),
        },
        "assumptions": {
            "span_days": span_days,
            "days_per_year": DAYS_PER_YEAR,
            "historical_measurement_sigma_mm_per_year": (
                args.historical_measurement_sigma
            ),
            "historical_change_sigma_scenarios_mm_per_year": scenarios,
            "historical_change_sigma_status": (
                "[H] sensitivity parameters, not calibrated observations"
            ),
            "historical_subsidence_threshold_mm_per_year": (
                args.historical_subsidence_threshold
            ),
            "historical_threshold_status": (
                "[H] PeatGuard operational threshold, not prescribed by Tay"
            ),
            "significance_sigma": args.significance_sigma,
            "nisar_uncertainty_is_total": False,
            "absolute_datum_compatible": False,
        },
        "native_historical_support_evaluation": {
            "reason": (
                "Statistics use the native Tay grid to avoid treating "
                "upsampled 80m pixels as independent historical samples."
            ),
            "overlap_cells": overlap_count,
            f"historical_expected_{int(span_days)}day_displacement_mm": _summary(
                expected_historical_displacement,
                native_overlap,
            ),
            f"nisar_observed_{int(span_days)}day_displacement_mm": _summary(
                displacement_on_historical,
                native_overlap,
            ),
            "nisar_partial_uncertainty_mm": _summary(
                uncertainty_on_historical,
                native_overlap,
            ),
            "spatial_pattern": _correlations(
                expected_historical_displacement[native_overlap],
                displacement_on_historical[native_overlap],
            ),
            "rmse_mm": rmse,
            "sign_agreement": sign_agreement,
            "significant_downward_cells": int(
                native_significant_downward.sum()
            ),
            "significant_upward_cells": int(native_significant_upward.sum()),
            "significant_fraction": significant_fraction,
        },
        "bayesian_sensitivity": scenario_results,
        "evidence_classes": class_counts,
        "short_span_directional_clusters": {
            "nisar_only_downward": _component_summary(evidence == 9),
            "nisar_only_upward": _component_summary(evidence == 10),
            "interpretation": (
                "Small isolated clusters remain commissioning flags, not "
                "subsidence or rebound claims."
            ),
        },
        "products": {
            name: str(path.resolve()) for name, path in products.items()
        },
        "appropriateness": {
            "historical_context_map": {
                "verdict": "appropriate",
                "confidence": "[E]",
            },
            "categorical_evidence_map": {
                "verdict": "appropriate for internal commissioning",
                "confidence": "[S]",
            },
            "current_2026_subsidence_rate": {
                "verdict": "not appropriate",
                "confidence": "[E]",
                "reason": (
                    f"A {int(span_days)}-day chain of NISAR acquisitions does not "
                    "support annual velocity, and the datum is not ground "
                    "controlled."
                ),
            },
            "bayesian_candidate_as_public_rate": {
                "verdict": "not appropriate",
                "confidence": "[E/S]",
                "reason": (
                    "The posterior is dominated by the historical prior and "
                    "assumes a short-span displacement represents the same "
                    "latent long-term rate."
                ),
            },
        },
        "limitations": [
            "Tay and NISAR observe different periods.",
            f"Seasonal peat shrink and swell can dominate a {int(span_days)}-day displacement.",
            "The NISAR reference is relative and not tied to stable GNSS.",
            "NISAR uncertainty omits troposphere, unwrapping, orbit, DEM, "
            "horizontal motion, and absolute-datum error.",
            "The 80 m fusion posting does not improve Tay's approximately "
            "180 m effective support.",
            "Tay CC BY-NC 4.0 requires permission for commercial use.",
        ],
    }
    report_path = args.output_dir / "historical_nisar_fusion_audit.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
