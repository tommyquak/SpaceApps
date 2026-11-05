#!/usr/bin/env python3
"""Audit and export compact QA products from a short NISAR GUNW stack."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import rasterio

from peatguard.export.cog import write_cog
from peatguard.timeseries.nisar_commissioning import (
    assess_datum_sensitivity,
    build_commissioning_quality_mask,
    los_to_vertical_displacement,
    propagate_ionosphere_uncertainty,
    rank_reference_patches,
    verify_file_manifest,
    summary_statistics,
    validate_tree_inversion,
)
from peatguard.timeseries.nisar_gunw import (
    ReferenceSelection,
    _mintpy_transform,
    _rasterize_geojson,
    assess_temporal_readiness,
)

NODATA = -9999.0


def _read_timeseries_final(path: Path) -> np.ndarray:
    with h5py.File(path, "r") as h5:
        return h5["timeseries"][-1].astype(np.float64)


def _read_aligned_uncertainty(
    paths: list[Path],
    shape: tuple[int, int],
    transform: rasterio.Affine,
    crs: str,
    *,
    expected_scene_tags: dict[str, str],
) -> np.ndarray:
    arrays: list[np.ndarray] = []
    for path in paths:
        with rasterio.open(path) as source:
            if source.shape != shape:
                raise ValueError(f"Uncertainty shape mismatch: {path}")
            if source.crs != rasterio.crs.CRS.from_string(crs):
                raise ValueError(f"Uncertainty CRS mismatch: {path}")
            if not source.transform.almost_equals(transform):
                raise ValueError(f"Uncertainty transform mismatch: {path}")
            tags = source.tags()
            units = tags.get("units")
            if units != "radians":
                raise ValueError(
                    "Ionosphere uncertainty must be tagged units=radians "
                    f"(got {units!r}): {path}"
                )
            # Source scene tags, where present, must match the values the
            # products will be labelled with; a silent mismatch would
            # publish tifs attributed to the wrong track or frame.
            for key, expected in expected_scene_tags.items():
                actual = tags.get(key)
                if actual is not None and actual != expected:
                    raise ValueError(
                        f"Uncertainty {key} tag {actual!r} disagrees with "
                        f"the requested value {expected!r}: {path}"
                    )
            data = source.read(1).astype(np.float64)
            if source.nodata is not None:
                data[data == source.nodata] = np.nan
            arrays.append(data)
    return np.stack(arrays)


def _masked_float32(data: np.ndarray, mask: np.ndarray) -> np.ndarray:
    output = np.full(data.shape, NODATA, dtype=np.float32)
    valid = mask & np.isfinite(data)
    output[valid] = data[valid].astype(np.float32)
    return output


def _write_product(
    data: np.ndarray,
    mask: np.ndarray,
    path: Path,
    *,
    crs: str,
    transform: rasterio.Affine,
    description: str,
    units: str,
    tags: dict[str, str],
) -> Path:
    output = write_cog(
        _masked_float32(data, mask),
        path,
        crs=crs,
        transform=transform,
        nodata=NODATA,
        band_names=[description],
        dtype="float32",
    )
    with rasterio.open(output, "r+") as dataset:
        dataset.update_tags(
            processing_status="commissioning_not_velocity",
            product_maturity="PROVISIONAL",
            source_product="NISAR L2 GUNW",
            source_sensor="NISAR L-SAR",
            units=units,
            not_ground_truth="true",
            public_displacement_release="not_recommended",
            **tags,
        )
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mintpy-dir", type=Path, required=True)
    parser.add_argument("--commissioning-dir", type=Path, required=True)
    parser.add_argument(
        "--ionosphere-uncertainty",
        type=Path,
        nargs="+",
        required=True,
    )
    parser.add_argument(
        "--pair-coherence-20m",
        type=Path,
        nargs="+",
        required=True,
    )
    parser.add_argument("--peat-geojson", type=Path, required=True)
    parser.add_argument("--canal-mask", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--minimum-coherence", type=float, default=0.50)
    parser.add_argument("--track-number", type=int, default=125)
    parser.add_argument("--frame-number", type=int, default=176)
    parser.add_argument("--orbit-direction", default="Ascending")
    parser.add_argument("--polarization", default="HH")
    args = parser.parse_args()

    inputs = args.mintpy_dir / "inputs"
    source_stack = inputs / "ifgramStack.h5"
    geometry_path = inputs / "geometryGeo.h5"
    manifest_path = args.commissioning_dir / "commissioning_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    support = manifest["reference_support"]
    replacements = {}
    if support.get("peat_geojson"):
        replacements[str(Path(support["peat_geojson"]).resolve())] = (
            args.peat_geojson.resolve()
        )
    if support.get("canal_mask"):
        replacements[str(Path(support["canal_mask"]).resolve())] = (
            args.canal_mask.resolve()
        )
    verified_sources = verify_file_manifest(
        manifest["source_files"],
        replacements=replacements,
    )
    verified_derived = verify_file_manifest(
        manifest.get("derived_files", {}),
    )
    reference = ReferenceSelection(**manifest["reference"])
    variants = {
        name: Path(path) for name, path in manifest["variants"].items()
    }
    inversions = {
        name: Path(paths["timeseries"])
        for name, paths in manifest["inversions"].items()
    }
    inversion_configuration: dict[str, dict[str, str]] = {}
    for name, path in inversions.items():
        with h5py.File(path, "r") as inversion:
            inversion_configuration[name] = {
                "weight_function": str(
                    inversion.attrs.get(
                        "mintpy.networkInversion.weightFunc", "missing"
                    )
                ),
                "observation_dataset": str(
                    inversion.attrs.get(
                        "mintpy.networkInversion.obsDatasetName", "missing"
                    )
                ),
                "mask_dataset": str(
                    inversion.attrs.get(
                        "mintpy.networkInversion.maskDataset", "missing"
                    )
                ),
                "mintpy_commit": str(
                    inversion.attrs.get("PEATGUARD_MINTPY_COMMIT", "missing")
                ),
            }

    with h5py.File(source_stack, "r") as stack:
        attrs = dict(stack.attrs)
        coherence = stack["coherence"][:].astype(np.float64)
        connected = stack["connectComponent"][:]
        wavelength_m = float(stack.attrs["WAVELENGTH"])
    with h5py.File(geometry_path, "r") as geometry:
        incidence = geometry["incidenceAngle"][:].astype(np.float64)
        land = geometry["waterMask"][:].astype(bool)

    transform = _mintpy_transform(attrs)
    crs = f"EPSG:{attrs['EPSG']}"
    quality = build_commissioning_quality_mask(
        coherence,
        connected,
        land,
        incidence,
        minimum_pair_coherence=args.minimum_coherence,
    )
    peat = _rasterize_geojson(
        args.peat_geojson, land.shape, transform, crs
    )
    quality_peat = quality & peat
    minimum_pair_coherence = np.min(coherence, axis=0)

    inversion_checks = {
        name: validate_tree_inversion(name, variants[name], inversions[name])
        for name in variants
    }
    final_vertical_mm = {
        name: los_to_vertical_displacement(
            _read_timeseries_final(path), incidence
        )
        * 1_000.0
        for name, path in inversions.items()
    }
    uncertainty = _read_aligned_uncertainty(
        args.ionosphere_uncertainty,
        land.shape,
        transform,
        crs,
        expected_scene_tags={
            "track_number": str(args.track_number),
            "frame_number": str(args.frame_number),
            "orbit_direction": args.orbit_direction,
        },
    )
    if uncertainty.shape[0] != coherence.shape[0]:
        raise ValueError(
            "Ionosphere uncertainty pair count differs from interferogram stack"
        )
    lower_sigma_mm, conservative_sigma_mm = (
        propagate_ionosphere_uncertainty(
            uncertainty,
            wavelength_m,
            incidence,
            reference=reference,
        )
    )

    strict_alternatives = rank_reference_patches(
        source_stack,
        geometry_path,
        peat_geojson=args.peat_geojson,
        canal_mask=args.canal_mask,
        minimum_separation_m=2_000.0,
    )
    stress_alternatives = rank_reference_patches(
        source_stack,
        geometry_path,
        peat_geojson=args.peat_geojson,
        canal_mask=args.canal_mask,
        minimum_separation_m=2_000.0,
        minimum_coherence=0.65,
        minimum_canal_distance_m=750.0,
    )
    alternatives_by_pixel = {
        (candidate.x, candidate.y): candidate
        for candidate in (
            reference,
            *strict_alternatives,
            *stress_alternatives,
        )
    }
    alternatives = tuple(alternatives_by_pixel.values())
    datum = assess_datum_sensitivity(
        variants["mission_full"],
        alternatives,
        minimum_separation_m=2_000.0,
    )
    full = final_vertical_mm["mission_full"]
    finite_lower = quality & np.isfinite(lower_sigma_mm)
    finite_conservative = quality & np.isfinite(conservative_sigma_mm)
    fraction_1sigma_lower = float(
        np.mean(np.abs(full[finite_lower]) >= lower_sigma_mm[finite_lower])
    )
    fraction_2sigma_lower = float(
        np.mean(
            np.abs(full[finite_lower]) >= 2.0 * lower_sigma_mm[finite_lower]
        )
    )
    fraction_2sigma_conservative = float(
        np.mean(
            np.abs(full[finite_conservative])
            >= 2.0 * conservative_sigma_mm[finite_conservative]
        )
    )

    readiness = assess_temporal_readiness(source_stack)
    coverage_gate = float(quality.sum() / land.sum()) >= 0.20
    algebra_gate = all(check.passed for check in inversion_checks.values())
    uncertainty_gate = fraction_2sigma_lower >= 0.50
    # With fewer than two evaluated candidates the only offset is the chosen
    # reference against itself (identically ~0 by construction), so the span
    # test carries no independent evidence and must not pass.
    datum_degenerate = datum.evaluated_count < 2
    datum_gate = (
        not datum_degenerate
        and datum.percentile_05_95_span_los_mm <= 10.0
    )
    public_displacement = (
        readiness.commissioning_displacement_ready
        and algebra_gate
        and coverage_gate
        and uncertainty_gate
        and datum_gate
    )

    correction_total = (
        final_vertical_mm["mission_full"] - final_vertical_mm["raw"]
    )
    correction_non_ionosphere = (
        final_vertical_mm["mission_full"]
        - final_vertical_mm["ionosphere"]
    )
    signal_to_uncertainty = np.divide(
        np.abs(full),
        conservative_sigma_mm,
        out=np.full(full.shape, np.nan, dtype=np.float64),
        where=np.isfinite(conservative_sigma_mm)
        & (conservative_sigma_mm > 0),
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    common_tags = {
        "reference_date": readiness.first_date,
        "secondary_date": readiness.last_date,
        "span_days": str(readiness.span_days),
        "track_number": str(args.track_number),
        "frame_number": str(args.frame_number),
        "orbit_direction": args.orbit_direction,
        "polarization": args.polarization,
        "pixel_posting_m": (
            f"{abs(transform.a):g}x{abs(transform.e):g}"
        ),
        "wavelength_m": str(wavelength_m),
        "minimum_pair_coherence": str(args.minimum_coherence),
        "relative_reference": (
            f"{reference.longitude:.8f},{reference.latitude:.8f}"
        ),
        "reference_status": reference.reference_status,
    }
    products = {
        "minimum_pair_coherence": _write_product(
            minimum_pair_coherence,
            land & np.all(connected > 0, axis=0),
            args.output_dir
            / "nisar_commissioning_minimum_pair_coherence_80m.tif",
            crs=crs,
            transform=transform,
            description="minimum_pair_coherence",
            units="1",
            tags={
                **common_tags,
                "public_layer_recommendation": "coherence_allowed_with_caveat",
            },
        ),
        "full_corrected_vertical_displacement": _write_product(
            full,
            quality,
            args.output_dir
            / (
                "nisar_commissioning_mission_full_relative_vertical_"
                f"displacement_{readiness.first_date}_{readiness.last_date}_80m.tif"
            ),
            crs=crs,
            transform=transform,
            description="relative_vertical_only_displacement",
            units="mm",
            tags={
                **common_tags,
                "corrections_applied": (
                    "mission_ionosphere,troposphere,solid_earth_tides"
                ),
                "vertical_projection": "los_divided_by_cos_pixel_incidence",
                "public_layer_recommendation": "do_not_publish",
            },
        ),
        "ionosphere_uncertainty_lower_bound": _write_product(
            lower_sigma_mm,
            quality,
            args.output_dir
            / (
                "nisar_commissioning_ionosphere_uncertainty_lower_bound_"
                "vertical_80m.tif"
            ),
            crs=crs,
            transform=transform,
            description="ionosphere_uncertainty_lower_bound_1sigma",
            units="mm",
            tags={
                **common_tags,
                "uncertainty_scope": (
                    "ionosphere_only_excludes_reference_troposphere_"
                    "unwrapping_orbit_dem_horizontal_motion"
                ),
            },
        ),
        "ionosphere_uncertainty_reference_conservative": _write_product(
            conservative_sigma_mm,
            quality,
            args.output_dir
            / (
                "nisar_commissioning_ionosphere_uncertainty_reference_"
                "conservative_vertical_80m.tif"
            ),
            crs=crs,
            transform=transform,
            description=(
                "ionosphere_uncertainty_with_reference_conservative_1sigma"
            ),
            units="mm",
            tags={
                **common_tags,
                "uncertainty_scope": (
                    "ionosphere_plus_reference_patch_excludes_troposphere_"
                    "unwrapping_orbit_dem_horizontal_motion"
                ),
            },
        ),
        "total_correction_sensitivity": _write_product(
            correction_total,
            quality,
            args.output_dir
            / "nisar_commissioning_total_correction_sensitivity_80m.tif",
            crs=crs,
            transform=transform,
            description="full_corrected_minus_raw_vertical_displacement",
            units="mm",
            tags=common_tags,
        ),
        "non_ionosphere_correction_sensitivity": _write_product(
            correction_non_ionosphere,
            quality,
            args.output_dir
            / (
                "nisar_commissioning_non_ionosphere_correction_"
                "sensitivity_80m.tif"
            ),
            crs=crs,
            transform=transform,
            description=(
                "full_corrected_minus_ionosphere_only_vertical_displacement"
            ),
            units="mm",
            tags=common_tags,
        ),
        "signal_to_ionosphere_uncertainty": _write_product(
            signal_to_uncertainty,
            quality,
            args.output_dir
            / (
                "nisar_commissioning_signal_to_ionosphere_uncertainty_80m.tif"
            ),
            crs=crs,
            transform=transform,
            description="absolute_displacement_over_ionosphere_sigma",
            units="1",
            tags={
                **common_tags,
                "interpretation": (
                    "diagnostic_only_not_total_signal_to_noise_ratio"
                ),
            },
        ),
    }

    coherence_20m_arrays: list[np.ndarray] = []
    coherence_20m_transform: rasterio.Affine | None = None
    coherence_20m_crs: str | None = None
    for path in args.pair_coherence_20m:
        with rasterio.open(path) as source:
            if coherence_20m_transform is None:
                coherence_20m_transform = source.transform
                coherence_20m_crs = str(source.crs)
                coherence_20m_shape = source.shape
            elif (
                source.shape != coherence_20m_shape
                or str(source.crs) != coherence_20m_crs
                or not source.transform.almost_equals(
                    coherence_20m_transform
                )
            ):
                raise ValueError(f"20 m coherence grid mismatch: {path}")
            data = source.read(1).astype(np.float64)
            if source.nodata is not None:
                data[data == source.nodata] = np.nan
            coherence_20m_arrays.append(data)
    if len(coherence_20m_arrays) != coherence.shape[0]:
        raise ValueError(
            "20 m coherence pair count differs from interferogram stack"
        )
    coherence_20m_stack = np.stack(coherence_20m_arrays)
    valid_coherence_20m = np.all(
        np.isfinite(coherence_20m_stack), axis=0
    )
    minimum_pair_coherence_20m = np.full(
        valid_coherence_20m.shape, np.nan, dtype=np.float64
    )
    minimum_pair_coherence_20m[valid_coherence_20m] = np.min(
        coherence_20m_stack[:, valid_coherence_20m], axis=0
    )
    if coherence_20m_transform is None or coherence_20m_crs is None:
        raise ValueError("No 20 m coherence rasters supplied")
    coherence_20m_posting = (
        f"{abs(coherence_20m_transform.a):g}x"
        f"{abs(coherence_20m_transform.e):g}"
    )
    products["minimum_pair_coherence_20m"] = _write_product(
        minimum_pair_coherence_20m,
        valid_coherence_20m,
        args.output_dir
        / "nisar_commissioning_minimum_pair_coherence_20m.tif",
        crs=coherence_20m_crs,
        transform=coherence_20m_transform,
        description="minimum_pair_coherence_mission_wrapped_grid",
        units="1",
        tags={
            **common_tags,
            "pixel_posting_m": coherence_20m_posting,
            "effective_support_m": coherence_20m_posting,
            "source_resolution_group": "wrappedInterferogram",
            "public_layer_recommendation": (
                "coherence_allowed_with_provisional_caveat"
            ),
        },
    )

    masks = {"quality_all": quality, "quality_peat": quality_peat}
    if public_displacement:
        release_decision = "public_displacement_allowed"
    elif algebra_gate and coverage_gate:
        release_decision = "coherence_and_coverage_only"
    else:
        release_decision = "no_public_release"
    report = {
        "status": "commissioning_not_velocity",
        "release_decision": release_decision,
        "inputs": {
            "source_stack": str(source_stack.resolve()),
            "geometry": str(geometry_path.resolve()),
            "commissioning_manifest": str(manifest_path.resolve()),
            "ionosphere_uncertainty": [
                str(path.resolve()) for path in args.ionosphere_uncertainty
            ],
            "pair_coherence_20m": [
                str(path.resolve()) for path in args.pair_coherence_20m
            ],
            "peat_extent": str(args.peat_geojson.resolve()),
            "canal_mask": str(args.canal_mask.resolve()),
        },
        "temporal_readiness": readiness.to_dict(),
        "tay_methodology_transfer": {
            "reference": (
                "Tay et al. (2025), Satellite radar advances carbon "
                "emissions accountability over tropical peat, "
                "doi:10.1038/s43247-025-02926-6"
            ),
            "adopted_now": [
                "L-band HH phase-preserving observations",
                "mission ionosphere, troposphere, and solid-Earth-tide "
                "correction variants retained separately",
                "pixel-wise line-of-sight to vertical sensitivity using "
                "local incidence angle",
                "coherence-derived uncertainty and quality screening",
                "MintPy inverse phase-variance weighting configured as "
                "networkInversion.weightFunc=var",
                "Tay long-term rate retained as an independent benchmark",
            ],
            "diagnostic_only_when_available": [
                "pending three-date fine GSLC phase-closure triangle",
                "pending mission GUNW versus fine GSLC phase and coherence "
                "checks",
                "the current two-pair tree cannot empirically distinguish "
                "weighted from unweighted redundant-network inversion",
            ],
            "not_yet_applicable": [
                "bandwidth-20 interferogram-network convergence",
                "dedicated bandwidth-4 per-epoch ionosphere solution for "
                "custom non-nearest GSLC pairs",
                "inverse phase-covariance weighted redundant SBAS inversion",
                "MintPy closure-phase bias correction",
                "stable GNSS-controlled absolute reference",
                "Tay-style neighbouring non-deforming-pixel troposphere "
                "correction",
                "multi-season long-term velocity",
            ],
            "network_geometry_ready": (
                readiness.tay_network_geometry_ready
            ),
            "reference_status": reference.reference_status,
            "policy": (
                "Transfer the method where sensor and temporal support allow; "
                "never fit NISAR values to the Tay raster."
            ),
        },
        "reference": reference.to_dict(),
        "provenance_verification": {
            "manifest": str(manifest_path.resolve()),
            "source_file_count": len(verified_sources),
            "derived_file_count": len(verified_derived),
            "all_manifest_files_verified": True,
            "mintpy_commit": manifest["mintpy_commit"],
        },
        "datum_sensitivity": datum.to_dict(),
        "datum_sensitivity_inputs": {
            "chosen_reference_included": True,
            "strict_alternative_count": len(strict_alternatives),
            "stress_alternative_count": len(stress_alternatives),
            "unique_evaluated_candidate_count": len(alternatives),
            "strict_minimum_coherence": 0.70,
            "stress_minimum_coherence": 0.65,
            "stress_minimum_canal_distance_m": 750.0,
            "minimum_spatial_separation_m": 2_000.0,
        },
        "inversion_validation": {
            name: check.to_dict()
            for name, check in inversion_checks.items()
        },
        "inversion_configuration": inversion_configuration,
        "quality": {
            "minimum_pair_coherence": args.minimum_coherence,
            "land_pixels": int(land.sum()),
            "quality_pixels": int(quality.sum()),
            "quality_fraction_of_land": float(quality.sum() / land.sum()),
            "peat_pixels": int(peat.sum()),
            "quality_peat_pixels": int(quality_peat.sum()),
            "quality_fraction_of_peat": float(
                quality_peat.sum() / peat.sum()
            ),
            "coherence_20m_valid_pixels": int(
                valid_coherence_20m.sum()
            ),
            "coherence_20m_fraction_ge_0_5": float(
                np.mean(
                    minimum_pair_coherence_20m[valid_coherence_20m] >= 0.5
                )
            ),
        },
        "statistics": {
            region: {
                "minimum_pair_coherence": summary_statistics(
                    minimum_pair_coherence, mask
                ),
                "raw_vertical_mm": summary_statistics(
                    final_vertical_mm["raw"], mask
                ),
                "ionosphere_corrected_vertical_mm": summary_statistics(
                    final_vertical_mm["ionosphere"], mask
                ),
                "mission_full_corrected_vertical_mm": summary_statistics(
                    full, mask
                ),
                "total_correction_sensitivity_vertical_mm": (
                    summary_statistics(correction_total, mask)
                ),
                "non_ionosphere_correction_sensitivity_vertical_mm": (
                    summary_statistics(correction_non_ionosphere, mask)
                ),
                "ionosphere_uncertainty_lower_bound_vertical_mm": (
                    summary_statistics(lower_sigma_mm, mask)
                ),
                "ionosphere_uncertainty_reference_conservative_vertical_mm": (
                    summary_statistics(conservative_sigma_mm, mask)
                ),
            }
            for region, mask in masks.items()
        },
        "significance": {
            "fraction_abs_displacement_ge_1sigma_ionosphere_lower_bound": (
                fraction_1sigma_lower
            ),
            "fraction_abs_displacement_ge_2sigma_ionosphere_lower_bound": (
                fraction_2sigma_lower
            ),
            "fraction_abs_displacement_ge_2sigma_ionosphere_reference_"
            "conservative": fraction_2sigma_conservative,
            "warning": (
                "Ionosphere uncertainty is not total uncertainty; omitted "
                "terms can only weaken the measurement claim."
            ),
        },
        "release_gates": {
            "temporal_commissioning_displacement": (
                readiness.commissioning_displacement_ready
            ),
            "temporal_annual_velocity": readiness.annual_velocity_ready,
            "temporal_tay_network_geometry": (
                readiness.tay_network_geometry_ready
            ),
            "inverse_phase_variance_weighting_configured": all(
                configuration["weight_function"] == "var"
                for configuration in inversion_configuration.values()
            ),
            "exact_inversion_algebra": algebra_gate,
            "quality_coverage_at_least_20_percent_of_land": coverage_gate,
            "at_least_50_percent_ge_2sigma_ionosphere_lower_bound": (
                uncertainty_gate
            ),
            "reference_05_95_span_at_most_10mm_los": datum_gate,
            "datum_gate_degenerate_single_candidate": datum_degenerate,
            "public_displacement_release": public_displacement,
            "public_coherence_release": algebra_gate and coverage_gate,
            "gate_threshold_status": (
                "conservative_internal_release_policy_not_a_published_"
                "universal_standard"
            ),
        },
        "products": {name: str(path) for name, path in products.items()},
        "limitations": [
            (
                f"{readiness.pair_count} pairs and {readiness.acquisition_count} "
                f"acquisitions over {readiness.span_days} days cannot support "
                "an annual velocity."
            ),
            "The reference patch is a relative datum, not surveyed stable ground.",
            "Vertical conversion assumes horizontal motion is negligible.",
            "Mission ionosphere uncertainty excludes other error sources.",
            "NISAR PROVISIONAL products are fully calibrated but only partially validated.",
        ],
    }
    report_path = args.output_dir / "nisar_commissioning_audit.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
