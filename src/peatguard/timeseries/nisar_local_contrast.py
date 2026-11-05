"""Native GUNW local contrasts without a presumed stationary reference.

[S] This is a relative-observable adaptation, not a replication of Tay et al.
(2025), doi:10.1038/s43247-025-02926-6. Subtract target and background on
fixed support BEFORE network inversion. A common phase datum cancels, but
background motion and spatially varying propagation errors do not.

[E] GUNW corrections are not applied by default:
https://nisar-docs.asf.alaska.edu/gunw/
[S] Product labels are additionally separated at invalid-sample boundaries,
following the safeguard suggested by H. Fattahi on 27 August 2026:
https://forum.earthdata.nasa.gov/viewtopic.php?t=8137
"""

from __future__ import annotations

from contextlib import ExitStack
from dataclasses import asdict
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform

import h5py
import numpy as np
from pyproj import CRS, Transformer

from peatguard.insar.nisar_gunw import inspect_gunw, inspect_gunw_grid
from peatguard.timeseries.local_contrast import extract_local_contrast, solve_contrast_network
from peatguard.timeseries.nisar_commissioning import phase_to_los_displacement


def decode_gunw_mask(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Decode the current uint32 GUNW mask, not the older uint8 schema.

    [E] Encoding is specified in the native product's mask description.
    [S] Conservatively reject water, anomalies and reserved bits. Bit 24
    denotes an interpolated ionosphere estimate, not an invalid SLC sample;
    retain it but disclose it separately. This is NOT a stationary-land mask.
    """
    mask = np.asarray(mask)
    if mask.dtype != np.dtype("uint32"):
        raise ValueError("Requires the uint32 GUNW mask schema; no implicit conversion")
    low = mask & 255
    valid = ((low != 255) & (low < 100) & ((low // 10) > 0) & ((low % 10) > 0)
             & (((mask >> 8) & 65535) == 0) & ((mask >> 25) == 0))
    return valid, ((mask >> 24) & 1).astype(bool)


def _window_values(dataset: h5py.Dataset, selection: tuple[slice, slice]) -> np.ndarray:
    values = dataset[selection].astype(np.float64)
    fill = dataset.attrs.get("_FillValue")
    if fill is not None:
        values[values == fill] = np.nan
    return values


def _array_hash(values: np.ndarray) -> str:
    digest = hashlib.sha256()
    digest.update(str(values.shape).encode())
    digest.update(values.dtype.str.encode())
    digest.update(np.ascontiguousarray(values).tobytes())
    return digest.hexdigest()


def _units(dataset: h5py.Dataset) -> str:
    value = dataset.attrs.get("units", "")
    return value.decode() if isinstance(value, bytes) else str(value)


def _lineage(handle: h5py.File) -> dict:
    base = "science/LSAR/GUNW/metadata/processingInformation"
    result = {}
    for role in ("Reference", "Secondary"):
        values = np.asarray(handle[f"{base}/inputs/l1{role}SlcGranules"][()]).reshape(-1)
        names = sorted(v.decode() if isinstance(v, bytes) else str(v) for v in values)
        if not names or any(not n.strip() for n in names) or len(set(names)) != len(names):
            raise ValueError("RSLC lineage must contain distinct, nonempty granule identities")
        result[role.lower() + "_rslc_granules"] = names
    bandwidths = {}
    for role in ("reference", "secondary", "unwrappedInterferogram"):
        bandwidth = float(handle[f"{base}/parameters/{role}/frequencyA/rangeBandwidth"][()])
        if not np.isfinite(bandwidth) or bandwidth <= 0:
            raise ValueError("Invalid native range bandwidth")
        bandwidths[role] = bandwidth
    result["range_bandwidth_hz"] = bandwidths
    return result


def audit_local_contrasts(
    gunw_paths: list[Path], sites_path: Path, *, target_radius_m: float = 500,
    background_inner_m: float = 1000, background_outer_m: float = 2000,
    min_coherence: float = 0.3, min_pixels: int = 10, min_fraction: float = 0.5,
    polarization: str = "HH",
) -> dict:
    """Read only small native-grid windows and return a diagnostic report.

    No spatial resampling, fitted ramp, external reference or Tay labels are
    used. Raw and ionosphere-subtracted variants share EXACTLY the same valid
    support. Troposphere/SET are deliberately left unapplied in both variants.
    Radii and quality choices are engineering parameters, not accuracy gates.
    """
    if not np.isfinite([target_radius_m, background_inner_m, background_outer_m,
                       min_coherence, min_fraction]).all():
        raise ValueError("Parameters must be finite")
    if not 0 < target_radius_m < background_inner_m < background_outer_m:
        raise ValueError("Require 0 < target radius < background inner < background outer")
    if not 0 < min_coherence <= 1 or not 0 < min_fraction <= 1 or min_pixels < 1:
        raise ValueError("Invalid quality thresholds")
    if not gunw_paths or polarization not in {"HH", "VV"}:
        raise ValueError("Supply at least one GUNW and a co-polarization")
    paths = sorted(Path(p).resolve() for p in gunw_paths)
    infos = [inspect_gunw(p) for p in paths]
    grids = [inspect_gunw_grid(p, "unwrappedInterferogram", polarization=polarization)
             for p in paths]
    first, grid = infos[0], grids[0]
    identity = (first.track_number, first.frame_number, first.orbit_direction,
                first.composite_release_id, first.product_version,
                first.product_specification_version)
    for info, other in zip(infos, grids):
        current = (info.track_number, info.frame_number, info.orbit_direction,
                   info.composite_release_id, info.product_version,
                   info.product_specification_version)
        if current != identity or not np.isclose(info.wavelength_m, first.wavelength_m,
                                                rtol=0, atol=1e-12):
            raise ValueError("Separate tracks, frames, releases, versions and radar frequencies")
        if (other.epsg != grid.epsg or other.shape != grid.shape
                or not np.allclose(other.x_coordinates, grid.x_coordinates, rtol=0, atol=1e-6)
                or not np.allclose(other.y_coordinates, grid.y_coordinates, rtol=0, atol=1e-6)):
            raise ValueError("Native GUNW grids differ; no implicit phase/component resampling")
    crs = CRS.from_epsg(grid.epsg)
    if not crs.is_projected or any(a.unit_conversion_factor != 1 for a in crs.axis_info):
        raise ValueError("Native grid must use projected metres")
    pairs = [(i.reference_date, i.secondary_date) for i in infos]
    # Also fails closed on duplicated/reversed pairs or disconnected periods.
    network = solve_contrast_network(pairs, np.zeros(len(pairs)))
    network = {k: network[k] for k in ("dates", "rank", "cycle_count", "gauge_date")}
    sites_bytes = Path(sites_path).read_bytes()
    geojson = json.loads(sites_bytes)
    if geojson.get("type") != "FeatureCollection" or "crs" in geojson:
        raise ValueError("Sites must be a WGS84 GeoJSON FeatureCollection without legacy CRS")
    transform = Transformer.from_crs(4326, grid.epsg, always_xy=True)
    rows = []
    source_records = []
    with ExitStack() as stack:
        handles = [stack.enter_context(h5py.File(p, "r")) for p in paths]
        groups = [h[g.group_path] for h, g in zip(handles, grids)]
        masks = [h[g.group_path.rsplit("/", 1)[0] + "/mask"] for h, g in zip(handles, grids)]
        source_dates = {}
        expected_bandwidth = None
        for info, group, mask, handle in zip(infos, groups, masks, handles):
            if mask.dtype != np.dtype("uint32") or mask.shape != grid.shape:
                raise ValueError("Unsupported native mask schema or grid")
            for name in ("unwrappedPhase", "coherenceMagnitude", "connectedComponents",
                         "ionospherePhaseScreen"):
                if group[name].shape != grid.shape:
                    raise ValueError(f"Native layer grid mismatch: {name}")
            for name in ("unwrappedPhase", "ionospherePhaseScreen"):
                if _units(group[name]) != "radians":
                    raise ValueError(f"Expected radians for {name}")
            lineage = _lineage(handle)
            for role, date in (("reference", info.reference_date), ("secondary", info.secondary_date)):
                names = lineage[role + "_rslc_granules"]
                if date in source_dates and source_dates[date] != names:
                    raise ValueError("Shared-date RSLC lineage differs between pairs")
                source_dates[date] = names
            for bandwidth in lineage["range_bandwidth_hz"].values():
                if expected_bandwidth is None:
                    expected_bandwidth = bandwidth
                if not np.isclose(bandwidth, expected_bandwidth, rtol=0, atol=1e-3):
                    raise ValueError("Separate processing with different range bandwidths")
            rec = asdict(info)
            rec["path"] = str(info.path)
            rec["bytes"] = info.path.stat().st_size
            rec["mtime_ns"] = info.path.stat().st_mtime_ns
            rec["lineage"] = lineage
            source_records.append(rec)

        seen_ids = set()
        for index, feature in enumerate(geojson["features"]):
            geometry = feature["geometry"]
            if geometry["type"] != "Point":
                raise ValueError("Site geometries must be points")
            lon, lat = geometry["coordinates"][:2]
            if not np.isfinite([lon, lat]).all() or not (-180 <= lon <= 180 and -90 <= lat <= 90):
                raise ValueError("Invalid longitude/latitude")
            props = feature.get("properties") or {}
            site_id = str(props.get("site_id", props.get("id", index + 1)))
            if site_id in seen_ids:
                raise ValueError("Site identifiers must be unique")
            seen_ids.add(site_id)
            x, y = transform.transform(lon, lat)
            radius = background_outer_m
            row = {"site": site_id, "longitude": lon, "latitude": lat,
                   "status": "outside_full_window", "available": False}
            rows.append(row)
            if (x-radius < grid.x_coordinates[0] or x+radius > grid.x_coordinates[-1]
                    or y+radius > grid.y_coordinates[0] or y-radius < grid.y_coordinates[-1]):
                continue
            c0, c1 = np.searchsorted(grid.x_coordinates, [x-radius, x+radius])
            r0, r1 = np.searchsorted(-grid.y_coordinates, [-y-radius, -y+radius])
            selection = (slice(r0, r1), slice(c0, c1))
            yy, xx = np.meshgrid(grid.y_coordinates[r0:r1], grid.x_coordinates[c0:c1],
                                 indexing="ij")
            distance2 = (xx-x)**2 + (yy-y)**2
            target = distance2 <= target_radius_m**2
            background = ((distance2 >= background_inner_m**2)
                          & (distance2 <= background_outer_m**2))
            raw, ion, coherence, components, interpolated, hashes = [], [], [], [], [], []
            for group, mask in zip(groups, masks):
                values = {name: _window_values(group[name], selection) for name in (
                    "unwrappedPhase", "coherenceMagnitude", "connectedComponents",
                    "ionospherePhaseScreen")}
                mask_values = mask[selection]
                valid, interp = decode_gunw_mask(mask_values)
                cc = values["connectedComponents"]
                valid &= np.isfinite(values["ionospherePhaseScreen"])
                valid &= (cc > 0) & (cc < 65535) & np.isfinite(cc)
                hashes.append({**{k: _array_hash(v) for k, v in values.items()},
                               "native_mask": _array_hash(mask_values)})
                raw.append(np.where(valid, values["unwrappedPhase"], np.nan))
                ion.append(values["ionospherePhaseScreen"])
                coherence.append(values["coherenceMagnitude"])
                components.append(np.where(valid, cc, 0))
                interpolated.append(interp)
            raw, ion, coherence, components = map(np.asarray, (raw, ion, coherence, components))
            kwargs = dict(min_coherence=min_coherence, min_pixels=min_pixels,
                          min_target_fraction=min_fraction, min_background_fraction=min_fraction)
            variants = {}
            for name, phase in (("raw", raw), ("ionosphere_subtracted", raw-ion)):
                result = extract_local_contrast(phase, coherence, components, target,
                                                background, **kwargs)
                if result["available"]:
                    solution = solve_contrast_network(pairs, result["contrast_rad"])
                    solution["epoch_relative_los_mm"] = (1000 * phase_to_los_displacement(
                        np.asarray(solution["epoch_contrast_rad"]), first.wavelength_m)).tolist()
                    result["network_solution"] = solution
                variants[name] = result
            row.update(status=variants["raw"]["status"], available=variants["raw"]["available"],
                       window_row_col=[int(r0), int(c0)], window_shape=list(target.shape),
                       requested_target_pixels=int(target.sum()),
                       requested_background_pixels=int(background.sum()),
                       consumed_window_sha256=hashes, variants=variants,
                       window_coordinates_sha256={"x": _array_hash(grid.x_coordinates[c0:c1]),
                                                  "y": _array_hash(grid.y_coordinates[r0:r1])},
                       ionosphere_interpolated_fraction_requested_regions=[
                           float(np.mean(v[target | background])) for v in interpolated])
            if row["available"]:
                result = variants["raw"]
                row["ionosphere_interpolated_fraction_selected_support"] = {
                    region: [float(np.mean(v.ravel()[result[region + "_support_flat_indices"]]))
                             for v in interpolated] for region in ("target", "background")}
                if any(result[key] != variants["ionosphere_subtracted"][key] for key in
                       ("target_support_flat_indices", "background_support_flat_indices")):
                    raise ValueError("Correction variants changed selected spatial support")
                raw_end = variants["raw"]["network_solution"]["epoch_relative_los_mm"][-1]
                ion_end = variants["ionosphere_subtracted"]["network_solution"]["epoch_relative_los_mm"][-1]
                row["endpoint_ionosphere_sensitivity_mm"] = ion_end - raw_end

    return {
        "schema_version": 1, "observable": "local_target_minus_background_los_change",
        "status": "experimental_relative_diagnostic_only", "absolute_datum": False,
        "release_ready": False, "annual_velocity_ready": False, "ground_accuracy_validated": False,
        "background_assumed_stationary": False, "sources": source_records,
        "sites_source": {"path": str(Path(sites_path).resolve()),
                         "sha256": hashlib.sha256(sites_bytes).hexdigest()},
        "native_grid": {"epsg": grid.epsg, "pixel_spacing_m": [grid.x_spacing, grid.y_spacing],
                        "transform": list(grid.transform), "shape": list(grid.shape),
                        "resampling": "none"},
        "implementation_sha256": {
            name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ("nisar_local_contrast.py", "local_contrast.py", "nisar_commissioning.py")},
        "runtime_versions": {"python": platform.python_version(),
                             **{name: version(name) for name in ("numpy", "scipy", "h5py", "pyproj")}},
        "parameters": {"target_radius_m": target_radius_m, "background_inner_m": background_inner_m,
                       "background_outer_m": background_outer_m, "minimum_coherence": min_coherence,
                       "minimum_pixels": min_pixels, "minimum_fraction": min_fraction,
                       "polarization": polarization, "spatial_connectivity": 4},
        "network": network, "date_pairs": pairs, "sites": rows,
        "available_sites": sum(r["available"] for r in rows), "total_sites": len(rows),
        "phase_to_los_scale_m_per_rad": -first.wavelength_m / (4*np.pi),
        "los_sign_convention": "positive_toward_sensor; target_minus_background; not_vertical",
        "restoration_attribution": False,
        "limitations": [
            "Relative LOS change only: background may move; shared regional motion is unobservable.",
            "Ionosphere-only subtraction is not full atmospheric or solid-Earth-tide correction.",
            "Spatial atmosphere, scattering changes and within-component unwrap errors may remain.",
            "Fixed common spatial support is conditional on these input pairs, not prospective validation.",
            "Local 4-connected relabeling is conservative and may reject paths outside the crop.",
            "A tree network has no independent measured closure check; exact fit is not accuracy.",
            "No observational covariance supplied: epoch uncertainty is unknown, not zero.",
            "No vertical projection, annualization, Tay calibration or restoration attribution.",
            "Input identity uses native metadata and hashes of consumed windows, not full-file hashes.",
        ],
    }
