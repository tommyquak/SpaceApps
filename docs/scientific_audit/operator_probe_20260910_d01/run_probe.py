"""Bounded, preselected D-01 GSLC aggregation experiment, without downloads.

Run from the repository root with PYTHONPATH=src. Inputs
are opened read-only. Outputs are written once to a fresh sibling data folder.
This is an operator diagnostic, not a deformation or Tay replication test.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.windows import Window, bounds, transform

from peatguard.insar.nisar_gslc import form_gslc_interferogram
from peatguard.insar.nisar_validation import (
    _reproject_fine_block,
    _valid_raster_values,
    compare_fine_to_mission,
    wrap_phase,
)


ROOT = Path(os.environ.get("DATA_ROOT", "."))  # monorepo data root holding grid_products/
HERE = Path(__file__).resolve().parent
INPUTS = ROOT / "grid_products/_deletable_work/nisar"


def crop_mission(source: Path, destination: Path, window: Window) -> None:
    with rasterio.open(source) as src:
        profile = {
            **src.profile,
            "width": int(window.width),
            "height": int(window.height),
            "transform": transform(window, src.transform),
        }
        with rasterio.open(destination, "w", **profile) as dst:
            dst.write(src.read(window=window))
            dst.update_tags(**src.tags())
            for band, description in enumerate(src.descriptions, 1):
                if description:
                    dst.set_band_description(band, description)


def common_support_metrics(products: dict, mission_phase: Path, mission_coh: Path) -> list:
    reports = []
    with (
        rasterio.open(products["wrapped_phase"]) as phase,
        rasterio.open(products["coherence"]) as coherence,
        rasterio.open(products["complex_cross_product"]) as cross,
        rasterio.open(mission_phase) as mission,
        rasterio.open(mission_coh) as mission_quality,
    ):
        window = Window(0, 0, mission.width, mission.height)
        unit_phase, unit_coh, unit_strength = _reproject_fine_block(
            phase, coherence, mission, window
        )
        weighted_phase, weighted_coh, weighted_strength = _reproject_fine_block(
            phase, coherence, mission, window, cross
        )
        mission_values = mission.read(1)
        mission_coherence = mission_quality.read(1)
        base = (
            np.isfinite(unit_phase) & np.isfinite(weighted_phase)
            & np.isfinite(unit_coh) & np.isfinite(weighted_coh)
            & _valid_raster_values(mission_values, mission.nodata)
            & _valid_raster_values(mission_coherence, mission_quality.nodata)
            & (unit_coh >= 0.25) & (weighted_coh >= 0.25)
            & (mission_coherence >= 0.25)
        )
        for threshold in (0.0, 0.25):
            common = base & (unit_strength >= threshold) & (weighted_strength >= threshold)
            row = {
                "minimum_coherence": 0.25,
                "minimum_strength_in_both_modes": threshold,
                "common_pixels": int(common.sum()),
                "modes": {},
            }
            for mode, values in (("unit_phasor", unit_phase), ("complex_weighted", weighted_phase)):
                metrics = {}
                for sign in (1, -1):
                    delta = wrap_phase(sign * values[common].astype(float) - mission_values[common])
                    offset = float(np.angle(np.exp(1j * delta).sum()))
                    residual = wrap_phase(delta - offset)
                    metrics["same_sign" if sign == 1 else "opposite_sign"] = {
                        "circular_offset_rad": offset,
                        "circular_mae_rad": float(np.abs(residual).mean()),
                        "circular_rmse_rad": float(np.sqrt((residual ** 2).mean())),
                        "resultant_length": float(np.abs(np.exp(1j * residual).mean())),
                    }
                row["modes"][mode] = metrics
            reports.append(row)
    return reports


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE / "data")
    output = parser.parse_args().output_dir.resolve()
    output.mkdir(exist_ok=False)
    site_source = ROOT / "deliverables/tay_full_grid/dam_sites.geojson"
    site = json.loads(site_source.read_text())["features"][0]
    lon, lat = site["geometry"]["coordinates"]
    pairs = (("20260622", "20260704"), ("20260704", "20260716"))
    summary = {
        "selection_rule": "first existing dam-site feature, before looking at phase/coherence",
        "site_source": str(site_source),
        "site": site,
        "interpretation": (
            "Local operator diagnostic only. No Tay fitting or deformation validation. "
            "Both operators follow GSLC geocoding, not mission radar-coordinate multilooking. "
            "The small crop cannot satisfy the production million-pixel support gate."
        ),
        "pairs": [],
    }
    for first, second in pairs:
        pair_output = output / f"{first}_{second}"
        pair_output.mkdir()
        reference, secondary = [
            next((INPUTS / "gslc_subset").glob(f"*_{day}T*_peatguard_aoi.h5"))
            for day in (first, second)
        ]
        phase_source = INPUTS / "gunw_products" / f"nisar_gunw_{first}_{second}_wrapped_phase_20m.tif"
        coherence_source = INPUTS / "gunw_products" / f"nisar_gunw_{first}_{second}_coherence_20m.tif"
        with rasterio.open(phase_source) as source:
            x, y = Transformer.from_crs(4326, source.crs, always_xy=True).transform(lon, lat)
            row, col = source.index(x, y)
            window = Window(col - 50, row - 50, 100, 100)
            if window.col_off < 0 or window.row_off < 0:
                raise ValueError("Preselected window lies outside source")
            bbox = Transformer.from_crs(source.crs, 4326, always_xy=True).transform_bounds(
                *bounds(window, source.transform), densify_pts=21
            )
            mission_bounds = bounds(window, source.transform)
        mission_phase = pair_output / "mission_phase_20m.tif"
        mission_coherence = pair_output / "mission_coherence_20m.tif"
        crop_mission(phase_source, mission_phase, window)
        crop_mission(coherence_source, mission_coherence, window)
        products = form_gslc_interferogram(
            reference, secondary, pair_output, bbox_wgs84=bbox,
            block_size=256, export_cross_product=True,
        )
        pair_report = {
            "pair": [first, second],
            "bbox_wgs84": bbox,
            "mission_window": [window.col_off, window.row_off, window.width, window.height],
            "mission_bounds": mission_bounds,
            "input_inventory_not_content_hashes": [
                {"path": str(path), "size_bytes": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns}
                for path in (reference, secondary, phase_source, coherence_source)
            ],
            "individually_gated": [],
        }
        for threshold in (0.0, 0.25):
            for mode in ("unit_phasor", "complex_weighted"):
                report = compare_fine_to_mission(
                    products["wrapped_phase"], products["coherence"],
                    mission_phase, mission_coherence,
                    output_json=pair_output / f"{mode}_strength_{threshold}.json",
                    aggregation=mode,
                    fine_cross_product_path=(products["complex_cross_product"] if mode == "complex_weighted" else None),
                    min_phasor_strength=threshold,
                )
                pair_report["individually_gated"].append({
                    "mode": mode, "minimum_strength": threshold,
                    "phase": report["phase"], "coherence": report["coherence"],
                    "release_gates": report["release_gates"],
                })
        pair_report["identical_pixel_support"] = common_support_metrics(
            products, mission_phase, mission_coherence
        )
        pair_report["output_raster_hashes"] = {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in pair_output.glob("*.tif")
        }
        summary["pairs"].append(pair_report)
    summary["output_data_bytes_before_summary"] = sum(p.stat().st_size for p in output.rglob("*") if p.is_file())
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"output": str(output), "pairs": [
        {"pair": pair["pair"], "identical_pixel_support": pair["identical_pixel_support"]}
        for pair in summary["pairs"]
    ], "bytes": summary["output_data_bytes_before_summary"]}, indent=2))


if __name__ == "__main__":
    main()
