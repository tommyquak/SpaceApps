# Sample data

Small rasters so the viewer runs out of the box (`make demo`). Everything
here is a derived, compact copy; the full-resolution products and the raw
NISAR granules stay outside the repository. All files are Cloud-Optimized
GeoTIFFs readable by TiTiler, QGIS or rasterio.

Area of interest: Sebangau and the former Mega Rice Project, Central
Kalimantan, Indonesia (114.21 to 114.83 E, 3.00 to 2.42 S).

| File | Size | What it is | Maturity |
| --- | --- | --- | --- |
| `nisar_coherence_80m.tif` | 3.2 MB | Minimum interferometric coherence across the four 12-day-chain NISAR GUNW pairs 4 Jul to 14 Sep 2026 (ascending track 125, frame 176, HH), 4x4 block mean of the 20 m product to 80 m | PROVISIONAL, quality layer, not displacement |
| `nisar_commissioning_minimum_pair_coherence_80m.tif` | 2.8 MB | Minimum pair coherence of the first commissioning GUNW pairs, 22 Jun to 21 Aug 2026, same track and frame, native 80 m GUNW grid | PROVISIONAL, quality layer |
| `nisar_A_disp60_hp10km_mm.tif` | 19 MB | 60-day relative line-of-sight displacement from the chained ascending GUNW unwrapped phase, quadratic ramp removed per pair on stable ground, 10 km spatial high-pass, mm, EPSG:32750, 80 m | EXPERIMENTAL, fails the release gate (see below) |
| `nisar_D_disp60_hp10km_mm.tif` | 19 MB | Same for the descending chain | EXPERIMENTAL |
| `subsidence_velocity.tif` | 6.5 MB | Tay et al. 2025 ALOS-2 L-band long-term rate, 12 Oct 2018 to 9 Aug 2024, mm/yr, negative = subsidence, EPSG:4326 (~15 m grid, valid only where Tay measured) | Published independent benchmark |
| `restoration_zones.tif` | 6.0 MB | Restoration priority 0 to 1: 0.75 x canal_risk + 0.25 x fire history, where canal_risk = 0.45 x canal proximity + 0.55 x normalised Tay subsidence | Screening layer [S] |
| `fire_history.tif` | 3.7 MB | Recency- and frequency-weighted burn history 2015 to 2024 from MODIS MCD64A1 burned area, 0 to 1, 2x2 block mean of the 20 m product | Screening layer [S] |

Nothing here is ground truth and nothing is a water-table proxy. The NISAR
layers carry their maturity in the GeoTIFF tags (`product_maturity`,
`processing_status`, `not_ground_truth`, `public_displacement_release`).

## Provenance and licences

- NISAR L2 GUNW (`NISAR_L2_GUNW`, provisional CRID): NASA/ISRO NISAR mission,
  distributed by the Alaska Satellite Facility DAAC through NASA Earthdata.
  NASA data policy: free and open, no restrictions on use; cite
  "NASA/ISRO NISAR L2 GUNW, provisional, ASF DAAC". Provisional products are
  calibrated but only partially validated and will be reprocessed by JPL.
- Tay et al. 2025 long-term rate: Tay, C. et al. (2025), Communications Earth
  & Environment, "Satellite radar advances carbon emissions accountability
  over tropical peat"; raster from the authors' DR-NTU repository under its
  own terms. Used only as an independent benchmark and as a published input
  to the restoration screen; never fitted or calibrated against.
- NASA FIRMS MODIS MCD64A1 burned area and active-fire detections: NASA
  LANCE/FIRMS, free and open (NASA Earth science data policy).
- Landsat Collection 2 Level-2 (used upstream for the burn-index epochs in
  `viz/`): USGS, public domain.
- OSM settlements and waterways (`src/peatguard/webdata/`): OpenStreetMap
  contributors, ODbL.

## Processing lineage

| Layer | Code | Documentation |
| --- | --- | --- |
| GUNW download and QA export | `scripts/export_nisar_gunw.py`, `src/peatguard/insar/nisar_gunw.py`, `cluster/full/nisar_pull.sbatch` | `docs/NISAR_INTEGRATION_BASELINE_PLAN.md` |
| Minimum pair coherence (20 m and 80 m) | `src/peatguard/timeseries/nisar_commissioning.py`, `scripts/audit_nisar_commissioning.py`, `cluster/full/analysis/nisar_stage2c.py` | `docs/scientific_audit/NISAR_SCIENTIFIC_AUDIT.md` |
| 60-day high-passed displacement (A and D) | `cluster/full/analysis/nisar_stage2c.py` (composite, ramp, 10 km high-pass), `cluster/full/analysis/nisar_stage2d.py` (site versus control) | `docs/scientific_audit/RECOVERY_DECISION_20260910.md`, `docs/scientific_audit/NISAR_REFERENCE_ENSEMBLE_PREREG_20260910.md` |
| Release gate the displacement fails | `src/peatguard/timeseries/nisar_commissioning.py`, `tests/test_nisar_commissioning.py` | `docs/scientific_audit/NISAR_SCIENTIFIC_AUDIT.md` (executive decision) |
| Tay benchmark and restoration priority | `scripts/build_tay_restoration.py`, `src/peatguard/analysis/restoration_zones.py`, `src/peatguard/analysis/risk_score.py` | `docs/scientific_audit/TAY_REPLICATION_DECISION_20260910.md` |
| Fire history | `src/peatguard/analysis/fire_history.py`, `src/peatguard/ingest/fire.py` | `viz/degraded_peat_fire_recurrence/README.md` |
| Downsampling for this bundle | 4x4 (coherence) and 2x2 (fire) block means with nodata excluded, rasterio, deflate, 512 px tiles, overviews | this file |

## Why the displacement layers are experimental

The chained GUNW phase is dominated by per-track long-wavelength error: the
two tracks anti-correlate and nominally stable ground moves 16 to 47 mm
between pairs. There is no mineral-ground reference within 60 km of the AOI,
the per-pair residual is +-20 to 40 mm and flips sign between pairs, and the
GUNW phase still contains ionosphere after the mission screen. After the 10 km
high-pass only local anomalies remain, which is what the site-versus-control
diagnostic needs; the layer is a pattern check, not a rate, and the sign
convention (positive = range increase) is reported as found [S]. A sustained
subsidence rate needs a year or more of 12-day acquisitions, a bandwidth-20
network, closure-phase bias correction and a GNSS-controlled datum [E].
