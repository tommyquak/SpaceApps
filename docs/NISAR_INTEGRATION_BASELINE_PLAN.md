# PeatGuard NISAR integration and historical-baseline plan

Status: implementation in progress, 2026-07-25.

## Storage boundary

All large, reproducible downloads and processing intermediates are contained
under `grid_products/_deletable_work/`. This includes full GUNW files, cropped
complex GSLC files, fine-resolution diagnostic rasters, MintPy workspaces, and
the ALOS-2/ALOS-1 raw ISCE2/MintPy pilots. Compact validated COGs, reports, and
web-ready copies remain outside this purgeable root. A size manifest and exact
cleanup target will be produced at handoff.

## Scientific position

PeatGuard will keep five distinct concepts separate:

1. **Current measurement:** NISAR L-band observations from June 2026 onward.
2. **Current geodetic datum:** an independently stable GNSS reference such as
   CPKY, used only to fix the additive displacement reference.
3. **Independent benchmark:** Tay et al. (2025) ALOS-2 products, used to test
   agreement without fitting or force-calibrating PeatGuard estimates to Tay.
4. **Reproducible near-historical segment:** open ALOS-2 PALSAR-2 ScanSAR
   Level 1.1 observations, processed independently through the Tay method and
   kept as a separate 2018–2025 segment.
5. **Deep historical comparison:** ALOS-1 PALSAR-1 observations from 2007–2010,
   processed as their own time series and never blended numerically with a 2026
   NISAR velocity as if both described the same period.

This separation is necessary because drainage, fire, land cover and canal
blocking can change peat-subsidence rates substantially between epochs [S].

NISAR products acquired from 2026-06-17 onward are currently **PROVISIONAL**:
calibrated but only partially validated [E]. Every PeatGuard NISAR product and
web layer must retain that maturity label.

## Live NISAR inventory selected for PeatGuard

Operational AOI: `114.21 E, 3.00 S, 114.83 E, 2.42 S`.

The preferred high-resolution stack is:

- Path 125, frame 176, ascending.
- Frequency A, 40 MHz (`4005` in GSLC/RSLC filenames).
- HH polarization for interferometry.
- Acquisitions: 2026-06-22, 2026-07-04 and 2026-07-16.
- Each GSLC footprint covers 100% of the operational AOI.
- Each GSLC is 20.18–20.75 GiB; all three total 61.1 GiB.
- The 40 MHz GSLC Frequency A layer is posted at 5 × 5 m [E].

The NISAR Data User Guide states that GSLC phase is geocoded and
topographically flattened, so cross-multiplying matching GSLCs yields a
phase-flattened interferogram [E]:

- <https://nisar-docs.asf.alaska.edu/gslc/>
- <https://nisar-docs.asf.alaska.edu/amplitude/>
- <https://nisar.asf.earthdatacloud.nasa.gov/NISAR-SAMPLE-DATA/DOCS/NISAR_D-102269_RevE_NASA_SDS_Product_Specification_L2_GSLC_Nov8_2024_w-sigs.pdf>

Mission GUNW provides:

- Wrapped interferogram and coherence at 20 m posting.
- Unwrapped phase, coherence, connected components and ionosphere screen at
  80 m posting.
- Nearest-neighbour 12-day pairs only.
- Atmospheric, ionospheric and solid-Earth correction layers that are supplied
  but not applied by default [E].

Reference: <https://nisar-docs.asf.alaska.edu/gunw/>

### Supported time-series loader

MintPy release 1.6.3 predates the operational NISAR loader. MintPy merged
native NISAR loading and auxiliary-stack support in pull request 1487 on
2026-04-22. PeatGuard pins upstream revision
`53030ef4ae390ff17e845a342e086993eca81115`, verified on 2026-07-23.
`prep_nisar.py` now builds:

- `ifgramStack.h5` from mission unwrapped phase/coherence/components.
- `ionStack.h5` from the supplied ionosphere phase screen.
- `tropoStack.h5` by DEM-aware interpolation of wet plus hydrostatic phase.
- `setStack.h5` by DEM-aware interpolation of solid-Earth-tide phase.
- `geometryGeo.h5` using the supplied NISAR radar-grid lookup tables.

The upstream routine prepares all four stacks, but the standard
`smallbaselineApp.py` correction sequence only consumes `ionStack.h5`
directly. PeatGuard will therefore test both ERA5/PySolid and explicit
inversion/subtraction of the mission-provided troposphere/SET stacks. The
correction routes remain separate outputs until their sign and magnitude
tests pass.

The required DEM will be the mission's open Modified Copernicus DEM for NISAR
(DOI `10.5067/NIDEM-1`), whose heights are converted from EGM2008 to the WGS84
ellipsoid. This matches the vertical reference needed for interpolation into
the GUNW radar-grid correction cubes [E].

## Tay et al. (2025) methodology transfer

PeatGuard will use the processing logic in Tay et al. (2025) wherever the
NISAR acquisition count and products support it, without using Tay's output
values as a calibration target [E]. Tay processed 92 ALOS-2 PALSAR-2 ScanSAR
L1.1 acquisitions from 2018-10-12 to 2024-08-09 with ISCE2 `alosStack`, then
used MintPy 1.5.3 for a weighted SBAS time series [E]. The paper's method
cannot be reproduced in full from three NISAR dates over 24 days [E].

| Tay method element | NISAR transfer | Current status |
|---|---|---|
| L-band HH phase observations | Use Frequency A HH GSLC/GUNW phase | implemented [E] |
| Ionosphere estimation on a bandwidth-4 network | Apply and audit NISAR's supplied screen for mission nearest-neighbour GUNW; for custom non-nearest GSLC pairs, derive per-epoch ionosphere from a dedicated bandwidth-4 dual-frequency or validated mission-supported route | implemented at 80 m for nearest pairs; custom network pending [E/S] |
| Bandwidth-20 interferogram network | Form GSLC pairs by acquisition-rank bandwidth and test convergence at 4, 8, 12 and 20 | three-date triangle only [E] |
| Inverse phase-variance/covariance weighting | Use MintPy `networkInversion.weightFunc=var` and preserve coherence plus effective looks | configured in commissioning inversion; effect requires a redundant network [E/S] |
| Closure-phase bias correction | Use MintPy closure correction only after the graph contains redundant cycles; retain uncorrected comparison | current triangle is diagnostic only [E] |
| Solid-Earth-tide correction | Compare the supplied NISAR SET phase with PySolid | mission stack implemented, independent comparison pending [E/S] |
| Spatially correlated troposphere correction from nearby non-deforming pixels | Implement as a separate candidate once a stable-control network exists; compare against mission screens, raw phase and ERA5 | not yet applicable [E/S] |
| Stable CPKY GNSS reference at 113.926 E, 2.220 S | Use CPKY only if its current stability is independently confirmed and it is connected to the AOI in the same unwrap component; otherwise require equivalent stable GNSS | both full GUNWs contain CPKY in connected component 1 with 5 × 5 median coherence 0.773 and 0.735; current fine crop excludes it, so a full-frame reference branch remains required [E/S] |
| Pixel-wise LOS-to-vertical projection | Divide by cosine of local incidence angle and retain LOS as the primary observable | implemented [E/S] |
| Mean phase-standard-deviation quality screen | Propagate coherence-derived phase variance through the inversion; evaluate Tay's 0.7 cm threshold rather than copying it across different looks and posting | pending redundant inversion [S] |
| Approximately 180 m multilooked science grid | Produce 20/40/80/180 m validation tiers; do not sacrifice the 5 m-posted diagnostic source or claim Tay-level accuracy at 5 m | planned and partly implemented [E/S] |

The temporal release audit encodes graph redundancy explicitly [E]. It records
cycle count, bridge-edge count, minimum acquisition degree, maximum connection
bandwidth, the largest fully populated bandwidth, and the fraction of
Tay-required pairs present. A long nearest-neighbour chain is connected but
has no independent closure information, so it cannot pass PeatGuard's annual
velocity gate. A Tay-network candidate must include every pair in the
bandwidth-20 network, rather than a single rank-20 edge, and still pass
correction, reference, coherence, uncertainty, closure and
independent-validation gates [E/S].

Mission GUNW is generated only for nearest-neighbour 12-day pairs [E]. Its
pair ionosphere screens therefore cannot simply be assigned to custom
bandwidth-20 GSLC interferograms. Before a Tay-compatible custom network is
inverted, PeatGuard must estimate per-epoch ionosphere on a dedicated
bandwidth-4 network using Frequency A/B or another validated mission-supported
method, then construct the correction for each longer pair [E/S].

Tay's custom tropospheric correction is especially important to test. The
paper reports that PyAPS/ERA5, GACOS and height-correlation corrections
worsened or failed to improve accuracy over its flat Central Kalimantan
study area, while a correction estimated from nearby non-deforming pixels
improved epoch-to-epoch agreement [E]. PeatGuard will therefore retain ERA5
as a comparison branch, not assume it is the preferred NISAR correction.
No tropospheric branch will become operational without stable-control or
field validation [E/S].

## Resolution architecture

PeatGuard will publish resolution as both **pixel posting** and **effective
support**. A 5 m output pixel does not imply an independent, reliable 5 m
deformation estimate.

| Product | Posting | Effective support | Role |
|---|---:|---:|---|
| GSLC complex amplitude/phase | 5 m | Native mode-dependent; posting is not resolution | Source only |
| Wrapped phase | 5 m | Inherits native mode-dependent support | Fine diagnostic |
| Local coherence, 5-pixel window | 5 m | about 25 m | Fine quality layer |
| Fine unwrapped displacement | target 20 m | 20–40 m | Experimental decision-support |
| Mission GUNW displacement | 80 m | about 80 m | Conservative validation layer |
| NISAR time-series velocity | initially 80 m | 80 m+ | Operational only after enough epochs |

The website must not show a two-epoch displacement as `velocity`, and must not
upsample an 80 m GUNW into a layer labelled 5 m.

## Implemented data path

### AOI staging

`peatguard.ingest.nisar_gslc`:

1. Searches ASF for path 125/frame 176/40 MHz GSLC products.
2. Resolves authenticated Earthdata URLs to short-lived cloud URLs.
3. Opens HDF5 through HTTP range requests.
4. Copies only AOI-intersecting 512 × 512 chunks.
5. Preserves HH, coordinate vectors, projection, identification metadata and
   mission quality masks.
6. Checkpoints completed rows in a `.part` HDF5 and resumes after interruption.
7. Renames to the final AOI file only after all rows are written.

This avoids downloading 61.1 GiB of full-frame data when only the PeatGuard
crop is needed.

### Fine interferogram

`peatguard.insar.nisar_gslc`:

1. Validates track, frame, CRS, shape and coordinate equality.
2. Reads complex HH data block-wise with a coherence halo.
3. Applies the mission GSLC valid-data mask from both dates. Positive values
   from 1 through each product's `numberOfSubSwaths` identify fully focused
   samples; 0 and 255 are invalid/partially focused or fill and are excluded,
   following Section 4.4 of the NISAR L2 GSLC specification [E].
4. Excludes every nonzero `inputDataExceptionMask` value. The current product
   metadata identifies value 2 as an LSAR qFSP-H1 sample-slip anomaly [E].
5. Forms `reference × conjugate(secondary)`.
6. Exports wrapped phase and normalized local coherence as tiled GeoTIFFs.
7. Uses nearest-neighbour phase overviews to avoid averaging across the
   `-pi/+pi` boundary.
8. Tags outputs as `wrapped_phase_not_displacement_or_velocity`.

Three GSLC epochs permit three fine-resolution pairs:

- 2026-06-22 to 2026-07-04
- 2026-07-04 to 2026-07-16
- 2026-06-22 to 2026-07-16

The 24-day third pair is not distributed as mission GUNW, but it provides a
phase-closure test for the two 12-day pairs [S].

All three subsets completed and passed full integrity checks [E]. Each contains
177,792,386 finite, nonzero HH complex samples on the same 5 m EPSG:32750
grid. The mission masks were read from row zero through row 12,863, and every
`inputDataExceptionMask` sample is zero.

## Mandatory real-data validation gates

No NISAR layer becomes decision-driving or deploys as a subsidence rate until
the applicable gates pass.

### Gate 1: Product integrity [E]

- Expected ASF byte count equals local byte count for every full GUNW.
- `h5py` opens each file and reads product type/version.
- GSLC subset source granule, AOI, coordinate spacing and row checkpoint are
  recorded.
- All pair grids match exactly before cross-multiplication.

### Gate 2: Geolocation [E]

- GSLC output EPSG and affine transform come directly from product coordinate
  vectors.
- Reproject canal intersections and stable landmarks from PlanetScope onto the
  5 m grid.
- Require no systematic displacement larger than one GSLC pixel before using
  the 5 m overlay.

### Gate 3: Phase sign and convention [E/S]

- Aggregate the 5 m complex interferogram to the matching GUNW 20 m wrapped
  grid using complex, not arithmetic phase, averaging.
- Compare both candidate signs against GUNW wrapped phase modulo `2 pi`.
- Select the documented convention once, then freeze it in metadata and tests.
- A constant phase reference may be removed for comparison. No multiplicative
  scaling or regression to Tay is allowed.

Result: the expected same-sign convention was selected for both 12-day pairs
[E]. After one constant circular offset, phase MAE was 0.379 and 0.399 rad and
the displacement-phase gate failed [E]. Fine phase remains diagnostic only.

### Gate 4: Coherence [E/S]

- Aggregate PeatGuard local coherence to 20 m and 80 m.
- Compare with mission GUNW coherence using valid-overlap bias, MAE, Pearson
  and Spearman correlation.
- Report coverage at thresholds 0.25, 0.35, 0.5 and 0.7.
- Inspect forest, degraded peat, canal-buffer and open-water strata separately.

Result: both supporting-coherence gates passed [E]. Pearson correlation with
mission GUNW was 0.913 and 0.902, absolute bias was 0.010 and 0.014, and MAE
was 0.055 and 0.057. The validated release is coherence only.

### Gate 5: Atmospheric and ionospheric correction [E]

- Apply the supplied GUNW ionosphere phase screen with the documented sign.
- Compare the supplied hydrostatic and wet troposphere screens against raw
  phase, ERA5 and, once stable controls exist, the Tay-style neighbouring
  non-deforming-pixel correction.
- Treat ERA5 as a comparison branch rather than a default truth because Tay
  et al. (2025) found global atmospheric-model corrections worsened its local
  field agreement [E].
- Compare the supplied solid-Earth-tide phase against PySolid.
- Keep raw and corrected products.
- Record correction magnitude and uncertainty. A downloaded sample had a mean
  ionosphere screen of about 3.6 radians, so this correction is material.

### Gate 6: Fine unwrapping [S]

- Start from mission 80 m connected components.
- Use coarse-to-fine unwrapping at 80 m, 40 m and target 20 m.
- Require low phase-residue density and no visible component-boundary steps.
- Test 24-day closure:
  `phase(22 Jun–4 Jul) + phase(4 Jul–16 Jul) - phase(22 Jun–16 Jul)`.
- Require at least 10,000 pixels passing the three-pair coherence mask, wrapped
  closure MAE at most 0.05 rad, RMSE at most 0.10 rad, resultant length at
  least 0.99, and at least 99% of pixels within 0.10 rad.
- Reject local 5/20 m areas whose closure exceeds the predeclared tolerance.

Result: native-grid closure passed over 114,216,083 coherence-qualified
pixels with MAE 0.000000064 rad and 100% within 0.01 rad [E]. This proves
phase algebra and grid integrity only. It does not override the failed
fine-to-mission displacement-phase gate.

### Gate 7: Displacement plausibility [S]

- Convert corrected phase using the NISAR wavelength and verified sign.
- Report LOS displacement first.
- Convert to vertical only where a vertical-only assumption is defensible and
  label it as an assumption.
- Reject isolated extreme values and report masking rules, rather than clipping
  them into plausible-looking rates.

### Gate 8: Time-series readiness [E/S]

- Three epochs are a commissioning stack, not an operational velocity series.
- Publish coherence as supporting. Publish pair displacement only if its
  separate correction, datum, uncertainty, and fine-to-mission gates pass.
- A connected nearest-neighbour chain is not sufficient. Require no graph
  bridges, at least two interferogram connections per acquisition, and
  redundant closure cycles before any trend candidate.
- Build and compare acquisition-rank bandwidths 4, 8, 12 and 20. Require
  convergence and test Tay's bandwidth-20 closure-bias route before calling
  the result method-compatible.
- Wait for at least 24 NISAR epochs and at least 365 days before even treating
  a redundant result as an annual-velocity candidate. Earlier short-span
  slopes must be labelled provisional diagnostics.

### Gate 9: Independent validation [E/S]

- Compare the eventual NISAR velocity spatial pattern to Tay et al. (2025)
  ALOS-2 on overlap.
- Report correlation, bias, MAE and shifted-grid sensitivity.
- Do not fit NISAR to Tay.
- Treat field dipwells/peat cameras as the later TRL 5–6 validation step.

## Alternate-source baseline decision

### Candidate A: Tay et al. (2025) ALOS-2 published peat motion

Availability:

- Open derived data: <https://doi.org/10.21979/N9/JGQE4O>
- Already present locally as `deliverables/tay_restoration/subsidence_tay.tif`.
- Same broad Central Kalimantan peat system.
- Dataset licence: CC BY-NC 4.0. A commercial PeatGuard product requires
  permission from the rights holder [E].

Source audit completed on 2026-07-23:

- The official 21,914,912-byte long-term-rate raster was downloaded directly
  from DR-NTU.
- Its MD5 is `32a196798d8b01cc3e9690c29e4804ac`, exactly matching the
  repository metadata [E].
- The source and PeatGuard derivative have identical dimensions, CRS,
  transform and bounds.
- Every one of PeatGuard's 101,979 retained pixels is exactly the source value
  multiplied by 1,000. Maximum and mean numerical difference are both zero.
  This is unit conversion, not calibration [E].
- PeatGuard excludes 144 nonzero source cells outside the project-defined
  `abs(rate) < 100 mm/yr` plausibility rule. Tay's repository README does not
  prescribe this threshold, so it must remain visible as a PeatGuard filter
  rather than being attributed to Tay [E].
- The source posting is about 158 by 208 m at the grid centre. It is not a
  few-metre deformation product and cannot independently validate 5 m NISAR
  motion [E].
- The repository calls the raster a linear long-term rate and records
  `[Units = m]`. Interpreting the slope as m/yr and converting to mm/yr is
  dimensionally consistent with the paper, but the README omits `/yr`; this
  metadata ambiguity is retained explicitly [S].

Machine-readable result:
`grid_products/external/nisar/baselines/alternate_baseline_audit.json`.

Strengths:

- L-band and independently processed.
- Field peat-camera comparison reported by the paper.
- Closest available bridge between historical ALOS-1 and current NISAR.

Limits:

- It is a derived product, not an acquisition stack PeatGuard controls.
- Its observation epoch differs from NISAR.
- It cannot serve as a fitted prior without destroying independence.

Verdict: **adopt now as the independent benchmark, never as a calibration
target or temporal extension of the NISAR network** [E]. It is also not a
substitute for a 2026 current rate because its averaging period is
2018-10-12 to 2024-08-09 [S].

### Candidate B: open ALOS-2 PALSAR-2 ScanSAR L1.1 complex data

NASA ASF now exposes JAXA's full-aperture ScanSAR Level 1.1 archive as the
public collection `C3315903479-ASF` [E]:

- Collection: <https://catalog.data.gov/dataset/alos-2-palsar-2-level-1-scansar-product>
- JAXA source and terms: <https://www.eorc.jaxa.jp/ALOS/en/dataset/alos_open_and_free_e.htm>
- The Level 1.1 complex I/Q samples preserve magnitude and phase and are
  explicitly intended for interferometric processing [E].
- The products are WBD dual-polarization ScanSAR, with separate HH and HV
  image files for each of five sub-swaths.

Live inventory checked on 2026-07-25 for the PeatGuard AOI:

- Descending Path 28 has 129 products from 2014-11-28 to 2025-03-07.
- The archive contains 91 of Tay's 92 acquisition dates. The only missing
  date is 2021-04-23.
- The archive records these products as frame 3700, while Tay et al. report
  frame 3650. The exact date agreement and footprint make this the likely
  matching acquisition series, but the CEOS leader and processed overlap
  must resolve this identifier difference before PeatGuard claims an exact
  reproduction [E/S].
- Three acquisitions postdate Tay's published endpoint: 2024-09-20,
  2024-10-04 and 2025-03-07.
- The three post-Tay dual-polarization ZIPs total 182.21 GB. Each granule is
  approximately 60.7 GB.

The ASF ZIPs support HTTP byte ranges and store every CEOS entry without ZIP
compression [E]. A verified 2024-08-09 granule contains five HH files ranging
from 5.26 to 7.61 GB, five equivalent HV files, and small leader, trailer and
metadata files. PeatGuard can therefore fetch only the HH sub-swath covering
the AOI instead of downloading each complete 60.7 GB archive. All selected
entries and processing intermediates must remain under
`grid_products/_deletable_work/baselines/alos2_scansar/`.

Proposed bounded pilots:

1. **Four-epoch transfer test:** 2024-08-09 plus the three post-Tay dates.
   Download only the required HH sub-swath and metadata, then confirm that
   ISCE2 `alosStack` can ingest the current full-aperture CEOS products.
   This is an ingestion, coherence and sign test, not a rate product [E].
2. **Thirteen-epoch endpoint pilot:** 2023-08-11 to 2025-03-07. The ten dates
   overlapping Tay provide a hold-out comparison and the three later dates
   test extension. Use a redundant network and transfer Tay's ionosphere,
   phase-covariance, closure, SET and troposphere tests where supported [S].
3. Treat the 154-day gap from 2024-10-04 to 2025-03-07 as a likely
   decorrelation risk. Do not connect or publish the March epoch unless its
   coherent network edges and uncertainty pass the same release gates [S].

Strengths:

- Same ALOS-2 L-band ScanSAR acquisition family, Path 28 and nearly identical
  date list as Tay.
- Raw phase is available, so PeatGuard can reproduce the method rather than
  infer deformation from a published rate raster.
- It provides the closest open phase-preserving observations to the first
  NISAR acquisitions currently found for this AOI.

Limits:

- Four new endpoint epochs cannot support a robust annual rate by themselves.
- The latest open Path 28 acquisition currently found is 2025-03-07, leaving
  a gap to NISAR in June 2026.
- Long temporal gaps may decorrelate over tropical peat.
- Even selective extraction still requires tens of gigabytes plus much larger
  ISCE2 intermediates.
- JAXA attribution and the distribution platform's terms must accompany any
  derived release. Commercial-use terms must be confirmed before a PeatGuard
  commercial product [E].

Verdict: **adopt as the highest-priority reproducible baseline pilot, ahead
of ALOS-1** [S]. Keep its time series separate from NISAR, use the published
Tay raster only as an independent hold-out, and do not label a four-epoch
result as a 2025 annual subsidence rate.

### Candidate C: ALOS-1 PALSAR-1 L1.1 SLC

ASF provides open L1.1 SLC with phase preserved [E]:

- Dataset DOI: <https://doi.org/10.5067/PAFUZQHBZF3A>
- Product guide: <https://docs.asf.alaska.edu/datasets/palsar/>

Direct peat precedent:

- Hoyt et al. (2020), *Nature Geoscience*, mapped 2007–2011 Southeast Asian
  peat subsidence at 90 m from ALOS PALSAR and found mean regional subsidence
  of 2.2 cm/yr [E].
- DOI: <https://doi.org/10.1038/s41561-020-0575-4>

Exact live ASF inventory:

#### Recommended pilot

- AOI: current high-priority canal pilot,
  `114.26–114.409976 E, 2.700072–2.57 S`.
- Path 422, frame 7130, FBD.
- 14 dates, 2007-07-09 to 2010-12-02.
- 100% pilot coverage.
- 18.27 GiB total L1.1 download.

#### Full operational AOI

The full 69 × 64 km operating area crosses adjacent paths/frames:

- Path 421/frame 7120: 12 dates, 15.74 GiB.
- Path 421/frame 7130: 12 dates, 15.61 GiB.
- Path 422/frame 7120: 14 dates, 18.23 GiB.
- Path 422/frame 7130: 14 dates, 18.27 GiB.
- Total: 52 granules, about 67.85 GiB before processing intermediates.

Strengths:

- Same L-band class as ALOS-2/NISAR.
- Multi-year time series with direct tropical-peat precedent.
- Free and reproducible from raw phase-preserving products.

Limits:

- Historical 2007–2010 rate is not current conditions.
- FBD resolution and 46-day cadence are coarser than NISAR.
- Processing and mosaicking four stacks substantially increases storage and
  quality-control burden.

Verdict: **adopt the one-stack canal pilot first**. Expand to four stacks only
if the pilot agrees sensibly with the published Tay spatial pattern and passes
coherence/atmosphere/unwrap checks [S].

### Candidate D: Hoyt et al. (2020) derived map

The paper is strong direct evidence and its data statement cites Zenodo DOI
`10.5281/zenodo.3694667`. The DOI, Zenodo record page and API returned 404 when
checked on 2026-07-23, so the cited derived raster is presently not a
dependable automated input [E].

Verdict: **use the paper as method precedent, not as a depended-on data
source while its repository is inaccessible**. Contact the authors or Zenodo
support later if their exact Central Kalimantan raster would materially reduce
processing [E/S].

### Candidate E: JERS-1

ASF has JERS-1 L-band archive coverage from the pre-/early-Mega Rice Project
period. Its value is historical land-state context.

Limits:

- Older orbit knowledge and product handling.
- Sparse/irregular stack over the AOI.
- A 1990s displacement estimate would not represent a current drainage regime.

Verdict: **do not use as a deformation baseline**. Consider only as a
pre-disturbance amplitude/context experiment after NISAR and ALOS-1 are stable
[S].

### Rejected as deformation baselines

- Sentinel-1 C-band velocity: independently shown by PeatGuard to be
  noise-dominated over this AOI.
- SRTM/Copernicus/TanDEM-X DEM differencing: metre-scale DEM errors cannot
  resolve centimetre-per-year peat motion over short spans [E].
- GEDI/ICESat-2 alone: sparse footprints and canopy/ground-return complications
  prevent a wall-to-wall peat velocity baseline [S].
- ALOS annual backscatter mosaics: amplitude products do not retain the phase
  needed for InSAR deformation [E].

## Baseline adoption sequence

1. Keep Tay as the current independent benchmark.
2. Resolve the ALOS-2 sub-swath covering the PeatGuard AOI from the CEOS
   metadata, then run the four-epoch selective-HH ingestion and coherence
   test.
3. If the transfer test passes, process the thirteen-epoch ALOS-2 endpoint
   pilot and compare its overlapping 2023–2024 rates and displacement pattern
   against Tay without fitting.
4. Publish no ALOS-2 extension unless the post-Tay epochs retain coherent,
   redundant connections and pass the correction and uncertainty gates.
5. Process ALOS-1 path 422/frame 7130 over the canal pilot as a separate
   2007–2010 historical layer.
6. Compare ALOS-1, ALOS-2 and later NISAR by spatial patterns and strata, not
   pixelwise cross-epoch fusion.
7. Keep each sensor-period segment as a separate layer and provenance record.

## Historical and NISAR fusion experiment

An uncertainty-aware fusion trial was completed on 2026-07-25 with:

- Tay et al. (2025) 2018–2024 long-term vertical rate.
- NISAR mission-corrected relative vertical displacement from 2026-06-22 to
  2026-07-16.
- NISAR's conservative ionosphere-plus-reference-patch uncertainty. This is
  still a partial uncertainty because it excludes troposphere, unwrapping,
  orbit, DEM, horizontal-motion and absolute-datum error.

The experiment never fitted a bias, scale or datum to Tay [E]. It produced:

1. A Gaussian inverse-variance sensitivity calculation in which Tay is an
   explicit historical prior and the annualized NISAR commissioning
   displacement is a diagnostic likelihood.
2. A categorical evidence map that preserves the two periods and labels
   historical subsidence, historical non-subsidence, current NISAR
   inconclusive coverage, and significant short-span directional motion
   separately.

Evaluation was performed on Tay's native support to avoid treating an
upsampled 80 m-posted map as independent 80 m historical measurements [E]:

- Common cells: 11,364.
- Expected 24-day displacement from the historical rate: spatial standard
  deviation 0.801 mm.
- Observed NISAR 24-day displacement: spatial standard deviation 12.605 mm.
- Pearson correlation 0.103; Spearman correlation 0.150.
- RMSE 17.342 mm; sign agreement 43.7%.
- No common cell exceeded the conservative positive or negative 2-sigma
  threshold.

The Bayesian sensitivity used the 3.6 mm/yr eight-site Tay rate RMSE as the
historical measurement term [E] and hypothetical between-period change terms
of 20, 50 and 100 mm/yr [H]. Under the reference 50 mm/yr scenario, the median
NISAR weight was only 0.421%, and even this is an upper bound because the
NISAR uncertainty is incomplete. The posterior remained almost identical to
Tay (Pearson 0.997). It therefore does not provide a current 2026 rate [E/S].

The categorical map is the retained interpretation [S]:

- Historical subsidence with current NISAR inconclusive: 5,943 pixels on the
  80 m-posted output grid.
- Historical non-subsiding/rebound with current NISAR inconclusive: 37,197
  pixels.
- Significant current directional motion with historical overlap: zero.
- Sixteen significant downward and one significant upward 80 m-posted pixels
  occur outside Tay coverage. They remain commissioning flags rather than
  subsidence claims because the uncertainty is partial and the clusters are
  extremely small.

Compact outputs are under
`grid_products/external/nisar/historical_fusion/`. The machine-readable audit
is `historical_nisar_fusion_audit.json`. The quantitative candidate rasters
are tagged `do_not_publish_as_current_subsidence`; only the categorical map is
appropriate for internal commissioning context [E/S].

## NISAR and PlanetScope complementarity

PlanetScope remains the geometry source for narrow canals. NISAR should improve
confidence and prioritisation, not blur a 3 m canal mask to radar resolution.

### Immediate integration

- Reproject Planet canal centerlines onto the NISAR 5 m grid.
- Compute NISAR coherence and wrapped-phase statistics in canal-distance bands
  (0–25, 25–100, 100–300 and 300–1000 m).
- Add a per-canal `nisar_evidence` score describing radar coherence, valid
  displacement coverage and observed cross-canal gradient.
- Preserve Planet geometry unchanged.

Implemented result [E]:

- Fine NISAR coherence contains 177,504,771 valid 5 m-posted pixels and uses
  about 25 m support.
- In mapped peat, the 0 to 100 m minus 300 to 1,000 m block-mean coherence
  difference is -0.02678, with a 95% spatial-bootstrap interval of
  -0.04445 to -0.01031 using 4 km blocks.
- This is an observability association, not evidence that canals caused the
  coherence difference.

### Later integration

- Use NISAR GCOV or carefully normalized GSLC amplitude as cloud-independent
  support for water/canal persistence.
- Flag Planet canal candidates unsupported by either radar darkness,
  interferometric response or network continuity for manual review.
- Detect post-blocking changes with Planet visual evidence plus NISAR
  displacement/coherence response.

### Prohibited shortcut

Do not feed provisional NISAR displacement into the Planet canal detector as
training truth. That would create circular evidence and conflate canal
location with drainage response.

## Website release gates

Current web layers:

- `NISAR minimum-pair coherence (20 m posting)` as supporting evidence.
- `NISAR fine minimum-pair coherence (5 m posting, ~25 m support)` as
  separately validated supporting evidence.
- No displacement or velocity, because the three-epoch, two-edge commissioning
  network has no redundancy and the propagated correction uncertainty is
  larger than the observed short-span signal [E].

Later, after applicable validation gates:

- `NISAR wrapped phase (5 m, diagnostic)` behind an experimental/provisional
  label.
- `NISAR corrected displacement (20 m experimental)`.
- `NISAR velocity` only after time-series readiness.
- `NISAR evidence` in dam-site/canal rationale.

Every NISAR layer should expose:

- Product maturity: PROVISIONAL.
- Acquisition pair/span.
- Pixel posting and effective support.
- LOS or assumed-vertical component.
- Correction status.
- Valid-data/coherence mask.
- A direct note that it is not ground truth.
