# PeatGuard NISAR scientific and release audit

Audit date: 2026-07-24
Study area: Sebangau and the former Mega Rice Project, Central Kalimantan
Release status: NISAR coherence support is releasable; NISAR displacement and
velocity are not releasable

## Confidence notation

- **[E] Established:** directly supported by mission documentation, a named
  paper, an exact data audit, or a reproducible PeatGuard test.
- **[S] Standard inference:** a conventional interpretation that follows from
  established methods but has not yet been independently validated for this
  exact PeatGuard product.
- **[H] Hypothesis:** a useful proposition that still needs direct testing.

## Executive decision

PeatGuard should integrate the current NISAR release as a **radar-observability
and quality layer**, not as a subsidence-rate replacement [E]. Two deployable
supporting products are available: minimum coherence on the mission GUNW
20 m grid, and independently formed GSLC coherence at 5 m posting with about
25 m support. Both span the two available 12-day pairs and are labelled
`PROVISIONAL` [E].

PeatGuard must not publish a NISAR velocity yet [E]. The selected stack has
three acquisitions, two nearest-neighbour interferograms, a 24-day span, and
no independent network cycle. The fully corrected relative displacement is
also less stable than the correction and datum sensitivity:

- Median fully corrected vertical-equivalent displacement over quality pixels:
  -7.80 mm over 24 days.
- Median lower-bound propagated ionosphere uncertainty: 44.00 mm.
- Median conservative reference-aware uncertainty: 54.00 mm.
- Relative-datum 5th to 95th percentile span: 15.15 mm LOS.
- Fraction of quality pixels exceeding two-sigma: below 0.02%.

These values fail the commissioning release gate by a wide margin [E]. They do
not prove that NISAR is unsuitable for peat subsidence. They show that the
current three-date commissioning network is too short and too weakly
constrained to support a PeatGuard subsidence rate [E].

The independent baseline remains the Tay et al. (2025) ALOS-2 L-band
long-term-rate product [E]. It is appropriate for spatial benchmarking, but
not for calibrating NISAR, extending the NISAR time series, representing a
current 2026 rate, or validating few-metre motion [E/S].

PeatGuard does, however, transfer Tay's processing methodology where
applicable [E/S]. The current implementation uses L-band HH phase,
coherence-derived inverse phase-variance weighting, explicit ionosphere,
troposphere and solid-Earth-tide variants, pixel-wise incidence-angle
projection, and an independent Tay comparison. The remaining Tay elements
need a much larger NISAR stack: bandwidth-20 network convergence,
closure-phase bias correction, a stable GNSS-controlled datum, Tay's local
non-deforming-pixel troposphere correction, and a multi-season velocity [E].

## Evidence base

### Literature map

The associated literature-map corpus retrieved 840 unique records and retained
305 core records, of which 254 had abstracts. It includes:

- 219 records tagged InSAR or radar subsidence.
- 160 tagged validation or uncertainty.
- 117 tagged ionosphere or atmosphere correction.
- 102 tagged time-series or network inversion.
- 70 tagged water table or hydrology.
- 38 tagged ALOS or PALSAR.
- 32 records geographically tagged Indonesia or Kalimantan.
- 29 tagged tropical peatlands.
- 6 tagged NISAR.

The corpus is a structured discovery set, not a formal PRISMA systematic
review [S]. Method decisions below rely on the named primary papers and mission
documentation, while the corpus is retained for traceability:

- `nisar_peatguard_bibliography.csv`
- `nisar_peatguard_literature_summary.json`
- `nisar_peatguard_literature_overview.png`
- `nisar_peatguard_theme_cooccurrence.png`

### Primary references

1. Tay et al. (2025), *Communications Earth & Environment*,
   [doi:10.1038/s43247-025-02926-6](https://doi.org/10.1038/s43247-025-02926-6).
   This is the independent ALOS-2 L-band benchmark, not a calibration target.
2. Hoyt et al. (2020), *Nature Geoscience*,
   [doi:10.1038/s41561-020-0575-4](https://doi.org/10.1038/s41561-020-0575-4).
   This provides direct ALOS-1 tropical-peat subsidence precedent.
3. Zhou et al. (2019), *Remote Sensing*,
   [doi:10.3390/rs11212592](https://doi.org/10.3390/rs11212592).
   This provides a Central Kalimantan ALOS-1 time-series and phase-closure
   precedent.
4. Zahro et al. (2025), *Scientific Reports*,
   [doi:10.1038/s41598-025-08390-8](https://doi.org/10.1038/s41598-025-08390-8).
   This supports L-band time-series evaluation of tropical peat restoration.
5. Umarhadi et al. (2021), *Land Degradation & Development*,
   [doi:10.1002/ldr.4061](https://doi.org/10.1002/ldr.4061).
   This directly compares C-band and L-band SBAS for peat subsidence.
6. Yunjun et al. (2019), *Computers & Geosciences*,
   [doi:10.1016/j.cageo.2019.104331](https://doi.org/10.1016/j.cageo.2019.104331).
   This is the named MintPy small-baseline inversion and noise-reduction basis.
7. Dadap et al. (2021), *AGU Advances*,
   [doi:10.1029/2020AV000321](https://doi.org/10.1029/2020AV000321).
   This is the optical canal-geometry reference. It is not the source of the
   Tay radar-subsidence product.
8. Hooijer et al. (2012), *Biogeosciences*,
   [doi:10.5194/bg-9-1053-2012](https://doi.org/10.5194/bg-9-1053-2012).
   This is the basis for interpreting long-term peat subsidence in relation to
   carbon loss, with explicit process caveats.

Mission sources:

- [NISAR data availability](https://nisar-docs.asf.alaska.edu/availability-overview/)
- [NISAR PROVISIONAL known issues](https://nisar-docs.asf.alaska.edu/provisional-known-issues/)
- [NISAR GUNW product guide](https://nisar-docs.asf.alaska.edu/gunw/)
- [NISAR GSLC product guide](https://nisar-docs.asf.alaska.edu/gslc/)
- [NISAR L2 GSLC Product Specification, Rev E, JPL D-102269](https://nisar.asf.earthdatacloud.nasa.gov/NISAR-SAMPLE-DATA/DOCS/NISAR_D-102269_RevE_NASA_SDS_Product_Specification_L2_GSLC_Nov8_2024_w-sigs.pdf)
- [NISAR phase conventions](https://nisar-docs.asf.alaska.edu/phase/)
- [ASF ALOS PALSAR archive guide](https://docs.asf.alaska.edu/datasets/palsar/)

## Sensor finding

PeatGuard's full Sentinel-1 ISCE2 and MintPy replication found that C-band
velocity is noise-dominated over this AOI [E]. The Sebangau run had an
approximately -0.03 Pearson correlation with the independent Tay ALOS-2
surface, and the method-faithful Block B replication remained near -0.04 even
with higher coherence [E]. The failure persisted after georeferencing,
low-canopy, degraded-peat, and near-canal checks [E].

This is consistent with the physical and published case for L-band over
vegetated tropical peat [E]. It does not mean Sentinel-1 is universally
incapable of peat monitoring. It means PeatGuard must not use the local
Sentinel-1 velocity as a deformation baseline for this intact-canopy operating
area [E].

## NISAR product maturity and resolution

NISAR observations acquired from 2026-06-17 onward are released at
`PROVISIONAL` maturity: calibrated and partially validated [E]. Known issues
include ionospheric effects, radio-frequency interference, SLC ripple,
Frequency B calibration limitations, incomplete products, and some
mode-transition artefacts [E]. PeatGuard therefore propagates maturity in the
catalog, UI, audit JSON, and layer explanation.

Resolution must be reported as both posting and effective support:

| Product | Pixel posting | Effective support | Valid role now |
|---|---:|---:|---|
| GSLC complex HH | 5 m | mode and processing dependent; posting is not resolution | source |
| GSLC cross-product wrapped phase | 5 m | inherits the native mode-dependent support | experimental |
| Local coherence, 5-pixel window | 5 m | about 25 m | validated supporting quality |
| Mission GUNW wrapped phase/coherence | 20 m | about 20 m | quality/diagnostic |
| Mission GUNW unwrapped phase/coherence | 80 m | about 80 m | conservative commissioning |
| NISAR velocity | not released | 80 m or coarser initially | future |

A 5 m GeoTIFF does not create an independently accurate 5 m subsidence
measurement [E]. The official NISAR Data User Guide establishes that the
selected 40 MHz Frequency A GSLC is posted on a 5 by 5 m grid, while effective
resolution remains mode and processing dependent [E]. The implemented route
to a few-metre-posted product is phase-preserving GSLC
cross-multiplication. It yields 5 m-posted wrapped phase and 5 m-posted local
coherence with roughly 25 m support [E]. The two coherence fields agree
strongly with mission GUNW, but the fine phase does not pass PeatGuard's
displacement-phase gate. A defensible unwrapped displacement should therefore
initially be solved and validated at about 20 to 40 m support [S]. PeatGuard
must not upsample the 80 m mission displacement and call it 5 m.

## Current NISAR commissioning stack

Selected acquisitions:

- Path 125, frame 176, ascending.
- Frequency A, 40 MHz.
- HH polarization.
- 2026-06-22, 2026-07-04, and 2026-07-16.
- Exact wavelength used: 0.2419632429378531 m.

The mission GUNW stack was ingested through the native MintPy NISAR loader at
upstream revision `53030ef4ae390ff17e845a342e086993eca81115` [E]. Source and
derived SHA-256 values were verified before analysis [E]. Exact inversion
round-trip tests passed for raw, ionosphere-corrected, and full
mission-corrected variants with maximum absolute errors below
0.000000059 m [E].

The fine GSLC path requires the mission valid-data mask in both acquisitions
[E]. Following Section 4.4 of the NISAR L2 GSLC Product Specification, only
positive subswath identifiers from 1 through the product's
`numberOfSubSwaths` contribute to phase or local coherence. Values 0 (invalid
or partially focused) and 255 (fill) are excluded [E]. The three live products
used here each report one subswath. Coherence is not spatially filled into a
centre pixel that the mission mask marks invalid.
The live products define `inputDataExceptionMask` value 0 as no anomaly and
value 2 as a NISAR LSAR qFSP-H1 sample slip. PeatGuard preserves this mask and
excludes every nonzero exception value from phase and coherence [E].

All three staged GSLC subsets share the same EPSG:32750 grid and
12,863 by 13,822 shape [E]. Across 533,377,158 complex samples, every real and
imaginary component is finite and no complex sample is zero. All three
exception masks contain only zero. The repaired 22 June valid-data mask was
verified from its first row rather than trusted from a legacy resume
checkpoint [E].

Quality coverage:

- 11,107,825 valid 20 m coherence pixels.
- 34.62% of valid 20 m pixels have minimum pair coherence at least 0.5.
- 210,506 quality pixels at 80 m, or 30.96% of mapped land.
- 75,093 quality mapped-peat pixels at 80 m, or 40.61% of mapped peat.

The high coherent-pixel count supports publishing the coherence layer [E].
Coherence is a quality proxy, not displacement, subsidence, hydrology, or
ground truth [E].

Fine-product validation:

- Three 5 m-posted GSLC pairs were formed for 22 June to 4 July, 4 July to
  16 July, and 22 June to 16 July [E].
- Native-grid closure passed over 114,216,083 pixels where all three pair
  coherences were at least 0.25. Closure MAE was 0.000000064 rad and every
  evaluated pixel was within 0.01 rad [E]. This is an internal algebra and
  grid-integrity result, not geophysical validation.
- Against matching mission 20 m GUNW, local coherence Pearson correlations
  were 0.913 and 0.902; biases were -0.010 and -0.014; MAEs were 0.055 and
  0.057 [E]. Both supporting-coherence gates passed.
- Fine wrapped phase selected the expected sign, but after one allowed
  constant circular offset its MAE was 0.379 and 0.399 rad, RMSE was 0.501
  and 0.523 rad, and resultant length was 0.886 and 0.876 [E]. Both
  displacement-phase gates failed.
- The likely distinction is pair-specific GUNW coregistration and
  geometry-offset flattening beyond direct same-grid GSLC
  cross-multiplication [S]. NISAR documentation confirms that mission
  interferogram formation refines geometry-derived coregistration with
  data-driven offsets [E].

The public fine product is consequently limited to minimum-pair coherence. It
contains 177,504,771 valid pixels, with 31.49% at or above coherence 0.5. It
stores coherence as UInt8 with a 1/254 scale, retaining about 0.004 precision
in a formally ordered 217 MB COG [E].

## Tay-method compatibility audit

Tay et al. (2025) used 92 ALOS-2 PALSAR-2 ScanSAR L1.1 images over nearly six
years, ISCE2 `alosStack`, an ionosphere-estimation bandwidth of 4, a
closure/inversion network expanded to bandwidth 20, inverse
phase-variance/covariance weighting, MintPy closure-phase correction,
solid-Earth-tide correction, a stable CPKY GNSS reference, pixel-wise
LOS-to-vertical projection, and a locally estimated spatial troposphere
correction [E].

The following transfer decision is enforced:

| Element | Decision |
|---|---|
| L-band HH observations | adopted [E] |
| Inverse phase-variance weighting | configured through MintPy `weightFunc=var`; its effect requires a redundant network [E/S] |
| Raw and ionosphere/troposphere/SET variants | adopted for mission nearest-neighbour GUNW pairs [E] |
| Tay bandwidth-4 per-epoch ionosphere solution | future custom GSLC networks must estimate per-epoch ionosphere from a dedicated bandwidth-4 dual-frequency or validated mission-supported route; nearest-pair screens must not be extrapolated silently [E/S] |
| Pixel-wise incidence-angle projection | adopted, with vertical-only caveat [E/S] |
| Three-pair phase closure | diagnostic gate only [E] |
| Bandwidth-20 convergence and closure-bias correction | wait for a redundant multi-season network [E] |
| CPKY stable reference | both full GUNWs contain CPKY in connected component 1 with 5 × 5 median coherence 0.773 and 0.735; current GNSS stability still requires confirmation and the present fine crop excludes the station [E/S] |
| Tay local troposphere correction | future comparison branch requiring stable non-deforming controls [E/S] |
| Tay 180 m multilooking | validation tier, not a reason to relabel a 5 m-posted diagnostic as 180 m or vice versa [E] |
| Tay raster values | independent benchmark only; never fitted into NISAR [E] |

The temporal-readiness code now rejects a long nearest-neighbour tree as a
velocity network even when it spans a year [E]. It records graph cycles,
bridge pairs, minimum acquisition degree, maximum acquisition-rank bandwidth,
the largest fully populated bandwidth, and the fraction of Tay-required pairs
present. Annual readiness requires at least 24 acquisitions, at least 365 days
and a bridge-free redundant graph. Separate Tay-network geometry requires every
pair in the bandwidth-20 network, not merely one rank-20 edge. Passing geometry
is necessary but not sufficient:
correction convergence, closure residuals, uncertainty, reference stability,
coherence and independent validation still gate release [E/S].

Tay's atmospheric result also changes PeatGuard's correction policy. The
paper found that PyAPS/ERA5, GACOS and height-correlation approaches worsened
or failed to improve field agreement in the flat study area, while its nearby
non-deforming-pixel correction improved epoch-to-epoch results [E].
Consequently, ERA5 is retained as a sensitivity branch rather than assumed to
be truth. The current three-date stack cannot identify stable non-deforming
time series robustly, so no Tay-style tropospheric correction is applied as
an operational claim yet [E/S].

## Corrections and uncertainty

NISAR GUNW supplies ionospheric, wet-tropospheric, hydrostatic-tropospheric,
and solid-Earth-tide correction layers, but does not apply them to the
interferometric phase by default [E]. PeatGuard prepared raw,
ionosphere-corrected, and fully corrected variants separately.

The correction sensitivity is material:

- Median total correction sensitivity: 24.17 mm vertical-equivalent.
- Median non-ionosphere correction sensitivity: 5.12 mm
  vertical-equivalent.
- Median lower-bound ionosphere uncertainty: 44.00 mm
  vertical-equivalent.

For a target peat signal near 10 mm/yr, correction status and uncertainty
cannot be treated as metadata decoration [E]. A future velocity release must
retain raw and corrected products, compare mission screens with ERA5,
PySolid and the Tay-style local-control route where applicable, and report
the differences [S].

## Alternate baseline assessment

### Tay et al. (2025) ALOS-2

The official DR-NTU long-term-rate raster was independently downloaded and
audited [E]:

- Dataset DOI: [10.21979/N9/JGQE4O](https://doi.org/10.21979/N9/JGQE4O).
- File size: 21,914,912 bytes.
- MD5: `32a196798d8b01cc3e9690c29e4804ac`, matching repository metadata.
- Period: 2018-10-12 to 2024-08-09.
- Source and PeatGuard derivative have identical grids.
- Every retained PeatGuard value is exactly source value times 1,000.
- No bias fit, regression, or force calibration was applied.

PeatGuard's derivative excludes 144 nonzero source cells at or beyond
`abs(rate) = 100 mm/yr`. That is a visible PeatGuard plausibility policy, not
a Tay threshold [E].

The source pixels are approximately 158 by 208 m at the grid centre, consistent
with the paper's roughly 180 m effective product [E]. It cannot validate
independent 5 m deformation. Tay et al. report field-comparison results for
their long-term rate, including a 0.36 cm/yr across-site RMSE and a maximum
absolute site difference of 0.59 cm/yr in the detailed validation [E]. Those
figures support regional benchmarking, not millimetre-scale truth at every
pixel.

Appropriateness verdict:

- Independent spatial benchmark: appropriate [E].
- Calibration target or fitted prior: not appropriate [E].
- Temporal extension of the NISAR network: not appropriate [E].
- Current 2026 rate substitute: not appropriate [S].
- Few-metre motion validator: not appropriate [E].
- Public non-commercial demonstration: conditionally appropriate under
  CC BY-NC 4.0 [E].
- Commercial product use: permission required [E].

### Open ALOS-2 ScanSAR L1.1 and fusion result

NASA ASF collection `C3315903479-ASF` now provides open, full-aperture,
phase-preserving ALOS-2 ScanSAR Level 1.1 products over PeatGuard [E]. A live
catalog audit found 129 descending Path 28 acquisitions from 2014-11-28 to
2025-03-07. The archive contains 91 of Tay's 92 dates and three post-Tay dates:
2024-09-20, 2024-10-04 and 2025-03-07 [E].

This is the highest-priority reproducible historical pilot because it closely
matches Tay's sensor, orbit and acquisition series [S]. ASF ZIP byte-range
reads were verified, allowing PeatGuard to retrieve a single 5.26 to 7.61 GB
HH sub-swath instead of each approximately 60.7 GB dual-polarization granule
[E]. The apparent ASF frame 3700 versus Tay frame 3650 identifier difference
must be resolved from CEOS metadata and processed overlap before claiming an
exact reproduction [E/S].

An uncertainty-aware Tay and NISAR fusion trial was also completed [E]:

- 11,364 native Tay-support cells overlap the quality-masked NISAR product.
- Historical expected 24-day displacement has 0.801 mm spatial standard
  deviation, while observed NISAR has 12.605 mm.
- Pearson correlation is 0.103, Spearman correlation is 0.150, RMSE is
  17.342 mm, and no common cell exceeds the conservative 2-sigma threshold.
- With a hypothetical 50 mm/yr between-period change sigma [H], NISAR's
  median inverse-variance weight is only 0.421%, itself an upper bound because
  several uncertainty terms are omitted.

The quantitative posterior is therefore a Tay-dominated historical map, not
a current 2026 rate [E/S]. PeatGuard retains only a categorical internal
evidence interpretation: historical subsidence with current NISAR
inconclusive, historical non-subsiding/rebound with current NISAR
inconclusive, or a significant short-span directional flag. Tay was not used
to fit the NISAR bias, scale, or datum [E].

### ALOS-1 PALSAR-1

Phase-preserving ALOS-1 L1.1 SLC remains the best reproducible deep-historical
stack for 2007–2010, behind the open ALOS-2 pilot for near-historical
continuity [S]. Hoyt et al. (2020) mapped 2007 to 2011 Southeast Asian peat
subsidence at 90 m [E], and Zhou et al. (2019) processed Central Kalimantan
ALOS-1 tracks with phase-closure and overlap checks [E].

The recommended PeatGuard pilot is path 422/frame 7130 over the high-priority
canal subarea: 14 dates from 2007-07-09 to 2010-12-02 and 18.27 GiB of L1.1
downloads [E]. Expanding immediately to all four overlapping stacks would
require about 67.85 GiB before processing intermediates [E].

Use ALOS-1 as a separate historical epoch, not as pseudo-observations in a
2026 NISAR inversion [E]. Compare spatial patterns and strata across epochs,
while allowing real land-cover, fire, drainage, and restoration change [S].

The pilot remains operationally blocked until the current `sbas` VM host key
is independently confirmed. PeatGuard will not bypass or delete a conflicting
SSH known-host entry without that confirmation.

### Other candidates

- Hoyt et al. derived Zenodo map: method precedent only for now. The cited DOI
  `10.5281/zenodo.3694667` returned 404 on 2026-07-23 [E].
- JERS-1: possible pre-disturbance amplitude context, not a current
  deformation baseline [S].
- Sentinel-1 velocity: rejected locally as noise-dominated [E].
- DEM differencing: rejected for centimetre-per-year motion because DEM errors
  are metre scale [E].
- GEDI or ICESat-2 alone: rejected as a wall-to-wall velocity source [S].
- Annual ALOS backscatter mosaics: rejected for deformation because they do
  not preserve interferometric phase [E].

## NISAR and PlanetScope result

PlanetScope and NISAR measure different quantities. PlanetScope supplies canal
geometry from native roughly 3 m optical imagery, exported in PeatGuard on a
roughly 5.6 m raster grid. NISAR supplies radar phase stability and eventually
ground-motion response. PeatGuard must preserve the Planet canal geometry and
attach NISAR evidence rather than blur the canal network to radar resolution
[E/S].

The fine coherence analysis covers 177,504,771 valid NISAR pixels,
7,056,300 canal pixels on the 5 m grid, and 47,336,312 mapped-peat pixels [E].
The product has 5 m posting and about 25 m effective support. Median minimum
coherence generally rises with distance from mapped canals:

| Distance to canal | All valid median | Mapped peat median |
|---|---:|---:|
| 0 to 25 m | 0.362 | 0.366 |
| 25 to 100 m | 0.382 | 0.390 |
| 100 to 300 m | 0.409 | 0.461 |
| 300 to 1,000 m | 0.425 | 0.496 |
| over 1,000 m | 0.488 | 0.508 |

A spatial block bootstrap comparing 0 to 100 m with 300 to 1,000 m gives:

- All valid: block mean difference -0.01366, 95% interval
  -0.02432 to -0.00294.
- Mapped peat: block mean difference -0.02678, 95% interval
  -0.04445 to -0.01031.
- SBG-NE mapped peat: block mean difference -0.03727, 95% interval
  -0.07741 to 0.00515.

The fine analysis preserves the earlier result that radar observability is
lower near canals overall, and now resolves a negative mapped-peat
association at the predeclared 4 km bootstrap-block scale [E]. It still does
not independently validate or improve Planet canal geometry, identify
causality, or prove drainage impact [E]. Regional results are heterogeneous.
Plausible contributors include open water, wet-edge scattering, vegetation
differences, and drainage-related surface change [H]. The analysis cannot
identify which mechanism dominates.

Operational decision:

- Keep PlanetScope as the canal-geometry source [E].
- Use NISAR coherence as a future per-reach confidence/manual-review feature
  [S].
- Do not use current NISAR displacement as Planet training truth [E].
- Do not infer canal causality from coherence alone [E].

## Scientific code-policy correction

The full code audit found and removed two legacy behaviours that were not
scientifically defensible [E]:

1. `run_analysis_stage()` shifted every InSAR velocity field so its
   high-coherence mean equalled a hard-coded -25 mm/yr literature target.
   A literature mean is not a geodetic datum. PeatGuard now preserves the
   observation-derived relative reference. Any future absolute adjustment
   requires an explicit stable control area, GNSS, or another independently
   justified reference [E/S].
2. The same stage automatically multiplied subsidence by a universal
   `0.5 tCO2 ha-1 yr-1 per mm yr-1` factor and emitted a carbon-loss raster.
   Hooijer et al. (2012) requires separation of oxidation, compaction, and
   consolidation, plus peat properties and long-duration observations. The
   automatic raster and basin total were removed [E].

The canal-proximity radius, linear decay, and 0.45/0.55
proximity/subsidence blend remain visible project decision heuristics [H].
Hooijer et al. (2012) supports the process interpretation linking drainage,
long-term subsidence, and carbon loss under stated conditions; it does not
prescribe those numerical risk parameters [E]. They require sensitivity
reporting now and field hydrological calibration before being interpreted as
a physical response model [S].

Regression tests enforce the no-forced-velocity and no-universal-carbon
conversion policies [E].

A 15-scenario sensitivity audit varied the project canal-influence radius
from 300 to 1,500 m and the proximity weight from 0.30 to 0.60 over
1,706,541 finite Tay-covered pixels [E]. Relative to the current
1,200 m / 0.45 proximity reference:

- Scenario-wide Spearman rank correlation remained between 0.805 and 0.990.
- Top-10% spatial Jaccard overlap ranged from 0.311 to 0.924.
- Top-5% spatial Jaccard overlap ranged from 0.281 to 0.913.
- The largest mean absolute score change was 0.119 on the 0 to 1 scale.

Broad priority structure is therefore moderately robust, but exact top-site
membership is parameter sensitive [E]. PeatGuard should use the surface to
shortlist field surveys, retain the parameter values in provenance, and avoid
claiming that a site's exact rank is physically calibrated [E/S]. The
machine-readable audit is
`grid_products/external/nisar/baselines/restoration_priority_sensitivity.json`.

## Carbon and restoration interpretation

InSAR measures displacement, not water-table depth and not carbon emissions
[E]. Hooijer et al. (2012) supports a relation between long-term subsidence and
carbon loss in drained tropical peat, but early subsidence can include
consolidation and compaction, while later rates more strongly reflect oxidation
[E]. Fire, land use, peat properties, drainage history, and time since drainage
also matter [E/S].

PeatGuard may use subsidence as one drainage-degradation indicator in a
restoration-priority model [S]. It must not present a short NISAR displacement
or a single universal conversion as verified carbon credits [E]. Carbon claims
require the applicable tropical-peatland methodology, field measurements, and
uncertainty accounting.

## Release matrix

| Product | Release | Reason |
|---|---|---|
| NISAR 20 m minimum-pair coherence | yes, supporting | adequate coverage; honest quality meaning [E] |
| NISAR 5 m-posted minimum-pair coherence | yes, supporting | both coherence-comparison gates passed; about 25 m support [E] |
| NISAR 80 m corrected displacement | no | uncertainty and datum sensitivity exceed signal [E] |
| NISAR short-span velocity | no | 3 dates, 24 days, no cycle, seasonal aliasing [E] |
| NISAR 5 m wrapped phase | experimental only | diagnostic, not unwrapped motion [E] |
| NISAR 20 to 40 m fine displacement | future experimental | requires coarse-to-fine unwrap and closure [S] |
| Tay ALOS-2 long-term rate | independent benchmark | suitable with licence and epoch caveats [E] |
| Historical/NISAR categorical evidence | internal only | preserves epochs but current NISAR is inconclusive [S] |
| Historical/NISAR Bayesian candidate rate | no | NISAR median weight is 0.421%; output remains historical-prior dominated [E/S] |
| Open ALOS-2 endpoint rate | pilot planned | phase-preserving bridge through 2025, subject to coherence and network gates [S] |
| ALOS-1 historical rate | pilot planned | separate historical epoch [S] |
| PlanetScope canal output | yes, geometry | native roughly 3 m imagery, deployed roughly 5.6 m grid [E] |

## Conditions for a future NISAR velocity

Before a NISAR velocity becomes decision-driving, PeatGuard should require:

1. A temporally redundant interferogram network with independent cycles [E].
2. A span long enough to separate seasonal motion and atmospheric structure
   from a long-term rate, preferably including multiple seasons [S].
3. Stable phase sign and correction conventions frozen in tests [E].
4. Per-pixel correction uncertainty and relative-datum sensitivity [E].
5. Fine-to-mission phase and coherence validation [E/S].
6. Closure, residue, connected-component, and geolocation checks [E/S].
7. Independent comparison with Tay without fitting [E].
8. Field dipwell, peat-camera, or GNSS validation for TRL 5 to 6 [S].

## Current limitations and requested inputs

1. Confirm the changed ED25519 host key for the `sbas` VM through the OVH
   console before the ALOS-1 pilot is started.
2. Obtain rights-holder permission before using the Tay CC BY-NC raster or
   Planet Education and Research imagery in a commercial product.
3. Provide or plan field reference points, preferably stable GNSS plus
   dipwell/peat-camera stations across land-cover and drainage strata.
4. Decide whether future public products should show LOS displacement only or
   also an explicitly assumed vertical-equivalent layer.
5. Keep all raw downloads and processing intermediates in
   `grid_products/_deletable_work/`; remove them only after the final hash and
   size manifest is reviewed.

## Addendum: independent re-audit (2026-08-12)

Three independent audit passes were run against commit c1426d9 with a clean
working tree: a code-correctness audit of the ten NISAR scripts, the
`src/peatguard` NISAR modules and the policy tests; a data audit that
recomputed published statistics directly from the on-disk rasters without
reading the generating code; and a claims audit of this document, the
integration plan, the operator console and the pipeline dashboard against the
project's scientific house rules. The full pipeline test suite passed
(106 tests).

### Independently confirmed [E]

- The fine-versus-mission validation statistics reproduce from the rasters to
  about one part in 1e10 (circular MAE 0.379219692 rad, RMSE 0.501095063 rad,
  resultant length 0.886430095, 8,480,500 pixels). The mask combination is
  unique: dropping any stated threshold moves the MAE off the published value.
- All eight source and three derived SHA-256 hashes in the v3 commissioning
  manifest match on disk, byte counts included.
- The three correction variants are distinct data with physically sensible
  ordering (per-pair phase std 3.50 to 1.17 to 1.12 rad). The correction sign
  was confirmed empirically on the real stacks: subtracting the screens
  shrinks the quality-mask phase std (4.05 to 1.11 rad); adding grows it
  (7.18). Corrections are applied exactly once, before referencing.
- Wavelength derives from the product `centerFrequency`
  (1,239,000,000 Hz, so 0.2419632429378531 m) on every live NISAR path, and a
  synthetic-product test confirms metadata overrides the config nominal. Note
  that L2 products expose grid-level `centerFrequency`;
  `processedCenterFrequency` is an RSLC-swath field. The one path where the
  config nominal is authoritative is the (unused) ISCE2 stripmapApp route in
  `timeseries/mintpy_prep.py`.
- The published 5 m coherence layer equals the pixelwise minimum of the two
  pair coherences to half a quantization step; its footprint fully contains
  the operational AOI; values lie in [0, 1]; nodata is declared.
- `phase_to_los_displacement` matches the pinned MintPy implementation; the
  tree inversion re-integrates independently at 1e-6 m tolerance; closure runs
  on the native grid behind hard same-grid assertions; the GSLC interferogram
  path enforces earlier-reference ordering, a correctly normalised coherence
  estimator, and identical coordinate vectors with no resampling.

### Blockers

- **B1. The commissioning gate is advisory [E].**
  `scripts/audit_nisar_commissioning.py` writes every product unconditionally;
  pass or fail exists only as JSON strings that no downstream consumer checks.
  This has already had an effect: `scripts/fuse_nisar_historical.py` consumed
  a displacement raster tagged `do_not_publish` and wrote
  `bayesian_candidate_rate_process50_diagnostic_80m.tif` into
  `grid_products/external/`. `scripts/publish_nisar_fine_coherence.py`
  enforces its gate correctly and is the model to copy.
- **B2. `tests/test_scientific_policy.py` pins strings, not behavior [E].**
  All six assertions are source-text matches over one function. A
  demonstrated reintroduction of both the forced -25 mm/yr shift and a
  universal carbon factor passes every assertion, and a shift placed in
  `velocity.py`, `subsidence_class.py` or `fusion.py` is invisible to them.
  The "Scientific code-policy correction" section above therefore overstates
  the protection: the behaviors were removed [E], but the tests do not
  enforce the policy [E].

### Major defects

- **M1.** The fusion computes `posterior_sigma_mm_per_year`
  (`analysis/historical_nisar.py`) and never writes or reports it; the
  blended rate ships with no uncertainty layer.
- **M2.** The fusion RMSE and sign agreement compare a relative-datum NISAR
  field with an absolute Tay field without common-mode removal. About 68% of
  the reported mean-square error is the arbitrary datum
  (sqrt(11.74^2 + 12.60^2) = 17.22 of the reported 17.34 mm), and the 0.437
  sign agreement is not meaningful under a -11.7 mm mean offset. The fusion
  comparison statistics quoted in "Alternate baseline assessment" above must
  not be read as disagreement magnitudes.
- **M3.** `np.nansum` variance accumulation in
  `timeseries/nisar_commissioning.py` lets missing ionosphere-uncertainty data
  lower a pixel's sigma, easing the 2-sigma release gate. Currently benign:
  measured NaN fraction is zero in both uncertainty rasters.
- **M4.** Per-pair variance and phase summation assume a chain network; the
  acyclicity check does not enforce chain topology. Correct for the current
  two-edge chain; wrong for any future star or branched tree.
- **M5.** `insar/nisar_gunw.py` builds overviews with average resampling for
  `wrapped_phase_20m`, averaging wrapped radians across the +/-pi
  discontinuity (the exact artifact `insar/nisar_gslc.py` documents avoiding).
  Full-resolution data is unaffected; overview-level reads are corrupted.

### Claims and labeling findings

- The operator console's hollow-pin texts are now false as worded: "no direct
  InSAR coverage" and "areas L-band radar could not measure (~11%)" describe
  the Tay et al. (2025) product's coverage, not the L-band sensor class, and
  are contradicted by the two NISAR L-band coherence layers in the same
  panel. Four strings (`select_dam_sites_v2.py` remark, map legend, map
  explainer, pin-evidence substring match) plus seeded task rows must change
  together, and the rewording must not imply NISAR displacement fills the
  gap: NISAR observes these areas but releases coherence only.
- The console layer labels omit PROVISIONAL (the pipeline dashboard labels
  carry it); the maturity badge is gated on the 20 m product alone and its
  text ("20 m coherence only") predates the 5 m layer.
- `nisar_24day_equivalent_rate_diagnostic_80m.tif` and
  `bayesian_candidate_rate_process50_diagnostic_80m.tif` carry `mm/year`
  units and sit in `grid_products/external/`, against the storage boundary in
  the integration plan. The Bayesian blend crosses the project's own
  no-blend rule and is retained only as evidence that blending is
  unjustified (median NISAR weight 0.42%); both belong under
  `_deletable_work/` with that sentence attached.
- `tay_historical_rate_on_nisar_80m.tif` holds Tay pixels but carries
  NISAR-scoped tags (`nisar_span_days`, `nisar_product_maturity`); prefix or
  remove them so tag-trusting readers cannot mistake it for a NISAR product.
- Confidence-tag corrections: policy verdicts tagged [E] in the integration
  plan are decisions, not evidence, and should be untagged or [S]; the
  "mean ionosphere screen of about 3.6 radians" figure has no shipped
  evidence file and must be sourced or dropped; Tay's 0.36 cm/yr field RMSE
  drives `historical_measurement_sigma_mm_per_year: 3.6` and needs a
  page-level pointer.
- Corrections to this document: the executive summary quotes the quality_all
  correction sensitivity (24.17 mm); the mapped-peat value from the same JSON
  is 54.12 mm and belongs alongside it (the omission runs conservative). The
  two 20 m pixel counts quoted in "Current NISAR commissioning stack"
  (11,107,825 and the fine-validation 11,090,185) use different masks; both
  are correct and should say so.
- Packaging: the 5 m coherence GeoTIFF is tiled and compressed but has no
  overviews, so it fails COG validation; rebuild through the COG driver. The
  troposphere screens live as radar-grid cubes under
  `science/LSAR/GUNW/metadata/radarGrid/`, not as geocoded layers; any tropo
  step needs an interpolation stage.

### Improvement plan with accuracy impact

| # | Improvement | Effect on subsidence accuracy | Priority |
|---|---|---|---|
| 1 | Enforced release gate: quarantine directory, promotion on pass; fusion refuses `do_not_publish` inputs; thresholds in a versioned policy file | none (integrity) | fix now |
| 2 | Behavioral policy tests: mean preservation on synthetic velocity; no carbon-unit outputs; wavelength-override test at a frequency other than 1239 MHz | none (integrity) | fix now |
| 3 | NISAR pair-network builder: dedupe reprocessed granules, all pairs within a baseline limit from GSLC, extend `insar/pairs.py`; cycles enable closure-bias correction and functioning inverse-variance weighting | direct gain [S] | next |
| 4 | Automated 12-day ingestion: ASF poll, AOI subset (HTTP range-read pilot), stack append, gate re-run | gain via record growth [E] | next |
| 5 | Descending track (p133, p61): cross-geometry consistency, vertical decomposition | gain [S] | next |
| 6 | Datum: multi-patch common-mode reference with reported spread; widen subset to include CPKY | gain for absolute rates; attacks the 15.15 mm datum span [S] | next |
| 7 | Correction-sign gate inside `prepare_phase_variant` (std must shrink); troposphere-only variant; per-pair correction magnitudes in the manifest | robustness and diagnosis | next |
| 8 | Fusion statistics: demean before RMSE or drop RMSE and sign agreement; write posterior sigma | honest reporting | fix now |
| 9 | M3, M4, M5 fixes: all-pairs-finite masking, BFS-path variance, nearest-neighbour overviews for wrapped phase | correct uncertainty and rendering | fix now |
| 10 | Packaging and copy: COG overviews; hollow-pin strings and seeded rows; PROVISIONAL in console labels; ungated maturity badge; relocate the two mm/year diagnostics | public integrity | fix now |

None of these makes the current three-acquisition stack releasable; the gate
fails for physical reasons (24-day span, ionosphere lower bound near 44 mm
against a -7.8 mm median displacement) [E]. Items 3 to 6 raise the accuracy
of the eventual velocity and shorten the path to the annual gate; the
remainder guarantee that the stated accuracy is real when it ships.
