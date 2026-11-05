# Tay replication verdict

## Conclusion

As checked on 10 September 2026, **PeatGuard has not independently replicated Tay et al. (2025), and the currently verified local inputs are insufficient to run a faithful reproduction.** [E]

It is **not established that replication is physically impossible**. The available record contains important input and method mismatches, so failure of this pipeline is not a controlled test of the published method. A properly matched ALOS-2 reproduction remains a conditional research possibility, not a promised outcome. [S]

| Meaning of replication | Current verdict |
|---|---|
| Load and use the published historical map as a benchmark | Already possible; this is reuse, not independent reproduction |
| Independently reproduce the historical ALOS-2 result | Not achieved; input coverage and processing equivalence remain unresolved |
| Reproduce the paper's field-validated accuracy | Not established; the complete processing and independent field comparison are missing |
| Transfer the method to NISAR in 2026 | A separate, unvalidated experiment; it cannot reproduce observations from 2018-2024 |

This audit is read-only with respect to pipeline code and science products. It includes a fresh local raster comparison and live public catalog checks, but no processing job or scene download. Remote storage and backups have not been independently inventoried.

## A concrete data mismatch

Tay reports 92 acquisitions on path 28/frame 3650, ending 9 August 2024, and uses CPKY as its stable reference. Its reference-grid acquisition is 19 June 2020. [E] [Tay et al., Methods](https://doi.org/10.1038/S43247-025-02926-6).

Our [campaign manifest](../../cluster/full/scenes_full.txt) has 93 frame-3700 acquisitions: 91 through the published endpoint plus 20 September and 4 October 2024. The date 23 April 2021 is absent. More dates overall therefore do not mean the published acquisition set was reproduced. [E]

Live NASA CMR queries for collection `C3315903479-ASF`, path 28 and the published interval returned **zero frame-3650 granules and 91 frame-3700 granules**. A separate query for the missing date, including a one-day margin over the AOI, returned none; the reference-date positive control returned a known granule. This establishes the result of these catalog queries, not absence from JAXA, author-held data or every other archive. [E]

The official record for `ALOS2327953700-200619-WBDR1.1__D`, the frame-3700 scene on the reference date, has polygon vertices below in longitude/latitude order:

```text
112.462, -5.393
115.585, -6.068
116.259, -2.878
113.149, -2.217
112.462, -5.393
```

CPKY at 113.926, -2.220 lies outside this polygon. At that longitude the polygon's northern edge is approximately -2.382. Independently, our saved grid extends only from -3.400 to -2.340 latitude and therefore excludes CPKY even without trusting the catalog polygon. [E]

The frame discrepancy can no longer be treated merely as an unimportant naming difference. Catalog geometry is not a substitute for CEOS leader coordinates and actual valid-pixel geolocation, but this is a concrete spatial-coverage conflict that must be reconciled before a repeat run. [S]

Sources: [official reference-scene metadata](https://cmr.earthdata.nasa.gov/search/concepts/G3484105879-ASF.umm_json), [frame-3650 query](https://cmr.earthdata.nasa.gov/search/granules.json?collection_concept_id=C3315903479-ASF&temporal=2018-10-12T00%3A00%3A00Z%2C2024-08-09T23%3A59%3A59Z&attribute%5B%5D=int%2CPATH_NUMBER%2C28&attribute%5B%5D=int%2CFRAME_NUMBER%2C3650&page_size=2000), [frame-3700 query](https://cmr.earthdata.nasa.gov/search/granules.json?collection_concept_id=C3315903479-ASF&temporal=2018-10-12T00%3A00%3A00Z%2C2024-08-09T23%3A59%3A59Z&attribute%5B%5D=int%2CPATH_NUMBER%2C28&attribute%5B%5D=int%2CFRAME_NUMBER%2C3700&page_size=2000), and [hardcoded processing bounds](../../cluster/full/build_stack2.py:20).

## Processing differences that invalidate a conclusive failure claim

The paper propagates interferogram phase uncertainty through the inversion to epoch uncertainties before averaging across dates. It describes time-series geocoding after inversion, pixel-wise vertical projection, and nearby stable-pixel atmospheric correction. [E] [Tay et al., Methods and Phase variance and standard deviation](https://doi.org/10.1038/S43247-025-02926-6).

Our inspected code differs materially:

- **Wrong uncertainty quantity:** [sigphi_tay.py:43](../../cluster/full/analysis/sigphi_tay.py:43) averages per-interferogram standard deviations directly. It does not propagate covariance through the actual inversion, account for retained-pair selection, or distinguish absent zero coherence from observed low coherence. The reported nominal/effective-look pass fractions are therefore not a faithful implementation of the paper's screen. [E]
- **Unestablished effective looks:** [sigma_phi.py:18](../../cluster/full/analysis/sigma_phi.py:18) infers effective looks from closure scatter in one central crop, assuming independent equal-noise triplet legs. Pair-dependent processing and scattering can affect that scatter. The later frame-wide value of 41 is model-dependent, not an independently measured information limit. [E/S]
- **Different reference:** the coverage/coherence-based [reference selector](../../cluster/full/analysis/pick_ref.py:25) is not a demonstrated stable GNSS datum. Coverage is necessary for processing but does not establish ground stability. [E/S]
- **Different processing order:** the current builder average-warps already geocoded unwrapped interferograms, with separately nearest-resampled component labels. Mixing component offsets during averaging is a risk requiring a controlled test; it has not been proved to explain all observed errors. [E/H] See [build_stack2.py](../../cluster/full/build_stack2.py:114).
- **Additional trend removal:** the workflow specifies a linear deramp. Its effect on real spatial deformation is unquantified for this reproduction. Removing a fitted plane can remove geophysical structure as well as noise. [E/S] See [mintpy2.sbatch](../../cluster/full/mintpy2.sbatch).
- **Atmospheric control selection:** [tay_tropo.py:33](../../cluster/full/analysis/tay_tropo.py:33) may correct a closure-adjusted time series while selecting near-zero-motion controls from uncorrected velocity under our arbitrary datum. Its poor control availability does not prove physically stable neighbors are absent. [E/S]
- **Coarse salvage is another estimator:** [coarse_ml.py:5](../../cluster/full/coarse_ml.py:5) averages already filtered 17x56 interferograms another 10x8 and uses filtered-phase consistency as a coherence proxy. Regridding this to 90 m posting does not create 90 m independent resolution or restore original SLC coherence. [E/S]

Consequently, neither the old low temporal coherence nor the approximate effective-look threshold proves that appropriate original ALOS-2 observations cannot reproduce the published result. Conversely, the newer higher temporal coherence does not prove successful recovery. [S]

The previous recovery discussion also compared approximately 32 mm seasonal amplitude with 134 mm residual RMS. The [seasonal exporter](../../cluster/full/analysis/seasonal_amp.py:48) projects amplitude vertically but leaves residual RMS in LOS. Those figures are not directly commensurate, and residuals include unmodelled signals as well as noise. They should not be used as a conclusive physical signal-to-noise test. [E/S]

## Fresh comparison of the latest local maps

The three latest saved velocity/coherence TIFFs were read directly, rather than substituting statistics from the earlier 1,650-pair run. Their whole-grid temporal-coherence median is 0.540, not the earlier 0.0975. No source manifest ties these TIFFs to an exact input network, reference or correction sequence. The directory also contains older seasonal outputs. [E]

For this check, both velocity maps were area-averaged onto the original published raster's native grid. The source mask requires finite temporal coherence in (0,1]; legitimate zero velocities in our maps are retained. Benchmark nonfinite cells are excluded, and its zero-valued cells are treated as background under the existing project convention. The benchmark has no embedded nodata declaration, so that background interpretation remains an assumption. At least 99.9% of the contributing source area must be valid, and both candidate maps must share the comparison pixels. No reference offset, sign change, amplitude fit, detrend or outcome-selected outlier clipping was applied. [E/S]

| Average source temporal-coherence threshold | Shared pixels | Original saved velocity r | Saved closure-corrected velocity r |
|---|---:|---:|---:|
| No additional threshold | 32,125 | 0.124 | 0.104 |
| 0.4 | 30,417 | 0.161 | 0.132 |
| 0.5 | 28,493 | 0.169 | 0.137 |
| 0.7 | 6,729 | 0.056 | 0.003 |

For the unthresholded comparison, 2,000 whole-group bootstrap draws with seed 7 give the following descriptive 95% intervals. The grouping uses fixed approximate local metric coordinates derived from longitude/latitude, with longitude scaled at -2.8 degrees, rather than outcome-based grouping. [E]

| Group size | Groups | Original r interval | Corrected r interval |
|---|---:|---:|---:|
| 3 km | 426 | 0.068 to 0.184 | 0.054 to 0.159 |
| 5 km | 197 | 0.046 to 0.199 | 0.038 to 0.167 |
| 10 km | 71 | 0.018 to 0.224 | 0.012 to 0.185 |

These maps show weak spatial correspondence with the published field. The saved correction does not improve the point estimate of correlation. This is a product diagnostic, not a fully matched accuracy test: our TIFFs lack verified units, observation intervals and incidence-angle geometry, and are likely LOS while the benchmark is vertical. [E/S]

Constant unit conversion, constant LOS scaling or a constant datum offset cannot change Pearson correlation. Those issues therefore cannot simply explain away the weak r. Spatially varying projection, processing support, data coverage and differing time intervals still prevent attribution of the failure to a single physical cause. [S]

Input hashes:

```text
velocity.tif
d9664f87c3b0e902dfd46160243e16d27066414f4542c4a958be57d440fdd75c
velocityClosureCor.tif
7b900f53331168887ff952a1cabc01b2b63b0718d0bf1dcba7d078ca0d4b1ed4
temporalCoherence.tif
36a88ce7921aad448074e046dd432f7897c1cbe2986d7e38b65b03c8bf0c4578
InSAR_LongTermRate_20181012-20240809.tif
223b8616781ad23ad92cb15ce78fcb36953ce2a2ba9e7300bc923cf3041aacb8
```

The first three are in the ALOS-2 block output directory (`$DATA_ROOT/grid_products/external/alos2_block_2018_2024`, local data outside this repository). The fourth is the original Tay raster (`$DATA_ROOT/grid_products/_deletable_work/baselines/tay_2025/InSAR_LongTermRate_20181012-20240809.tif`, local data outside this repository). The benchmark README labels the slope product's unit as `m`; interpreting it as m/year remains a documented dimensional interpretation rather than a complete embedded unit declaration.

## What prevents a decisive reproduction run now

The local inventory contains the final maps and a separate old 20-epoch pilot, but not the full 93-date stack, the bandwidth-4/60 time series, native ALOS-2 SLCs, or the surviving full-campaign filtered complex pairs. The September 6 stack-loss statement is contemporaneously recorded in commit `aa2daa8`, not merely in the recent report. It does not prove that no remote backup exists. [E]

The public data repository supplies derived radar outputs. The paper's code statement identifies upstream packages rather than a complete run manifest; its field camera data have separate access requirements. Exact reproduction and independent accuracy verification therefore have additional evidence requirements beyond obtaining a final raster. [E] [Tay et al., Data and Code availability](https://doi.org/10.1038/S43247-025-02926-6), [DR-NTU dataset](https://doi.org/10.21979/N9/JGQE4O).

## Smallest justified next step

For the specific objective of replicating Tay, **resolve the source scene identity before trying another correction algorithm**. [S]

1. Obtain the authors' exact granule manifest and reconcile frame 3650 against our frame-3700 CEOS leader, footprint and valid pixels. Ask the relevant archive about the missing frame/date. Do not interpret a negative ASF query as proof of unavailability elsewhere.
2. Confirm a matched input covering a validation target and CPKY or a separately justified connected stable reference. Inventory remote raw/intermediate data and backups before planning reacquisition.
3. Reproduce a small fixed set of interferograms with documented native support, reference and ionospheric handling. Compare stage outputs if author intermediates are available. A few-date experiment can reject processing mistakes, but cannot establish a six-year rate.
4. Only after input/processing equivalence passes, run the required historical time span on a bounded spatial target, preserving the support needed by long-wavelength corrections. Implement actual epoch covariance propagation and a common-support vertical comparison. Avoid tuning parameters to maximize agreement with Tay.
5. Reserve independent evidence for accuracy verification. Reproducing an author output and reproducing its ground-validation claim are separate tests.

No author or archive messages have been sent. No code has been changed to implement these steps. The correct present claim is **“replication not achieved; matched-data feasibility unresolved”**, not **“Tay has been replicated”** or **“Tay is impossible to replicate.”**
