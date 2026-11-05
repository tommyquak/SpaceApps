# PeatGuard recovery decision

## Decision

[S] Continue NISAR commissioning through a bounded phase-comparison experiment, while retaining historical ALOS-1 and published Tay products as separate, dated evidence layers. Do not restart the full ALOS-2 campaign, release a current annual subsidence map, or expand the machine-learning model on the present evidence. The immediate objective is a reproducible, quality-controlled local displacement contrast, not a complete carbon-accounting product.

The evidence cut-off is 10 September 2026. [E] denotes a documented method, inspected code behavior, or recorded measurement, with its stated scope. [S] denotes an inference or engineering recommendation. [H] denotes a hypothesis that still needs an experiment. An [E] measurement from this repository is not thereby independently validated ground motion.

Three issues have been conflated: whether the workflow completes correctly, whether radar phase contains recoverable displacement information at the target sites, and whether a displacement change can be attributed to restoration. Software repairs can address the first issue. They can expose, but cannot manufacture, the evidence needed for the other two.

There is a promising new diagnostic. At 31 of 60 dam sites with adequate shared coverage, the existing filtered descending and ascending NISAR fields correlate at **r = 0.724**, with 27 of 31 site medians sharing the same sign. This remains positive under three spatial-group bootstrap choices. These are retrospective, filtered cross-geometry comparisons, not a measurement-accuracy result. The analysis explicitly keeps `release_ready=false` and `annual_velocity_ready=false`.

## Current evidence and its limits

### Historical ALOS-2 campaign

[E] The September 4 failure was not simply a missing script. The large retained network had no connected region shared by every interferogram. MintPy's implementation of the integer-ambiguity correction therefore had no common anchor on which to solve. Its Yunjun, Fattahi and Amelung (2019) correction and the Zheng et al. (2022) multilook-bias correction address different errors. A successfully generated Zheng product does not establish that integer ambiguities have been repaired.[^1][^2]

[E] Subsequent repository records include effective-look phase-quality calculations, temporal-baseline coherence tests, wider-bandwidth inversion, tropospheric tests, seasonal fits, and Zheng variants. They do not establish a usable annual rate field at dam sites. The analysis inventory reports mean spatial coherence around 0.22 at 14 days and 0.07 at 400 days. It also records a seasonal fit with approximately 134 mm residual RMS against 32 mm amplitude across the analyzed population. Selecting a quieter subset afterward is exploratory selection, not held-out validation.

[E] These later analyses supersede the old README's claim that the phase-quality calculation had not been performed. The full SLC stack was subsequently lost, according to the September 6 project history; filtered, heavily multilooked interferograms remain. That distinction matters for software choice. A program requiring the original complex SLC observations cannot recover them from a filtered phase cube.

[S] An earlier audit's implication that a mineral-ground anchor is intrinsically useless was too strong. An anchor outside peat can be appropriate. The real questions are whether the selected reference is valid, whether correction transfers through genuinely observed connected components to the target sites, and whether those sites retain useful phase information. Requiring the anchor itself to lie on peat would confuse datum selection with target coverage.

Evidence: [analysis inventory](../../cluster/full/analysis/README.md), [existing commissioning audit](../../docs/scientific_audit/NISAR_SCIENTIFIC_AUDIT.md), and the repository history through `9b1d53f`. No live cluster inventory or historical full-stack inversion was repeated for this decision.

### ALOS-1 baseline and transfer learning

[E] The historical four-frame mosaic is a real additional asset, not a proposed future rescue. Its recorded inter-frame LOS-rate shifts are 0, +13.18, +21.92 and +6.42 mm/year. The reference-frame datum shift is +15.91 mm/year using a location attributed to Zahro et al. (2025). These shifts are material relative to the peat rates under discussion. A precise spatial-bootstrap interval on a regional mean does not include an unmeasured reference-motion or mosaic-offset systematic error.

[S] Preserve this layer as a **2007-2009 relative historical rate**, with datum sensitivity and coverage attached. A reference location validated in another study or era is not automatically stable in 2007-2009 or 2026. Vertical projection also retains the assumption that horizontal motion is negligible. Do not subtract the historical map from a 2026 short-window contrast and call the difference restoration benefit.

[E] Experiment 0 already tested the proposed covariate rescue. Covariates-only R² falls from 0.721 under random splitting to 0.297 at 5 km and -0.018 at 20 km for ALOS-1; leave-one-path-out R² is -0.246. Tay-label covariates-only R² falls from 0.439 to 0.164 and -0.010 under the same random/5 km/20 km progression. The recorded verdict is `STRATUM_MEANS_ONLY`, but even its requested three distinct strata were not obtained. Descriptive strata therefore remain descriptive; the verdict is not proof that strata transfer to unobserved areas.

[S] More model complexity is not the next experiment. These results are consistent with learning spatial structure that does not transfer at the required scale. Spatially structured validation is specifically needed when observations are spatially dependent, as discussed by Roberts et al. (2017).[^3]

Evidence: [mosaic record](../../cluster/full/palsar1/results/mosaic.json) and [Experiment 0 metrics](../../cluster/full/palsar1/results/experiment0_metrics.json).

### Published benchmarks

[E] Tay et al. (2025) establish that a carefully referenced L-band time series can support peat-subsidence analysis. Their reported dataset has 92 acquisitions ending 9 August 2024, whereas this campaign describes 93 ending 4 October 2024. Their methods identify path 28/frame 3650; the campaign describes frame 3700. These are explicit reconciliation items, not proof of a geometry error. The study's bandwidth-4 ionosphere solution, wider inversion network and stable-reference treatment must not be collapsed into a single parameter called bandwidth. Tay remains an external benchmark and must never be used to force our output scale or datum to match.[^4]

[E] Zahro, Izumi et al. (2025) provide an especially relevant counterweight to assuming that L-band alone solves coverage. They report coherent coverage of 76.12% in their earlier ALOS period and 18.3% in their later ALOS-2 period. They explicitly discuss unequal observation periods and reject relaxed coherence thresholds as a reliability solution in dense canal-block areas.[^5] [S] This supports keeping an explicit unobservable category and studying target-specific support, rather than interpreting every blank pixel as a software defect.

## Fresh NISAR site comparison

The new diagnostic reads the two existing September 5 high-pass-filtered rasters and the 60-site GeoJSON. It does not fit one direction to the other, apply another ramp, annualize the values, or use Tay. Ascending pixels are resampled onto the descending grid with nearest-neighbor interpolation, and medians use only the same finite, unmasked pixel support in both fields.

[E] The actual sampling windows are 1,200 by 1,200 metres, not exactly one kilometre. A nominal 500 m half-window rounds outward on the 80 m grid, producing 15 by 15 pixels. A site needs at least half of that window in common valid support, with a minimum of ten valid pixels. Windows extending outside the descending raster are rejected. Valid numeric zeros are retained.

| Quantity | Recorded result | Interpretation |
|---|---:|---|
| Paired sites | 31 / 60 | Coverage-qualified subset, not all proposed dams |
| Pearson correlation | 0.724 | Agreement of filtered patterns |
| Same-sign site medians | 27 / 31 | Consistency diagnostic |
| Median descending minus ascending | +2.25 mm | Cross-field offset, not bias against truth |
| RMS descending-ascending difference | 7.55 mm | Pairwise disagreement, not ground-accuracy RMSE |
| 3 km grouping | 16 groups; r interval 0.470 to 0.860 | Descriptive percentile bootstrap |
| 5 km grouping | 14 groups; r interval 0.455 to 0.868 | Descriptive percentile bootstrap |
| 10 km grouping | 11 groups; r interval 0.458 to 0.880 | Descriptive percentile bootstrap |

The bootstrap resamples whole spatial groups and retains the pairing and original within-group site counts. It uses 2,000 draws and seed 7. The minimum of eight groups is an engineering screening choice, not a theorem establishing independence. Leave-one-group-out correlation remains positive across all tested groupings; the lowest recorded value is approximately 0.634.

[S] This is better evidence of a repeatable filtered spatial pattern than a bootstrap treating neighboring sites as independent replicates. It is still vulnerable to shared spatial artifacts, filter selection on the same observations, differing acquisition windows, component offsets and differences in viewing geometry. In particular, the input rasters do not carry enough verified acquisition/correction metadata to establish that they measure the same time interval. High-pass filtering removes long wavelengths and cannot recover absolute or frame-wide motion. Positive agreement does not make this an annual-rate product.

[E] A separately regenerated baseline table finds ALOS-1 samples at 53 sites, Tay samples at 51, seasonal-amplitude samples at 48, descending NISAR at 60 and ascending NISAR at 37. These independent-layer counts use a different sampler and do not contradict the stricter 31 common-support sites. The original deliverable was not overwritten. The seasonal screen is labeled exploratory, NISAR is labeled a provisional contrast, and the historical datum limitation is explicit.

Artifacts: [full contrast diagnostic with input hashes](../../docs/scientific_audit/nisar_site_contrasts_20260910.json), [regenerated site table](../../docs/scientific_audit/site_baseline_20260910/dam_sites_measured_baseline.csv), and [reproducible checker](../../scripts/validate_nisar_site_contrasts.py).

## Data and software options

### NISAR availability and observability

[E] ASF currently lists calibrated but provisional NISAR products released on 20 July 2026, using CRID P05023 for acquisitions from 17 June onward. Provisional means only limited-site validation is complete. The official guide records a permanent acquisition gap from 27 July at 22:03:25 UTC to 10 August at 00:55:27 UTC. It describes validated reprocessing as a future Q4 activity, not a completed guarantee. Mixing lower-CRID beta data with provisional products needs explicit handling.[^6]

[E] Mission GUNW provides nearest-neighbor interferograms, with unwrapped phase at 80 m posting and a complex wrapped layer at 20 m posting derived from a multilooked radar-coordinate interferogram. Ionosphere and other correction layers are supplied but not automatically removed from the phase. The layer's posting is not its independent resolution.[^7] [S] A nearest-neighbor chain does not itself supply redundant observed triangles; a missing interval must be declared rather than silently treated as zero displacement.

[E] GSLC preserves complex phase and is therefore a candidate source for custom observed connections.[^8] [S] Product availability makes a pilot feasible in principle, but actual same-frame coverage, masks, acquisition continuity and processing versions must be inventoried before choosing dates. This report did not download new scenes or establish that a desired bridge across the instrument gap is coherent.

### Candidate code

| Candidate | Source and inspected revision | Fit to this problem | Decision |
|---|---|---|---|
| MintPy selector plus existing integer correction | Repository selector; Yunjun et al. (2019) | Can identify a smaller complete network with a usable common anchor; does not create missing target support | Keep as bounded historical diagnostic |
| Dolphin | `c2f7c24a055f3de2be13e5fda717baa03cbabad2` | Phase linking from coregistered complex observations; existing covariance-level API does not make incomplete filtered phase pairs equivalent to SLCs | Preferred phase-linking pilot after NISAR input comparison passes |
| SPURT | `4be99dec066e41567b7adfae12fe3d8b03a4a5de` | Spatial-temporal unwrapping; SLC or single-reference stack interfaces, not a drop-in arbitrary 1,650-pair reader | Compare with conventional unwrapping on the same small pilot |
| Sweets | `8e51d7c2973dbc181c3ca43b4e9968f65c3dde42` | Advertises NISAR GSLC input, but inspected development dependencies make reproducibility nontrivial | Implementation reference, not an unpinned production dependency |
| CorPhU | `bb204af5c6374d6259292b5d76a556d0d46ecec5` | Triplet-local correction avoids one whole-stack intersection; legacy ISCE pair-directory and wrapped-data interface needs an adapter | Secondary salvage experiment only |
| MiaplPy | `0d5f4b0cf0de46658c31b850a25c32a5258b773d` | Another complex-stack phase-linking route; not recovery of discarded original observations | No full campaign rebuild solely to change tool |

The source links and access dates are recorded below.[^9][^10][^11][^12][^13] Dolphin, SPURT and Sweets offer BSD-3-Clause or Apache-2.0 licensing in their inspected sources; CorPhU and MiaplPy use GPL-family licensing. Preserve upstream notices and review distribution implications before incorporating code into a distributed product. No candidate was installed or represented as a validated dependency in this patch.

[S] Two implementation traps deserve explicit tests. First, use the NISAR product's center frequency to derive wavelength, not a Sentinel-1 default or rounded L-band constant. Second, SPURT-adapter component labels or interpolated ambiguity fields are not evidence of observed phase support. A controlled trial should disable ambiguity interpolation and report the original quality mask and temporal connectivity separately.[^9][^10]

[E] Phase linking can produce a single phase history whose reconstructed edges close algebraically. Heimpel, Hajnsek and Frey (2026) discuss this rank-one phase representation and temporal-quality limitations.[^14] [S] Near-zero closure in those reconstructed edges cannot validate their physical correctness. Hold out actual measured pair information or later acquisitions, rather than testing the solver only against edges derived from its own solution.

## The next bounded experiment

### 1. Resolve the existing GSLC/GUNW discrepancy

[E] This project already formed three fine GSLC interferograms for 22 June, 4 July and 16 July 2026 on ascending path 125/frame 176. Native closure was approximately 6.4e-8 rad, but that was an algebra/grid check. Against matching GUNW, coherence correlations were 0.913 and 0.902 while phase MAE was 0.379 and 0.399 rad after one permitted circular offset. Both comparisons failed the existing displacement-phase gate. Repeating the claim that fine GSLC is a new, untested route would miss this work.

[E] The comparison uses different aggregation operators. Our validator reconstructs unit phasors from the exported fine phase and averages their sine/cosine components. The mission implementation forms complex cross-products and multilooks those before output.[^17] These estimators differ: complex samples `[9+0j, 1j]` have a unit-phasor mean phase of 0.7854 rad but an amplitude-preserving complex mean phase of 0.1107 rad. [H] How much of the measured 0.38-0.40 rad error this explains is not yet known. Mission GUNW also performs multilooking before geocoding, so matching amplitude weighting alone does not guarantee equivalent processing support.[^7]

[H] Pair-specific registration and flattening differences remain another plausible cause; ASF documents geometry-derived coregistration refined using data-driven offsets.[^15] [E] A separate grid-validation defect was confirmed and patched: relative tolerance in coordinate comparison could admit a five-metre shift at large UTM coordinates. The check now uses only an absolute tolerance, with shifted-grid regression tests. This code defect does not prove that the three staged acquisitions were actually shifted.

[S] Use one small fixed crop containing stable mineral controls and a peat target, with the same two acquisition pairs already staged. Record exact source IDs, masks, frequency, polarization, sign, grid and correction flags. Compare the current unit-phasor aggregation, amplitude-preserving complex aggregation and, if available, matching mission radar-coordinate processing support. Keep the footprint and quality mask fixed and permit only one constant circular offset. Do not fit a spatial surface to make the phase comparison pass.

Retain the existing engineering comparison gates: phase MAE at most 0.10 rad, RMSE at most 0.15 rad, and resultant length at least 0.98. These are project acceptance thresholds, not universal accuracy bounds from a paper. If known processing operators are not expected to produce equivalent phase, document that incompatibility and redesign an independently justified comparison before commissioning the custom product; do not relax the thresholds in response to the current residuals.

If an operator change explains the error, freeze it and test a second crop and a pair not used to select the change. If the residual remains structured, inspect pair offsets, flattening and correction layers. If those cannot be reconciled on a small crop, use mission GUNW for commissioning and leave the fine displacement path unreleased. The full-frame rerun is conditional on this result.

### 2. Test additional information, not just another solver

[S] Only after phase-input compatibility passes, choose a compact sequence with actual redundant observed connections and a stable anchor that links to the peat component. Include the instrument-gap interval only if real observations and usable phase connections exist; otherwise split the time series into separately referenced connected periods. Missing intervals are not zero-motion intervals.

Compare the same masked data and same network through conventional unwrapping and a pinned Dolphin/SPURT pilot. Keep component support, per-pair residuals, source-product versions and unwrapped ambiguity decisions in the output. Reserve at least one measured connection or later acquisition from parameter selection. The pass question is whether added observations resolve ambiguity and predict held-out phase better, not whether the reconstructed solution's triangles close.

[H] The useful product may be a stable local displacement contrast at a subset of sites, not a wall-to-wall rate. Count failed sites as unavailable and report the missingness. Do not obtain apparent coverage by filling disconnected components or lowering quality thresholds until a map looks complete.

### 3. Test the 31-site pattern prospectively

[S] Freeze the current site definitions and proposed filter choices before inspecting additional observations. Record shared acquisition windows where possible, use independently selected stable controls for nuisance fitting, and reserve different controls for checking residual motion. Evaluate no high-pass, 5 km, 10 km and 20 km high-pass variants as a prespecified sensitivity analysis, not a contest won by the highest same-data correlation.

Report spatial-group confidence intervals, leave-group-out sensitivity, component linkage and mask coverage. Compare consistent observables across viewing directions with explicit LOS/vertical assumptions. If the same local pattern survives new observations and independent controls, advance it to an experimental displacement product. An annualized short-window contrast is still not a measured annual trend; trend, seasonality and systematic-error behavior require longer observations.

### 4. Keep restoration attribution separate

[S] Maintain separate columns for historical ALOS-1 rate, published Tay rate, exploratory seasonal diagnostics and provisional NISAR contrasts. Do not fuse their magnitudes into a calibrated current rate. Operational dam priorities can retain supported canal and historical-risk evidence without claiming that the new radar layer has measured intervention impact.

[S] For intervention evaluation, plan repeated pre/post measurements with multiple comparison sites and known intervention dates. Underwood's (1994) Beyond-BACI design addresses why temporal change at one treatment location is insufficient to isolate an intervention effect.[^16] Site displacement, water level, canal condition and carbon consequences are distinct observations. Field reference verification and water-level monitoring would add evidence; neither can be inferred merely because two filtered satellite fields correlate.

## Code and operational boundaries

The associated edits focus on the active scientific workflow and evidence extraction, not an exhaustive certification of the web application or every historical scratch script. The selector now respects the current retained network by default, validates reference and sample support, requires a provenance report for pair export, and refuses to overwrite output files. A diagnostic anchor mask may still be emitted when the current reference is invalid, allowing explicit reference reselection without pretending the pair list is ready.

The MintPy changes add current-network and product provenance checks and stop the normal path before expensive inversion when integer ambiguity status is unresolved. An explicit exploratory override is not a science-success result. The Zheng path distinguishes its inversion bandwidth from its closure-loop parameter and rejects incompatible downstream networks before doing costly work. Comparison should consume a verified current product rather than whichever familiar filename happens to exist.

The NISAR edits isolate each run's outputs, preserve stage failures, and prevent old rasters from being copied into a newly dated success directory. Stable-ground fitting must fail when stable support is insufficient; it must not silently fit away the target peat field. Multiple frames or CRIDs must be explicitly separated or selected before accumulating a time sequence. Duplicate reprocessed pairs are rejected; selection among their individual product versions remains a manual inventory decision.

These changes do not repair existing science rasters retrospectively, supply missing SLC observations, or resolve the failed GSLC phase comparison. No cluster job, large download, deployment, commit or push is part of the recorded local experiment. Legacy products remain available for diagnosis but cannot acquire new provenance merely from a new success marker. The next remote pilot must use a fresh isolated work directory and include the new helper modules with the edited scripts.

[E] The complete configured local Python test suite passes: **237 passed, one skipped**. The skipped seed-metadata integration test needs a cluster smoke-test artifact tree that is not present locally. Eleven warnings remain in unrelated existing tests, including deprecation and synthetic-data warnings. Shell syntax and whitespace checks pass. These checks establish local code behavior, not successful execution of the remote processing environment or validation of the radar products.

## Reproduction

From the repository root, the local contrast diagnostic can be repeated using a new output filename:

```bash
.venv/bin/python3 scripts/validate_nisar_site_contrasts.py \
  --descending grid_products/external/nisar/20260905/nisar_D_disp60_hp10km_mm.tif \
  --ascending grid_products/external/nisar/20260905/nisar_A_disp60_hp10km_mm.tif \
  --sites deliverables/tay_full_grid/dam_sites.geojson \
  --block-km 3 5 10 --draws 2000 \
  --output /private/tmp/peatguard_nisar_contrasts_recheck.json
```

The report stores source-raster and site-file SHA-256 values, actual spatial support, excluded-site status, bootstrap group sizes and explicit scientific limitations. The code does not interpret a successful execution as release authorization. The historical baseline table can also be regenerated into a new directory with `scripts/site_baseline_table.py --output-dir ...`; its independent-layer samples are not substitutes for paired common-support validation.

## Sources

Primary publications and official software/product documentation were checked for this decision on 10 September 2026. Software revisions are audit anchors, not claims of production compatibility.

[^1]: Yunjun, Z., Fattahi, H. and Amelung, F. (2019). [Small baseline InSAR time series analysis: Unwrapping error correction and noise reduction](https://doi.org/10.1016/j.cageo.2019.104331). Computers & Geosciences 133, 104331. Implementation behavior additionally checked against the locally staged MintPy source.
[^2]: Zheng, Y., Fattahi, H., Agram, P., Simons, M. and Rosen, P. (2022). [On Closure Phase and Systematic Bias in Multilooked SAR Interferometry](https://authors.library.caltech.edu/records/ksehb-esr90). IEEE Transactions on Geoscience and Remote Sensing 60, 5226611. DOI 10.1109/TGRS.2022.3167648.
[^3]: Roberts, D. R. et al. (2017). [Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure](https://florianhartig.github.io/publications/roberts2017cross/). Ecography 40, 913-929. DOI 10.1111/ecog.02881. Coauthor publication record.
[^4]: Tay et al. (2025). [Satellite radar advances carbon emissions accountability over tropical peat](https://www.nature.com/articles/s43247-025-02926-6). Communications Earth & Environment 6, 971. DOI 10.1038/s43247-025-02926-6. Published 26 November 2025.
[^5]: Zahro, Q., Izumi, Y. et al. (2025). [Tropical peatland restoration evaluation using time-series L-band InSAR analysis](https://www.nature.com/articles/s41598-025-08390-8). Scientific Reports 15, 23749. Published 3 July 2025.
[^6]: Alaska Satellite Facility. [Available NISAR Data](https://nisar-docs.asf.alaska.edu/availability-overview/). Live release status, provisional CRID, permanent gap and planned validated reprocessing. See also [Provisional known issues](https://nisar-docs.asf.alaska.edu/provisional-known-issues/) before a new inventory.
[^7]: Alaska Satellite Facility. [Geocoded Unwrapped Interferogram](https://nisar-docs.asf.alaska.edu/gunw/). Product formation, posting and unapplied correction layers.
[^8]: Alaska Satellite Facility. [Geocoded Single Look Complex](https://nisar-docs.asf.alaska.edu/gslc/). Phase-preserving complex input product.
[^9]: ISCE framework. [Dolphin README at the inspected revision](https://github.com/isce-framework/dolphin/blob/c2f7c24a055f3de2be13e5fda717baa03cbabad2/README.md); [phase-linking core](https://github.com/isce-framework/dolphin/blob/c2f7c24a055f3de2be13e5fda717baa03cbabad2/src/dolphin/phase_link/_core.py); [displacement configuration](https://github.com/isce-framework/dolphin/blob/c2f7c24a055f3de2be13e5fda717baa03cbabad2/src/dolphin/workflows/config/_displacement.py).
[^10]: ISCE framework. [SPURT SLC-stack interface](https://github.com/isce-framework/spurt/blob/4be99dec066e41567b7adfae12fe3d8b03a4a5de/src/spurt/io/_slc_stack.py); [Dolphin spatial-temporal unwrap adapter](https://github.com/isce-framework/dolphin/blob/c2f7c24a055f3de2be13e5fda717baa03cbabad2/src/dolphin/unwrap/_unwrap_3d.py).
[^11]: ISCE framework. [Sweets README](https://github.com/isce-framework/sweets/blob/8e51d7c2973dbc181c3ca43b4e9968f65c3dde42/README.md) and [dependency declarations](https://github.com/isce-framework/sweets/blob/8e51d7c2973dbc181c3ca43b4e9968f65c3dde42/pyproject.toml).
[^12]: Benoit et al. (2020). [CorPhU paper](https://doi.org/10.1093/gji/ggaa120). Geophysical Journal International 221, 1959; [published source README](https://github.com/braimbow/CorPhU/blob/bb204af5c6374d6259292b5d76a556d0d46ecec5/README.md) and [triplet implementation](https://github.com/braimbow/CorPhU/blob/bb204af5c6374d6259292b5d76a556d0d46ecec5/triplet.py).
[^13]: InSAR Laboratory. [MiaplPy source at inspected revision](https://github.com/insarlab/MiaplPy/tree/0d5f4b0cf0de46658c31b850a25c32a5258b773d). Complex-stack workflow and licensing.
[^14]: Heimpel, M., Hajnsek, I. and Frey, O. (2026). [Heuristic Quality Coefficients for Interferometric Phase Linking](https://arxiv.org/html/2510.24512v2). ISPRS Journal of Photogrammetry and Remote Sensing 237, 1-21. DOI 10.1016/j.isprsjprs.2026.04.015. Revised manuscript 21 April 2026.
[^15]: Alaska Satellite Facility. [Range-Doppler Wrapped Interferogram](https://nisar-docs.asf.alaska.edu/rifg/). Geometry and data-driven coregistration, interferogram formation.
[^16]: Underwood, A. J. (1994). [On Beyond BACI: Sampling Designs that Might Reliably Detect Environmental Disturbances](https://esajournals.onlinelibrary.wiley.com/doi/10.2307/1942110). Ecological Applications 4, 3-15.
[^17]: ISCE3. [Complex cross-multiplication and multilooking implementation](https://github.com/isce-framework/isce3/blob/develop/cxx/isce3/signal/Crossmul.cpp), inspected 10 September 2026. Local comparison implementation: [nisar_validation.py](../../src/peatguard/insar/nisar_validation.py:112). The upstream link is a moving development branch, not a pinned production dependency.
