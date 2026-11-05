# A route without additional JAXA imagery or CPKY data

## Decision

[S] Continue with **local relative LOS change from retained native NISAR GUNW observations**, not an exact historical Tay reproduction. No additional JAXA scene or CPKY observation is required by the implemented path. A missing externally stable reference prevents identifying common regional motion, but does not prevent measuring a target-minus-background contrast within valid connected observations.

[E] The new implementation has been run on the existing local data. Fixed-support checks pass at 21 of 60 dam sites over 60 days and 40 of 60 over the shorter 24-day interval. These are counts of supported diagnostic observations, not independently validated displacement measurements. Neither result authorizes annual subsidence, carbon-loss estimates or restoration-effect claims.

This document supersedes the earlier recommendation to make access to the authors' exact JAXA manifest or CPKY a prerequisite for all further work. Those remain relevant to exact reproduction, but they are not prerequisites for this alternative observable. The previous [Tay replication audit](../../docs/scientific_audit/TAY_REPLICATION_VERDICT_20260910.md) remains valid about what has not been replicated.

Confidence: [E] documented behavior or measured output; [S] engineering inference; [H] untested physical hypothesis. Code tests and support checks do not establish ground accuracy.

## What the reference-free calculation actually means

For each measured interferogram, compute:

```text
contrast = mean(target phase) - mean(local background phase)
```

[S] If both areas share a component offset `b`, then `(target + b) - (background + b)` cancels `b`. We do not need to know whether that background is motionless. We do need valid connectivity, a fixed spatial definition and checks for remaining spatial errors.

[S] Inverting these observed contrasts gives changes relative to the background, with the first epoch set to zero. This temporal zero is a mathematical convention, not a ground-stability claim. If target and background both subside equally, the contrast is zero. The method cannot recover that shared motion. Spatial atmosphere, scattering changes and internal unwrapping errors may still contaminate the contrast.

Target-reference covariance matters for InSAR uncertainty, as discussed by Fattahi and Amelung (2015), *InSAR bias and uncertainty due to the systematic and stochastic tropospheric delay*. [E] [Author institutional record](https://scholarship.miami.edu/esploro/outputs/journalArticle/InSAR-bias-and-uncertainty-due-to/991031576363602976). The implementation consequently does not turn pixel scatter or a perfect network fit into a measurement-accuracy estimate.

Zhang et al. (2024), *Automated Reference Points Selection for InSAR Time Series Analysis on Segmented Wetlands*, is a useful precedent for region-specific, coherence/connectivity-aware references. Its application was wetland water-depth change, not validation of our peat-subsidence observable. [E] [Paper](https://doi.org/10.1109/LGRS.2024.3390568), [author code](https://github.com/geopaulzhang/auto-ref-point). We implemented the simpler explicitly relative contrast rather than importing its algorithm or assuming its selected references are stationary.

## Implemented changes

The new [native GUNW command](../../scripts/audit_nisar_local_contrasts.py) bypasses the unavailable external reference and reads small windows directly from retained native products. It does not rerun the old ALOS-2 chain or consume the high-pass-filtered September composites.

The [native reader](../../src/peatguard/timeseries/nisar_local_contrast.py) and [contrast/network core](../../src/peatguard/timeseries/local_contrast.py) implement:

- Matching track, frame, release, product version, frequency, native grid, shared-date RSLC identities and processing bandwidth. Ambiguous combinations are rejected, not silently resampled or deduplicated.
- The actual retained uint32 mask encoding, rejecting invalid subswaths, water, anomaly flags, reserved bits and component fill values. Valid zero phase is retained.
- Four-connected spatial relabeling within each original component after validity screening. A product label alone is not accepted as proof of connectivity through valid samples.
- Fixed target and background pixels across every included pair. Selection uses target support and geometric tie-breaking, not phase values, desired signs or agreement with Tay.
- Raw and ionosphere-subtracted variants on identical selected support. The output records where ionosphere estimates were interpolated. Neither variant is called fully atmosphere-corrected.
- Actual network inversion and, when an observation covariance is supplied to the core, propagation through its linear operator. Without covariance, uncertainty is null, not zero. Actual observed cycle closure is separate from fitted residuals.
- Source identities, exact selected pixel indices, coordinate/window hashes, implementation hashes and runtime versions. The source inventory does not claim full-file content hashes for the multi-gigabyte inputs.
- Explicit `absolute_datum=false`, `background_assumed_stationary=false`, `ground_accuracy_validated=false`, `annual_velocity_ready=false` and `release_ready=false`.

[E] NASA's NISAR expert acknowledged a provisional-product case where a component crossed invalid subswath samples and suggested masking/relabeling as a mitigation. That motivates the additional connectivity screen; it is not proof all within-component errors are removed. [Official discussion, Fattahi response of August 27](https://forum.earthdata.nasa.gov/viewtopic.php?t=8137).

[E] Mission GUNW supplies ionosphere and external correction layers without applying them automatically. The new diagnostic subtracts only its named ionosphere layer and leaves troposphere/solid Earth tides unapplied. [ASF GUNW documentation](https://nisar-docs.asf.alaska.edu/gunw/).

## Real-data results

The input is ascending track 125/frame 176, frequency A/HH, CRID P05023, 40 MHz. The three actual observed pairs are June 22 to July 4, July 4 to July 16 and July 16 to August 21, 2026. The last edge spans 36 days; no missing acquisition or zero-motion edge was fabricated. [E]

Parameters were chosen before the native phase results: 500 m target radius, 1-2 km comparison annulus, coherence at least 0.3 in every pair, at least 10 pixels and 50% coverage in each requested region, with valid spatial connectivity. These are engineering screens, not published accuracy thresholds. No threshold was lowered to obtain these counts. [E/S]

| Observation interval and geometry | Supported sites | Observed redundant cycles |
|---|---:|---:|
| June 22 to August 21, 1-2 km background | 21 / 60 | 0 |
| June 22 to July 16, 1-2 km background | 40 / 60 | 0 |
| June 22 to August 21, 1-1.5 km background | 21 / 60 | 0 |

[E] Shorter-window coverage is not a replacement for the failed longer-window sites. Each window has its own fixed support; the two estimates must not be spliced into one time series. Likewise, the narrower comparison geometry supports a slightly different set: 20 sites are shared with the primary 60-day run.

[E] On the primary 21 sites, ionosphere subtraction changes the endpoint contrast by a median absolute 0.35 mm and maximum 1.05 mm. This is sensitivity to the supplied correction, not ionosphere uncertainty. Between 94.6% and 100% of selected target pixels carry the ionosphere-interpolation flag in the individual pairs, so small correction differences must not be read as independently demonstrated sub-millimetre accuracy.

[E] Changing the comparison annulus and local connectivity window gives a median absolute endpoint difference of 1.29 mm, maximum 3.66 mm, across the 20 shared sites; two change sign. This is geometry/reference sensitivity, not a confidence interval or evidence that the better-looking geometry is correct. Both definitions remain recorded.

[S] The present use is to identify where a defined relative observation is supported and where further validation is needed. Do not rank restoration outcomes by these provisional millimetre values. The network is a tree, so exact inversion fit provides no independent error check.

Authoritative measured reports with source-lineage and support provenance:

- [60-day primary](../../docs/scientific_audit/nisar_native_local_contrasts_60day_20260910_v2.json)
- [24-day diagnostic](../../docs/scientific_audit/nisar_native_local_contrasts_24day_20260910_v2.json)
- [60-day geometry sensitivity](../../docs/scientific_audit/nisar_native_local_contrasts_60day_inner_ring_20260910_v2.json)

The earlier non-v2 JSON is retained as the initial local run, before provenance hardening. Its primary count is unchanged; use v2 for reproducibility.

## Fine-resolution route: a real repair test, not another assumption

[E] The GSLC exporter now optionally preserves complex cross-product amplitudes. The validator can explicitly compare amplitude-preserving aggregation against its existing unit-phasor diagnostic, with date/grid/source/sign checks. The default legacy mode is preserved for compatibility.

[E] A fixed 2 km crop at the first existing dam site, D-01, was selected before inspecting phase/coherence. On identical comparison pixels, the two 12-day pairs' phase MAE falls from 0.374 to 0.305 rad and from 0.369 to 0.302 rad, approximately 18%. This shows that discarding amplitudes contributed to this local mismatch. Both results still fail the existing phase thresholds. No new resolution or displacement-accuracy claim is justified.

[S] Keep the custom fine displacement route unreleased. Use the native GUNW diagnostic for current local commissioning. Its independence from the custom fine operator removes that operator mismatch as an implementation dependency, not as a general radar-error problem. [Full reproducible operator experiment](../../docs/scientific_audit/operator_probe_20260910_d01/README.md).

## Next steps that do not depend on JAXA or CPKY

1. Freeze the local geometry and explicit relative observable. Keep failed sites visible. Preserve separate historical Tay/ALOS layers without using their values to tune or calibrate these contrasts.
2. Add full correction sensitivity using the already retained mission/MintPy correction products, only after verifying their native support, phase convention and source lineage. Treat ionosphere-interpolated coverage explicitly, not as measured correction certainty.
3. Obtain a genuinely redundant, independently formed NISAR connection on mission-compatible processing support, or another independently justified check. Another nearest-neighbor observation alone extends the record but does not create a closure loop. Do not use algebraic closure from a fitted phase history as validation.
4. Test on a withheld observation/window before changing product status. Repeated local comparisons can test reproducibility; neither this nor cross-track agreement establishes absolute ground accuracy or restoration causality.
5. If valid phase support remains unavailable at a site, retain the canal/historical-risk evidence and label the current displacement observation unavailable. Do not fill it with ML predictions presented as measurements.

[S] Absolute regional subsidence cannot be recovered by software from an unanchored relative observation without added assumptions or evidence. This alternative deliberately retains the identifiable local difference instead of assuming the unavailable datum is zero. [H] Whether those differences become a reliable peat-motion monitoring product requires the remaining correction and independent-observation checks.

NISAR P05023 is calibrated but provisional, with remaining documented artifacts and incomplete broader validation. That status is preserved in the native source records. [E] [ASF known issues](https://nisar-docs.asf.alaska.edu/provisional-known-issues/).

## Reproduce the primary local check

From the monorepo data root (`$DATA_ROOT`), use a new output path:

```bash
.venv/bin/python3 scripts/audit_nisar_local_contrasts.py \
  --gunw grid_products/_deletable_work/nisar/gunw/*125_A_176*.h5 \
  --sites deliverables/tay_full_grid/dam_sites.geojson \
  --output /private/tmp/peatguard_native_relative_recheck.json
```

The reader rejects incompatible files if future additions make that glob ambiguous. The report records every actual input and selected support. Do not run against files being downloaded or modified.

[E] Verification after the implementation: **321 tests passed, one skipped**, with 11 existing warnings. The skipped integration test requires an absent cluster seed artifact. CLI help and whitespace checks passed. Only local bounded reads and small derived artifacts were used. No original data were overwritten, no new imagery downloaded, and no remote job, deployment, commit or push was performed.
