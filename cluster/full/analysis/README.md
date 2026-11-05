# Analysis scripts behind the reported numbers

Every headline claim in the parent README traces to a script here. They were
scratch-only until 2026-09-05, which meant the results could not be reproduced
from the repository. All are read-only against the MintPy work directory.

## Interpretation update (2026-09-10)

The records below are historical diagnostics, not a release certificate.
Stage 2c/2d agreement is retrospective filtered-pattern consistency, not a
validated subsidence signal or intervention effect. The new common-support
site check uses spatial-group resampling: 31/60 sites, r=0.724, with positive
descriptive intervals at 3/5/10 km grouping. It still has no verified original
acquisition/correction lineage and cannot establish annual rates.
See the [recovery decision](../../../docs/scientific_audit/RECOVERY_DECISION_20260910.md)
and its reproducible diagnostic. The current stage scripts also require
`nisar_safety.py` alongside them, honor `NISAR_OUT`, and fail on ambiguous
product groups, broken chains or insufficient stable fitting support.

## No-external-reference route (2026-09-10)

[S] The unavailable JAXA/CPKY inputs are no longer prerequisites for all work.
The new native-GUNW local-contrast path measures target-minus-background LOS
change without assuming the background is stationary. It checks product masks,
connected valid islands, fixed spatial support and source lineage. Local tests
find support at 21/60 sites over 60 days and 40/60 over 24 days. Neither is
validated ground motion; both observed networks have zero redundant cycles.
See the [implemented alternative and measured checks](../../../docs/scientific_audit/NO_EXTERNAL_REFERENCE_ROUTE_20260910.md).

[E/S] Interpretation corrections to historical rows below: `sigphi_tay.py`
averages interferogram uncertainty rather than propagating it to epochs, so
its pass fractions do not implement Tay's actual screen. `tay_tropo.py` has
reference/control-selection mismatches, so its null is not a faithful method
test. `seasonal_amp.py` compares vertical amplitude with LOS residual scatter;
the cited 134/32 ratio is not a commensurate physical signal-to-noise bound.

| script | produces |
|---|---|
| `dam_area_quality.py` | temporal coherence at the 25 candidate dam sites vs the SE hills vs the whole AOI. Source of "ZERO px above 0.5 at the dam sites, median 0.096". |
| `refpoint.py` | per-pixel connectComponent coverage; showed MintPy's chosen reference was valid in only 5.0% of 1650 interferograms. |
| `pick_ref.py` | selects a reference on coverage with an edge margin and 5x5 uniformity. Chose (1013,1601) at 73.9%. |
| `sigma_phi.py` | effective independent looks from closure statistics (41 vs 952 nominal) and Tay's Eq 4 sigma-phi. |
| `closure_strat.py` | closure quality stratified by temporal baseline (7% bad under 30 d, 75% beyond 400 d). |
| `ion_mag.py` | per-date ionosphere magnitude and gradient vs snaphu's pi/px limit. |
| `ion_chain.py` | per-PAIR ion estimates before ion_ls, showing 247 of 362 are physically sane. |
| `remaining_bl.py` | temporal-baseline distribution of completed vs remaining pairs. |

**Caveat carried by all of them:** they read the MintPy work directory directly.
Several refuse to run while a job is active, because reading a progressively
written HDF5 produced two wrong conclusions in this project. Keep that guard.

## Method-audit measurements (2026-09-05)

| script | job | what it measured |
|---|---|---|
| `coh_tbase_tay.py` | 821087 | mean spatial coherence by temporal baseline, Tay strip vs whole frame (0.22 at 14 d, 0.07 at 400 d) |
| `sigphi_tay.py` | 821118 | Tay's sigma-phi <= 0.7 cm rule on our stack at nominal (76.5%) and effective (0.06%) looks; cross-tab with tempCoh |
| `zheng_verdict.py` | login | is the Zheng correction an offset or a field (12.5% offset; std 15.9 mm/yr) |
| `tay_prepost.py` | login | block-bootstrapped r vs Tay before and after Zheng on the same pixels |
| `overlap.py` | login | overlap of tempCoh-gated pixels with Tay's data (0 px at 0.7); centroids |
| `neff.py` | 821127 | speckle ENL cross-check; found no single-date amplitude product, abandoned |
| `summarize_bw.py` | 821123 | summary of the bandwidth-60 inversion against the strip and Tay's mask |

All read `~/full/mintpy_tay` and `~/full/tay/InSAR_LongTermRate.tif`; the Tay raster declares no NoData so exact zeros are treated as padding.
| `tay_tropo.py` | 821146 | Tay's tropospheric neighbour-stack, faithfully implemented; null on this frame (317 qualifying neighbour px) |
| `nisar_stage1.py` | 821165 | NISAR provisional GUNW inventory and 12-day coherence over the frame and the Tay strip |
| `nisar_stage2.py` | 821251 | NISAR coherence at the 60 dam sites, per-pixel NISAR/ALOS-2 ratio; its hardcoded reference fell outside the descending frame (defect 24) |
| `nisar_stage2b.py` | 821252 | 60-day composites on both NISAR tracks with a data-chosen reference, and the two-track agreement |
| `nisar_stage2c.py` | 821264 | per-pair quadratic deramp on stable ground, ion-uncertainty mask, 10 km high-pass; retrospective filtered two-track agreement r = +0.61, not validated displacement |
| `nisar_stage2d.py` | 2d | site-level D vs A agreement and site-versus-control contrast with bootstrap intervals |
| `era_coh.py` | 821312 | 2016-17 pair coherence vs every 2019-21 pair at the same baseline, strip / dam sites / frame |
| `elnino.py` | 821338 | 2019 El Nino response of the bw60 (and bw20) series vs Tay's InSAR_MaxSubsidence-2019ElNino.tif, block bootstrap; r = +0.22 to +0.28 |
| `elnino_box.py` | 821340 | box estimator (window mean minus pre-mean) for the El Nino comparison, per-site table, OLS slope; r = +0.28 to +0.36 |
| `season_coh.py` | 821353 | per-pair strip/dam/frame coherence with calendar months; month matrix and dry-to-dry annual inventory (annual pairs at the floor) |
| `seasonal_amp.py` | 821387 | annual sinusoid plus trend per pixel on the bw60 series; amplitude is noise-dominated (residual 134 mm RMS vs 32 mm amplitude), component 3 episodic-only on ALOS-2 |
| `seasamp_strat.py` | login | amplitude and phase by fit-residual quartile: the clean half of the footprint carries a 35 mm annual cycle at 2.5x its noise expectation, January minimum |
| `zheng_apply_bw60.py` | 821409 | applies a Zheng bias series (BIAS) to a MintPy series (BW) and tests the rate vs Tay; the bandwidth-4 bias moves the mean, not the scatter |
| `bw4_regional.py` | bw4reg | block-bootstrapped regional mean rates from the bandwidth-4 series after the full Zheng bias: strip vs Tay, dam zone, stable ground as the null |
