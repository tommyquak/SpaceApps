# Full Tay-fidelity baseline campaign (NUS SoC Slurm)

93 scenes path 28 / frame 3700, 2018-10-12 to 2024-10-04, bandwidth-20 network
= 1650 pairs, ISCE2 **alosStack** (not pairwise alos2App), multilooked **17 range
x 56 azimuth** per Tay et al. 2025, geocoded to a 3 arcsec grid.

## Current decision (2026-09-10)

The historical campaign below is not a validated subsidence product. Later
work includes bandwidth/bias/seasonality tests, an ALOS-1 historical mosaic,
NISAR local-contrast diagnostics and a failed spatial-transfer ML experiment.
Read the [recovery decision](../../docs/scientific_audit/RECOVERY_DECISION_20260910.md)
for the current evidence, inspected upstream code and bounded next experiment.

[S] Priority: reconcile the already-observed fine GSLC versus mission GUNW
phase discrepancy on a small crop, then test new observed connections and
held-out acquisitions. Do not restart the full ALOS-2 campaign or treat
zero closure in reconstructed phase-linking edges as independent validation.
The September 6 history records loss of the full ALOS-2 SLC stack, so methods
requiring complex SLCs are not drop-in repairs for the surviving filtered pairs.

[E] A new local, retrospective check finds r=0.724 at 31/60 dam sites on shared
1.2 km sampling support, with positive spatial-group bootstrap intervals at
3/5/10 km grouping. That is filtered-pattern consistency, not ground accuracy,
annual velocity or restoration attribution. Its source rasters lack verified
acquisition/correction lineage. The JSON report explicitly blocks release.

### Workflow safeguards and compatibility

Deploy `workflow_provenance.py` together with `mintpy2.sbatch`,
`apply_zheng.sh` and `compare_tay.py` when a future remote trial is authorized.
New records bind stage products to the retained network, reference and file
versions. Large phase cubes are not fully content-hashed. Legacy caches without
records fail closed; they are not retroactively certified. Use a new isolated
work directory, or explicit `REBUILD=1`, which preserves the previous work
directory before rebuilding. No remote deployment is part of this local edit.

`PHASE_CLOSURE_PAIRS` must name an existing valid pair list; unknown pairs and
bad syntax stop the run, while a keep-all list is accepted. The actual MintPy
retained flags must match it. An unrepaired or unverifiable integer-ambiguity
state stops the normal path before inversion. `ALLOW_UNREPAIRED=1` is an
explicit diagnostic override and cannot produce normal workflow success.

`ZHENG_BW` (default 20) specifies the retained inversion bandwidth;
`ZHENG_NL` (default 20) specifies the closure-loop parameter. The retained
network and required source connections must satisfy the pinned estimator
before inversion. A shorter exploratory network may use `APPLY_ZHENG=0`,
which is explicitly not the full baseline workflow. A successful workflow
marker establishes execution/provenance checks only, not scientific validity.

NISAR stage scripts require `analysis/nisar_safety.py` alongside them. The
pull runner now writes directly to a fresh run-specific `NISAR_OUT`; failed
stages cannot publish a provisional completion marker or reuse old global
rasters. Select a single group with `NISAR_A_TRACK`, `NISAR_A_FRAME`,
`NISAR_D_TRACK`, `NISAR_D_FRAME` and `NISAR_CRID` if the cache contains several.
`NISAR_{A,D}_START/END` can select a contiguous acquisition interval. Source
manifests give actual dates; legacy `disp60` names do not guarantee 60 days.
Missing temporal links and insufficient stable-ground fitting support stop the
analysis. Existing scheduler resubmission policy is unchanged.

## Archived state (2026-09-04)

The chain runs end to end and the product is **NOT scientifically valid**. Do not
quote rates from it.

- 1650/1650 pairs, zero failures; temporal coherence median 0.0975, 7,449 px above temporal coherence 0.7 (NOT Tay's sigma-phi criterion, which we have not computed)
- **All 7,449 of those pixels sit on mineral hills in the SE.** At the 25
  candidate dam sites, ZERO pixels exceed 0.5 and the median is 0.096. The
  product is inverted relative to what the deliverable needs.
- **PROVISIONAL** (measured on the PRE-closure-correction product): benchmark vs Tay, Pearson r = **-0.04** (block-bootstrap CI -0.153..0.053), our
  rates ~150x larger in magnitude. This is a **same-sensor benchmark, not
  independent validation** - both are ALOS-2.
- Tay's published raster covers only **6.63%** of our grid (141,697 px, ~1,210
  km2 of 18,243). Count NONZERO cells; their GeoTIFF declares no NoData so
  gdalinfo VALID_PERCENT is misleading.
- ~68% of closure triplets still carry unwrap errors (PROVISIONAL, pre-correction).
  phase_closure runs but produces nothing: 'number of common regions: 0' after
  62 minutes, because maskConnComp is empty by construction (max per-pixel
  coverage 75.3%). The workflow now reports this as INCOMPLETE rather than as a
  clean skip.

### Phase-closure recovery preflight

MintPy main still requires one connected region shared by every retained
interferogram for the Yunjun et al. (2019) phase-closure correction. Lowering
`connCompMinArea` cannot repair an entirely empty intersection. Before another
inversion, run the read-only selector after no job is using the HDF5 stack:

```bash
python ~/full/select_phase_closure_network.py \
  ~/full/mintpy_tay/inputs/ifgramStack.h5 \
  --min-bandwidth 2 --max-bandwidth 20 --min-region-pixels 2500 \
  --write-pairs ~/full/phase_closure_pairs.txt \
  --json ~/full/phase_closure_network.json
```

The selector streams `connectComponent`, finds the largest tested complete
sequential bandwidth with a common region and closure triangles, and samples
integer closure there. It respects `dropIfgram` by default; `--include-dropped`
explicitly opts previously dropped pairs back into the diagnostic. It does not
edit the input stack. Pair exports require `--json`, and every requested output
path must be new. It validates the reference component and finite sampled phase
before emitting a pair list; default sampling is 100 reproducible random pixels.
The JSON records source metadata/file-version fingerprints and a pair-list hash,
not a hash of the full phase cube.

Exit 0 means the diagnostic completed, not that the peat measurement is valid.
Exit 2 means no complete sequential candidate passed the tested structural
criteria. Exit 3 means the reference, phase samples or inputs failed validation.
An optional `--write-region-mask` anchor mask and JSON can still be emitted for an
invalid reference, without a pair list, to support explicit anchor reselection.
An optional `--aoi-bbox W S E N --min-aoi-pixels N` filters anchor regions, but
mineral-ground anchors are allowed: correction transfer to actual target
components must be tested separately. Treat the sampled bad
triplet fraction as a diagnostic only: it does not by itself prove MintPy's
"majority rightness" assumption. If no complete sequential network is selected,
do not fabricate a relaxed common mask. A topology-aware targeted-pair trial is
still possible; otherwise re-unwrapping or a different spatial-temporal
unwrapper is required.

[H] A secondary historical salvage option is CorPhU from Benoit et al. (2020),
doi:10.1093/gji/ggaa120. Unlike MintPy's implementation, it builds a mask for
each interferogram triplet, so it does not require one region shared by the
whole stack. Published GPL-3 source survives at
`github.com/braimbow/CorPhU`, but it is not a drop-in dependency: the 2019 code
expects raw ISCE pair directories plus wrapped `.int` files, uses deprecated
HDF5 APIs, and does not read a MintPy stack. Port it behind an experimental
HDF5 adapter and validate closure before and after; do not place it directly in
the production path.

Historical audit record, see git log f2c46a5, 6444ef1, 2f198f0. At that time open: sigma-phi
metric with effective looks, Tay's troposphere method, ion decorrelated-area
masking, ion ablation, stable reference datum, time-matched comparison,
environment pinning. (Regression tests now EXIST and pass: tests/test_defect_regressions.py, 11 static checks, one per defect class.)

Historical alosStack flow: seed.sbatch -> wave1.sbatch (per-date resample)
-> ion.sbatch (per-pair split-spectrum, 80x448) -> barrier.sbatch (ion_ls, per-date
grids) -> main.sbatch (array over date blocks: form, mosaic, diff, look_coherence,
ion_correct, filt, unwrap_snaphu, geocode) -> mintpy2.sbatch (build_stack2 warp to
the 3 arcsec grid, then smallbaselineApp with var weighting, topographic residual
OFF and solid Earth tides ON per Tay et al. 2025, no weather tropo since they
report PyAPS and GACOS both worsened accuracy, then closure_phase_bias.py per
Zheng et al. 2022 applied PRE-deramp).

LEGACY, not part of this path: controller.sbatch, fanout2.sbatch, gate.sbatch and
rehearsal.sbatch belong to the superseded pairwise alos2App campaign.
rehearsal.sbatch refuses to run without ALLOW_SUPERSEDED_REHEARSAL=1.

Rationale and the pilot failure this fixes: ../PILOT_RESULT.md.

## Ultra-audit v2 (2026-09-01, pre-gate)
25-agent audit confirmed 18 findings; all applied before the gate opened:
fail-closed squeue snapshots with two-tick vanish confirmation; successor-at-top
singleton controller chain; retry-safe mintpy submission (singleton on mintpy2);
atomic scene extraction (.part rename + .complete markers, stale-partial reaping);
window accounting counts in-flight downloads (true 21-scene cap); infra aborts
release QUEUED instead of charging science attempts (INFRA cap 6); geometric
retry also fires on the silent cull-fallback (rc=0 "use offset without culling"),
first-attempt logs preserved; real snaphu conncomp geocoded (nearest interp
globally, labels never blended) and loaded into the stack (pilot pairs fall back
to coherence>0); slantRangeDistance dataset for topographicResidual;
closure_phase_bias CLI corrected (-a, -c local) behind a sequential-network
completeness gate that never silently lowers --nl.

## Ultra-audit round 3 (applied)
velocityClosureCor now fit with the merged smallbaselineApp.cfg (same excludeDate
as velocity.h5, isolating the closure correction); finale clears the scene store
and cancels downloads before CAMPAIGN_DONE; scancel fences released jids;
submit-hold-mark-release pattern closes the unwatched-pack window (hold-release
smoke-tested); pack publishes are whole-file atomic (tmp+rename, --no-requeue;
zombie residual documented in fanout2 header); unique-pair DONE counting; NPAIRS
derived from pairs_full.txt. Operational: xcni rack excluded (job_container
launch failures; LAUNCHFAIL protection held with zero attempts charged).
