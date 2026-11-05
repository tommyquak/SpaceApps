# Pilot campaign: documented deviations from the Tay method

Recorded 2026-09-01 per the scientific-conformity audit. Each deviation is a
deliberate pilot-scope decision, not an oversight, and each names its exit.

1. Pairwise alos2App instead of alosStack single-reference coregistration.
   Reason: the 500 GB home quota cannot hold alosStack's stack workspace;
   pairwise runs on node-local NVMe with only small products in home.
   Effect: per-pair independent geometry, adequate for multilooked SBAS.
   Exit: full alosStack run on provisioned scratch for the final baseline.

2. Network bandwidth 4 (70 pairs) instead of Tay's bandwidth 20.
   Reason: pilot scope; bandwidth 20 over 20 scenes is 190 pairs.
   Effect: closure-bias correction will not converge at pilot scale; the
   pilot reports displacements and network diagnostics, not corrected rates.
   Exit: bandwidth escalation on the full 90+ scene archive.

3. Ionosphere correction at ISCE2 alos2App defaults (doIon on, 80 looks).
   Reason: defaults are the documented alos2App configuration; no evidence
   yet to prefer alternatives on this stack.
   Exit: compare against the mission-style bandwidth-4 estimation when the
   full baseline runs.

4. No water body mask supplied to alos2App.
   Reason: sardem NASA_WATER source failed (retired endpoint) and the mask
   is optional for interferogram formation.
   Effect: ionosphere filtering and coherence statistics include water.
   Exit: derive a mask from the Copernicus WBM layer before the full run.

5. Frame 3700 only; the 3700 vs 3650 framing question vs Tay is unresolved.
   Reason: frame 3700 contains Tay's exact first and last scene dates and
   fully covers the AOI (verified empirically from pair-test geocoding).
   Exit: confirm Tay's frame set when comparing products; never calibrate.

6. Multilook posting at alos2App WBD defaults (5 range x 28 azimuth looks).
   Reason: defaults accepted for the pilot; posting recorded in filenames.
   Exit: re-look to the Tay-comparable ~90 m grid for the full baseline.

7. Swath selection F2-F4, validated empirically after processing: the
   pair-test geocoded footprint spans 112.997-115.585 E, -5.931 to -2.334,
   fully containing the AOI (114.21-114.83 E, -3.0 to -2.42).

8. Ionosphere screens observed zeroed (empirical, 2026-09-01): the ion.geo
   layer of a sampled 42-day pair is entirely zeros although all runIon steps
   executed with about 97 percent burst synchronization. This is alos2App QA
   rejecting and zeroing the estimate, so affected pairs are effectively
   ionosphere-uncorrected. Quantify the zero-screen fraction across all 70
   pairs at stack-assembly time and carry it in the pilot error budget.

9. Water body root cause upgrade (2026-09-01): the missing water body file was
   not cosmetic. alos2App auto-downloads SRTM SWBD tiles from the retired
   LP DAAC endpoint, every fetch fails, and the resulting empty mask zeroed
   the unwrapped phase of all 70 pairs at the unwrap step (radar and geocoded
   products both 0 percent valid; coherence and wrapped phase unaffected).
   Fixed by supplying a synthetic all-land water body on the DEM grid; water
   pixels are handled by coherence masking in the inversion instead. Full
   pair reprocess required; per-pair checkpoints made the restart one
   command.

10. Three pairs dropped for deterministic matching failure (2026-09-01):
    2018-10-12_2018-12-21 (70 d), 2019-03-15_2019-06-07 (84 d) and
    2019-07-05_2019-09-27 (84 d, spanning the 2019 fire season) fail
    reproducibly at runSlcOffset with "Bad match at level 1": surface change
    over these long baselines defeats cross-correlation tie-point matching.
    Network attrition 70 to 67 pairs; every acquisition retains at least
    three links, so connectivity and redundancy are preserved.
