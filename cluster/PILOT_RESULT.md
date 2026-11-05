# ALOS-2 ScanSAR pilot: result and diagnosis (2026-09-01)

## Outcome
MintPy attempt four (job 786802) completed: 67-pair weighted SBAS inversion over the AOI,
20 epochs 2018-10-12 to 2019-11-08, products velocity/temporalCoherence/avgSpatialCoh
(.h5 + .tif) in ~/pilot/mintpy. The processing chain is validated end to end.

The velocity field is NOT a usable subsidence product and must not be released:
median LOS -91 mm/yr (vertical ~ -113 mm/yr), ~30x the Tay et al. 2025 benchmark and
beyond physical peat subsidence rates (Hooijer et al. 2012).

## Measured diagnosis
Per-pair median LOS displacement on a common coherence>=0.45 mask, by temporal baseline:

| dt (d) | pairs | med disp (mm) | apparent mm/yr |
|-------:|------:|--------------:|---------------:|
| 14 | 10 | -0.5 | -14 |
| 28 | 11 | -11.4 | -148 |
| 42 | 17 | -12.2 | -106 |
| 56 | 9 | -9.2 | -60 |
| 70 | 10 | -86.5 | -451 |
| 84 | 10 | +96.1 | +418 |

1. 14 d pairs nearly clean: the interferometric machinery is sound.
2. 28-56 d pairs carry a near-constant ~-11 mm per-interferogram offset independent of
   span: the closure / fading-signal phase bias of multilooked InSAR over vegetated,
   moisture-cycling ground [E: Ansari et al. 2021, IEEE TGRS; Zheng et al. 2022, IEEE
   TGRS; S: attribution here]. Divided by short spans it inflates to -60..-150 mm/yr;
   at mean span ~45 d it reproduces the map median (~-90 mm/yr) almost exactly.
3. 70-84 d pairs approach +-1 L-band fringe (119 mm) with sign flips: snaphu cycle
   slips at coherence ~0.27, injecting the heavy tails (std ~160 mm/yr).

Corroboration: velocity magnitude falls monotonically with temporal coherence
(-190 mm/yr at tcoh 0.25-0.4 vs -53 mm/yr at tcoh>=0.9); a planar fit explains only 7%
of variance (not a ramp or reference artifact); 2018-2019 is solar minimum, so
ionosphere is a secondary suspect.

## Implication
The pilot measured the exact failure mode PILOT_DEVIATIONS.md predicted (closure-bias
correction cannot converge at 20 scenes). Tay et al. 2025 avoid it via archive depth:
long pairs dilute a constant per-igram bias, and the Zheng 2022 correction removes the
rest. The full 92+ scene baseline with bandwidth-20 and closure correction is therefore
required for a physical velocity, not optional refinement. Never force-calibrate to Tay.

## Cheap follow-ups on the existing stack (no reprocessing)
- Enable mintpy.unwrapError.method=phase_closure and re-invert (targets 70-84 d slips).
- Invert the 14 d chain alone as a minimally biased cross-check.
