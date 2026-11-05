# Strip test: matched-input check of the Tay processing order on frame 3700

Cluster-only experiment, 2026-09-10 (jobs 835439, 835440, 835583, 835442, 835443).
It answers one question raised by the
[Tay replication verdict](../../../docs/scientific_audit/TAY_REPLICATION_VERDICT_20260910.md):
on the part of Tay et al. (2025) area that frame 3700 does cover, do our confirmed
processing differences (geocode-then-average before inversion, linear deramp,
coherence-picked reference) explain the lack of agreement with the published rate,
or is the input itself the limit?

No scene was downloaded and no SLC exists any more; the test uses only the surviving
1,650 filtered 17x56 interferograms in radar coordinates (`~/full/pairs_tay`), the
ionosphere-stage geometry (`pairs_ion`) and the geometry of the original MintPy run.

## Design

- Window: rows 0..498, columns 361..770 of the 17x56 radar grid (409 x 498 px,
  158,214 px inside the strip 113.8..114.5 E, -2.90..-2.20 S), lat -3.042..-2.368,
  lon 113.701..114.575, incidence 39.8..45.0 deg. This is the northern edge of the
  scene; CPKY is 17 km further north.
- Stage 1 `strip_prep.py`: interpolate the 80x448 ionosphere-stage lat/lon to 17x56,
  derive the window, sample geometry (incidence, azimuth, height, slant range) from the
  original geocoded geometry into a radar-coordinate `geometryRadar.h5`.
- Stage 2 `strip_unwrap.py` (array): crop each filtered interferogram, re-unwrap once
  with snaphu (SMOOTH, MCF init, NCORRLOOKS 41 = measured effective looks). 1,650/1,650,
  zero failures. Both variants share these crops, so the unwrapper cannot explain a
  B versus A' difference.
- Variant B `strip_stackB.py` + `strip_runB.sbatch` (Tay order): radar-coordinate
  ifgramStack with the dates, baselines and attributes of the original run; reference
  = interior pixel with the highest proxy coherence among the top 0.1 percent by
  connected-component presence (no interior pixel is inside a component in more than
  76 percent of pairs, so MintPy's reference mask check is switched off); var weights,
  SET on, no deramp, no tropo, no topographic residual; invert, then geocode
  (3 arcsec); `compare_tay.py` on the geocoded velocity.
- Variant A' `strip_geoA.py` + `strip_geoA.sbatch` + `strip_runA.sbatch` (our order):
  pad each crop to the frame, `geocode.py` with the track XML, average onto the 3 arcsec
  grid with `build_stack2.py`, MintPy with maxCoherence reference and linear deramp,
  `compare_tay.py`.

Pre-registered reading (written before any result): B at or above r = +0.30 with the
95 percent block-bootstrap interval clear of zero means order, deramp and reference were
material; both variants below +0.20 means the limit is the data.

## Results

See `results/`. Variant B: r = +0.06 [-0.08, +0.20] on 49,610 comparable 90 m pixels;
tempCoh > 0.7 leaves 1,815 px at r = -0.31 [-0.54, +0.27]. Our std 28.5 mm/yr against
Tay's 10.9; demeaned RMSE 29.8 mm/yr. Inside the strip the median pixel has 7,546 of
15,010 closure triplets (50 percent) with a non-zero integer ambiguity; temporal
coherence median 0.262, 8.8 percent at or above 0.5.

Caveats: same-sensor benchmark, not independent validation; Tay's datum (CPKY) is
outside the scene, so only correlation and demeaned residuals are interpretable; the
comparison used the cluster's pre-provenance `compare_tay.py` (the provenance-checked
version in this tree has not been deployed to the cluster).

## Reproduce on the cluster

```bash
cd ~/full/strip_test
sbatch strip_prep.sbatch                      # window.json, geometryRadar.h5
sbatch --array=0-7 strip_unwrap.sbatch        # ROOT/pairs/<pair>/{unw,coh,cc}.npy
sbatch strip_runB.sbatch                      # mintpyB/, compare_tay
sbatch --array=0-7 strip_geoA.sbatch          # ROOT/pairsA/<pair>/*.aoi90.tif
sbatch --dependency=afterok:<geoA> strip_runA.sbatch
```

Paths are set in `strip_common.py` (PAIRS, ION, REFSTACK, GEOGEO, ROOT, OUT).
Scratch use is about 5 GB.
