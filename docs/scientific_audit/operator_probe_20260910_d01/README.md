# Fixed-site GSLC aggregation probe, September 10, 2026

## Result

[E, local measurement] Preserving complex amplitudes reduced disagreement with mission wrapped phase in both tested pairs. It did not repair the phase gate or validate displacement.

On exactly the same selected pixels, with both aggregation strengths at least 0.25 and both coherence measures at least 0.25:

| Pair | Shared pixels | Unit-phasor MAE | Complex-weighted MAE | Unit-phasor RMSE | Complex-weighted RMSE |
|---|---:|---:|---:|---:|---:|
| June 22 to July 4 | 6,757 | 0.373801 rad | 0.305430 rad | 0.522061 rad | 0.441939 rad |
| July 4 to July 16 | 7,094 | 0.368622 rad | 0.301976 rad | 0.505575 rad | 0.422373 rad |

MAE falls by about 18% on this common support. Complex-weighted resultant lengths are 0.9143 and 0.9194. Both pairs prefer the expected same-sign convention.

[E, repository policy] The existing phase gates require MAE <=0.10 rad, RMSE <=0.15 rad and resultant >=0.98. Both modes fail these numerical gates here. The 2 km crop also intentionally cannot satisfy the million-pixel spatial-support gate. These are engineering gates, not published uncertainty bounds.

## Selection and comparability

The first feature of the existing `deliverables/tay_full_grid/dam_sites.geojson` was fixed before inspecting phase or coherence: D-01 at longitude 114.68309869375906, latitude -2.683670460424211. A 100 by 100 pixel window on the existing 20 m mission grid defines the same approximately 2 km crop for both pairs. No search for a favorable location was performed.

Only bounded reads of existing local June 22, July 4 and July 16 GSLC subsets and matching mission rasters were used. Original inputs were not modified. Generated data and reports occupy approximately 5.7 MB.

The modes share source quality masks and fine coherence. Their aggregation-strength masks can differ because amplitude weighting changes the resultant. Separate reports retain each mode's ordinary selection. The table instead uses the intersection, recomputing one constant circular offset for each mode on those identical pixels. Results with zero minimum strength are also recorded; those improve from 0.4682 to 0.3641 rad and 0.4465 to 0.3454 rad respectively. No scale, planar fit, phase-error clipping, Tay fitting or physical calibration was applied.

## Interpretation

[S] The amplitude-discarding operator contributes to local disagreement. This bounded result does not establish its importance elsewhere, statistical significance across independent sites, or a geophysical signal.

[E, code behavior] Both tested operators aggregate already geocoded GSLC cross-products. Neither reproduces mission radar-coordinate multilooking followed by geocoding. The arithmetic mean of local fine coherence also is not a reconstruction of mission coherence.

[S] Retain the complex export for future operator tests, but keep fine displacement gated. The remaining mismatch needs separately predeclared tests of product processing compatibility, geocoding support and filtering. This experiment provides no historical Tay replication, absolute reference, annual velocity or independent accuracy validation.

## Reproduction and artifacts

`run_probe.py` contains the complete bounded experiment. Run from the repository root with `PYTHONPATH=src` and the repository Python environment. Supply `--output-dir` pointing to a new directory; existing directories are refused.

`data/summary.json` records fixed windows, input file sizes and modification times, individual and common-support metrics, and SHA-256 hashes of generated raster crops. Large source HDF5 files were not hashed in full, so the input inventory is explicitly not a content-hash guarantee. Each pair folder also contains the four full validator reports and five raster products. The source granule and ordered acquisition-date tags are checked in complex-weighted mode.
