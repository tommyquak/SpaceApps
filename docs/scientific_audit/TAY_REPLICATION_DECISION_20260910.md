# Tay replication feasibility: decision

Written 2026-09-10 after the [replication verdict](TAY_REPLICATION_VERDICT_20260910.md),
the [recovery decision](RECOVERY_DECISION_20260910.md) and the
[route without JAXA or CPKY](NO_EXTERNAL_REFERENCE_ROUTE_20260910.md), all at main 989aa5e.
The verdict's findings were treated as claims to verify, not conclusions. Tags: [E] documented
or measured, [S] inference, [H] hypothesis. Nothing here tunes a result toward Tay et al. (2025),
and no NISAR agreement is counted as historical replication.

## Decision

[S] **Exact replication of Tay et al. (2025) is not feasible with data we hold or can fetch
without JAXA, and our processing is no longer the suspect.** The matched-input strip test shows
that correcting the confirmed processing differences (inversion order, deramp, reference) on the
half of Tay's area that frame 3700 covers does not produce agreement with the published rate.
The remaining unresolved item is the input scene itself (path 28 frame 3650, JAXA-provided) and
Tay's datum (CPKY GNSS), neither of which exists at ASF or anywhere in our storage.

What is reproducible now: the chain (alosStack 17x56, bandwidth-20, split-spectrum ionosphere,
MintPy var weights, SET, Zheng) on frame 3700, including a Tay-order inversion over the overlap
strip. What is not: Tay's input and datum, and therefore Tay's map or its accuracy claim.

Consequences:
1. Keep the correct claim: "replication not achieved; the published input is unavailable to us".
   Do not claim Tay is impossible to replicate, and do not spend further cluster time on frame
   3700 processing variants for the purpose of matching Tay.
2. Tay's raster stays an external benchmark layer; ALOS-1 2007-2009 stays a separate relative
   historical layer; NISAR local contrasts stay the forward observable, per the recovery decision.
3. The only actions that can still change the historical verdict are external and need approval
   before sending: a JAXA request for path 28 frame 3650 (92 dates), an author request for the
   granule manifest, and a request for the CPKY GNSS series. These were not sent.

## 1. Frame 3650 versus frame 3700 (reconciled)

[E] Tay Methods: 92 ScanSAR L1.1 images, path 28 frame 3650, 2018-10-12 to 2024-08-09, reference
image 2020-06-19, reference point at the CPKY GNSS (-2.220, 113.926). The acknowledgement thanks
JAXA "for the acquisition and provision of the ALOS-2 data": the scenes came from JAXA, not ASF.

[E] ASF holds one ALOS-2 collection (C3315903479-ASF). On path 28 it lists frames 3700, 3750 and
3800 only. Frame numbers step by 50 per 2.5 degrees of latitude with 3.9 degree scenes, so frame
3650 spans about -3.6 to +0.3 and contains CPKY and all of Tay's Fig. 3; frame 3700 (catalog
-6.09 to -2.21) contains roughly the southern half of Tay's Fig. 3.

[E] Verified image coverage, not catalog footprint: the processed 17x56 geometry ends at -2.368 S
at 113.81 E, 17 km south of CPKY; the catalog polygon's northern edge at CPKY's longitude is
-2.382. The northern limit was the scene, not our processing bbox.

[E] Dates: our stack holds 93 dates, 91 inside Tay's interval (Tay: 92; 2021-04-23 is absent from
ASF on path 28 for every frame) plus 2024-09-20 and 2024-10-04 after Tay's last date.

[S] Tay's six southern camera sites (about -2.35 to -2.45) straddle the real scene edge; their
exact coordinates are in the DR-NTU per-site time-series files, which the NTU server did not
serve today (HTTP 503 on every retry).

## 2. Recoverable raw data and intermediate products (inventoried)

[E] Cluster, read-only inventory 2026-09-10:
- Raw ALOS-2 CEOS: none. `~/full/ceos` holds two dangling symlinks into the deleted
  `/mnt/scratch/t/tommy/full/ceos`; the move log records 23 dates (384 GB) moved there on 3 Sep
  and removed on 6 Sep. No IMG/LED/zip for ALOS-2 exists under home or scratch.
- Coregistered SLC stack: none (deleted 6 Sep).
- Surviving: 1,650 filtered wrapped interferograms at 17x56 in radar coordinates (26 GB) with
  geocoded unwrapped, coherence and component products; the MintPy run (98 GB); split-spectrum
  ionosphere products with the reference track geometry (1.1 GB); 170c, 170d and 85b variants;
  the PALSAR-1 zips for paths 421 and 422.
[E] Remote: the Vultr snapshot of 2026-05-19 holds the ALOS-1 PALSAR stack (14 SLCs, 79
interferograms), not ALOS-2. The OVH VM holds Sentinel-1 work already copied to the Mac and now
presents a changed SSH host key (not trusted until checked in the OVH console). GCS holds
Sentinel-1 only. The Mac holds final rasters only.
[S] The only sources for ALOS-2 scenes are ASF (frame 3700, about 16.7 GB per date as a range-fetched
crop, 1.5 TB for the campaign) and JAXA (frame 3650).

## 3. Does the chain reproduce Tay's method? (separated by evidence class)

Confirmed defects, read in code:
1. `mintpy.deramp = linear` (Tay: none).
2. Geocode each interferogram, average onto 90 m, then invert (Tay: invert in radar
   coordinates, geocode afterwards).
3. Reference by MintPy maxCoherence (Tay: fixed pixel at CPKY, outside our scene).
4. Uncertainty: per-interferogram sigma-phi averaged (Tay: per-pixel propagation of
   (1 - g^2) / (2 N g^2) through the inversion to each date, mask on mean sigma > 0.7 cm).
5. Zheng closure-bias correction at one bandwidth level (Tay: per-pixel levels).
6. Neighbour-stack troposphere run on our datum: 0.01 percent eligible neighbours, null effect.

Matched: 17x56 looks, bandwidth-20 network, inverse-variance weights, ionosphere from
bandwidth-4 with decorrelated areas masked, SET on, no topographic residual, ISCE 2.6 alosStack
and MintPy.

Unresolved assumptions: effective looks (measured coherence floor gives 41 to 75; Tay's mask
assumes nominal 952); whether Tay's frame-3650 scene unwraps better over the same ground.

Missing inputs: frame 3650 scenes, CPKY GNSS series, Tay's exact manifest, Tay's per-site
time series (NTU 503), any SLC stack.

## 4. The matched-input experiment (defects 1 to 3 tested together)

Design and scripts: [cluster/full/strip_test](../../cluster/full/strip_test/README.md).
Only surviving products were used; no scene was downloaded; about 5 GB of scratch and four short
jobs. The reading was pre-registered before any result: variant B at or above r = +0.30 with the
95 percent block-bootstrap interval clear of zero means the defects were material; both variants
below +0.20 means the data are the limit.

- Window: rows 0 to 498, columns 361 to 770 of the 17x56 grid (409 x 498 px, 158,214 inside the
  Tay overlap strip), lat -3.042 to -2.368, lon 113.701 to 114.575, incidence 39.8 to 45.0 deg.
- Re-unwrapping of all 1,650 crops: zero failures. Both variants share these crops.
- Variant B (Tay order: radar-coordinate stack, interior reference, var weights, SET, no
  deramp, invert then geocode). The first attempt's reference selector fell through to pixel
  (0, 0) because no interior pixel lies inside a connected component in 80 percent of pairs
  (in-component fraction p50 0.45, p90 0.65, max 0.76); it was cancelled and rerun with the
  reference at the highest proxy coherence among the top 0.1 percent by component presence
  (row 181, column 116, present in 76 percent of pairs; MintPy's reference mask check off).
- [E] **B result** (job 835583): r = +0.06 [-0.08, +0.20] on 49,610 comparable 90 m pixels;
  tempCoh > 0.3: +0.06 [-0.09, +0.21] (47,457 px); > 0.5: +0.04 [-0.13, +0.22] (27,846 px);
  > 0.7: -0.31 [-0.54, +0.27] (1,815 px). Our std 28.5 mm/yr against Tay's 10.9; demeaned RMSE
  29.8 mm/yr. Inside the strip the median pixel has 7,546 of 15,010 closure triplets (50 percent)
  with a non-zero integer ambiguity; temporal coherence median 0.262, 8.8 percent at or above
  0.5, 0.40 percent at or above 0.7.
- [E] **A' result** (our order on the same crops, job 835443): pending at the time of writing; the
  B result alone settles the question the experiment was designed for (see the reading below). The
  A' line will be appended to `cluster/full/strip_test/results/` when the job reports.

[S] Reading: B is below the +0.20 line with an interval that includes zero. Correcting order,
deramp and reference does not recover Tay's pattern on this scene. Half of all closure triplets
carry unwrapping errors even in a small re-unwrapped window, which is a property of the phase
data at the measured effective looks, not of the inversion order.

Caveats: same-sensor benchmark, not independent validation; correlation and demeaned residuals
are the only interpretable measures because Tay's datum is outside the scene; the comparison used
the cluster's pre-provenance `compare_tay.py` (the provenance-checked version in this tree has not
been deployed to the cluster).

## 5. Reconciliation with the 10 Sep record

- The verdict's "remote storage and backups have not been independently inventoried" is now
  closed (section 2).
- The verdict's candidate explanations (order, deramp, reference) are now tested on matched
  input and found not material on this strip; sigma-phi propagation and per-pixel Zheng levels
  remain untested but cannot create agreement where the unwrapped phase does not carry it.
- The route document's decision to drop JAXA and CPKY as prerequisites for the forward NISAR
  observable is unaffected; they remain prerequisites for exact historical reproduction only.
- The non-v2 file `nisar_native_local_contrasts_60day_20260910.json`, described as retained, is
  untracked in the working tree.

## 6. What would change this decision

1. JAXA supplies path 28 frame 3650 for Tay's 92 dates (request not sent; needs approval).
2. The authors supply the granule manifest and CPKY series (request not sent; needs approval).
3. Tay's per-site time series (DR-NTU, nine files) become downloadable: a site-level comparison
   of our strip series at their southern cameras inside frame 3700 would then be run without any
   new processing. A retry loop is watching the server.

## 7. Same-day addendum: what NISAR can do right away (reference-ensemble test)

The user's route after this decision, peat motion relative to evidence-selected mineral-ground
references with the regional signal preserved, was pre-registered and run twice on the retained
GUNW chains (see NISAR_REFERENCE_ENSEMBLE_PREREG_20260910.md, Results 1 and 2). Outcome [E]: the
peat plain has no mineral ground within about 60 km of the dam sites; references 90 to 220 km
apart disagree by tens to hundreds of mm over 60 days; and the per-pair residual after the
product's own corrections is a smooth +-20 to 40 mm pattern over the AOI that reverses sign
between consecutive pairs on the quiet dawn track (worse on the dusk track). A 60-day chain
therefore cannot resolve peat motion against mineral ground, and an interim reading of lowering
was withdrawn the same day. [S] The forward observable NISAR can deliver now is the reference-free
local contrast; a datum-referenced rate needs a year or more of 12-day pairs, an in-AOI residual
model, or a GNSS anchor, none of which a short window supplies.

## 8. Addendum 2026-09-12: Tay's manifest and validation sites, from the DR-NTU files

The nine per-site time-series files (doi 10.21979/N9/JGQE4O, CC BY-NC, downloaded 2026-09-12 after
two days of HTTP 503) each carry all 92 epochs. Dates saved as
`tay_2025_acquisition_dates_from_dataset.txt`.
- [E] Tay's 92 dates = our 91 common dates + **2021-04-23**. Our 93 = those 91 + 2024-09-20 and
  2024-10-04. The 2021-04-23 acquisition exists in Tay's stack and is absent from ASF on path 28 at
  every frame (CMR checks of 2026-09-10). This is direct evidence that the published input was not
  the ASF frame-3700 collection, consistent with the paper's JAXA acknowledgement.
- [E] Site locations (OpenStreetMap Nominatim, 2026-09-12): Kalampangan -2.279, 114.011; KHDTK
  Tumbang Nusa -2.368, 114.139; University of Palangka Raya -2.217, 113.895. Our verified frame-3700
  data edge runs from -2.368 at 113.81 E to about -2.49 at 114.45 E, so Kalampangan lies about
  13 km, Tumbang Nusa about 5 km and UPR about 20 km outside it; Hampangen and Misik are north of
  Palangkaraya. None of Tay's validation sites is inside our scene, so a site-level comparison of
  our strip series is not possible on frame 3700, and the site files add no test we can run.
- [E] Per-site mean phase standard deviation in Tay's files: 0.33 to 0.52 cm at the plantation,
  cropland-Misik and RePeat sites, 0.72 to 0.74 cm at Kalampangan and South, 0.82 (KHDTK) and
  0.99 cm (UPR) at the forest sites, i.e. above the 0.7 cm mask Tay applies to the raster. Their
  displacement ranges are -4.6 to -10.9 cm over the period, with mean sigma near 0 only at the
  reference epoch.
