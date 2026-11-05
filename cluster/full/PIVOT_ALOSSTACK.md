# alosStack pivot: implementation specification (v5, audit rounds 1-6 applied)

Authentic alosStack (share/isce2/alosStack, the Tay et al. 2025 code path) on
93 dates, bandwidth-20 (1650 pairs), ion network nsub_ion=4 (362 pairs).
Pairwise campaign stopped; 229 pairs banked as cross-check. Looks: nrlks1=1,
nalks1=14 (looks1), nrlks2=5, nalks2=2 (looks2 = 5x28 total). Reference date
181123 (staged). Every rule below cites the audited upstream contract.

## 1. Crop window (frame-level, applied per swath)
- ONE absolute azimuth window [T0, T0+len*dt] for all swaths of every date;
  per-swath range windows derived from absolute [R0, R0+width*dr]. Offsets are
  exact integers on the resampled common grid (resample_common_grid snaps
  origins to nrlks1/nalks1 multiples).
- Divisibility (REQUIRED, not conservative): lineOff % 28 == 0, cropLen % 28
  == 0, sampOff % 5 == 0, cropWidth % 5 == 0, measured on the resampled
  full-frame grid, so ml1 (1x14) and ml2 (5x2) rasters crop on exact look
  bins and every downstream int() truncation is exact.
- AOI must lie strictly inside every swath's resampled azimuth extent and
  keep the s2/s3 and s3/s4 range-overlap strips alive (swathMosaic phase
  estimation needs them; the phaseDiffFixed snap list is only a fallback).
- There is NO per-date SLC mosaic in alosStack: crop-after-date-mosaic was
  evaluated and rejected. SLCs stay per swath; mosaicking happens per pair on
  interferograms only.
- Assertions: after metadata rewrite + mosaic_parameter re-run, frame dims
  must equal (cropWidth, cropLen/14) and track fields match; after the first
  pair mosaic, binary size must equal frame dims x 8 B. Both run in a cropped
  3-date rehearsal before the fleet.

## 2. Metadata rewrite (ONE target)
Rewritten to cropped geometry: ONLY dates_resampled/<ref>/<ref>.track.xml and
f1_3700/<ref>.frame.xml, by loadTrack, mutating per swath: numberOfSamples=
width_s, numberOfLines=cropLen, startingRange += sampOff_s*rangePixelSize,
sensingStart_s = T0 exactly, i.e. += lineOff_s*azimuthLineInterval with per-swath integer lineOff_s = (T0 - sensingStart_s)/dt (frame-level lineOff % 28); saveTrack
with cwd inside dates_resampled/<ref>; then RE-RUN mosaic_parameter.py on the
cropped tree (never hand-write frame/track mosaic fields). Untouchable: prf,
rangePixelSize, rangeSamplingRate, azimuthLineInterval, azimuthPixelSize, wavelength, range/
azimuth bandwidth, doppler and fmrate polynomials, orbit, operationMode,
frame/swath numbers, burst parameters. The dates/<d> track/frame xmls (ALL
dates) persist VERBATIM full-frame (ion_filt and unwrap_snaphu read original
parameters; pair_up copies them as-is). Per-swath .slc.xml/.vrt are
regenerated at unpack with cropped dims. No .pck files exist in these trees.

## 3. Subband regeneration (wave-3 ion jobs only) and the demodulation ramp
Only ion_subband consumes subbands (~362 pairs); cmd_2 scope does not. At ion
staging, per swath: render cropped full-band .slc.xml AND .slc.vrt
(SourceFilename relative, ByteOrder LSB, ImageOffset 0, LineOffset 8*width);
get (wavelengths, bw, centers, prefixes) from StackPulic.subbandParameters of
the cropped reference track; call contrib.alos2proc rg_filter(slc, 2,
[lower, upper], bw, centers, 257, 2048, 0.1, 0, 0.0) with cwd = swath dir.
PHYSICS (audit r2 finding 11): upstream demodulates each date on its own
original grid, so authentic subband phase carries -2pi*bc*(off_ref-off_sec).
Wave 1 therefore ALSO persists each date's applied range/azimuth offset
polynomials (tiny text), and each regenerated subband SLC is multiplied by
exp(-i*2pi*bc_k*(off_d(az, rg) + sampOff_s)) at cropped 0-based (i, k), az =
lineOff_s + i + 1, rg = sampOff_s + k + 1 (Amendment A.7: rg_filter with
zero_cf=0, offset=0 demodulates sample j of the file it is given by
exp(-i*2pi*bc*j); the authentic subband was demodulated on the original grid
at j = rg - 1 + off_d, the regenerated one on the cropped file at j = k, so
the ratio carries the per-swath constant sampOff_s in addition to off_d; the
constant cancels in every interferogram). DECISION GATE on the wave-0 smoke
output: plain regeneration is EXPECTED to fail vs the authentic subbands;
the ramp-corrected variant must pass masked mean |dphi| < 1e-2 rad at 1x14
looks AND match the lower-minus-upper combination to the same tolerance;
evaluated per SLC the sampOff_s term is mandatory (evaluated on subband
interferograms it cancels).
DECISION (2026-09-03, superseding the gate): the fallback is now the primary
path. Regeneration existed only because persisting subbands did not fit the
~425 GB home quota; the crop store moved to /mnt/scratch (54 TB free) after
that quota killed SEED and a wave-1 lane, so wave 1 persists each date's
AUTHENTIC lower/upper subbands (~5.8 GB/date, ~540 GB for 93 dates) exactly as
resample_common_grid produced them. This removes the ramp reconstruction
(Amendment A.7), the sec-3 A/B gate, and round-14 finding 20 (rg_filter
zero-pads outside the file it is given, so a subband regenerated from the crop
carries a 128-sample transient at every interior crop edge that no ramp can
restore and the whole-crop gate cannot see). It also removes the repeated
regeneration work: dates recur across ion blocks at stride 20, so regeneration
would have run about twice per date. ION blocks now unpack authentic subbands
like any other crop. Set PERSIST_SUBBANDS=0 in wave1.sbatch to revert to
regeneration. The text above is kept for the physics record.

## 4. Wave-1 persist manifest (per date d)
A. dates_resampled/d/f1_3700/s{2,3,4}/d.slc cropped, complex int16, sidecar
   JSON {width_s, length, sampOff_s, lineOff, scale, offset_polys} = single
   source of truth (never infer dims from file size), plus d.slc.ovf: the
   samples above the int16 range (about 1 - pctl of the swath, the
   brightest scatterers) stored exactly as complex64 with their flat
   indices (Amendment A.11; tens to hundreds of KB per swath).
B. dates_resampled/d/insar/d_1rlks_14alks_rg_rect.off cropped at the ml1
   window PLUS its rewritten .off.xml (float, cropped ml1 dims; imageMath in
   diff_interferogram loads via the xml). Values are window-relative; valid
   under pure crop. Reference date
   has none (no consumer expects it). rg_rect.off is a SYMLINK to rg.off when
   rect_range_offset found an identity affine (fewer than 50 culled DEM
   offsets, or a tiny land ratio): crop it (dereferencing the link, regular
   file of ml1_width x ml1_length x 4 B verified) BEFORE discarding
   rg.off/az.off (Amendment A.8).
C. dates/d/d.track.xml + f1_3700/d.frame.xml VERBATIM full-frame.
D. baseline/baseline_center.txt only (.rmg grids dropped: no consumer).
   Written ONCE by an authoritative compute_baseline over the complete
   persisted dates/ tree, run serially at the BARRIER.
Once: E. cropped-rewritten resampled ref track/frame xmls (sec 2);
F. cropped ref insar/ geometry ml1 {lat,lon,hgt,wbd,los} + affine_transform
   .txt (ml2 geometry regenerated by look_geom at stage time, or persisted);
G. dates/<ref> xmls + f1_3700/mosaic/swath_offset.txt verbatim
   (resample_common_grid reads swath_offset.txt unconditionally).
Budget: ~2.8 GB/date int16 crops -> ~260 GB + MBs of metadata.

## 5. Stages
SEED (T0, one job, full frame, ref already staged): read_data(ref+181012),
estimate_swath_offset, estimate_frame_offset, estimate_slc_offset(neighbor),
resample_common_grid(both), mosaic_parameter (dates/ and dates_resampled/
variants), rdr2geo; PUBLISH ref bundle (gates every wave-1 job: swath_offset
.txt + ref track loadable); then form+mosaic the ref-181012 pair at looks1,
radar_dem_offset -> affine_transform.txt, rect_range_offset(neighbor),
look_geom. ~2-4 h; the fanout gate opens ~2 h in.
WAVE 1 (per date, full-frame-then-crop, packed 3-4 dates per node lane set):
read_data -virtual, mosaic_parameter(date), estimate_slc_offset (ALOS2_FORCE_GEOM_OFFSET=1 armed when
|date - 181123| > 90 days, per the measured ~3-month ampcor failure), resample_common_grid -subband on
node-local, geo2rdr vs full-frame ref geometry, rect_range_offset with the
full-frame affine (each lane WAITS for affine_transform.txt to appear in
the published bundle before its rect step; earlier steps proceed), THEN
crop per sec 1-2 and persist per sec 4. CEOS
logistics: home holds 23 staged dates; pin the reference CEOS; process
staged tranche first (unlocks seed + first ion pairs), download
chronologically from the first unstaged date, delete each date's CEOS after
its crops land (net -13.5 GB/date); tranche cap ~12-14 CEOS beside the crop
store. Wall ~7-9 h, download-bound (1.14 TB remaining at 55 MB/s).
ION PASS (pre-barrier, 362 pairs only): date-block jobs (blocks of ~24
dates at stride 20, i.e. 5 blocks over 93 dates at the nsub_ion=4 window) stage crops, regenerate ramp-corrected subbands
ONCE per staged date, run per pair: cmd_2 scope (form/mosaic/diff, staging each date's persisted cropped rg_rect.off) +
ion_subband/ion_unwrap/ion_filt + the ionpair1 check-prep; persist ONLY
pairs_ion/<pair>/ion/ion_cal products plus each pair dir's <ref>.track.xml
copy (KBs; ion_ls loads it from the first pair dir). ~1-2 h behind wave-1 tranches.
BARRIER: automated ion QC on the persisted ion_cal products (per-pair
filt_ion std and residual thresholds; any pair failing goes into -exc_pair
automatically and is logged ION_EXCLUDED for review; connectivity of the
remaining ion network is asserted before release), then ion_ls (serial, minutes); its output dates_ion/ (per-date .ion
grids + xmls, small) is persisted to home.
MAIN PASS (post-barrier, all 1650 pairs): 7 date-block jobs (B=32, stride
12; each pair assigned to its first containing block, ~200-430 pairs/block;
each date staged <=3x total ~0.6 TB home reads). Each 128-core job runs GNU
parallel at 48-64 lanes, OMP_NUM_THREADS=1-2, per pair in ONE lane:
form/mosaic/diff (ml1), look_coherence (produces the ml2 diff that
ion_correct renames and corrects), ion_correct (folds the upstream 14 h
serial loop into lanes), filt, unwrap_snaphu, geocode. NO looks1 pair product
ever persists to home (655 GB float / ~330 GB int16 would breach quota).
Persist only looks2 geocoded products MintPy ingests. Wall ~2-4 h.
FINALE INTERFACE (audit r4 finding 1): build_stack2 discovers pairs as
~/full/pairs/<YYYY-MM-DD_YYYY-MM-DD>/ dirs with a DONE marker and fanout2
file names. MAIN PASS lanes therefore persist per pair into that layout:
gdalwarp each geocoded raster onto the aoi90 grid with the exact names
filt_*.unw.geo.aoi90.tif, *.cor.geo.aoi90.tif,
filt_*.unw.conncomp.geo.aoi90.tif (dates translated YYMMDD ->
YYYY-MM-DD), the reference *los.geo.aoi90.tif into one pair dir, DONE
touched last. BANK RELOCATION: before the MAIN pass launches, mv
~/full/pairs -> ~/full/pairs_pairwise_banked (cross-check material) and
recreate an empty ~/full/pairs; find_pairs in build_stack2 now discovers the
MAIN output root ONLY (pilot and pairwise banks dropped from the corrected
stack by decision, recorded here); MAIN lanes refuse any pair dir already
containing an isce.log. BPERP BRIDGE: build_stack2.bperp_of gains a fallback that
reads baseline/baseline_center.txt (manifest D) when no isce.log exists;
without this the stack silently loads bperp=0 and corrupts
topographicResidual. The bperp fallback is a small audited change to the
rehearsal-validated builder and is re-verified in the cropped 3-date
rehearsal. FINALE: mintpy2 path on that stack.

SEED metadata commands (implementing sec 2, added 2026-09-02): crop_pack
rewrite <dates_resampled/<ref>> <ref_date> mutates ONLY the per-swath
numberOfSamples/numberOfLines/startingRange/sensingStart of the resampled
reference track+frame xmls (saveTrack with cwd inside the date dir), refusing
a tree whose dims no longer match ref_geometry.json (already rewritten, or
not the tree the geometry came from) and asserting sensingStart_s == T0 for
every swath; the caller then re-runs mosaic_parameter.py, and crop_pack
verify <dates_resampled/<ref>> <ref_date> asserts the sec-1 contract on the
result (frame ml1 == cropWidth x cropLen/nalks1, frame and track origins ==
the window origin, per-swath dims == the window). tests/test_seed_metadata.py
runs the whole chain on real ISCE products rebuilt from the smoke gate tree
(cluster only; skips without the gate tree or ISCE_STACK) and checks that a
second rewrite and a mismatched window are both refused.

## 5b. Block-job staging protocol (every ION and MAIN block job)
Stage into the node-local tree before pair_up: the once-bundle E+F+G, the
manifest-C verbatim dates/<d> xmls for every date the block touches, AND the
reference date directory skeleton even when 181123 is outside the block's
date range (stackDateStatistics raises without it and enumerates frames and
swaths from it). Unpack int16 crops to complex64 via the sidecar scale,
restoring the overflow samples from the .ovf side file (unpack_slc refuses a
missing or foreign side file);
render .slc.xml/.vrt for every unpacked SLC and the rg_rect.off xmls BEFORE
pair_up; mirror the same bundle into the ion pair tree for ION blocks; MAIN blocks
additionally stage dates_ion/ at the node-local stack root (MAIN lanes invoke ion_correct.py with cwd = the pair dir and
-ion_dir ../../dates_ion; the script prepends ../ after chdir insar, i.e.
../../../dates_ion resolved from insar/, so dates_ion/ sits at the
node-local stack root; ion_correct raises if absent).
DEM/water staging rule (r1 finding 1): every job that stages
elevation.dem.wgs84, elevation_3s.dem.wgs84 or waterBody.wbd to node-local
runs fixImageXml.py -i <file> -f on each, fails fast on a missing xml
sidecar, and never '|| true's the copy.

## 6. Geocoding config (audit r2 finding 2; r3 finding 1)
'dem for geocoding' = elevation_3s.dem.wgs84 (3 arcsec, ~90 m, matching
5x28 looks: geocode.py takes the output posting from the DEM it is given,
so the 1-arcsec DEM would emit 9x larger products). The interactive build
on 2026-09-02 FAILED inside ISCE's stitcher (a DEM xml loaded as a plain
Image: 'the attribute corresponding to the key reference is not present')
and left no product; build_dem3.py now derives it from the 1-arcsec
elevation.dem.wgs84 by keeping every third point (the grid is
pixel-is-point, so the kept points stay on the same geodetic locations)
with a 3x3 mean to suppress aliasing, and renders the xml through
isceobj.createDemImage. Verified 2026-09-02: 7201x8401, ISCE loads it with
delta 0.000833333 and origin 111.0/0.0, the 3x3 mean reproduces to 6e-5 m,
AOI heights 0..557 m (mean 53).
BBOX MECHANISM: create_cmds.py line 1262 is BUGGY for bbox (emits a literal
'-bbox /'), so wave scripts call geocode.py DIRECTLY with
-bbox -3.4/-2.34/113.8/115.2 (s/n/w/e slash-joined), bypassing create_cmds
for the geocode step; no third env patch. Per pair geocode filt unw + cor
only (skip the msk duplicate: bit-identical under the all-land water body);
ONCE geocode the reference los (MintPy geometry) and, nearest-interpolated,
the snaphu conncomp per pair (build_stack2 prefers real conncomp).
~100 MB/pair -> ~165 GB worst case, inside quota; geocode <1 min/pair.

## 7. Queue accounting (MaxJobs 16, MaxSubmit 32)
seed(1) + wave-1 tranche jobs(4-8) + ion blocks(4-7) + ion_ls(1) + main
blocks(7) + finale(1) <= 32 across two submission rounds, chained with
--dependency=afterok; no churning submitter needed.

## 8. ETA
smoke + spec audits today; seed+wave1 overnight (download-bound ~7-9 h); ion
pass + barrier + main pass ~4-6 h; finale ~2-4 h. Corrected velocity
expected ~Sep 3-4.

## Environment patches recorded
isce_runSlcOffset_geom_override.patch,
isce_alosstack_estimate_slc_offset_geom_override.patch (.orig backups).

## Operational note (2026-09-02 afternoon)
The stopped pairwise controller chain revived once: scancel raced the
successor-at-top spawn, and the survivor kept downloading scenes for ~3 h.
Killed via the designed switch (touch CAMPAIGN_DONE retires every
generation) plus scancel of live jids. Lesson for the runbook: stopping the
chain means CAMPAIGN_DONE first, scancel second. Side effect: 23 CEOS dates
now staged (wave-1 tranche 1 grows), home at ~427 GB.

## Amendment A (audit round 10, 2026-09-02): crop_pack.py contract
Round 10 audited the first wave-1 code (crop_pack.py) and confirmed four
defects against secs 1-4; the fixes below are now part of the contract.
1. Per-swath azimuth start. Each resampled swath starts k_s*14 lines after
   the frame origin (resample_common_grid snaps sensingStart to
   sensingStartMinimum + round(offsetAzimuth)*14*dt; mosaic_parameter defines
   the frame origin as the minimum). The frame-level lineOff (% 28) therefore
   maps to a swath-local lineOff_s = lineOff - 14*k_s, which is a multiple of
   14, NOT necessarily of 28. crop_pack derives and stores lineOff_s per
   swath in crop_window.json ("swaths") and in every sidecar ("lineOff";
   the frame-level value is "frame_lineOff"). Assertions: k_s >= 0 and
   k_s % 14 == 0; 0 <= lineOff_s and lineOff_s + cropLen <= numberOfLines_s.
   The sec-2 rewrite sets sensingStart_s = T0 for every swath, which is only
   consistent with the binaries when the crop used lineOff_s.
2. Window validation (sec 1 MUSTs now enforced at derivation): the azimuth
   window is clipped to the extent common to all swaths (snapped inward to
   28; recorded as "az_clipped"); the range window is widened outward (on
   the 5-sample quantum) so any s/s+1 overlap strip it touches is included
   whole ("rg_widened"); a swath entirely outside the range window is an
   error, as is a strip not whole after cropping; pack refuses windows
   outside the input and sizes that disagree with the geometry.
3. Applied offset polynomials (sec 3). The persisted polynomials are the ones
   resample_common_grid actually applied, not a refit of raw cull.off:
   secondary dates: each cull.off point is mapped to the resampled swath
   grid (x1 = x*dr_o/dr_r + (R0_o - R0_r)/dr_r, y1 likewise with dt and
   sensingStart), the TOTAL offset (x + dx) - x1 is fitted with
   OffsetField.getFitPolynomials(2, 2, maxOrder=True) (Fortran offsetpoly,
   weights 1 + 1/snr, azimuth mean shift), i.e. the identical call; the
   reference date: the analytic grid-snap polynomials of the
   idate == dateIndexReference branch (rangeOrder 1/azimuthOrder 0 and
   0/1). Serialisation is the Poly2D convention {coeffs, meanAzimuth,
   normAzimuth, meanRange, normRange}; evaluation coordinates are the 1-based
   line/sample of the RESAMPLED SWATH grid as resamp_slc.f90 evaluates them
   (az = lineOff_s + i + 1, rg = sampOff_s + k + 1 for 0-based cropped
   (i, k)). "orig" values are the offset-ADJUSTED original swath parameters
   (swath_offset.txt matching column accumulated as resample_common_grid
   does), "res" values the PRE-crop resampled parameters. No sidecar may ship
   without polynomials; pack fails otherwise.
4. Reference geometry file (manifest item H). crop_pack.py geometry writes
   ~/full/ref_geometry.json {per swath orig/res parameters, lines/samples
   from the frame origin, frame origin, frame ml1 dims, subband parameters}
   from the raw and the pre-crop resampled reference tracks. SEED order is
   fixed: geometry -> derive -> pack --reference (per swath) -> sec-2
   rewrite -> mosaic_parameter re-run. build_geometry refuses a tree that
   already looks rewritten. Single-frame stacks only (frame_offset.txt
   accumulation is not replicated; frame 3700 only).
Sidecar keys (sec 4A, superseding): {swath, width, length, sampOff, lineOff
(swath-local), frame_lineOff, frame_sampOff, in_width, in_length, scale,
dtype, pctl, peak, max_amp, nonzero_fraction, overflow_count, overflow_file,
overflow_layout, source, date, offset_polys{rg, az, fit: isce|analytic
(numpy only with CROP_PACK_NUMPY_FIT=1, recorded), n, source, eval}}.

Round 11 (2026-09-02, Fable 5.1; 3 lenses, 9 findings, 5 verified 2-of-2,
4 minor applied unverified) added:
5. Per-swath range windows carry NO quantum. resample_common_grid snaps
   range origins to nrlks1 = 1 multiples and numberOfSamples_s =
   round(N*dr_o/dr_r), so samples_from_origin and width_s are arbitrary
   integers; only the frame-level sampOff_frame and cropWidth are on the
   5-sample grid (sec 1), the range analogue of item 1. crop_pack asserts
   instead that min/max of the per-swath frame extents equal the frame
   window exactly (what mosaic_parameter's re-run and swathMosaic's
   outWidth depend on). The earlier per-swath % 5 check aborted derive for
   ~99 percent of real geometries (frame 3700, 181123 reference, from the
   smoke run's swath_offset.txt: s3 at 4100 and s4 at 9412 samples from the
   frame origin, swath ends 4440/9770/15206, so 9412 % 5 = 2 would have
   aborted derive; the tests/selftest.py fixture uses synthetic 7623/15206).
6. Azimuth line offsets are recovered in the time domain: sensingStart_s =
   t_min + timedelta(seconds=k*14*dt) is microsecond-quantised (timedelta
   and the ISCE xml %f field), so the residual is tested as
   |dt - k*ali| <= 1.5 us, not as 1e-3 lines (0.5e-6*prf reaches 1.7e-3
   lines at ALOS-2 ScanSAR PRFs; frame 3700 s2 k = 42 has a 1.5e-3-line
   residual). The k >= 0 and k % 14 == 0 checks remain the integrity test.
   Range offsets keep the 1e-3-sample tolerance (float64 arithmetic only).
7. Ramp restoration constant (sec 3 text amended): multiply the regenerated
   subband by exp(-i*2pi*bc_k*(off_d(az, rg) + sampOff_s)); the sidecar
   (rg polynomial, swath-local sampOff/lineOff, eval statement) is
   sufficient; the azimuth polynomial is not needed for the range ramp and
   bc_k is the same normalised value for every swath.
8. crop_pack cropoff is the crop for every BIL manifest raster: rg_rect.off
   (float32 1-band), ml1 lat/lon/hgt (double 1-band), wbd (byte), los
   (float32 2-band), band by band with an explicit dtype; a dangling
   rg_rect.off link is refused; sec 4B discard order fixed (crop before
   discarding rg.off/az.off). render uses the matching FLOAT/DOUBLE/BYTE.
9. pack cross-checks provenance: date and swath come from the SLC path
   (.../<date>/f1_*/s<n>/<date>.slc, or --date), a --cull path must carry
   the same <date> and s<n> components, --reference is accepted only for
   the reference date and --cull only for secondaries, and the sidecar
   records "date" so wave-3 asserts sidecar date == job date before
   applying the ramp.
10. build_geometry imports StackPulic from $ISCE_STACK/alosStack (or the
   directory of resample_common_grid.py on PATH) when it is not on
   PYTHONPATH; the SEED job needs ISCE_STACK exported.

Golden replay (2026-09-02, after round 11): tests/test_golden_replay.py
rebuilds the 181123 reference geometry from the smoke run's printed
swath_offset table and minima and checks crop_pack against the polynomials
resamp_slc actually applied (tests/golden/: authentic 181012 cull.off for
swaths 2-4, the six printed reference polynomials, the 6-digit applied
polynomials and all six full-precision offsetpoly.f90 terms). reference_polys
matches all six analytic polynomials (including the microsecond-quantised
constants 4.12325 and -2.61835 and the stretch terms -0.27393/-0.314959);
fit_offset_polys reproduces the applied secondary polynomials over the whole
resampled swath grid (allowance 1e-4 px) with the azimuth mean shift pinned
to the cull file's own azmin rather than the stored value (round 13 finding
15): with the ISCE OffsetField fit on the cluster the reproduction is exact
(0 px on all six polynomials against the full-precision offsetpoly prints),
and the numpy replica agrees to <= 1.1e-5 px. That residual on the s2 azimuth
terms is offsetpoly.f90's own conditioning (singular-value ratio 7e-11), not
a transform difference; the s2 grid is the reference grid so its transform is
the identity. tests/selftest.py is the round-11 self-test on frame-3700-like
numbers (fixtures in a temp directory).

Round 12 (2026-09-02, Fable 5.1; 4 lenses, 8 findings, 6 confirmed 2-of-2,
1 split applied, 1 rejected) added:
11. Lossless overflow. The int16 quantum is P99.995(|z|)/32000, so by
   construction ~5e-5 of every swath's samples (the brightest, most coherent
   scatterers, up to 12-40x the percentile on real SAR) exceed the int16
   range. Clipping I and Q independently moved their phase toward an odd
   multiple of 45 deg (up to 70 deg interferometric error per pixel on the
   reproduction, cell errors to 158 deg after 5x28 looks on real bursts),
   and amplitude saturation, while phase-preserving, removes their weight
   from the multilooked cell, so neither is acceptable. pack_slc stores the
   overflow samples exactly (complex64 with int64 flat indices) in
   <out>.ovf, records overflow_count/overflow_file/pctl/peak/max_amp in the
   sidecar, and unpack_slc restores them, refusing a missing side file, a
   size that disagrees with the count, unsorted or out-of-range indices, or
   a side file whose pixels are not at the int16 limit in the pack. Packing
   refuses non-finite samples. Selftest: a speckle pair with 60 point
   targets at 5-50x the limit and different phases per date must come back
   exactly on the targets and match the float pair's 5x28 multilooked
   interferogram to < 1e-2 rad in every cell (< 1e-3 rad in target cells);
   the old clip failed this at 0.7 rad per pixel and 2.9 rad per cell.
12. The --cull provenance check is positional (.../<date>/f<i>_<frame>/
   s<n>/cull.off by path position, the upstream layout), not membership of
   the path components, which a per-date job cwd or a <date>-named lane
   directory above dates/ satisfied for a foreign cull.off.
13. render validates before writing: cropoff names map through ISCE_DTYPE
   (float32/float64/int8 were passed verbatim and renderVRT raised KeyError
   after the .xml was on disk), the scheme must be BIL/BIP/BSQ, and the
   binary size must equal width x length x bands x itemsize.
14. Golden replay dims are derived, not hardcoded: the resampled swath
   extents are round(N_ref * ali_s / ali_2) x round(N_ref_samples) from the
   181123 raw swaths (172356 x 4440, 172355 x 5670, 172355 x 5794, as
   resamp_slc printed for every date), asserted against the printed values;
   the earlier constants used the 181012 raw s3 width and the s2 line
   count. The reconstructed swath origins are asserted at 0/4100/9412 and
   the frame width at 15206.

Round 13 (2026-09-02, Fable 5.1 lenses, Opus 5 verification; 3 findings,
all 3 confirmed) added:
15. The golden replay no longer borrows the azimuth mean shift from the unit
   under test: it recomputes upstream's azmin (min y1 of the cull points, as
   Offset.getFitPolynomials does), asserts the stored meanAzimuth equals it,
   and evaluates the upstream polynomial with it. A sidecar recording
   meanAzimuth 0 while the fit was shifted displaced eval_poly by 125 lines
   on s3 and 151 on s4 yet printed GOLDEN_REPLAY_OK with unchanged numbers.
   The upstream polynomial now also uses the full-precision range and azimuth
   quadratic terms offsetpoly.f90 prints (cmd1.out 5409/5422, 6808/6821,
   8014/8027), so the grid allowance drops from ~1e-3 px to 1e-4 px: the ISCE
   fit then reproduces all six applied polynomials exactly (0 px) and the
   numpy replica to 1.1e-5 px.
16. The selftest's foreign-.ovf refusal reached only the size check (the
   copied side file had a different overflow count), so unpack_slc's int16
   limit ownership check was never exercised. The suite now covers three
   distinct refusals: missing, wrong-sized, and a right-sized side file whose
   sorted in-range indices point at pixels below the int16 limit, each
   asserted on the message it must produce. With the ownership check removed
   the third case is accepted, which is what it now catches.
17. render writes the BARE BASENAME with cwd = the binary's directory (as
   upstream create_xml does), not the path it was called with. Consumers that
   open the data through the xml follow the stored name: look_geom loads
   <date>_1rlks_14alks.los.xml and hands the image to mroipac Looks, which
   clones it and calls createImage() without resetting filename
   (Looks.py:47-50). Reproduced on the cluster ISCE: an xml rendered under the
   persist tree and copied to a node-local tree made Looks open the PERSIST
   .vrt and multilook the wrong data silently (and would fail outright if that
   tree were absent or renamed); with the basename it reads the staged copy.
   Selftest renders from another directory and asserts the xml file_name, the
   relative VRT SourceFilename, the sidecar files beside the binary, and that
   cwd is restored.

Round 14 (2026-09-02, Opus 5; 4 lenses, 8 findings, 1 confirmed 2-of-2,
2 split, 5 rejected) closed the loop. Applied:
18. pack_slc publishes the sidecar atomically (tmp + os.replace) like both
   binaries. A lane killed by the SLURM time limit mid-json.dump left two
   complete files beside a truncated sidecar, i.e. a directory that looks
   finished to a presence-based check while the 0.93 GB crop is unreadable
   (the sidecar is sec 4A's single source of truth). check_run.py now gates
   DONE on content, and tests/test_check_run.py breaks a good date eight ways
   to prove the gate catches each one.
Recorded, not fixed (neither breaks a run; the audit loop stops here by
decision, so they are logged for whoever picks the work up):
19. Amendment A.2's window validation is untested. The selftest fixture makes
   every adjustment a no-op, and deleting the azimuth clip, the range
   widening and the frame-edge clip from crop_pack leaves the whole suite
   green. On the real frame the azimuth clip IS load-bearing: the AOI reaches
   the frame's north edge, so without it lineOff stays 0 and derive aborts
   with "swath 2: azimuth window -42+59472 outside [0, 172356)". SEED's
   derive output is the live check (expect az_clipped=true, lineOff 56,
   lineOff_s 14/42/56).
20. rg_filter zero-pads outside the file it is given, so regenerating a
   subband from the CROPPED SLC leaves a filter transient of hnfilter = 128
   samples at every crop range edge that is not an original swath edge. No
   phase ramp restores it and sec 3's whole-crop mean cannot see it (beyond
   128 columns the reproduction agrees to 2.1e-13 rad). Decide before the ION
   blocks run: mask the 128-column bands in the ion estimate, or widen the
   ion-side crop by more than hnfilter.

## Implementation and storage (2026-09-03)
Stage scripts, each ending in a content validity check that gates its DONE
marker (check_run.py; tests/test_check_run.py breaks a good run eight ways to
prove the gate catches them):
  seed.sbatch     SEED, validated: bundle E/F/G + reference crops
  wave1.sbatch    one date per job; releases the date's CEOS only after the
                  check passes. Resamples ONE band: with the ion pass
                  regenerating subbands from the crops, -subband tripled the
                  most expensive step for output discarded with the node-local
                  tree (34 min -> about 12 min per date).
  ion.sbatch      ion block; regenerates subbands + the ramp (ion_regen.py)
                  because ion_subband.py reads them but its rg_filter call is
                  commented out upstream
  barrier.sbatch  ion QC (ion_qc.py: per-pair thresholds AND a network
                  connectivity gate, since ion_ls solves over that network),
                  compute_baseline, ion_ls -> dates_ion
  main.sbatch     main block; geocode.py directly with -bbox, output named for
                  build_stack2's find_pairs (which regrids, so no aoi90
                  template is needed)
  make_blocks.py  reproduces the certified layouts: ion 362 pairs / 5 blocks of
                  24 at stride 20, main 1650 pairs / 7 blocks of 32 at stride
                  12, every pair in exactly one block (asserted)
  build_dem3.py   the 3-arcsec geocoding DEM the stitcher failed to build
  dl_selective.py fetches only the CEOS members we keep: every member is
                  stored, so 16.3 GiB of each 56.6 GiB granule (28.8 percent).
                  NOT yet verified byte-identical; the whole-zip downloader
                  remains the default until it is.

STORAGE (measured 2026-09-03 by writing until it failed; no tool reports these
limits, `df` shows the whole 167/153 TB filesystems and the RPC quota service
is down). home ~425 GB, /mnt/scratch/t/tommy ~512 GB. Both were hit in one
night: home full killed SEED at the first pack, killed a wave-1 lane mid-
resample when its log write failed, and broke SSH key auth for 2.5 hours
because both hops read ~/.ssh/authorized_keys from it; scratch full then killed
the next pack. Plan: scratch = crop store (270 GB at 93 dates) + a 10-date CEOS
working set; home = the bulk CEOS buffer now, the ~165 GB of geocoded MAIN
output later; the download window is capped at 30 dates.

## Ionosphere result (2026-09-03, from the fleet's own ion pass)
filt_ion std, the per-pair uncertainty of the ionospheric phase estimate,
measured over 326 pairs: 0.0994 (2018), 0.0997 (2019), 0.0998 (2020), 0.0997
(2021), 0.1492 (2022), 0.4209 (2023), 0.1822 (2024) rad. Four years pinned at
0.0995 then a rise peaking in 2023, tracking solar cycle 25 from minimum toward
maximum. This is signal, not failure, and it sets the barrier's QC policy: a
fixed 0.35 rad cut, or a robust median+3*MAD cut (0.135, the MAD being 0.0067
because 2018-2021 dominate), would delete the recent dates for having an active
ionosphere, which is where the correction matters most. ion_ls reads each
pair's .std and weights the inversion by it, so noisy-but-valid pairs are
down-weighted rather than discarded. STD_MAX therefore defaults to 0.8 (nothing
observed exceeds 0.53) and the structural checks (non-finite, all-zero, flat
field) plus the network-connectivity gate carry the QC.

## Audit certification (2026-09-02)
Nine adversarial rounds (8, 12, 12, 8, 5, 4, 2, 1, 0 confirmed findings; 52
total, every one applied and re-verified) certified this spec plus
build_stack2.py, stack_smoke.sbatch, and the two environment patches clean
for wave-script implementation. Highlights the loop caught before any code
ran at scale: the 655 GB ion-barrier product trap, the 4-5 TB full-frame
geocode default, the subband demodulation-ramp physics, the finale interface
mismatch, the baseline-file column parse that silently zeroed every
perpendicular baseline, and the geometric-validation call that could never
have run.
