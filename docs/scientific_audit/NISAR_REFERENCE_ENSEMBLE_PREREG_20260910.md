# Pre-registration: NISAR peat motion relative to a mineral-ground reference ensemble

Written 2026-09-10 before any displacement value was inspected for this design. Tags: [E] documented
or measured, [S] inference, [H] hypothesis. This is a proposed route to test, not a demonstrated
capability. The annulus diagnostic in
[NO_EXTERNAL_REFERENCE_ROUTE_20260910.md](NO_EXTERNAL_REFERENCE_ROUTE_20260910.md) is kept as a
comparison, not as the product.

## Question

Over the retained native NISAR GUNW chains (June to September 2026), does peat in the Sebangau and
ex-Mega Rice Project area move relative to independently selected mineral ground, with the regional
signal preserved (no annulus, no high-pass)? The potential output is "LOS change of peat relative to
the mineral-ground reference ensemble, assuming that ensemble is stable over the window", projected
to vertical under a stated no-horizontal-motion assumption. Short-window lowering is not an annual
trend, irreversible peat loss or a restoration effect; seasonal rebound is documented by Tay et al.
(2025) [E].

## Inputs (all already on the cluster; no download)

- GUNW P05023, frequency A, HH, 40 MHz, 80 m native grids. Descending track 61 frame 90 (UTM 49S):
  2026-06-30 to 07-12, 07-12 to 07-24, 07-24 to 08-17, 08-17 to 08-29 (four pairs, 60 days).
  Ascending track 125 frame 176 (UTM 50S): 07-04 to 07-16, 07-16 to 08-21, 08-21 to 09-02
  (three pairs, 60 days). Both frames span lon 112.3 to 115.1, lat -4.57 to -1.60 in common.
  Pass times from the product: ascending 22:07 UTC = 05:07 WIB, descending 10:45 UTC = 17:45 WIB.
- Corrections: [E] the product's unwrapped phase still contains the ionosphere (20 km low-pass of
  the phase correlates +0.90 with the ionospherePhaseScreen on the 07-16 to 08-21 ascending pair,
  slope 0.64; subtracting it lowers the low-pass std from 1.31 to 0.89 rad), and the ASF GUNW page
  states that corrections are not removed from the data. The metadata flags
  `*CorrectionApplied = True` are therefore not taken as evidence of application. We subtract the
  ionosphere screen and the hydrostatic, wet and solid-Earth-tide phases from the metadata cube,
  interpolated to each pixel's DEM height. Variants kept: raw, ionosphere only, all corrections.
- Evidence layers for reference selection: Copernicus 3 arcsec DEM (`~/pilot/dem`), PeatGrids
  thickness mean and standard deviation, ESA WorldCover 2021 (tiles S03E111, S03E114, lat -3 to 0
  only), Macrostrat geologic map unit at each candidate centroid (API, recorded verbatim).
- Targets: the 60 dam sites (500 m radius), PeatGrids peat (thickness at least 50 cm) inside
  lon 113.7 to 115.3, lat -3.45 to -2.2, and the Tay 2025 raster footprint (values not used).

## Reference ensemble selection (before any phase is read)

Evidence rules applied on a 200 m grid over the common frame area:
1. PeatGrids mean thickness 0 cm over at least 95 percent of the region and its 1 km buffer.
2. DEM median height at least 30 m, or at least 20 m above the minimum within 5 km (mineral
   terrace or hill; peat domes here do not reach that relief [S]).
3. WorldCover class not water (80), herbaceous wetland (90) or mangrove (95) where tiles exist;
   regions outside tile coverage are allowed and flagged `worldcover: unavailable`.
4. Contiguous qualifying area at least 4 km2; candidates at least 15 km apart; at most eight,
   chosen greedily by area with a spread rule (the next candidate maximises its minimum distance
   to those already chosen). Palangkaraya's built-up area is included only if it satisfies rules
   1 to 3 like any other candidate.
5. Geology: Macrostrat unit, lithology and age recorded per centroid. Pre-Cenozoic units count as
   strong mineral evidence; "Cenozoic sedimentary" units rely on rules 1 to 3.
6. One candidate is held out from the ensemble as a stability control, selected by a fixed seeded
   draw (seed 7) among the qualifying candidates, before any phase is read. The list is frozen with
   a SHA-256 of the candidate file.

High coherence is not a selection criterion. Connectivity is checked afterwards per pair: a
reference is usable for a pair only if at least 50 percent of its pixels are valid (mask decode as
in `nisar_local_contrast.py`) and share the target's dominant connected component after 4-connected
relabelling within valid samples.

### Amendment A1 (2026-09-10, after the evidence-only selection ran, before any phase was read)

The first selection returned a 7,258 km2 qualifying upland as one candidate. A reference must be
compact so that the ensemble spread measures atmosphere and reference motion at datum scale, so a
maximum area is added: a qualifying region larger than 25 km2 is reduced to the 5 km by 5 km block
inside it that contains the most qualifying cells (ties broken by higher median height). This is an
evidence-only rule; no displacement, coherence or phase was consulted.

### Amendment A2 (2026-09-10, before any phase result was inspected)

Natural Earth 10 m land polygons show candidates 2 and 4 of the first frozen list 0.7 km and
0.1 km outside the coastline: offshore blobs in the Copernicus DEM, not hills. Rule added: a
candidate must lie on land at least 1 km inside the coastline (land polygons rasterised to the
200 m grid and eroded by 5 cells). Selection rerun; the earlier list (SHA-256 4bea9d20...) is
superseded. The two ensemble jobs started on the earlier list were cancelled unread.

### Superseded first selection (frozen 2026-09-10, candidates SHA-256 4bea9d20c237761a..., before any phase read)

| Rank | Lon, lat | Area km2 | DEM median m | Prominence m | WorldCover | Macrostrat unit | Role |
|---|---|---:|---:|---:|---|---|---|
| 1 | 112.382, -1.900 | 24.9 | 106 | 41 | tree cover/shrub | Mesozoic volcanic and intrusive rocks Mesozoic to Mesozoic | ensemble |
| 2 | 114.890, -4.071 | 5.0 | 51 | 51 | unavailable | no map unit at point  | ensemble |
| 3 | 115.069, -1.664 | 24.9 | 103 | 25 | tree cover/shrub | Cenozoic sedimentary rocks Cenozoic to Cenozoic | ensemble |
| 4 | 112.854, -3.272 | 5.2 | 46 | 46 | unavailable | no map unit at point  | ensemble |
| 5 | 113.731, -1.684 | 24.9 | 98 | 41 | tree cover/shrub | Cenozoic sedimentary rocks Cenozoic to Cenozoic | ensemble |
| 6 | 113.997, -3.354 | 4.4 | 46 | 46 | unavailable | Cenozoic sedimentary rocks Cenozoic to Cenozoic | ensemble |
| 7 | 114.373, -1.640 | 4.5 | 99 | 40 | tree cover/shrub | Cenozoic sedimentary rocks Cenozoic to Cenozoic | ensemble |
| 8 | 114.626, -3.708 | 17.7 | 47 | 47 | unavailable | Mesozoic sedimentary and volcanic rocks Middle Triassic to Early Cretaceous | holdout_control |

Nine large regions were reduced to 5 km blocks under amendment A1. Candidates 2 and 4 have no
Macrostrat map unit at the point (possibly offshore or unmapped); their usability is decided by the
per-pair validity screen, not by hand.

### Current selection (frozen 2026-09-10, candidates SHA-256 aaad8ea69e7727b2..., all on land per Natural Earth, before any phase read)

| Rank | Lon, lat | Area km2 | DEM median m | Prominence m | WorldCover | Macrostrat unit | Role |
|---|---|---:|---:|---:|---|---|---|
| 1 | 112.382, -1.900 | 24.9 | 106 | 41 | tree cover/shrub | Mesozoic volcanic and intrusive rocks (Mesozoic to Mesozoic) | ensemble |
| 2 | 115.022, -3.687 | 24.9 | 841 | 720 | unavailable (lat < -3) | Mesozoic sedimentary and volcanic rocks (Middle Triassic to Early Cretaceous) | ensemble |
| 3 | 114.605, -1.651 | 23.9 | 102 | 32 | tree cover/shrub | Cenozoic sedimentary rocks (Cenozoic to Cenozoic) | ensemble |
| 4 | 113.522, -1.985 | 9.0 | 68 | 11 | tree cover/shrub | Cenozoic sedimentary rocks (Cenozoic to Cenozoic) | ensemble |
| 5 | 115.077, -2.070 | 16.1 | 75 | 16 | tree cover/shrub | Cenozoic sedimentary rocks (Cenozoic to Cenozoic) | ensemble |
| 6 | 114.036, -1.652 | 24.9 | 113 | 26 | tree cover/shrub | Cenozoic sedimentary rocks (Cenozoic to Cenozoic) | ensemble |
| 7 | 115.069, -1.664 | 24.9 | 103 | 25 | tree cover/shrub | Cenozoic sedimentary rocks (Cenozoic to Cenozoic) | ensemble |
| 8 | 114.629, -3.698 | 14.5 | 46 | 46 | unavailable (lat < -3) | Mesozoic sedimentary and volcanic rocks (Middle Triassic to Early Cretaceous) | holdout_control |

Two candidates sit on pre-Cenozoic units (Mesozoic volcanic and sedimentary rocks, the western
upland and the Meratus foothills); the holdout control is also on a Mesozoic unit. The five
northern candidates are on the Cenozoic sedimentary map unit and rely on the relief and
no-peat rules. WorldCover is unavailable south of 3 degrees S.

## Computation

For each pair and each usable reference r: phase relative to r = corrected phase minus the mean
corrected phase over r. Target statistics per pair: median over the target's valid pixels. Chain
along the nearest-neighbour tree (earliest epoch zero). Ensemble: median across references; spread =
interquartile range and max minus min across references; leave-one-reference-out. LOS to vertical
by 1/cos(incidence) from the metadata cube, stating the no-horizontal-motion assumption [S].
Sign: the product's phase sign is checked against the wrapped layer; LOS lengthening (motion away
from the sensor) is reported as negative vertical.

## Pre-registered reading

- Ensemble usable: at the final epoch the inter-reference IQR is below 10 mm LOS and no reference
  differs from the ensemble median by more than 15 mm. Otherwise: "reference ensemble unstable or
  contaminated; peat sign unresolved".
- Stability control: the held-out mineral region, treated as a target, lies within IQR + 5 mm of
  zero. Otherwise the stable-mineral assumption fails for this window.
- Detection: the peat median (Tay footprint and PeatGrids peat separately) relative to the ensemble
  median at the final epoch has the same sign in every leave-one-out subset, magnitude above twice
  the IQR and above 10 mm LOS. Then: "relative lowering (or rise) of peat against mineral ground
  over the window, conditional on reference stability". Per site: same rule with the site's own
  support; sites failing connectivity are reported unavailable, not filled.
- Cross-track: the descending and ascending peat medians over their overlapping interval (07-04 to
  08-29 approximately) must agree in sign; disagreement leaves the sign unresolved.
- A 60-day result is never annualised. No value is compared with Tay's rate for tuning.

## Failure modes named in advance

Residual troposphere after the model correction (wet delay over swamp forest is not in a
weather model at 80 m), reference regions in a different connected component from the peat,
unwrapping errors inside a component, ionosphere screen errors where it was interpolated, and
genuine mineral-ground motion (fluvial terraces can settle). The ensemble spread and the held-out
control are the only internal checks; neither proves ground accuracy.

## Result 1 (jobs 835844 D, 835845 A, 2026-09-10; list SHA-256 aaad8ea6...)

[E] Read against the gates above, unchanged:
- Descending 61/90: six of eight candidates had zero valid pixels in every pair (they lie outside
  the swath: the frames are skewed parallelograms and the candidate area was defined from the
  bounding-box intersection, a selection error). The remaining candidate (rank 4, terrace NW of
  Palangkaraya) shared the peat AOI's dominant component in all four pairs (1.5 million shared
  pixels) but was only 45 and 47 percent valid in the first two pairs, below the 50 percent rule,
  so no chain exists. 13 of 60 dam sites shared a component with it. Ensemble: zero usable
  references. Reading: unresolved.
- Ascending 125/176: only the Meratus block (rank 2) and the holdout coastal hill (rank 8), 40 km
  apart, were valid (81 to 99 percent). They share a component with 241,551 peat pixels in the
  eastern AOI, with 520 Tay-footprint pixels and with none of the 60 dam sites. Ensemble of one:
  IQR undefined. Holdout relative to the Meratus reference at the final epoch: raw -24.8 mm,
  ionosphere-corrected -19.9 mm, all corrections -35.3 mm LOS. The stability gate (IQR + 5 mm)
  fails by a wide margin: two mineral hills 40 km apart differ by 20 to 35 mm over 60 days.
  Connected peat relative to the same reference: -13.1 / -12.6 / -21.2 mm. Not interpretable while
  the control fails.
- Corrections: subtracting the ionosphere screen cut the 20 km low-pass phase std from 1.1 to 28.4
  rad down to 0.5 to 3.2 rad on every pair. The model troposphere plus tide layers lowered it
  further on two pairs and raised it on three; over 40 km between the two hills they widened the
  mineral-mineral difference from 20 to 35 mm.

[S] Reading: the residual long-wavelength error after the product's own corrections is 10 to
30 mm LOS over tens of km per 60-day chain, the same size as the expected dry-season peat signal.
A mineral reference therefore has to be within a few km of the peat target, or the ensemble must
bracket the target closely enough that residual gradients cancel. Distant uplands cannot serve as
the datum for this window. The design is not refuted; this candidate geometry is.

### Amendment A3 (before the second pass)
1. Candidate area = pixels valid in the first pair of both tracks (product mask, coverage only),
   eroded 1 km, instead of the bounding-box intersection.
2. Add a built-up candidate class: WorldCover built-up (50) regions of at least 1 km2 with
   PeatGrids 0 cm or with an independent stability record (the Palangkaraya InaCORS-checked point of
   Zahro et al. 2025, -2.2135, 113.9204, within 3 km), any relief. Tay's own datum (CPKY) lies in
   Palangkaraya; this class is the closest available analogue and is evidence-based.
3. Maximum candidates raised to 12 and minimum separation lowered to 10 km so that references
   can bracket the peat; the greedy spread rule now maximises the minimum distance to the peat
   AOI centroid last, not first (nearest candidates first).
4. Gates unchanged.

### Amendment A4 (before the second pass is read)
Selection inside the common footprint of both tracks leaves four candidates, all west or north-west
of the peat (Palangkaraya built-up terrace at 90 km from the AOI centroid, three western terraces at
150 to 214 km): the peat plain has no mineral ground closer, and the eastern uplands lie outside the
descending swath. Two additions, both evidence-only:
1. Per-track candidate lists (each track's own footprint) so the ascending ensemble can bracket the
   peat from the east (Meratus block) and south (coastal hill) as well as the north-west.
2. A mineral-plane datum variant: per pair, a plane fitted through the usable ensemble references'
   mean phases at their centroids (at least three), with targets read relative to that surface. Only
   mineral ground constrains the surface; the peat's regional signal is not filtered. Reported beside
   the constant-datum ensemble median. Gate for this variant: the held-out mineral control, read
   against the plane, lies within 10 mm LOS of zero at the final epoch; otherwise the residual is not
   planar at this scale and the variant is not usable.

## Result 2 (second pass, jobs 835873/835874 common list, 835875 ascending list, 835876 descending list)

[E] Gates applied as written. LOS mm, positive toward the sensor; 60-day chains; "ion" = ionosphere
subtracted; "all" = ionosphere + model troposphere + solid-Earth tide subtracted.

Constant-datum ensemble (median of references):
- Ascending, six-reference list bracketing the peat: reference IQR 10.3 mm and range 81 mm (ion),
  29 and 43 mm (all): gate FAILS. Descending, four western references: IQR 79 to 106 mm, range 174
  to 258 mm: FAILS. Common list on descending: IQR 81 to 92 mm: FAILS. Common list on ascending,
  "all": IQR 3.2 mm, range 6.5 mm, holdout -6.6 mm: gates pass, but only two references were usable
  and only Palangkaraya shares a component with the peat, so the target value (Tay footprint
  -19.3 mm) is a one-reference reading, not an ensemble result.
- Sites: no dam site shares a component with any reference in every pair on ascending; 13 of 60 on
  descending (Palangkaraya only).

Mineral-plane datum (plane through the usable references per pair):
- Ascending list (6, 5, 5 usable references per pair): holdout control read against the plane:
  raw +7.4, ion +2.5, all -11.9 mm at the final epoch (ion series 0, +7.8, +10.9, +2.5). Gate
  passes for raw and ion, fails for all.
- Under ion: Tay footprint 0, -0.6, -7.0, -27.2 mm at 4 Jul, 16 Jul, 21 Aug, 2 Sep; whole-AOI
  peat 0, +0.5, -1.0, -7.3 mm; 35 of 60 dam sites carry a full chain, final median -42.6 mm
  (p25 -66.2, p75 -19.5), 32 of 35 negative. Under all: Tay footprint -22.6, peat +0.6, sites
  median -20.5, 30 of 35 negative. At 36.7 deg incidence, -27 mm LOS is -34 mm vertical under the
  no-horizontal-motion assumption.
- Descending list (all references west of the peat): holdout -13.1 (ion), +17.6 (all): gate FAILS.
  The descending pairs also carry 3 to 6 times the ascending long-wavelength residual after the
  ionosphere subtraction (2.3 to 3.2 rad versus 0.5 to 0.9 rad at 20 km on the first two pairs),
  consistent with its 17:45 WIB dusk pass. Not usable for this test at present.
- Cross-track gate: cannot be evaluated; the descending run does not contradict, it is unusable.

[S] Reading. The pre-registered detection (usable constant-datum ensemble, sign agreement across
tracks) is NOT met. What survives is narrower: on the dawn ascending track, with a planar datum
constrained by five to six mineral references that bracket the peat and predicted a held-out
mineral region to within 2.5 mm at the end of the chain (10.9 mm mid-chain), the northern peat
(Tay footprint) and most connected dam sites lowered relative to mineral ground between 4 July and
2 September 2026 by roughly 20 to 45 mm LOS. Three cautions: one held-out control is a single
sample of the noise floor; 20 of the 27 mm at the Tay footprint come from the last 12-day pair
(21 Aug to 2 Sep), so a pair-specific residual cannot be excluded until the next pairs arrive; and
this is a dry-season, two-month lowering, which Tay et al. (2025) show can rebound. It is not an
annual rate, not irreversible loss and not a restoration effect. The wording that fits the
evidence is "lowering relative to the mineral-plane datum, assuming the bracketing mineral ground
is stable", with the descending track still to be brought to a usable state.

What would strengthen or overturn it: the 14 September and later ascending pairs (persistence or
rebound of the last-pair drop); a second held-out mineral region on the eastern side; the
descending track after a better ionosphere treatment or on quieter dates; and, for absolute
context, any GNSS in Palangkaraya (CPKY or InaCORS) over the same window.

### Result 2, corrected reading after inspecting the fields (same day, before any further run)

[E] The cumulative ascending plane-datum field (fig_ascending_plane_datum_cumulative.png) is a
smooth swing of about 100 mm across the AOI: -40 to -60 mm in the north and east, +20 to +40 mm
over the Sebangau peat in the south-west. The per-pair fields (fig_ascending_plane_datum_per_pair.png)
show the pattern reversing sign between pairs: east-half minus west-half means of -30, +12 and
-27 mm for 4 to 16 July, 16 July to 21 August and 21 August to 2 September. Per-pair increments
at the 35 connected dam sites (ion): -13.5, +1.4, -30.5 mm; at the held-out mineral region:
+7.8, +3.2, -8.4 mm.

[S] A pattern that is smooth over 100 km, reverses sign from one pair to the next and returns
with the same shape is residual long-wavelength propagation error (ionosphere, troposphere), not
ground motion: motion accumulates in the same place with the same sign. The mineral-plane datum
fitted through references 90 to 220 km away cannot represent this curvature inside the AOI, and
the single held-out region passed its gate because it happens to lie where the residual is small.
The figures of -27 mm (Tay footprint) and -43 mm (site median) in Result 2 are therefore
positions on a residual gradient, and the interim statement that the northern peat and most
connected sites "lowered relative to mineral ground" is withdrawn. No peat motion is detected or
excluded by this pass.

What this establishes [S]: per-pair residual amplitude after the product's corrections is
+-20 to 40 mm over the AOI on the quiet dawn track (worse on the dusk track). A 60-day chain of
three to four pairs cannot beat that. Two routes remain, neither of which is a short-window
result: (a) accumulate the 12-day series over a year or more so that the sign-changing residual
averages while motion accumulates (rate precision of order 10 mm/yr relative to mineral ground
if residuals are independent between pairs [S]); (b) estimate the per-pair residual surface
inside the AOI from something other than peat, which the geometry does not provide (no mineral
ground within 60 km). GNSS in Palangkaraya would give the absolute anchor but not the in-AOI
residual. The annulus diagnostic remains the only reference-free local observable, and it is
not a subsidence product.
