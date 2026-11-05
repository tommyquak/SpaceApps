# Closure compensation and subsidence baseline plan

Status: plan approved-for-review, 2026-08-31. Nothing in the pipeline has been
changed; this document folds three research strands into one sequence.

Provenance: (1) the four-epoch NISAR stack state after the 2026-08-21 pull
(3 mission pairs, 60-day span, zero closure cycles, release gate unchanged at
relative-displacement commissioning); (2) a 15-agent literature review of
closure-cycle compensation (7 themes, 108 citations checked by adversarial
verifiers, 3 minor mismatches corrected, 0 fabricated sources); (3) the
sensor-fusion baseline memo of 2026-08-31, whose headline (free ALOS-2 ScanSAR
L1.1 at ASF) was independently re-verified by CMR query and an authenticated
8 MB ranged download (HTTP 206, valid CEOS L1.1 volume file, 56.6 GiB granule).

## 1. The verdict

The zero-closure problem is a limitation of the mission GUNW products, not of
the data PeatGuard holds. Every research theme converged on the same fact: we
hold the full 5 m dual-pol GSLCs for all four epochs, so the three missing
non-nearest pairs (0622-0716, 0704-0821, 0622-0821) can be formed locally,
turning the 3-pair chain into the complete 4-node graph: 6 pairs, 3 independent
closure triangles. That enables published closure-based unwrap screening
(Yunjun et al. 2019) and a first-order fading-signal bias diagnostic (Zheng et
al. 2022, the same correction family Tay et al. 2025 used) immediately [E].

The literature is unanimous that a nearest-neighbour chain is the worst-case
network for fading-signal bias over moisture-dynamic ground like peat: the bias
lives in short multilooked pairs, decays with temporal baseline, and accumulates
through a chain sum without cancellation (Ansari, De Zan, Parizzi 2021;
Maghsoudi et al. 2022) [E]. The long pairs are therefore low-bias anchors of
the inversion, not optional extras.

What no published method rescues: an annualized velocity from 60 days of data,
or a converged bias correction from 3 triplets. Published tropical peat
shrink-swell is roughly 0.5 to 6 cm over a seasonal cycle versus ~0.33 cm of
secular subsidence in 60 days at 2 cm/yr, so the 60-day result stays "bounded
surface motion with stated uncertainty", never a rate [E]. The multi-year
redundant network that closure correction actually needs comes from the second
strand: the ALOS-2 ScanSAR L1.1 archive, free at ASF, 267 granules over the AOI,
2014-2025, the same sensor/mode/level as the Tay benchmark.

## 2. NISAR track: the compensation stack (all free, laptop-scale)

- N0. Complete-graph reconstruction. Form the three non-nearest fine
  interferograms from the staged GSLCs (existing script
  `nisar_gslc_interferogram.py`), compute the three multilooked closure
  triangles, run closure-based unwrap screening and the Zheng-style bias
  diagnostic as a mask/flag, not a correction [E]. Effort: days; scripts exist
  for pair formation and closure validation.
- N1. Phase linking on the full 4x4 coherence matrix as the primary estimator:
  EMI (Ansari, De Zan, Bamler 2018) or EVD fallback (CAESAR, Fornaro et al.
  2015), rooted in Guarnieri and Tebaldini 2008. This consumes all 6 pairs with
  variance-optimal weighting, bridging the 0.28-coherence 36-day fire-season
  pair through the two intermediate epochs, and yields a per-pixel
  temporal-coherence quality metric (SqueeSAR tradition, Ferretti et al. 2011)
  [E]. Software route: dolphin (Staniewicz et al. 2024, JOSS, the algorithmic
  basis of the future DISP-NISAR product) accepts geocoded SLC stacks; MiaplPy
  documents radar-coordinate input only, so treat GSLC ingest there as an
  unpublished adaptation. Caveat stated plainly: no published demonstration
  exists at N=4; frame the result as the seed mini-stack of a growing archive,
  not a defensible velocity [S].
- N2. Polarimetric lift (critic addition): our GSLCs are DHDH (HH+HV), so
  polarimetric coherence optimization (Navarro-Sanchez, Lopez-Sanchez,
  Vicente-Guijalba 2010) can raise coherence of the 36-day pair before phase
  linking at zero data cost [E method, S applicability here]. A delta-k
  split-bandwidth check (Bamler and Eineder 2005) gives a closure-free
  unwrap-error flag per pair [E/S].
- N3. External corrections per pair, replacing redundancy with physics:
  mission ionosphere screen or split-spectrum (equatorial L-band validation
  3-6 mm/yr, Marconato et al. 2024), solid Earth tide and ocean tide loading
  (deterministic), datum tied to stable non-peat reference plus the CPKY GNSS
  station, exactly Tay's published practice. Weather-model troposphere is NOT
  applied: Tay et al. 2025 found PyAPS and GACOS worsened accuracy over this
  same flat lowland; the negative result is citable [E].
- N4. Multi-geometry redundancy: run the descending 133/090 and 061/089 chains
  and the ascending 024/176 frames through the same MintPy chain independently,
  project to vertical, and report cross-geometry residual histograms as the
  reliability statistic (Wright et al. 2004; Fuhrmann and Garthwaite 2019;
  overlap-residual precedents: Weiss et al. 2020, 3.1-3.7 mm/yr overlap std;
  Li et al. 2021, 90 percent of overlap pixels within 9 mm; NISAR's own Solid
  Earth requirements are written as asc+desc vector components) [E]. Several
  of these GUNWs are already on disk.
- N5. Honest uncertainty machinery: chain-propagated per-epoch sigma from
  coherence (Rodriguez and Martin 1992 CRLB, MintPy variance weighting) plus an
  atmospheric variance floor (Emardson et al. 2003); adopt KFTS (Dalaison and
  Jolivet 2020) as the forward engine that updates cheaply per new 12-day
  epoch, with process noise grown across the 36-day gap; ISBAS-style per-pixel
  subsets (Sowter et al. 2013; tropical peat precedent Marshall et al. 2018)
  preserve coverage where the 36-day pair decorrelates, flagged per pixel [E].
- N6. Peat-specific defensibility, the way published peat studies argued
  reliability without closure correction (Hoyt et al. 2020; Zahro and Izumi
  2025; Alshammari et al. ISBAS tradition): empirical noise floor from stable
  non-peat control areas, GNSS-checked reference (CPKY sits inside the full
  GUNW footprint), rainfall plausibility checks [E].

Release-gate effect: N0-N2 make closure diagnostics and a redundant inversion
computable, so the audit's "redundant_network_ready" can flip for the custom
GSLC branch; "annual_velocity_ready" stays false on epoch count, by design.

## 3. ALOS-2 track: the baseline (data free, compute rented)

- A1. Pilot: ~20 descending ScanSAR L1.1 scenes through ISCE2 alosStack (Tay's
  app; full-aperture ScanSAR InSAR per Liang and Fielding 2017) plus MintPy
  weighted SBAS with Tay's parameters. Requires a rented VM (the old sbas VM is
  gone; its IP now presents a foreign host key): ~16 vCPU / 30 GB RAM / 2-3 TB
  volume, ~US$50-100 per month for 1-2 months. Granules are 30-60 GiB each;
  download on the VM, keep zips off the working volume (topsStack lesson).
- A2. Full baseline: 2018-2025 (extendable to 2014), both geometries, with the
  converged closure-bias correction that a 90+ scene archive supports (Tay's
  bandwidth-20 network; bias worth up to 30 cm false cumulative subsidence
  over 5.5 years on scrub/forest per the peat-specific theme) [E].
- A3. Validation: compare, never calibrate, against Tay's public DR-NTU product
  and Zahro/Izumi 2025 (1.74 cm/yr 2007-09 vs 0.29 cm/yr 2017-23 over PHU
  Kahayan-Sebangau). Comparison harness already exists from the Sentinel-1
  falsification work.
- A4. Absolute anchors (Arm C of the fusion memo): request CPKY RINEX from BIG
  (one email); Vernimmen et al. 2019 open lidar DTM of the ex-MRP; ICESat-2
  ground photons over open degraded peat (Michaelides et al. 2021 precedent,
  ~1.7 cm crossover uncertainty). Critic addition for the field roadmap:
  subsidence poles (Hooijer et al. 2012 practice) are village-maintainable and
  feed the TRL 5-6 ground-truth stage; active L-band transponders (Mahapatra
  et al. 2014) are the published hardware route if absolute tie points are
  ever funded (passive L-band corner reflectors are impractically large).
- A5. Counterfactual model (Arm B): drainage-age / peat-depth / water-table
  regression trained on A2 (Hooijer 2012; Umarhadi 2022 drainage-age decay;
  Ma et al. 2022 water-table relationship), projecting business-as-usual into
  the NISAR era. NISAR-measured minus modeled counterfactual is the MRV
  quantity. The ASF archive currently ends 2025-03-27, so A5 and N4 jointly
  span the ~15-month bridge to the NISAR era.

## 4. How the tracks join

The ALOS-2 baseline hands NISAR a same-band, same-ground pre-restoration
trajectory; NISAR replicates the Tay method forward from 2026 with closure
handled per Section 2 until the archive matures. Agreement between our A2
product and Tay's (compared, not calibrated) upgrades the TRL 4 claim from
"cross-validated against independent satellite data" to "independently
reproduced the benchmark method on open data". Restoration impact is then
measured as NISAR trajectory minus the A5 counterfactual, with A4 anchors
fixing the absolute datum.

## 5. Sequence and spend

| Step | What | Cost | Wall clock |
|---|---|---|---|
| N0-N2 | Custom pairs, closure diagnostics, phase linking, HV lift | $0 | ~1 week |
| N3-N6 | Corrections, multi-geometry, uncertainty, noise floor | $0 | ~1-2 weeks |
| A1 | 20-scene ScanSAR pilot on rented VM | ~$50-100/mo | 2-6 weeks |
| A4 | BIG email, DTM + ICESat-2 downloads | $0 | background |
| A2-A3 | Full 2018-2025 baseline + Tay comparison | VM months | ~1-2 months |
| A5 | Counterfactual model | $0 | after A2 |

## 6. Integrity rules carried forward

No velocity claim before a multi-season NISAR record; no force calibration to
Tay at any stage; PROVISIONAL labels retained on all NISAR-era products; every
released layer states posting versus effective support; the 60-day NISAR
result is bounded motion, not a rate; Tay's product license is CC BY-NC, so
commercial use of the comparison requires permission.

## 7. Citation corrections from adversarial verification

- Yunjun, Fattahi, Amelung 2019: correct subtitle is "Unwrapping error
  correction and noise reduction" (Computers and Geosciences 133, 104331).
- Li et al. 2021 (adjacent-track InSAR): last author is Zechang Lu, not Zhong Lu.
- Staniewicz et al. 2024 (dolphin): JOSS 9(103), 6997; description verified.
