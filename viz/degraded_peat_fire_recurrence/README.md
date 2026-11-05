# Suspected degraded peat and fire recurrence, PeatGuard dam-site AOI

Screening products for the PeatGuard dam-site AOI (114.20 to 114.83 E, 2.41 to 3.01 S, ex-Mega Rice
Project, Central Kalimantan). The one message of the animation: fire keeps coming back to the same
drained, degraded peat. All numbers below are copied from `stats.json`, which `build.py` writes.

Confidence tags: [E] established (published dataset or finding), [S] standard inference (our own
screening statistics, methods after a named source), [H] hypothesis.

## Files

| File | What it is |
|---|---|
| `map_suspected_degraded_peat.png` | 3000 px wide, 200 dpi. Main map: suspected degraded peat (single red, #c8102e at alpha 0.28, dark outline), intact peat tan, non-peat hillshade grey, rivers, 60 dam sites (red if burning this week), 8 settlement labels, scale bar, north arrow. Three secondary panels show which criterion fires (canal-drained, burn-scarred, measured sinking). Canal lines are drawn only in panel (a). |
| `map_fire_recurrence_2001_2024.png` | Fire return count (years with fire, 0 to 24) per 500 m cell, sequential light orange to dark red, degraded-peat outline, sites and labels, colour bar, contrast numbers, and a bar chart of MODIS C6.1 detections per year with El Nino years (2002, 2006, 2009, 2015, 2019, 2023) marked. |
| `fire_years_vs_canal_distance.png` | Mean fire years per work-grid cell 2001 to 2024 against distance to the nearest mapped canal, bands 0 to 0.5, 0.5 to 1, 1 to 2, 2 to 3, 3 to 5, over 5 km, with n cells per band. |
| `fires_return_to_degraded_peat_2001_2026_SLIDE.mp4` | 1080 x 1080 map-only version for a slide: year (about 90 px) with El Nino tag, one headline line ("N percent of degraded peat burned in 3 or more years so far"), 0 / 5 / 10 km scale bar, four-swatch legend. The 0.63 deg wide box is trimmed by 0.015 deg on each lon edge to fill the square. |
| `fires_return_to_degraded_peat_2001_2026.mp4` | 1920 x 1080 backup with the side panel: title, year, headline, sensor and count lines, per-year bar chart, fire-years-vs-canal-distance inset, legend, sources. |
| `fires_return_to_degraded_peat_2001_2026.gif` | 720 x 720, from the square frames, 3 fps sampling (one sub-frame per year, so no flicker; this keeps the file small). |
| `stats.json` | Every number used: areas, class counts, recurrence contrasts (three grids), sensitivity to the peat cutoff, canal-distance bands, per-year detections and processing counts, per-site fire years, animation frame metadata. |
| `fires_viirs/viirs_YYYY.json` | Clipped, filtered, deduplicated VIIRS detections 2012 to 2024 (cache; delete to force a re-clip). |
| `frames/` | Rendered sub-frames (4 per year, panel and square versions) and ffmpeg concat lists. |
| `build.py` | The single script that produces everything. |

Raw NASA FIRMS Indonesia yearly CSVs are cached in `grid_products/external/firms_viirs/raw/` (about 450 MB).

## Definitions

- Peat: PeatGrids mean thickness >= 50 cm [E for the dataset: Widyastuti et al. 2022, Environ. Res. Lett.;
  S for the 50 cm cutoff, which follows the Indonesian regulatory definition of peatland]. Water
  (`water_mask.npy` and `rivers.npy`) is removed from every land class. PeatGrids is a 1 km global
  product (EPSG:4326, NaN nodata) resampled bilinearly to the 0.0006 deg work grid. Over this box the
  thickness is 55 to 628 cm (median 171 cm), so 100 percent of the land is peat under the 50 cm cutoff
  and the non-peat class is empty (0 km2). At a 150 cm cutoff 63 percent of the box would be peat; the
  sensitivity of the contrast to the cutoff is in `stats.json` (`recurrence_sensitivity_to_peat_cutoff_500m`).
- Suspected degraded peat [S, a screening layer, not ground-truthed]: peat AND at least one of
  (a) within 1 km of a mapped canal (binary dilation of `canals.npy`, the Dadap et al. 2021 canal map
  max-resampled to the work grid, by a disk of radius 15 cells, about 1 km);
  (b) burn-scarred: Landsat normalised burn ratio (Jul to Oct median) < 0.45 in at least 2 of the 11
  valid epochs 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2022, 2023, 2024, 2025, NaN ignored (2010 has
  37 percent gaps; 2026 is excluded because its Jul to Oct window was incomplete at build time; including
  it changes the degraded area by 2 km2). The 0.45 cutoff follows the NBR burn-severity scheme of Key and
  Benson 2006 (FIREMON Landscape Assessment) [S];
  (c) measured sinking: Tay et al. 2025 ALOS-2 L-band long-term rate (2018-10-12 to 2024-08-09) < -5 mm/yr
  where covered; the Tay product covers only 10.9 percent of the box.
  Intact peat = peat meeting no criterion. The three criterion masks are kept separately (panels in the
  static map; areas below).
- Fire return count [E for the detections, S for the gridding]: number of years 2001 to 2024 with at
  least one active-fire detection per 0.0045 deg (about 500 m) cell. 2001 to 2011: NASA FIRMS MODIS C6.1
  (1 km footprint); each detection marks every 500 m cell whose centre lies within 500 m of the
  detection. 2012 to 2024: NASA FIRMS VIIRS 375 m, S-NPP for all years plus NOAA-20 (JPSS-1) from 2018;
  each detection marks the 500 m cell containing it. VIIRS filter: `type == 0` (presumed vegetation
  fire), confidence nominal or high (low dropped), NOAA-20 detections dropped when an S-NPP detection
  lies within 375 m and 1 hour (14,147 of 16,372 NOAA-20 detections in the box were duplicates of S-NPP
  and were dropped; 2,225 kept). Two reference grids are also computed: a MODIS-only 500 m grid for all
  years, and the original MODIS-only 1 km grid (0.009 deg, 15 x 15 work cells).
- Class per fire cell: majority of degraded / intact / non-peat work cells by cell centre; cells with no
  land are excluded.

## Areas (work grid, 0.0006 deg, cell 0.00445 km2)

| Layer | km2 | share of box |
|---|---|---|
| Box | 4674 | 100 percent |
| Water (rivers, open water) | 148 | 3.2 percent |
| Peat (>= 50 cm) | 4526 | 96.8 percent |
| Suspected degraded peat (any criterion) | 4322 | 92.5 percent of the box, 95.5 percent of the peat |
| Intact peat | 204 | 4.4 percent |
| Non-peat | 0 | 0 |
| (a) within 1 km of a canal | 4289 | 91.8 percent |
| (b) burn-scarred | 1403 | 30.0 percent |
| (c) measured sinking (Tay 2025, 10.9 percent covered) | 94 | 2.0 percent |
| (a) only | 2869 | |
| all three | 44 | |

Plain statement: 92.5 percent of the box is classed as suspected degraded peat. That is why the whole
box reads red in the maps and the animation. The driver is criterion (a): the Dadap 2021 canal mask,
after max-resampling, puts 91.8 percent of the land within 1 km of a canal. The intact class is small
(204 km2, 804 cells of 500 m), so the degraded-versus-intact contrast rests on a small reference.

## Headline recurrence contrast

500 m cells, majority class, 2001 to 2024, MODIS to 2011 and VIIRS from 2012 (the headline grid):

| Class | cells | >= 1 fire year | >= 3 fire years | >= 5 fire years | mean fire years | max |
|---|---|---|---|---|---|---|
| Suspected degraded peat | 17,956 | 78.3 percent | 40.5 percent | 15.9 percent | 2.35 | 12 |
| Intact peat | 804 | 73.5 percent | 31.8 percent | 10.7 percent | 1.95 | 10 |
| Non-peat | 0 | n/a | n/a | n/a | n/a | n/a |

MODIS-only 1 km grid, 2001 to 2024 (sensor-consistent): degraded 4,503 cells, 78.5 / 43.4 / 17.7 percent,
mean 2.48; intact 187 cells, 72.2 / 32.1 / 14.4 percent, mean 2.07.
MODIS-only 500 m grid: degraded mean 2.10, intact mean 1.87; >= 3 years 35.4 versus 29.9 percent.

The contrast is real but modest: suspected degraded peat has about 0.4 more fire years per cell and a
9 point higher share of cells burned in 3 or more years than intact peat. Fire is widespread on both
classes (three quarters of intact cells also burned at least once). With the intact class this small and
the screening layer this broad, the product shows that nearly the whole box is drained peat that burns
repeatedly; it does not isolate a clean unburned control [S].

Running counter in the animation (share of degraded 500 m cells burned in 3 or more years so far):
0 percent through 2002, 8 percent by 2006, 14 percent by 2009, 32 percent by 2015, 38 percent by 2019,
41 percent by 2023 and at the final 2026 frame; burned at least once: 78 percent by 2024.

## Fire years against distance to the nearest canal (`fire_years_vs_canal_distance.png`)

Work-grid cells, water excluded, fire years from the 500 m mixed-sensor grid, distance from a
Euclidean distance transform of `canals.npy`:

| Distance band (km) | n cells | mean fire years | >= 1 | >= 3 |
|---|---|---|---|---|
| 0 to 0.5 | 869,870 | 2.40 | 80 percent | 41 percent |
| 0.5 to 1 | 92,063 | 2.52 | 78 percent | 44 percent |
| 1 to 2 | 40,703 | 2.56 | 86 percent | 45 percent |
| 2 to 3 | 11,721 | 1.06 | 54 percent | 13 percent |
| 3 to 5 | 2,358 | 0.36 | 27 percent | 5 percent |
| over 5 | 0 | n/a | | |

Within 2 km of a canal the recurrence is flat at about 2.4 to 2.6 fire years per cell; beyond 2 km it
falls to 1.1 and then 0.4. No land in the box is more than 3.6 km from a mapped canal, so the far bands
are small. The canal mask is generous (max-resampled), so distances are lower bounds and the 0 to 0.5
band absorbs 85 percent of the land. Association of fire with canals and drainage in this landscape:
Hoscilo et al. 2011, Int. J. Wildland Fire 20(4), 578-588 [E]; the band numbers above are our own
screening statistics [S].

## Sources

- PeatGrids peat thickness: Widyastuti et al. 2022, Environ. Res. Lett. [E]. File
  `grid_products/external/peatgrids/PEATGRIDS_THICKNESS_CM_MEAN.tif`.
- Canals: Dadap et al. 2021 canal map, resampled with max to the work grid (`canals.npy`) [E for the
  map, S for the 1 km buffer].
- Landsat NBR epochs (`nbr_YYYY.npy`), thresholds after Key and Benson 2006 [S].
- Subsidence: Tay et al. 2025, ALOS-2 L-band long-term rate (`tay_rate_mmyr.npy`) [E], used as an
  independent benchmark, not force-calibrated.
- Fires: NASA FIRMS MODIS C6.1 (`fires_YYYY.json`, 2001 to 2024, clipped to the AOI) [E]; NASA FIRMS
  VIIRS 375 m S-NPP (2012 to 2024) and NOAA-20 / JPSS-1 (2018 to 2024) Indonesia yearly archives from
  `https://firms.modaps.eosdis.nasa.gov/data/country/<viirs-snpp|viirs-jpss1>/YYYY/<name>_YYYY_Indonesia.csv`
  (the NOAA-20 folder is `viirs-jpss1`, listed in the FIRMS `yearly_summary_files.txt`) [E]. This
  week's cells: `fires_2026_week.json` (VIIRS S-NPP, NOAA-20 and MODIS, 1 to 8 Oct 2026, 0.05 deg bins).
- Canal-fire association: Hoscilo et al. 2011, Int. J. Wildland Fire [E].
- Hillshade and DEM: the AOI products in `deliverables/aoi_3d_degradation/`.

## Caveats the reader must know

- Active-fire detections are thermal anomalies. They include agricultural and land-clearing burns, and a
  detection does not measure burned area. The MODIS footprint is about 1 km, VIIRS about 375 m.
- 2001 to 2011 is MODIS only; 2012 to 2024 is VIIRS. VIIRS sees many more small fires (2019: 11,278
  VIIRS versus 2,132 MODIS detections in the box), so per-year cell counts and the fire return count are
  not sensor-consistent across 2011/2012. The MODIS-only grids in `stats.json` are the consistent series;
  the contrast is similar on all three grids.
- The 2026 frame mixes VIIRS and MODIS for one week in 0.05 deg bins; the points are placed at random
  inside each bin (seeded), so positions in that frame are approximate and it is not comparable with the
  yearly frames.
- The burn patches in the animation are glow fields for visualisation: detections rasterised with their
  footprint, jittered by 300 m (seeded), Gaussian-smoothed (sigma 0.9 km for MODIS years, 0.5 km for
  VIIRS years, times a seeded 0.8 to 1.3 factor per year), normalised to one isolated detection's peak,
  drawn from 0.75 detection equivalents upward in three bands, with the outer threshold raised whenever
  the drawn area would exceed 3 times the true footprint area. Per-year footprint and drawn areas are in
  `stats.json` (`animation.frames`). They are not burned-area polygons.
- The suspected-degraded layer is a screening layer. Criterion (a) alone classes 92 percent of the box, so
  the layer says little about which parts of the box are worse; panels (b) and (c) carry that detail.
- PeatGrids is a 1 km model-based product; the 50 cm cutoff leaves no non-peat land in this box.
- Dam sites: all 60 sit in cells classed degraded; fire years per site cell (500 m) are in `stats.json`
  (`site_fire_years_500m`, range 1 to 9).
- The GIF samples one sub-frame per year (no flicker) to stay at about 8 MB; the MP4s flicker at 12 fps
  with 4 seeded sub-frames per year. Each year is held 11 frames (0.92 s), the 2026 frame 36 frames
  (3 s), total 25 s.

## Regeneration

From the repo root:

```
.venv/bin/python3 deliverables/degraded_peat_fire_recurrence/build.py
```

Add `--skip-download` to use only the cached FIRMS CSVs. Requirements: rasterio, numpy, scipy,
pandas, scikit-image, matplotlib, Pillow in `.venv`, and ffmpeg at `/opt/homebrew/bin/ffmpeg`. The
first run downloads about 450 MB of FIRMS Indonesia CSVs (public, no key) and writes the clipped cache
in `fires_viirs/`. A full run takes about four minutes.
