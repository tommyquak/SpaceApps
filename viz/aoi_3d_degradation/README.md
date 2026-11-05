# PeatGuard AOI 3D degradation animation, 1990 to 2026

Off-screen 3D flyover of the dam-site AOI in the ex-Mega Rice Project / Sebangau
area, Central Kalimantan (lon 114.20 to 114.83 E, lat -3.01 to -2.41, WGS84).
Twelve epochs (1990, 1995, 2000, 2005, 2010, 2015, 2020, 2022, 2023, 2024, 2025,
2026) are rendered as a continuous slow camera orbit over a terrain surface draped
with a Landsat burn-index class map, canals, rivers, fire detections and, from 2018
onward, a modelled surface lowering derived from the Tay et al. 2025 ALOS-2 rate.

## Deliverables

| File | Content |
| --- | --- |
| `PeatGuard_AOI_3D_degradation_1990_2026.mp4` | 1920x1080, 12 fps, H.264 High profile level 4.0, yuv420p, CRF 18, faststart. 654 frames, 54.5 s, 53.6 MB. Plays in PowerPoint on Mac and Windows. |
| `PeatGuard_AOI_3D_degradation_1990_2026.gif` | 960x540, 5 fps (every second frame of the 10 fps target; 10 fps came out at 51.3 MB, over the 40 MB cap), palettegen/paletteuse, 273 frames, 54.6 s, 33.3 MB. |
| `PeatGuard_AOI_3D_2026_still.png` | Last frame (`frames/2026_035.png`), 1920x1080. |
| `contact_epochs.png` | 3x4 contact sheet, frame 018 of each epoch. |
| `frames/YYYY_NNN.png` | 36 frames per epoch, 432 in total. `frames/list.txt` is the ffmpeg concat list. |

Timing: each epoch plays its 36 frames at 12 fps (3.0 s), then holds its last frame
for 1.5 s so the year and caption can be read; the final 2026 frame holds for 2.0 s.
12 x 3.0 + 11 x 1.5 + 2.0 = 54.5 s. The camera parameter t runs 0 to 1 across the whole
sequence with the end of each epoch segment exclusive, so consecutive epochs chain
without a duplicated or skipped camera pose (step between epochs = step within an
epoch = 0.0023).

## What the animation shows

- Terrain: Copernicus DEM relief (WGS84 ellipsoid heights, minimum removed, relative
  relief only), vertically exaggerated and lit with a Horn hillshade. Nothing in the
  scene is a hydrological gradient: this DEM cannot resolve the peat-dome drainage
  direction (gate-tested earlier in the project, see the memory note
  `project_dem_peat_gradient_unusable`). Caveat: the DEM is a 2011 to 2015 TanDEM-X
  product, so faint straight canal relief is visible in the terrain at every epoch,
  including 1990 and 1995 when the canals did not yet exist.
- Drape: a per-epoch Landsat NBR class image (dense forest, degraded/regrowth,
  sparse/cleared, bare/recent burn), fixed thresholds across all epochs, so colour
  changes between epochs are changes in the index, not in the ramp.
- Canals (Dadap map) in light blue from the 2000 epoch onward (the map is a single
  snapshot; the MRP canals were dug 1996 to 1998). OSM rivers in dark blue at all
  epochs.
- Fires: NASA FIRMS MODIS detections for the epoch year, smoothed into translucent
  orange-to-red patches (Gaussian width 2 km times a seeded 0.8 to 1.3 per year,
  300 m seeded positional jitter). No fires are drawn before 2001 (no FIRMS record)
  and the fire line says so. For 2026 the input is the live 7-day feed aggregated to
  0.05 deg cells; the dam-site markers at sites with a FIRMS detection within 7 days
  turn red.
- Sinking (2020, 2022 to 2026 epochs): the modelled lowering since 2018 is added to
  the surface height with an extra exaggeration and tinted magenta (lowered) / teal
  (raised) with a true-metre colour bar (+-0.30 m, the p99 of the 2026 field). No data
  outside the Tay footprint (10.9 percent of the box is covered; 35 of 60 dam sites fall
  on a valid Tay pixel).
- Dam sites (green spheres, 60 sites from `apps/web/public/live/dam_sites.geojson`)
  and seven OSM settlement labels.

## Exaggeration factors (both printed on screen)

- `VEXAG = 40`: terrain height x40. Printed top-right on every frame as
  "terrain height x40".
- `SINK_EXAG = 25`: the modelled lowering is multiplied by 25 on top of the terrain
  x40 (so x1000 in total relative to true metres). Printed as "sinking x25 on top
  of terrain x40" on every frame of an epoch that has sinking (2018 onward).

## Sinking is modelled, not measured per epoch

Tay et al. 2025 measured one long-term ALOS-2 L-band rate over 2018-10-12 to
2024-08-09 [E]. The per-epoch lowering shown here is
`rate_mm_per_yr * (year - 2018) / 1000` metres [S], i.e. a linear extrapolation of
that single rate to each epoch, including 2025 and 2026 which lie outside the
measurement window. It is labelled "modelled lowering since 2018 (from the measured
2018 to 2024 rate, Tay et al. 2025)" on screen. It is not an observation of what
happened in any given year, it carries the Tay rate's own uncertainty, and it says
nothing about the water table (InSAR gives subsidence, not a water-table proxy). No
sinking is drawn for epochs before 2018. Rate statistics on the work grid: p5 / p50 /
p95 = -13.1 / +2.2 / +22.7 mm/yr (negative = sinking), 38.7 percent of valid pixels
sinking; 2026 modelled lowering p5 = -0.10 m, minimum -1.06 m.

## NBR classes and thresholds

NBR = (NIR - SWIR2) / (NIR + SWIR2) from Landsat Collection 2 Level-2 surface
reflectance (bands nir08, swir22, QA_PIXEL cloud/shadow/cirrus masked), Jul to Oct
median composite per epoch (window widened to the neighbouring years where fewer than
three scenes exist), resampled to the 0.0006 deg work grid. Classes, fixed for every
epoch:

| Class | NBR |
| --- | --- |
| Dense forest | NBR >= 0.45 |
| Degraded / regrowth | 0.25 <= NBR < 0.45 |
| Sparse / cleared | 0.05 <= NBR < 0.25 |
| Bare / recent burn | NBR < 0.05 |
| Persistent water | NBR < -0.10 in at least half of the valid epochs |

The cut points follow the spirit of the Key and Benson 2006 FIREMON Landscape
Assessment NBR/dNBR severity scheme [S]; they are a standard-inference adaptation
for an absolute-NBR peat-swamp drape, not a calibrated burn-severity product, and the
cumulative "burn scar" (NBR < 0.05 on land at any epoch) is a proxy, not a validated
burned-area map [H]. The 2010 epoch has cloud and Landsat-7 SLC-off gaps over 37
percent of the box; in the drape those gaps are filled from the 2005 and 2015
composites, and the class percentages in the caption refer to the clear fraction only.

Class fractions used in the captions come from `burn_drapes_stats.json`.

## Data sources

| Layer | Source | Citation |
| --- | --- | --- |
| Terrain | Copernicus DEM GLO-30 as supplied with the NISAR products (`grid_products/_deletable_work/nisar/dem/nisar_dem_wgs84_ellipsoid_peatguard.tif`), WGS84 ellipsoid heights | European Space Agency, Copernicus DEM GLO-30, doi 10.5270/ESA-c5d3d65 (TanDEM-X derived) [E] |
| Burn index | Landsat Collection 2 Level-2 surface reflectance (Landsat 5 TM, 7 ETM+, 8 and 9 OLI) via Microsoft Planetary Computer, `scripts/compute_nbr_epochs.py`, 0.0003 deg composites in `grid_products/external/nbr/nbr_YYYY.tif` | USGS Landsat Collection 2 Level-2 Science Products; NBR after Key, C.H. and Benson, N.C. 2006, Landscape Assessment (LA), in FIREMON: Fire Effects Monitoring and Inventory System, USDA Forest Service RMRS-GTR-164-CD [E for the index, S for the thresholds] |
| Fires 2001 to 2024 | NASA FIRMS MODIS Collection 6.1 yearly Indonesia archive CSVs, clipped to the box (`fires_YYYY.json`, counts in `fires_counts.json`) | Giglio, L., Schroeder, W., Justice, C.O. 2016, The Collection 6 MODIS active fire detection algorithm and fire products, Remote Sensing of Environment 178, 31-41 [E] |
| Fires 2026 | `https://peatguard.com/api/fires` (NASA FIRMS last 7 days to 2026-10-08: VIIRS SNPP, VIIRS NOAA-20, MODIS), aggregated to 0.05 deg cells (`fires_2026_week.json`) | NASA FIRMS; VIIRS algorithm Schroeder, W. et al. 2014, Remote Sensing of Environment 143, 85-96 [E]. No annual FIRMS file for 2025 at render time. |
| Sinking | Tay et al. 2025 ALOS-2 ScanSAR long-term LOS rate, 2018-10-12 to 2024-08-09 (`InSAR_LongTermRate_20181012-20240809.tif`), bilinear to the work grid | Tay, C., Jovani-Sancho, A.J., Yulianti, N., Evans, C.D., Callaghan, N., Jaya, A., Salman, R., Zheng, Y., Susilo, Dohong, A., Yun, S.-H. 2025, Satellite radar advances carbon emissions accountability over tropical peat, Communications Earth and Environment 6, 971, doi 10.1038/s43247-025-02926-6 [E]. Independent benchmark; not calibrated to PeatGuard products. |
| Canals | Dadap canal map, 5 m 0/1 raster (`grid_products/external/dadap_canals_sebangau.tif`) aggregated to the 65 m grid | Dadap, N.C., Hoyt, A.M., Cobb, A.R., Oner, D., Kozinski, M., Fua, P.V., Rao, K., Harvey, C.F., Konings, A.G. 2021, Drainage canals in Southeast Asian peatlands increase carbon emissions, AGU Advances 2, e2020AV000321, doi 10.1029/2020AV000321 [E]. Optical canal map; not to be conflated with Tay 2025. |
| Rivers, settlements | OpenStreetMap waterways (`waterway` in river, stream) and named places (`grid_products/external/osm_*_sebangau.geojson`) | OpenStreetMap contributors, ODbL [E] |
| Dam sites | PeatGuard ranked canal-blocking sites, `apps/web/public/live/dam_sites.geojson` (60 sites) | PeatGuard internal |
| Captions | `captions.json`, numbers from `burn_drapes_stats.json` and `fires_counts.json` | |

Work grid for every `.npy` layer: EPSG:4326, 0.0006 deg (about 65 m), 1000 rows x
1050 columns, origin top-left (114.20, -2.41); sidecar JSON
`{"shape":[1000,1050],"lon0":114.2,"lat0":-2.41,"step":0.0006}`.

## How to regenerate

All commands run from the repo root with `.venv/bin/python3` (rasterio 1.5, numpy,
scipy, matplotlib, PIL, pyvista). The layer files (`dem.npy`, `hillshade.npy`,
`nbr_YYYY.npy`, `drape_YYYY.png`, `canals.npy`, `rivers.npy`, `water_mask.npy`,
`tay_rate_mmyr.npy`, `lowering_YYYY.npy`, `fires_YYYY.json`, `sites.json`,
`settlements.json`) are already in this directory.

Render one epoch (36 frames, camera segment t0 to t1, end exclusive):

```
.venv/bin/python3 deliverables/aoi_3d_degradation/render.py --epoch 2015 \
    --t0 0.4167 --t1 0.5000 --frames 36 \
    --caption "El Nino fires: 2,923 FIRMS detections; forest down to 74 percent of the box, 22 percent degraded" \
    --outdir deliverables/aoi_3d_degradation/frames
```

Epoch k of 12 uses `--t0 k/12 --t1 (k+1)/12`; captions are in `captions.json`.
Other options: `--start-index N` (first frame number). Frame cost is about 0.1 s each
after a few seconds of scene build.

Build the concat list (36 frames per epoch at 1/12 s, last frame repeated 18 times
for the 1.5 s hold, 24 times for the final 2 s hold, then the last file once more
without a duration, as the concat demuxer requires):

```
cd deliverables/aoi_3d_degradation
.venv/bin/python3 - <<'EOF'
import os
epochs=[1990,1995,2000,2005,2010,2015,2020,2022,2023,2024,2025,2026]
fps=12; L=[]
for e in epochs:
    for i in range(36):
        L.append("file '%s/frames/%d_%03d.png'\nduration %.6f" % (os.getcwd(), e, i, 1/fps))
    for j in range(int(round((2.0 if e==epochs[-1] else 1.5)*fps))):
        L.append("file '%s/frames/%d_035.png'\nduration %.6f" % (os.getcwd(), e, 1/fps))
L.append("file '%s/frames/2026_035.png'" % os.getcwd())
open("frames/list.txt","w").write("\n".join(L)+"\n")
EOF
```

Encode the MP4:

```
/opt/homebrew/bin/ffmpeg -y -f concat -safe 0 -i frames/list.txt \
  -vf "fps=12,scale=1920:1080:flags=lanczos,format=yuv420p" \
  -c:v libx264 -preset slow -crf 18 -profile:v high -level 4.0 -pix_fmt yuv420p \
  -movflags +faststart PeatGuard_AOI_3D_degradation_1990_2026.mp4
```

Encode the GIF (5 fps keeps it under 40 MB; `fps=10` gives 51.3 MB):

```
/opt/homebrew/bin/ffmpeg -y -f concat -safe 0 -i frames/list.txt \
  -vf "fps=5,scale=960:-2:flags=lanczos,split[a][b];[a]palettegen=max_colors=256:stats_mode=diff[p];[b][p]paletteuse=dither=bayer:bayer_scale=5:diff_mode=rectangle" \
  PeatGuard_AOI_3D_degradation_1990_2026.gif
```

Still and contact sheet: copy `frames/2026_035.png` to
`PeatGuard_AOI_3D_2026_still.png`; the contact sheet tiles frame 018 of each epoch at
640x360 in a 3x4 grid with PIL.

## Known limitations of this render

- The DEM carries the modern canal relief at every epoch (see above).
- The "Canals (Dadap map)" legend row is shown in the 1990 and 1995 frames although
  no canals are drawn in them.
- In the 2024 to 2026 epochs the near corner of the terrain reaches the bottom edge
  of the frame and the Marabahan label sits over the attribution footer.
- Epoch changes are hard cuts of the drape under a continuous camera (the camera
  itself does not jump); the 2015 to 2020 cut is the most visible because the 2015
  fire glow and degraded drape disappear at once. (v1 to v3 only; v4 crossfades,
  see below.)
- 37 percent of the 2010 drape is gap-filled from the 2005 and 2015 composites, so
  that epoch is partly a blend of its neighbours.
- The 1990 caption wraps with "1998)" alone on its second line.

## Slide version (8 Oct 2026)
PeatGuard_AOI_3D_degradation_1990_2026_SLIDE.mp4: same flyover with no caption, legend or
attribution and a flat black background, for insertion into a PowerPoint slide that carries its
own title and legend (render.py --minimal --bg "#000000"). Put the source line on the slide:
Terrain Copernicus DEM (x40 vertical exaggeration); burn index Landsat 1990 to 2026 (NBR classes after
Key and Benson 2006); fires NASA FIRMS MODIS 2001 to 2024, VIIRS+MODIS for the 2026 week; sinking
modelled from the Tay et al. 2025 ALOS-2 rate 2018 to 2024 (x25 on top of the terrain exaggeration).

## Yearly slide version (8 Oct 2026, v3)
PeatGuard_AOI_3D_fires_yearly_1990_2026_SLIDE.mp4 (64 s): 1990, 1995, 2000 as the intro, then EVERY year
2001 to 2024 and the 2026 week. Land cover for a year is the nearest Landsat composite (the composites
exist only for 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2022 to 2026). Canals are not drawn. Fires are
NASA FIRMS MODIS detections smoothed into organic patches; ground that burned in earlier years stays as
a dark scar (faint for one fire year, dark for three or more), so repeat burning accumulates on screen.
Regenerate: render.py --epoch YYYY --minimal --bg "#000000" --no-canals ... then ffmpeg concat of
frames_yearly/list.txt. (The v3 frames in `frames_yearly/` were made with the v3 render.py, which is
no longer in the tree; the v4 render.py below replaces it.)

## Yearly slide version v4 (8 Oct 2026): degraded peat, char and live 3D fire

Deliverables:

| File | Content |
| --- | --- |
| `PeatGuard_AOI_3D_fires_yearly_1990_2026_SLIDE_v4.mp4` | 1920x1080, 12 fps, H.264 High 4.0, yuv420p, CRF 18, faststart. The judges' demo clip, for a half-width PowerPoint placeholder. |
| `PeatGuard_AOI_3D_fire_proof_2015.mp4` | 10 s, fixed camera, 2015 alone: the flame and smoke animation on its own. |
| `frames_v4/YYYY_NNN.png` | moving frames then static hold frames of each year; `YYYY_pre_NN.png` are the 8 pre-roll frames at the previous year's last pose used for the crossfade; `frames_v4/xfade/` holds the blended crossfade frames and `frames_v4/list.txt` the concat list (both written by assemble.py). |

What is on the terrain, and nothing else (render.py v4, slide mode):

1. A hillshaded forest-green base, one colour (#2f6b45 times the Horn hillshade, about #244f35
   to #3f8458), reading as intact peat forest. No Landsat land-cover classes are drawn any more;
   the drape_YYYY.png files and the NBR thresholds above belong to v1 to v3.
2. Pre-existing degraded peat as of each year in one muted red (#b3261e at alpha 0.55), defined
   by what the satellites observed before that year and never inferred from canals: the whole
   box is peat (PeatGrids depth >= 50 cm everywhere in the box) AND any of (a) Landsat vegetation
   loss, NBR < 0.45 in any composite strictly before the year (NaN gaps such as the 2010 cloud and
   SLC-off gaps do not count; persistent water excluded), (b) a FIRMS fire footprint (sigma 0.4
   km per detection) in any earlier year, (c) from 2018 only, when the measurement starts, the
   radar-verified sinking zone, Tay et al. 2025 ALOS-2 long-term rate below -5 mm/yr [E for the
   rate], cleaned with a binary opening then closing (diamond of radius 2 cells, about 5 cells
   across) and blobs under 40 cells dropped (0.74 percent of the box; NaN outside the 10.9
   percent Tay coverage is never verified). The union is opened with a 3-cell diamond so it reads
   as zones, not speckle. The red grows through the sequence: none in 1990 (nothing observed
   before the first composite), 9 percent of the box in 2001, 54 percent in 2005, 79 percent in
   2015, 86 percent in 2026; every flame therefore lands either on red (degraded before it
   burned) or on green. The radar-verified subset is additionally outlined by a dashed cyan line
   (#9fe8ff, 3 px, dashes 0.5 km on 0.3 km off cut by arc length, lifted 80 m above the surface
   in the x40 frame). The inferred canal-proximity criterion is deliberately not drawn; nothing
   on screen calls an inferred area verified.
   A headline under the year (40 px, from 2001) states "N percent of this year's fire is on peat
   already degraded before it burned", N from `deliverables/degraded_peat_fire_recurrence/
   fire_on_preexisting_degraded.json` (field frac_pre_degraded, computed on that deliverable's
   0.0045 deg grid; 22 percent in 2001 rising to 88 in 2015 and 97 to 100 from 2016); for 2026
   the line reads "of the fire in the FIRMS last 7 days" and the fraction (86 percent) is
   computed here the same way against the pre-2026 red on the work grid. Recomputing the
   JSON years on the work grid gives numbers within a few points (2015: 90 versus 88).
3. Burned before: black at alpha 0.15 after one earlier FIRMS fire year rising linearly to 0.75 at
   eight or more, on top of the red, accumulating from 2002. A detection leaves a footprint of
   sigma 0.4 km (about one MODIS pixel). Over 2001 to 2024, 82 percent of the box burned in at
   least one year and 36 percent in four or more, so the single-year mark is kept faint and the
   repeat burners carry the layer.
4. This year's fire as live 3D flames: the fine detection density (sigma 0.6 km times a seeded
   0.8 to 1.3 per year) is gated, no flame below 0.5 detections/km2 and full height and opacity at
   2.0/km2, plus a small flame (alpha 0.45, sigma 0.25 km) right at a lone detection so quiet years
   show a few distinct fires. A flame mesh (a copy of the terrain grid) has z = terrain + 50 m +
   600 m x alpha x flicker, coloured #d7301f edge, #ff8c1a mid, #ffe066 core by density, per-point
   opacity = alpha, blanked (NaN) below alpha 0.05 so it never z-fights the terrain. flicker is a
   seeded noise field (sigma 3 cells) that crawls 0.8 and 1.4 cells per frame and swings between
   two fields every 9 frames, clipped to 0.7 to 1.3; it also scales opacity and shifts the colour,
   so the flames rise, fall and crawl at 12 fps. A grey-blue (#8e9bb0) smoke mesh sits 0.9 km above
   the flame base at 0.35 x flame alpha, as a plume of 8 puffs 4 cells apart drifting north-west 2
   cells per frame (new puffs fade in at the source). The coarse fused field (sigma 2 km) is baked
   into the texture only as a low ground glow at a quarter of its v3 alpha, and the flame field at
   half strength as that year's ground mark. In 2015 about a third of the box carries flames.
   Both meshes are updated in place each frame (points and rgba scalars), not rebuilt.
5. OSM rivers in thin blue #3f7fbf. Dam sites as white spheres (radius 0.55 km) with a thin dark
   rim (#101418, a 0.70 km sphere drawn back-faces only); seven settlement labels; north marker.
   No canals are drawn in the final cut (`canals_osm_1px.npy` and `canals_planet_1px.npy`, the
   1 px OSM and skeletonised PlanetScope rasterisations tried earlier, stay on disk unused;
   `canals.npy` is the 19 percent max-resampled Dadap blob, also unused).
6. Imagery base for 2026 only: `planet_nicfi_2026_09.npy` (1000 x 1050 x 3 uint8 RGB on the work
   grid, Planet NICFI monthly mosaic September 2026, 4.77 m source at zoom 13, about 19 m per
   pixel) replaces the green base in the 2026 segment, with the hillshade at half strength so
   the imagery keeps its colour and the red degraded layer at alpha 0.40 so the imagery shows
   through; legend footer "Imagery: Planet NICFI mosaic, Sep 2026". The historical years keep the
   green base because the imagery is a 2026 snapshot. The green 2024 frame dissolves into the
   imagery 2026 frame over 12 frames (the 2024 hold is extended to 12 to carry it).

Crossfades blend only the rendered scene: render.py writes the scene without the HUD
(`<frame>.raw.png`) for every hold and pre-roll frame, and assemble.py / cut_two_act.py blend
those and compose the slide-mode HUD once per frame, the year switching at the midpoint with a
four-frame text fade, so the numerals, headline and legend never double-expose.

Also: from 2018 the modelled lowering (Tay rate) is still applied to z with SINK_EXAG (geometry
only, no tint), and the note bottom-right says "terrain height x40; sinking x25".

Legend (bottom-left, every frame, panel rgba 10,14,20 at 0.72 with 10 px radius, DejaVu Sans 26 px
labels, 20 px title at 70 percent white, 28 x 18 swatches): Degraded before this year (Landsat
vegetation loss, earlier fire, or radar sinking); Verified by radar sinking (Tay et al. 2025) as
the dashed cyan line; Burned
before (black, darker = more years) from 2002; Burning this year from 2001 (2026: Burning now,
FIRMS, last 7 days); Dam site. Year top-left at 110 px bold. Type is sized for a clip at half the
width of a slide.

Camera: the same 35 degree oblique orbit, azimuth -40 to +40 degrees, at 158 km with a 22
degree view angle (less near-corner magnification than 132 km / 30 degrees) and the projection
shifted up by 0.28 of the viewport (VTK window centre), so the box fills 87 percent of the width at
the extreme azimuths and its near corner no longer leaves the frame.

Camera motion. Every camera parameter of every frame comes from one function of the global
output frame index g (`render.camera_t_at(g)`, used by `set_camera`): position, focal point, up
and view angle are fixed apart from the azimuth, which runs -40 to +40 degrees as
`0.25 u + 0.75 smoothstep(u)` with u = g / (frames - 1). The orbit never stops: a year's
"hold" is only its data layer persisting while the camera continues, and the dissolves run under
the moving camera. (The earlier v4 builds held the camera during holds and dissolves and derived
poses from per-year t ranges, which is what read as janky.) Pure smoothstep was rejected because it
starts at zero velocity, so the first frame steps grow several-fold frame to frame; the 25 percent
linear share keeps every consecutive azimuth step within about 2 percent of its neighbours while
the orbit still eases in and out. assemble.py and cut_two_act.py assert that the final list has
exactly one entry per global frame and that the azimuth is strictly monotonic with neighbouring
steps within 10 percent, and print the maximum step (full clip: 0.139 deg per frame at the middle
of the sequence, neighbouring steps within 2.2 percent; the mean is 0.101. Two-act cut: 0.268 deg,
within 4.3 percent).

No jitter between years. Diagnosis of the v3 pop: (1) the drape jumped to the nearest Landsat
epoch, (2) the fire glow and scar fields change seed per year, (3) the concat list repeated the
last frame for the hold then cut to the next year. In v4 the texture is year-independent except
for the accumulating red, char and the fire fields; the flicker noise has one seed and runs on the
global frame index, so the same frame gives the same flame state in two years; and every year
boundary is a pure data dissolve: render.py draws year Y+1's "pre-roll" at the global frames of
year Y's last xfade frames (same camera, same flicker), and assemble.py blends the raw scene
pairs in PIL with weights 1/(n+1) .. n/(n+1) and composes the HUD once per frame (year switching
at the midpoint with a four-frame text fade).

Schedule (frames per year, all under the moving camera): 1990, 1995, 2000 at 36 (3.0 s); 2001 to
2024 at 26 (2.17 s), the last 8 of each dissolving into the next year (12 into 2026); 2026 at 60
(5.0 s). 792 frames = 66.0 s.

Regenerate (from the repo root, about 3 minutes with 7 processes):

```
.venv/bin/python3 deliverables/aoi_3d_degradation/render.py --all --jobs 7 \
    --minimal --bg "#000000" --outdir deliverables/aoi_3d_degradation/frames_v4
.venv/bin/python3 deliverables/aoi_3d_degradation/assemble.py \
    --frames deliverables/aoi_3d_degradation/frames_v4 \
    --out deliverables/aoi_3d_degradation/PeatGuard_AOI_3D_fires_yearly_1990_2026_SLIDE_v4.mp4
```

### Two-act cut: PeatGuard_AOI_3D_degradation_then_fire_SLIDE.mp4 (about 34 s)

`cut_two_act.py` builds a second clip from the same layers, 411 frames = 34.25 s, rendered
directly at 12 fps (no frame reuse from the full clip). Act 1, degradation, 1990 to 2005 (global
frames 0 to 84): the green base with the red pre-existing degraded layer growing frame by frame,
no flames and no char even where fires occurred, headline "Peat degraded: N percent of the box"
climbing. The red is interpolated per frame between the masks of consecutive key years (1990,
1995, 2000, 2001 to 2006; cells new in the later mask fade in with the fraction of the segment
elapsed), so the Landsat-epoch and yearly steps both grow smoothly; the last 4 act 1 frames
dissolve into 2006. Act 2, fire, 2006 to 2026 (`render.py --spec twoact`): red as of the year,
black char from the 2001 onward fires so 2006 opens with the 2001 to 2005 burn history, flames,
fire-on-degraded headline, full legend; 13 frames per year, 23 for 2015 and 2019 (1.2 s longer),
21 for 2024 (12-frame dissolve into the imagery 2026), 51 for 2026 (3 s hold), 4-frame data
dissolves elsewhere. The camera is one eased function of the global frame index over the whole
cut and never stops. Work frames in `frames_two_act/`:

```
.venv/bin/python3 deliverables/aoi_3d_degradation/cut_two_act.py        # act 1, act 2, assemble
.venv/bin/python3 deliverables/aoi_3d_degradation/cut_two_act.py --skip-act1 --skip-act2   # assemble only
```

### Transparent, separable layers: layers_transparent/ (the land as a floating object)

`layers_transparent.py` renders the two-act schedule (411 frames, 34.25 s, the same global camera
function, lighting and frame count) four times with a transparent background (pyvista
`screenshot(transparent_background=True)`, RGBA frames), so the passes stack pixel for pixel in
PowerPoint: `land` (green hillshade, rivers, the NICFI imagery in the 2026 segment; no red, char or
fire), `peat` (the same terrain mesh textured with an RGBA texture whose alpha is the red degraded
layer plus the black char and 0 elsewhere, so only the red and char float in 3D), `fire` (flame
and smoke meshes only; fully transparent in act 1) and `all` (the full composite). Nothing but
geometry: no dam sites, labels, north mark, legend, text or year, and no cyan outline, so
`all` = `land` over `peat` over `fire` exactly. Year dissolves are premultiplied RGBA blends of the
two years rendered at the same global frame. Each pass is encoded twice at 12 fps with alpha:
`<pass>_prores4444.mov` (prores_ks profile 4444, yuva444p10le, StraightAlpha; PowerPoint on Mac)
and `<pass>_hevc_alpha.mov` (hevc_videotoolbox, alpha_quality 0.9, hvc1, bgra input,
AlphaChannelMode PremultipliedAlpha as read back through AVFoundation; PowerPoint 365 on Mac and
Windows 10 and later), plus `<pass>_2026_last.png` (the last 2026 frame, transparent PNG) and
`<pass>_2026_last_check.png` (the same over a checkerboard). ffmpeg cannot decode Apple's
auxiliary HEVC alpha layer, so its decode of the HEVC file shows no alpha; the ProRes file and the
PNGs are the reference. `report.json` holds frame counts, sizes and alpha statistics.

```
.venv/bin/python3 deliverables/aoi_3d_degradation/layers_transparent.py --test   # 2 s of 2015, alpha check
.venv/bin/python3 deliverables/aoi_3d_degradation/layers_transparent.py          # full set, 7 processes
```

One year only: `render.py --year 2015 --minimal --bg "#000000" --outdir .../frames_v4`. The
proof clip: `render.py --epoch 2015 --t0 0.64 --t1 0.64 --frames 120 --minimal --bg "#000000"
--outdir <dir>` then `ffmpeg -framerate 12 -i <dir>/2015_%03d.png` with the same x264 settings as
assemble.py. `--no-canals` is accepted for compatibility but canals are never drawn in v4. The
annotated HUD mode (no `--minimal`) still works with the v4 layers and the same legend.

Caveats: the verified zone is small and fragmented (46 blobs), so the cyan outline reads as a set
of small dashed rings rather than one boundary; that is the data at this threshold. The red
layer is time-correct only to the resolution of the Landsat composites (five-yearly to 2015,
then 2020, 2022 to 2025) and the yearly FIRMS record; by 2015 it covers 79 percent of the box
because one NBR < 0.45 epoch or one fire year anywhere is enough, so late years are mostly red
with the green base showing only in the never-degraded remnants. The Landsat vegetation-loss
criterion uses a fixed absolute NBR cut [S], not a calibrated degradation product. The 2026 fire
input is the 7-day feed, not a full year.
