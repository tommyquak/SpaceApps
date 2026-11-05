# Figure and animation code

The two folders here hold the code that produced the demo videos and the
static figures. They run on prepared work-grid layers that are too large for
the repository, so every script reads its inputs below a `DATA_ROOT`
directory (default: two levels above the script, which matches the original
monorepo layout). Set it explicitly:

```bash
export DATA_ROOT=/path/to/PeatGuard        # holds deliverables/ and grid_products/
export FFMPEG=/opt/homebrew/bin/ffmpeg     # optional; defaults to ffmpeg on PATH
```

Python dependencies for this folder: `pip install -e '.[viz]'` (pyvista with a
working off-screen VTK, Pillow, scikit-image) plus ffmpeg.

## `aoi_3d_degradation/`

Off-screen 3D flyover of the dam-site AOI (114.20 to 114.83 E, 3.01 to 2.41 S)
from 1990 to 2026. `render.py` renders frames, `assemble.py` concatenates them
into the yearly fire video, `cut_two_act.py` builds the two-act
"degradation, then fire" cut, and `layers_transparent.py` exports the drape
layers with alpha for compositing. `captions.json` holds the per-epoch
captions. `README.md` in the folder documents every layer, threshold and
exaggeration factor.

Inputs read from `$DATA_ROOT/deliverables/aoi_3d_degradation/`:

| Input | Source | Produced by |
| --- | --- | --- |
| `dem.npy`, `dem.json`, `hillshade.npy` | Copernicus GLO-30 DEM | work-grid preparation (not in this repository) |
| `nbr_YYYY.npy`, `burned_cum_YYYY.npy` | Landsat Collection 2 Level-2 NBR composites | `scripts/compute_nbr_epochs.py` |
| `canals.npy`, `canals_osm_1px.npy`, `canals_planet_1px.npy` | Dadap et al. 2021 canal map, OSM, PlanetScope detection | `src/peatguard/analysis/canal_detect_planet_v2.py` |
| `rivers.npy`, `water_mask.npy` | OSM waterways, water mask | `src/peatguard/analysis/water_mask.py` |
| `fires_YYYY.json`, `fires_2026_week.json` | NASA FIRMS MODIS C6.1 and the 7-day feed | `src/peatguard/ingest/fire.py` |
| `tay_rate_mmyr.npy`, `lowering_YYYY.npy` | Tay et al. 2025 ALOS-2 long-term rate | `scripts/build_tay_restoration.py` |
| `planet_nicfi_2026_09.npy` | Planet NICFI mosaic (non-commercial tier) | `src/peatguard/ingest/planetscope.py` |
| `sites.json`, `settlements.json` | PeatGuard dam sites, OSM settlements | `scripts/select_dam_sites_v2.py` |

Outputs: `frames_v4/`, `frames_two_act/`, `layers_transparent/` and the MP4s
listed in the root README, written next to the inputs.

## `degraded_peat_fire_recurrence/`

`build.py` is the single script that produces the suspected-degraded-peat
screening map, the 2001 to 2024 fire-return-count map, the fire-years versus
canal-distance chart, the `fires_return_to_degraded_peat_2001_2026` videos and
`stats.json`. The committed `stats.json`, `fire_on_preexisting_degraded.json`
and `fire_years_vs_canal_distance.json` are the numbers quoted in the README
and used by `aoi_3d_degradation/render.py` for the per-year
"fire on pre-degraded peat" headline.

Inputs: the `aoi_3d_degradation` layers above, raw NASA FIRMS VIIRS country
CSVs cached in `$DATA_ROOT/grid_products/external/firms_viirs/raw/`, and the
PeatGrids mean thickness raster
`$DATA_ROOT/grid_products/external/peatgrids/PEATGRIDS_THICKNESS_CM_MEAN.tif`.
Definitions, thresholds and the confidence tags for every number are in the
folder README.
