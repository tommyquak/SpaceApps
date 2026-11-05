"""FastAPI dashboard application for PeatGuard COG product visualization.

Serves a Leaflet-based web map with TiTiler-powered dynamic tile endpoints
for each pipeline product. Products are read from a local directory of
Cloud-Optimized GeoTIFFs (``PEATGUARD_WEBDATA_DIR``, default ``data/sample``
in this repository) or, when a file is absent locally, from a public GCS
bucket.

Usage:
    PEATGUARD_WEBDATA_DIR=data/sample uvicorn peatguard.dashboard.app:app --host 0.0.0.0 --port 8080
"""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Annotated
from urllib.parse import quote, urlparse

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, JSONResponse, Response
from fastapi.templating import Jinja2Templates
from rio_tiler.errors import PointOutsideBounds, TileOutsideBounds
from titiler.core.factory import TilerFactory

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

GCS_BUCKET = os.environ.get("PEATGUARD_STORAGE__GCS_BUCKET", "peatguard-data")
GCS_PREFIX = os.environ.get("PEATGUARD_GCS_PREFIX", "products")

# Public HTTPS URL pattern for COGs in GCS.  TiTiler reads COGs via HTTP
# range requests, so the bucket must either be public or accessed through
# signed URLs.  For simplicity we use the public storage.googleapis.com
# endpoint; swap to signed URLs if the bucket is private.
GCS_BASE_URL = f"https://storage.googleapis.com/{GCS_BUCKET}/{GCS_PREFIX}"

# Local-mode override: a directory of COGs is served in preference to the GCS
# bucket so the dashboard works offline against committed sample data.
# PEATGUARD_WEBDATA_DIR names the directory; it defaults to the package's
# webdata/ folder. Within it a product <filename> resolves to either
# <filename> itself or DEMO_products_<filename> (the legacy naming that
# mirrors the GCS layout).
_PACKAGE_WEBDATA_DIR = Path(__file__).resolve().parent.parent / "webdata"
_WEBDATA_DIR = Path(
    os.environ.get("PEATGUARD_WEBDATA_DIR", str(_PACKAGE_WEBDATA_DIR))
).expanduser().resolve()


def _local_candidates(directory: Path, filename: str) -> tuple[Path, Path]:
    return directory / filename, directory / f"DEMO_products_{filename}"


def resolve_cog_url(filename: str, aoi: str | None = None) -> str:
    """Return either an absolute local path (if a DEMO file exists) or the GCS URL.

    TiTiler accepts both file paths and HTTPS URLs as the ``url=`` parameter,
    so this lets the same product catalog work locally and in cloud mode
    without changing the frontend.

    When ``aoi`` names a region with a pre-clipped COG (webdata/aoi_clips/
    <AOI>/...), that small clip is served instead of the full-grid raster:
    focused views never read, decode, or ship a pixel outside the AOI, and
    TiTiler's per-tile window reads hit a file an order of magnitude smaller.
    """
    if aoi:
        for clip in _local_candidates(_WEBDATA_DIR / "aoi_clips" / aoi, filename):
            if clip.exists():
                return str(clip)
    for candidate in _local_candidates(_WEBDATA_DIR, filename):
        if candidate.exists():
            return str(candidate)
    return f"{GCS_BASE_URL}/{filename}"


def validated_cog_path(
    url: Annotated[str, Query(description="Dataset URL")],
) -> str:
    """Allow TiTiler to read only PeatGuard-owned COG locations.

    TiTiler's default dependency accepts an arbitrary local path or URL. The
    public map endpoints must not become a local-file reader or an SSRF proxy.
    """
    parsed = urlparse(url)
    if parsed.scheme:
        # Per-segment validation: the HTTP client normalises dot segments
        # after this check runs, so ".." (or a percent-encoded variant) in
        # the path could escape the bucket/prefix that a plain startswith
        # comparison appears to enforce. Reject any empty or dot segment and
        # any percent escape, then require the leading segments to equal the
        # bucket and prefix exactly.
        segments = parsed.path.split("/")
        expected_prefix = [GCS_BUCKET, *GCS_PREFIX.split("/")]
        if (
            parsed.scheme == "https"
            and parsed.hostname == "storage.googleapis.com"
            and len(segments) > len(expected_prefix) + 1
            and segments[0] == ""
            and all(
                seg and seg not in {".", ".."} and "%" not in seg
                for seg in segments[1:]
            )
            and segments[1 : len(expected_prefix) + 1] == expected_prefix
        ):
            return url
        raise HTTPException(
            status_code=403,
            detail="COG URL is outside the PeatGuard data catalog",
        )

    candidate = Path(url)
    try:
        resolved = candidate.resolve(strict=True)
        webdata_root = _WEBDATA_DIR.resolve(strict=True)
        resolved.relative_to(webdata_root)
    except (FileNotFoundError, OSError, RuntimeError, ValueError):
        raise HTTPException(
            status_code=403,
            detail="COG path is outside the PeatGuard data catalog",
        ) from None
    if not resolved.is_file() or resolved.suffix.lower() not in {".tif", ".tiff"}:
        raise HTTPException(
            status_code=403,
            detail="COG path is not an approved raster",
        )
    return str(resolved)


# Map opens centered on the Tay full grid (~69x64 km, Sebangau / ex-Mega Rice)
AOI_CENTER_LAT = -2.62
AOI_CENTER_LON = 114.58

# Discrete TiTiler colormap for subsidence_class.tif, URL-encoded for use as
# a tile query parameter. Class codes follow analysis/subsidence_class.py
# (0=NoData, 1=Severe, 2=Active Drying, 3=Moderate Drying, 4=Stable,
# 5=Rebound/Noise, 6=Water); values missing from the map render transparent.
_SUBSIDENCE_CLASS_CMAP = quote(
    json.dumps(
        {
            "1": "#a50026",
            "2": "#f46d43",
            "3": "#fee08b",
            "4": "#1a9850",
            "5": "#4575b4",
            "6": "#969696",
        }
    )
)

# Product definitions: filename -> display metadata
PRODUCTS = {
    "subsidence_velocity.tif": {
        "label": "Land Subsidence - Tay 2025 (mm/yr)",
        # rdbu (not reversed) runs red -> white -> blue from low to high, so
        # negative velocity = subsidence = red and uplift = blue, matching
        # the legend below.
        "colormap": "rdbu",
        "rescale": "-60,20",
        "description": (
            "Long-term land-motion rate from Tay et al. (2025), Comms Earth & "
            "Environment -- published ALOS-2 L-band InSAR. Red = subsidence, "
            "blue = uplift. Shown only where Tay measured (~11% of the grid)."
        ),
        "legend_type": "diverging",
        "legend_labels": ["-60 mm/yr (subsidence)", "0", "+20 mm/yr (uplift)"],
        "legend_colors": ["#d73027", "#ffffbf", "#4575b4"],
        "default_on": False,
        "tile_extra": "&nodata=-9999&tilesize=512",
    },
    "co2_emission_tay.tif": {
        "label": "Implied CO2 Emission - Tay 2025 model (t CO2/ha/yr)",
        "colormap": "ylorrd",
        "rescale": "0,40",
        "description": (
            "Peat oxidation emission implied by the Tay et al. (2025) 2018-2024 "
            "subsidence rate, converted with Eq. 5 of that paper (after "
            "Couwenberg and Hooijer): subsidence (cm/yr) x 0.08 g/cm3 dry bulk "
            "density x 55 percent carbon x 3.67, i.e. 16.15 t CO2/ha/yr per "
            "cm/yr. All subsidence is attributed to oxidation (appropriate for "
            "long-drained peat); uplift is set to zero. A baseline emission "
            "estimate, not a fire emission, not a NISAR measurement and not a "
            "credit volume."
        ),
        "legend_type": "sequential",
        "legend_labels": ["0", "20 t CO2/ha/yr", "40+ t CO2/ha/yr"],
        "legend_colors": ["#ffffb2", "#fd8d3c", "#bd0026"],
        "default_on": False,
        "tile_extra": "&nodata=-9999&tilesize=512",
        "public_role": "direct",
        "measurement_status": "derived_from_published_benchmark",
    },
    "elnino_max_2019.tif": {
        "label": "2019 El Nino Max Drawdown - Tay 2025 (cm)",
        # rdbu (not reversed): negative = subsidence/drawdown renders red.
        "colormap": "rdbu",
        "rescale": "-10,10",
        "description": (
            "Maximum peat-surface subsidence during the 2019 El Nino dry "
            "season (15 May - 16 Nov 2019) relative to the pre-El-Nino mean, "
            "from Tay et al. (2025) ALOS-2 L-band InSAR (DR-NTU "
            "10.21979/N9/JGQE4O). Strong negative values mark peat that "
            "draws down hardest under drought -- a drought-vulnerability "
            "indicator. Shown only where Tay measured; decorrelated areas "
            "are masked."
        ),
        "legend_type": "diverging",
        "legend_labels": ["-10 cm (max drawdown)", "0", "+10 cm"],
        "legend_colors": ["#d73027", "#ffffbf", "#4575b4"],
        "default_on": False,
        "tile_extra": "&nodata=-9999&tilesize=512",
    },
    "subsidence_class.tif": {
        "label": "Subsidence Classification",
        # The named colormap/rescale keys are kept only for the legacy
        # Leaflet template; /api/v1 tile URLs skip both (flags below) and
        # render through the explicit per-class colormap in tile_extra.
        # A continuous colormap over a categorical raster interpolates
        # between class codes and cannot match discrete legend swatches.
        "colormap": "rdylgn_r",
        "rescale": "0,4",
        "description": (
            "Severity classes: 1=Severe (< -50 mm/yr), 2=Active Drying "
            "(-50 to -20), 3=Moderate Drying (-20 to -5), 4=Stable "
            "(-5 to +5), 5=Rebound/Noise (> +5), 6=Water."
        ),
        "legend_type": "categorical",
        "legend_labels": [
            "Severe",
            "Active Drying",
            "Moderate Drying",
            "Stable",
            "Rebound/Noise",
            "Water",
        ],
        "legend_colors": [
            "#a50026",
            "#f46d43",
            "#fee08b",
            "#1a9850",
            "#4575b4",
            "#969696",
        ],
        "default_on": False,
        "tile_skip_named_colormap": True,
        "tile_skip_rescale": True,
        # nodata=0 masks the NoData class code so it never renders.
        "tile_extra": f"&colormap={_SUBSIDENCE_CLASS_CMAP}&nodata=0&tilesize=512",
    },
    "canal_mask.tif": {
        "label": "Canal Network",
        "colormap": "blues",
        "rescale": "0.5,1",
        "description": "Detected drainage canals from VV backscatter.",
        "legend_type": "binary",
        "legend_labels": ["No Canal", "Canal"],
        "legend_colors": ["transparent", "#08519c"],
        "default_on": False,
        # Binary mask: tell TiTiler that value 0 is nodata so all "no canal"
        # pixels render transparent. Value 1 gets coloured by `blues` clamped
        # at the top of the rescale range (dark blue).  tilesize=512 for
        # ArcGIS-grade resolution on retina displays.
        "tile_extra": "&nodata=0&tilesize=512",
    },
    "canal_mask_planet.tif": {
        "label": "Canal Network (PlanetScope-derived, 5.6 m output)",
        "colormap": "blues",
        # The mask is stored as 255 with AVERAGE-resampled overviews, so a
        # zoomed-out pixel encodes percent canal coverage. Rescaling 8..160
        # renders any overview pixel with >=~3% canal coverage -- thin canal
        # lines stay visible at every zoom (nearest-resampled binary
        # overviews dropped them until you zoomed in).
        "rescale": "8,160",
        "description": (
            "Drainage canals detected from native ~3 m PlanetScope SuperDove 8-band "
            "surface reflectance (v2.1 detector inside the four AOIs: NDWI "
            "+ RivaMap multiscale-singularity ridge response, hysteresis "
            "linking, width-based wide-river rejection, directional gap "
            "closing and collinear/geodesic centerline bridging; 8-scene "
            "2024 dry-season composite over Sebangau Central-East). The "
            "deployed raster is posted at approximately 5.6 m. "
            "Methodologically aligned with Dadap et al. (2021), AGU "
            "Advances. Non-commercial / research use only (Planet Education "
            "and Research licence)."
        ),
        "legend_type": "binary",
        "legend_labels": ["No Canal", "Canal (PlanetScope-derived)"],
        "legend_colors": ["transparent", "#0eb4e8"],
        "default_on": True,
        "tile_extra": "&nodata=0&tilesize=512",
        "posting_m": 5.6,
        "effective_support_m": 5.6,
    },
    "canal_risk.tif": {
        "label": "Canal Risk Score",
        "colormap": "rdylgn_r",
        "rescale": "0,1",
        "description": "Combined subsidence + canal proximity risk (0=low, 1=high).",
        "legend_type": "sequential",
        "legend_labels": ["0 (Low)", "0.5", "1.0 (High)"],
        "legend_colors": ["#1a9850", "#ffffbf", "#d73027"],
        "default_on": False,
        # nodata=-9999 masks the surrounding pipeline-fill area (re-cropped
        # to the actual AOI but kept defensive).  tilesize=512 for retina.
        "tile_extra": "&nodata=-9999&tilesize=512",
    },
    "fire_history.tif": {
        "label": "Fire History",
        "colormap": "ylorrd",
        "rescale": "0,1",
        "description": (
            "Recency- and frequency-weighted burn history from MODIS "
            "MCD64A1 burned-area (2015-2024). 0 = never burned, "
            "1 = recent and/or repeatedly burned."
        ),
        "legend_type": "sequential",
        "legend_labels": ["0 (no burn)", "0.5", "1.0 (recent/repeated)"],
        "legend_colors": ["#ffffb2", "#fd8d3c", "#bd0026"],
        "default_on": False,
        "tile_extra": "&tilesize=512",
    },
    "restoration_zones.tif": {
        "label": "Restoration Priority",
        "colormap": "rdylgn_r",
        "rescale": "0,1",
        "description": (
            "Unified restoration-priority surface: additive fusion of "
            "canal_risk (subsidence + canal proximity) and fire history. "
            "The headline layer for deciding where to act."
        ),
        "legend_type": "sequential",
        "legend_labels": ["0 (Low)", "0.5", "1.0 (Critical)"],
        "legend_colors": ["#1a9850", "#ffffbf", "#d73027"],
        "default_on": True,
        "tile_extra": "&nodata=-9999&tilesize=512",
    },
    "thermal_stress.tif": {
        "label": "Thermal Stress (ECOSTRESS)",
        "colormap": "rdylgn_r",
        "rescale": "0,1",
        "description": (
            "Diurnal land-surface-temperature range from ECOSTRESS day + "
            "night overpasses (a thermal-inertia surface-moisture proxy: "
            "wet peat swings little, drained peat swings hard). Diagnostic "
            "layer -- shown for context, not fused into restoration "
            "priority pending stronger validation."
        ),
        "legend_type": "sequential",
        "legend_labels": ["0 (wet/stable)", "0.5", "1.0 (dry/stressed)"],
        "legend_colors": ["#1a9850", "#ffffbf", "#d73027"],
        "default_on": False,
        "tile_extra": "&nodata=-9999&tilesize=512",
    },
    "vv_median_db.tif": {
        "label": "VV Backscatter (dB)",
        "colormap": "greys",
        "rescale": "-25,-5",
        "description": "Median VV backscatter composite in decibels.",
        "legend_type": "sequential",
        "legend_labels": ["-25 dB", "-15 dB", "-5 dB"],
        "legend_colors": ["#000000", "#808080", "#ffffff"],
        "default_on": False,
    },
    "water_mask.tif": {
        "label": "Water Mask",
        "colormap": "blues",
        "rescale": "0,1",
        "description": "Water bodies detected from low VV backscatter.",
        "legend_type": "binary",
        "legend_labels": ["Land", "Water"],
        "legend_colors": ["transparent", "#2171b5"],
        "default_on": False,
        # The mask COG is written with nodata=255 (water_mask.py); pass it
        # through so nodata pixels are masked rather than rendered as a
        # class value.
        "tile_extra": "&nodata=255&tilesize=512",
    },
    "coherence_median.tif": {
        "label": "Temporal Coherence",
        "colormap": "viridis",
        "rescale": "0,1",
        "description": "Median InSAR coherence (0=decorrelated, 1=stable).",
        "legend_type": "sequential",
        "legend_labels": ["0 (Low)", "0.5", "1.0 (High)"],
        "legend_colors": ["#440154", "#21918c", "#fde725"],
        "default_on": False,
    },
    "nisar_coherence_80m.tif": {
        "label": "NISAR Minimum Pair Coherence (80 m sample, Provisional)",
        "colormap": "viridis",
        "rescale": "0,1",
        "description": (
            "4x4 block mean of the 20 m minimum mission GUNW coherence across "
            "the four 12-day-chain pairs from 4 Jul to 14 Sep 2026 (ascending "
            "track 125, frame 176), downsampled to 80 m for the open-source "
            "sample bundle. A radar measurement-quality layer, not "
            "displacement, subsidence, or velocity. NISAR product maturity "
            "is PROVISIONAL."
        ),
        "legend_type": "sequential",
        "legend_labels": [
            "0 (decorrelated)",
            "0.5 (commissioning threshold)",
            "1.0 (high coherence)",
        ],
        "legend_colors": ["#440154", "#21918c", "#fde725"],
        "default_on": True,
        "tile_extra": "&nodata=-9999&tilesize=512",
        "maturity": "PROVISIONAL",
        "measurement_status": "supporting_quality_only",
        "public_role": "supporting",
        "posting_m": 80,
        "effective_support_m": 80,
        "acquisition_start": "2026-07-04",
        "acquisition_end": "2026-09-14",
        "span_days": 72,
        "pair_count": 4,
        "sensor": "NISAR L-SAR",
        "track": 125,
        "frame": 176,
        "polarization": "HH",
        "ground_truth": False,
    },
    "nisar_commissioning_minimum_pair_coherence_80m.tif": {
        "label": "NISAR Commissioning Pair Coherence (22 Jun to 21 Aug 2026, 80 m)",
        "colormap": "viridis",
        "rescale": "0,1",
        "description": (
            "Minimum coherence of the first commissioning GUNW pairs "
            "(22 Jun to 21 Aug 2026, ascending track 125, frame 176) at 80 m. "
            "Measurement-quality layer only; PROVISIONAL maturity."
        ),
        "legend_type": "sequential",
        "legend_labels": ["0 (decorrelated)", "0.5", "1.0 (high coherence)"],
        "legend_colors": ["#440154", "#21918c", "#fde725"],
        "default_on": False,
        "tile_extra": "&nodata=-9999&tilesize=512",
        "maturity": "PROVISIONAL",
        "measurement_status": "supporting_quality_only",
        "public_role": "supporting",
        "posting_m": 80,
        "effective_support_m": 80,
        "acquisition_start": "2026-06-22",
        "acquisition_end": "2026-08-21",
        "span_days": 60,
        "sensor": "NISAR L-SAR",
        "track": 125,
        "frame": 176,
        "polarization": "HH",
        "ground_truth": False,
    },
    "nisar_A_disp60_hp10km_mm.tif": {
        "label": "NISAR 60-day Relative LOS Displacement, Ascending (mm, high-pass 10 km, Experimental)",
        "colormap": "rdbu",
        "rescale": "-30,30",
        "description": (
            "Relative line-of-sight displacement over 60 days from the NISAR "
            "L2 GUNW ascending chain, after a 10 km spatial high-pass that "
            "removes long-wavelength ionosphere and troposphere. Red = motion "
            "away from the satellite. EXPERIMENTAL: the per-pair residual is "
            "+-20 to 40 mm and flips sign between pairs, there is no stable "
            "mineral-ground reference within 60 km, and this is NOT a "
            "subsidence rate, water-table proxy, or ground truth."
        ),
        "legend_type": "diverging",
        "legend_labels": ["-30 mm (away)", "0", "+30 mm (toward)"],
        "legend_colors": ["#d73027", "#ffffbf", "#4575b4"],
        "default_on": False,
        "tile_extra": "&tilesize=512",
        "maturity": "EXPERIMENTAL",
        "measurement_status": "relative_displacement_uncalibrated",
        "public_role": "diagnostic",
        "posting_m": 80,
        "effective_support_m": 80,
        "span_days": 60,
        "sensor": "NISAR L-SAR",
        "track": 125,
        "frame": 176,
        "polarization": "HH",
        "ground_truth": False,
    },
    "nisar_D_disp60_hp10km_mm.tif": {
        "label": "NISAR 60-day Relative LOS Displacement, Descending (mm, high-pass 10 km, Experimental)",
        "colormap": "rdbu",
        "rescale": "-30,30",
        "description": (
            "Relative line-of-sight displacement over 60 days from the NISAR "
            "L2 GUNW descending chain, after a 10 km spatial high-pass. "
            "Red = motion away from the satellite. EXPERIMENTAL, see the "
            "ascending layer: not a subsidence rate, not a water-table proxy, "
            "not ground truth."
        ),
        "legend_type": "diverging",
        "legend_labels": ["-30 mm (away)", "0", "+30 mm (toward)"],
        "legend_colors": ["#d73027", "#ffffbf", "#4575b4"],
        "default_on": False,
        "tile_extra": "&tilesize=512",
        "maturity": "EXPERIMENTAL",
        "measurement_status": "relative_displacement_uncalibrated",
        "public_role": "diagnostic",
        "posting_m": 80,
        "effective_support_m": 80,
        "span_days": 60,
        "sensor": "NISAR L-SAR",
        "polarization": "HH",
        "ground_truth": False,
    },
    "nisar_coherence_20m.tif": {
        "label": "NISAR Minimum Pair Coherence (20 m, Provisional)",
        "colormap": "viridis",
        "rescale": "0,1",
        "description": (
            "Minimum mission GUNW coherence across the four 12-day-chain "
            "pairs from 4 Jul to 14 Sep 2026 (ascending track 125, frame 176) "
            "on the native 20 m wrapped-interferogram grid. Updated 1 Oct 2026. "
            "This is a radar measurement-quality layer, not "
            "displacement, subsidence, or velocity. NISAR product maturity "
            "is PROVISIONAL: fully calibrated and partially validated."
        ),
        "legend_type": "sequential",
        "legend_labels": [
            "0 (decorrelated)",
            "0.5 (commissioning threshold)",
            "1.0 (high coherence)",
        ],
        "legend_colors": ["#440154", "#21918c", "#fde725"],
        "default_on": False,
        "tile_extra": "&nodata=-9999&tilesize=512",
        "maturity": "PROVISIONAL",
        "measurement_status": "supporting_quality_only",
        "public_role": "supporting",
        "posting_m": 20,
        "effective_support_m": 20,
        "acquisition_start": "2026-07-04",
        "acquisition_end": "2026-09-14",
        "span_days": 72,
        "pair_count": 4,
        "sensor": "NISAR L-SAR",
        "track": 125,
        "frame": 176,
        "polarization": "HH",
        "ground_truth": False,
    },
    "nisar_coherence_fine_5m.tif": {
        "label": (
            "NISAR Fine Minimum Pair Coherence "
            "(5 m posting, ~25 m support, Provisional)"
        ),
        "colormap": "viridis",
        "rescale": "0,1",
        "description": (
            "Minimum local coherence across the 22 Jun-4 Jul and "
            "4 Jul-16 Jul 2026 NISAR GSLC pairs. The raster is posted at "
            "5 m, but each estimate uses a 5 by 5 source-pixel window and "
            "therefore has approximately 25 m support. Both pairwise "
            "coherence fields passed comparison with mission GUNW "
            "(Pearson 0.913 and 0.902); the corresponding fine phase did "
            "not pass the displacement gate. This is a PROVISIONAL radar "
            "measurement-quality layer, not displacement, subsidence, "
            "velocity, or ground truth."
        ),
        "legend_type": "sequential",
        "legend_labels": [
            "0 (decorrelated)",
            "0.5 (quality threshold)",
            "1.0 (high coherence)",
        ],
        "legend_colors": ["#440154", "#21918c", "#fde725"],
        "default_on": False,
        # The compact COG stores coherence as UInt8 0..254 with a band scale
        # of 1/254. TiTiler must apply the scale before the 0..1 rescale.
        "tile_extra": "&unscale=true&tilesize=512",
        "point_unscale": True,
        "maturity": "PROVISIONAL",
        "measurement_status": "supporting_quality_only",
        "public_role": "supporting",
        "posting_m": 5,
        "effective_support_m": 25,
        "acquisition_start": "2026-06-22",
        "acquisition_end": "2026-07-16",
        "span_days": 24,
        "pair_count": 2,
        "sensor": "NISAR L-SAR",
        "track": 125,
        "frame": 176,
        "polarization": "HH",
        "ground_truth": False,
    },
}

# ---------------------------------------------------------------------------
# Application
# ---------------------------------------------------------------------------

app = FastAPI(
    title="PeatGuard Dashboard",
    description="Web map viewer for PeatGuard NISAR and peatland COG products",
    version="0.1.0",
)

# Local development accepts arbitrary origins so Expo Go can connect from a
# phone. A configured production allow-list must stay narrow; previously the
# unconditional appended "*" silently defeated PEATGUARD_API__CORS_ORIGINS.
_configured_cors_origins = os.environ.get(
    "PEATGUARD_API__CORS_ORIGINS", ""
).strip()
_cors_origins = (
    [
        origin.strip()
        for origin in _configured_cors_origins.split(",")
        if origin.strip()
    ]
    if _configured_cors_origins
    else ["*"]
)
app.add_middleware(
    CORSMiddleware,
    allow_origins=_cors_origins,
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.exception_handler(TileOutsideBounds)
async def tile_outside_bounds(
    request: Request,
    error: TileOutsideBounds,
) -> Response:
    """Return an empty tile when a map request misses a COG footprint."""
    del request, error
    return Response(status_code=204)


@app.exception_handler(PointOutsideBounds)
async def point_outside_bounds(
    request: Request,
    error: PointOutsideBounds,
) -> JSONResponse:
    """Treat a pixel query outside a COG as ordinary missing data."""
    del request, error
    return JSONResponse(
        {"detail": "Point is outside the raster footprint"},
        status_code=404,
    )


# Mount TiTiler tile endpoints at /cog
cog_tiler = TilerFactory(
    router_prefix="/cog",
    path_dependency=validated_cog_path,
)
app.include_router(cog_tiler.router, prefix="/cog", tags=["COG Tiles"])

# Jinja2 templates directory
TEMPLATES_DIR = Path(__file__).parent / "templates"
templates = Jinja2Templates(directory=str(TEMPLATES_DIR))


@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    """Serve the main dashboard page."""
    return templates.TemplateResponse(
        request,
        "index.html",
        {
            "products": PRODUCTS,
            "gcs_base_url": GCS_BASE_URL,
            "center_lat": AOI_CENTER_LAT,
            "center_lon": AOI_CENTER_LON,
        },
    )


@app.get("/health")
async def health():
    """Health check endpoint."""
    return {"status": "ok"}


# Products that are not present in the local data directory fall back to the
# GCS bucket. In local or demo mode the bucket is usually unreachable, so such
# products are hidden from the catalogue unless explicitly requested.
INCLUDE_REMOTE_PRODUCTS = os.environ.get(
    "PEATGUARD_INCLUDE_REMOTE_PRODUCTS", ""
).strip().lower() in {"1", "true", "yes"}


@app.get("/api/products")
async def list_products():
    """Return the list of available products and their tile URLs."""
    product_list = []
    for filename, meta in PRODUCTS.items():
        cog_url = resolve_cog_url(filename)
        if cog_url.startswith(GCS_BASE_URL) and not INCLUDE_REMOTE_PRODUCTS:
            continue
        tile_url = "/cog/tiles/WebMercatorQuad/{{z}}/{{x}}/{{y}}.png"
        tile_url += f"?url={cog_url}"
        if not meta.get("tile_skip_rescale"):
            tile_url += f"&rescale={meta['rescale']}"
        if not meta.get("tile_skip_named_colormap"):
            tile_url += f"&colormap_name={meta['colormap']}"
        if extra := meta.get("tile_extra"):
            tile_url += extra
        product_list.append(
            {
                "filename": filename,
                "label": meta["label"],
                "tile_url": tile_url,
                "cog_url": cog_url,
                "description": meta["description"],
                "legend_type": meta["legend_type"],
                "legend_labels": meta["legend_labels"],
                "legend_colors": meta["legend_colors"],
                "default_on": meta["default_on"],
                **{
                    key: meta[key]
                    for key in (
                        "maturity",
                        "measurement_status",
                        "public_role",
                        "posting_m",
                        "effective_support_m",
                        "acquisition_start",
                        "acquisition_end",
                        "span_days",
                        "pair_count",
                        "sensor",
                        "track",
                        "frame",
                        "polarization",
                        "ground_truth",
                    )
                    if key in meta
                },
            }
        )
    return {"products": product_list, "center": [AOI_CENTER_LAT, AOI_CENTER_LON]}
