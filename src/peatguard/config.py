"""Pipeline configuration management.

Loads settings from YAML config files with environment variable overrides.
All pipeline parameters are centralized here to eliminate hardcoded values
scattered across processing scripts.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Literal, Optional

import yaml
from pydantic import BaseModel, Field, field_validator, model_validator

logger = logging.getLogger(__name__)


class AOIConfig(BaseModel):
    """Area of interest definition."""

    bbox: list[float] = Field(
        description="Bounding box as [west, south, east, north] in geographic coordinates"
    )
    epsg: int = Field(default=32750, description="Target CRS EPSG code")

    @field_validator("bbox")
    @classmethod
    def validate_bbox(cls, v: list[float]) -> list[float]:
        if len(v) != 4:
            raise ValueError("bbox must have exactly 4 values: [west, south, east, north]")
        west, south, east, north = v
        if west >= east:
            raise ValueError(f"west ({west}) must be less than east ({east})")
        if south >= north:
            raise ValueError(f"south ({south}) must be less than north ({north})")
        return v


class Sentinel1Config(BaseModel):
    """Sentinel-1 acquisition parameters."""

    platform: str = "SENTINEL-1"
    polarization: str = "VV"
    subswath: str = "IW3"
    processing_level_slc: str = "SLC"
    processing_level_grd: str = "GRD_HD"
    beam_mode: str = "IW"
    min_temporal_gap_days: int = 14
    max_temporal_baseline_days: int = 48


class NISARConfig(BaseModel):
    """NISAR L-band acquisition parameters.

    NISAR uses L-band (24cm wavelength) which penetrates tropical forest
    canopy, enabling subsidence measurement under intact peat forest where
    C-band Sentinel-1 loses coherence.

    BETA products have been available since February 2026. PROVISIONAL
    products have been available since July 2026.
    """

    platform: str = "NISAR"
    polarization: str = "HH"
    beam_mode: str = "LSAR"
    processing_level: str = "GUNW"
    processing_level_gunw: str = "GUNW"
    processing_level_gslc: str = "GSLC"
    processing_level_rslc: str = "RSLC"
    min_temporal_gap_days: int = 12  # NISAR repeat cycle is 12 days
    max_temporal_baseline_days: int = 48
    # Nominal Frequency A value for configuration-only paths. Product metadata
    # must override this value whenever a NISAR HDF5 is available.
    wavelength_m: float = 0.2419632429378531
    orbit_height_m: float = 747000.0  # ~747 km orbit altitude
    range_pixel_size_m: float = 7.5  # approximate, varies by mode
    azimuth_pixel_size_m: float = 6.0  # approximate, varies by mode
    incidence_deg: float = 33.9  # typical mid-swath incidence angle


class ProcessingConfig(BaseModel):
    """InSAR and backscatter processing parameters."""

    resolution_m: float = 10.0
    coherence_threshold: float = 0.3
    wavelength_m: float = 0.056
    incidence_deg: float = 37.0
    dem_source: str = "SRTM"
    speckle_filter: str = "lee_sigma"
    speckle_window_size: int = 7
    goldstein_alpha: float = 0.8
    snaphu_cost_mode: str = "DEFO"
    do_ion: bool = True  # Ionospheric correction via split-spectrum (recommended for equatorial)


class ClassificationConfig(BaseModel):
    """Subsidence classification thresholds in mm/yr.

    Natural peat accumulation is ~1 mm/yr. Any subsidence below
    -5 mm/yr on peat indicates net carbon loss, so "stable" is
    restricted to the +/-5 mm/yr natural variability envelope.
    """

    severe_threshold: float = -50.0
    active_drying_threshold: float = -20.0
    moderate_drying_threshold: float = -5.0
    stable_threshold: float = 5.0


class WaterMaskConfig(BaseModel):
    """Water body detection parameters.

    Controls VV backscatter-based water masking to prevent false
    subsidence detections over water bodies. Water shows specular
    reflection (very low backscatter) that can mimic low coherence.
    """

    enabled: bool = True
    vv_threshold_db: float = Field(
        default=-18.0,
        description=(
            "VV backscatter threshold in dB. Pixels below this are water. "
            "-18 dB catches the river and major water bodies."
        ),
    )
    min_water_size_pixels: int = Field(
        default=25,
        description=(
            "Minimum connected component size (pixels) to retain as water. "
            "At 10m resolution, 25 pixels = 2500 sq m."
        ),
    )
    closing_radius: int = Field(
        default=3,
        description="Morphological closing radius to fill gaps in water bodies.",
    )
    exclude_canals: bool = Field(
        default=False,
        description=(
            "Remove canal pixels from water mask to avoid interfering "
            "with canal_detect.py results. Disabled by default: the canal "
            "mask is too broad (19.8% of scene) and removes all water."
        ),
    )


class PeatMaskConfig(BaseModel):
    """Peat extent mapping and masking parameters.

    Downloads peat polygons from the WRI/GFW Indonesia Peatlands dataset
    and rasterizes them to produce peat-masked subsidence products.
    Peat depth is estimated from distance-to-edge using the tropical
    peatland dome model (Page et al., 2006).
    """

    enabled: bool = Field(
        default=True,
        description="Enable peat extent mapping and masking.",
    )
    shallow_edge_m: float = Field(
        default=500.0,
        description=(
            "Distance from peat edge for shallow peat classification (meters). "
            "Pixels within this distance of the peat boundary are classified shallow."
        ),
    )
    moderate_edge_m: float = Field(
        default=1500.0,
        description=(
            "Distance from peat edge for moderate peat classification (meters). "
            "Pixels between shallow_edge_m and this distance are moderate; "
            "beyond this distance are classified as deep peat."
        ),
    )


class RiskScoreConfig(BaseModel):
    """Transparent project-priority parameters.

    Hooijer et al. (2012) supports the process link between drainage,
    long-term subsidence, and carbon loss after accounting for contributing
    processes. It does not prescribe this model's linear distance response,
    influence radius, or component weights.
    """

    proximity_weight: float = Field(
        default=0.45,
        description="Weight for canal proximity component (indirect proxy).",
    )
    subsidence_weight: float = Field(
        default=0.55,
        description="Weight for subsidence velocity component (direct measurement).",
    )
    max_influence_m: float = Field(
        default=600.0,
        description=(
            "Project-assumed canal influence radius in meters. Treat as a "
            "sensitivity parameter until calibrated with field hydrology."
        ),
    )
    severe_velocity_mm_yr: float = Field(
        default=-40.0,
        description=(
            "Velocity at which subsidence risk saturates (mm/yr). "
            "This is a project severity threshold, not a literature calibration."
        ),
    )
    coherence_threshold: float = Field(
        default=0.7,
        description=(
            "Risk-specific coherence cutoff: pixels below this value are "
            "excluded from the canal_risk product. Stricter than "
            "processing.coherence_threshold (which gates SBAS inversion "
            "quality), because the risk map is a decision-support product "
            "that should not show confident reds on unreliable phase."
        ),
    )


class FusionConfig(BaseModel):
    """Multi-sensor fusion parameters for combining C-band and L-band velocity.

    Fusion is only performed when both Sentinel-1 (C-band) and NISAR
    (L-band) velocity products are available. The coherence-weighted
    averaging strategy naturally favors the sensor with better signal
    quality at each pixel.
    """

    enabled: bool = Field(
        default=False,
        description=(
            "Enable multi-sensor fusion. Only takes effect when both "
            "C-band and L-band velocity products exist."
        ),
    )
    c_band_weight_boost: float = Field(
        default=1.0,
        description=(
            "Multiplicative boost for C-band coherence weight in dual-coverage "
            "areas. 1.0 = pure coherence weighting. Values > 1.0 increase "
            "C-band influence, which may be desirable because C-band has "
            "higher deformation sensitivity (shorter wavelength) over cleared land."
        ),
    )
    min_coherence: float = Field(
        default=0.3,
        description=(
            "Minimum temporal coherence for a sensor to contribute to the "
            "fused velocity. Pixels below this threshold are treated as "
            "invalid for that sensor. Matches the pipeline coherence_threshold."
        ),
    )


class ForestChangeConfig(BaseModel):
    """Hansen Global Forest Change data integration parameters.

    Downloads Hansen et al. (2013) "lossyear" and "treecover2000" tiles
    to provide temporal context for subsidence analysis. The clearing year
    directly relates to subsidence rate via Hooijer et al. (2012): recently
    cleared peatland subsides faster than areas cleared decades ago.
    """

    enabled: bool = Field(
        default=True,
        description="Enable Hansen Global Forest Change data download and analysis.",
    )
    reference_year: int = Field(
        default=2024,
        description=(
            "Reference year for computing time-since-clearing. "
            "Should match the InSAR observation period."
        ),
    )
    hansen_version: str = Field(
        default="GFC-2023-v1.11",
        description="Hansen GFC dataset version string for download URL construction.",
    )


class ValidationConfig(BaseModel):
    """Cross-validation settings for independent verification of InSAR velocity.

    Computes spatial correlations between subsidence velocity and
    independent datasets (backscatter, coherence, canal distance, NDVI)
    to provide evidence that the InSAR signal reflects real ground motion.
    """

    enabled: bool = Field(
        default=True,
        description="Enable cross-validation after analysis stage.",
    )
    attempt_ndvi_download: bool = Field(
        default=True,
        description=(
            "Attempt to download Sentinel-2 NDVI composite from Planetary Computer. "
            "Requires pystac-client and planetary-computer packages. If unavailable, "
            "validation proceeds with SAR-only products."
        ),
    )
    ndvi_max_cloud_cover: int = Field(
        default=30,
        description="Maximum cloud cover percentage for Sentinel-2 scene selection.",
    )
    ndvi_max_scenes: int = Field(
        default=20,
        description="Maximum number of Sentinel-2 scenes for the NDVI composite.",
    )
    coherence_threshold: float = Field(
        default=0.5,
        description=(
            "Coherence threshold for zonal analysis. Pixels above this "
            "are classified as 'cleared'; below as 'forested'."
        ),
    )


class ExportConfig(BaseModel):
    """Output format settings for ArcGIS Pro compatibility."""

    cog_blocksize: int = 512
    cog_compression: str = "DEFLATE"
    overview_levels: list[int] = Field(default=[2, 4, 8, 16])
    nodata: float = -9999.0


class StorageConfig(BaseModel):
    """Local and cloud storage paths."""

    output_dir: Path = Path("./output")
    scratch_dir: Path = Path("./scratch")
    catalog_db: Path = Path("./catalog.db")
    gcs_bucket: str = ""
    gcs_prefix: str = "peatguard-data"

    @field_validator("output_dir", "scratch_dir", "catalog_db", mode="after")
    @classmethod
    def resolve_paths(cls, v: Path) -> Path:
        """Ensure all paths are absolute to avoid cwd-relative issues."""
        return v.resolve()


class MintPyNetworkConfig(BaseModel):
    """MintPy network modification settings."""

    coherence_based: bool = True
    min_coherence: float = 0.3
    start_date: str = Field(
        default="20240101",
        description=(
            "MintPy network.startDate filter in YYYYMMDD form. Empty string "
            "disables the filter. Set to restrict the SBAS network "
            "to pairs whose reference date is on or after this cutoff -- "
            "useful for isolating a single year from a multi-year ingest to "
            "reduce temporal decorrelation in C-band."
        ),
    )


class MintPyConfig(BaseModel):
    """MintPy time-series processing settings."""

    reference_lalo: list[float] = Field(
        default_factory=list,
        description=(
            "Fixed reference point as [lat, lon] for velocity baseline. "
            "Must lie inside aoi.bbox. Empty list (default) falls back to "
            "MintPy coherence-based auto-selection; a site opts in only "
            "with a surveyed in-AOI stable point."
        ),
    )
    reference_min_coherence: float = Field(
        default=0.7,
        description="Minimum coherence for reference point (auto-select fallback).",
    )
    reference_yx: list[int] = Field(
        default=[92, 502],
        description=(
            "Fixed reference point as [row, col] pixel coordinates in the "
            "radar grid. Takes priority over reference_lalo. Use this when "
            "a lat/lon pin is unreliable because lat.rdr/lon.rdr lookup is "
            "off (synthetic placeholder, etc). Empty list = fall through to "
            "reference_lalo or auto-select. [92, 502] is the high-coherence "
            "pixel validated in the 2026-04-19 run."
        ),
    )
    tropospheric_correction: str = "pyaps"
    unwrap_error_correction: str = "phase_closure"
    era5_hour: int = Field(
        default=22,
        ge=0,
        le=23,
        description=(
            "UTC hour for the ERA5 pre-download when the prepared "
            "ifgramStack carries no CENTER_LINE_UTC attribute. PyAPS "
            "selects the GRIB hour from CENTER_LINE_UTC, so this must "
            "match the acquisition time. 22 is the 2026 Sentinel-1 "
            "descending pass over Kalimantan (~22:00 UTC)."
        ),
    )
    network_modification: MintPyNetworkConfig = Field(default_factory=MintPyNetworkConfig)

    @field_validator("reference_lalo")
    @classmethod
    def validate_reference_lalo(cls, v: list[float]) -> list[float]:
        if len(v) == 0:
            return v
        if len(v) != 2:
            raise ValueError("reference_lalo must be [lat, lon] (2 values) or [] for auto")
        lat, lon = v
        if not (-90 <= lat <= 90):
            raise ValueError(f"reference latitude {lat} out of range [-90, 90]")
        if not (-180 <= lon <= 180):
            raise ValueError(f"reference longitude {lon} out of range [-180, 180]")
        return v


class PeatGuardConfig(BaseModel):
    """Root configuration for the PeatGuard pipeline."""

    sensor: Literal["sentinel1", "nisar"] = Field(
        default="sentinel1",
        description=(
            "Active sensor selection. 'sentinel1' uses C-band TOPS mode "
            "(topsApp). 'nisar' uses L-band stripmap mode (stripmapApp)."
        ),
    )
    aoi: AOIConfig
    sentinel1: Sentinel1Config = Field(default_factory=Sentinel1Config)
    nisar: NISARConfig = Field(default_factory=NISARConfig)
    processing: ProcessingConfig = Field(default_factory=ProcessingConfig)
    classification: ClassificationConfig = Field(default_factory=ClassificationConfig)
    water_mask: WaterMaskConfig = Field(default_factory=WaterMaskConfig)
    peat_mask: PeatMaskConfig = Field(default_factory=PeatMaskConfig)
    risk_score: RiskScoreConfig = Field(default_factory=RiskScoreConfig)
    forest_change: ForestChangeConfig = Field(default_factory=ForestChangeConfig)
    fusion: FusionConfig = Field(default_factory=FusionConfig)
    validation: ValidationConfig = Field(default_factory=ValidationConfig)
    export: ExportConfig = Field(default_factory=ExportConfig)
    storage: StorageConfig = Field(default_factory=StorageConfig)
    mintpy: MintPyConfig = Field(default_factory=MintPyConfig)

    @property
    def is_nisar(self) -> bool:
        """Return True if the active sensor is NISAR."""
        return self.sensor == "nisar"

    @property
    def active_sensor_config(self) -> Sentinel1Config | NISARConfig:
        """Return the config object for the active sensor."""
        return self.nisar if self.is_nisar else self.sentinel1

    @model_validator(mode="after")
    def validate_reference_lalo_in_aoi(self) -> "PeatGuardConfig":
        """Reject a velocity reference point outside the AOI.

        MintPy references the whole velocity field to this point; a point
        outside aoi.bbox lies outside the clipped exported products and
        cannot be a valid baseline for them.
        """
        ref = self.mintpy.reference_lalo
        if ref and len(ref) == 2:
            lat, lon = ref
            west, south, east, north = self.aoi.bbox
            if not (south <= lat <= north and west <= lon <= east):
                raise ValueError(
                    f"mintpy.reference_lalo [{lat}, {lon}] is outside "
                    f"aoi.bbox {self.aoi.bbox}; use an in-AOI point or [] "
                    "for MintPy coherence-based auto-selection"
                )
        return self


def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge override dict into base dict."""
    merged = base.copy()
    for key, value in override.items():
        if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def _resolve_field_path(parts: list[str]) -> bool:
    """Return True if parts names a field reachable from PeatGuardConfig.

    Each intermediate part must be a nested BaseModel field; the final
    part must be a declared field on that model.
    """
    model: type[BaseModel] = PeatGuardConfig
    for depth, part in enumerate(parts):
        field = model.model_fields.get(part)
        if field is None:
            return False
        if depth == len(parts) - 1:
            return True
        annotation = field.annotation
        if isinstance(annotation, type) and issubclass(annotation, BaseModel):
            model = annotation
        else:
            return False
    return False


def _apply_env_overrides(data: dict) -> dict:
    """Apply environment variable overrides.

    Environment variables use double underscore as the section separator,
    with nested sections chained by further double underscores:
        PEATGUARD_<SECTION>__<KEY>=value
        PEATGUARD_<SECTION>__<SUBSECTION>__<KEY>=value

    For example: PEATGUARD_STORAGE__GCS_BUCKET=my-bucket or
    PEATGUARD_MINTPY__NETWORK_MODIFICATION__MIN_COHERENCE=0.4.
    Single underscores within section or key names are preserved.

    Keys that do not resolve to a declared config field are ignored with
    a warning, so a typo cannot silently leave a parameter unchanged.
    Keys without a double underscore (e.g. PEATGUARD_ALLOW_CBAND_FUSION)
    are not config overrides and are skipped.
    """
    env_prefix = "PEATGUARD_"
    unmapped: list[str] = []
    for env_key, env_value in os.environ.items():
        if not env_key.startswith(env_prefix):
            continue
        remainder = env_key[len(env_prefix):]
        if "__" not in remainder:
            continue
        parts = [part.lower() for part in remainder.split("__")]
        if not _resolve_field_path(parts):
            unmapped.append(env_key)
            continue
        target = data
        for part in parts[:-1]:
            node = target.get(part)
            if not isinstance(node, dict):
                node = {}
                target[part] = node
            target = node
        target[parts[-1]] = env_value
    if unmapped:
        logger.warning(
            "Ignoring PEATGUARD_ environment overrides that do not match "
            "any known config field: %s",
            ", ".join(sorted(unmapped)),
        )
    return data


def load_config(
    config_path: Optional[Path] = None,
    override_path: Optional[Path] = None,
) -> PeatGuardConfig:
    """Load pipeline configuration from YAML files.

    Args:
        config_path: Path to the base config YAML. Defaults to config/default.yaml
            relative to the package root.
        override_path: Optional path to an override YAML that is merged on top
            of the base config.

    Returns:
        Validated PeatGuardConfig instance.
    """
    if config_path is None:
        # Look for bundled config inside the package first, then fall back to source tree
        pkg_config = Path(__file__).parent / "data" / "default.yaml"
        src_config = Path(__file__).parent.parent.parent / "config" / "default.yaml"
        config_path = pkg_config if pkg_config.exists() else src_config

    with open(config_path) as f:
        data = yaml.safe_load(f)

    if override_path is not None:
        with open(override_path) as f:
            override_data = yaml.safe_load(f)
        if override_data:
            data = _deep_merge(data, override_data)

    data = _apply_env_overrides(data)
    return PeatGuardConfig(**data)
