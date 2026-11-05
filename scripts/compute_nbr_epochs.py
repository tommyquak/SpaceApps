"""Multi-temporal NBR predictor (non-GEE path).

Builds dry-season Landsat Collection-2 L2 surface-reflectance NBR composites for
the paper's epochs (1990-2022) over the Sebangau AOI, via Microsoft Planetary
Computer (open, anonymous signing). Composite-based stand-in for the paper's
LandTrendr-fitted NBR (which needs Earth Engine). Outputs nbr_<year>.tif.

  NBR = (NIR - SWIR2) / (NIR + SWIR2)   [Landsat common bands nir08, swir22]
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import planetary_computer as pc
import pystac_client
import odc.stac
import rioxarray  # noqa: F401  (registers .rio accessor)

CAT = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)
W, S, E, N = 114.20, -3.10, 114.83, -2.38   # AOI + small buffer
RES = 0.0003                                  # ~33 m in EPSG:4326
EPOCHS = [1990, 1995, 2000, 2005, 2010, 2015, 2020, 2022]
OUT = Path("grid_products/external/nbr")
OUT.mkdir(parents=True, exist_ok=True)


def qa_clear(qa):
    # Landsat C2 L2 QA_PIXEL: bit0 fill,1 dilated cloud,2 cirrus,3 cloud,4 shadow
    bad = (qa & 1) | (qa & (1 << 1)) | (qa & (1 << 2)) | (qa & (1 << 3)) | (qa & (1 << 4))
    return bad == 0


def epoch_nbr(yr):
    def _search(a, b):
        return list(CAT.search(collections=["landsat-c2-l2"], bbox=[W, S, E, N],
                    datetime=f"{a}/{b}", query={"eo:cloud_cover": {"lt": 80}}).items())
    items = _search(f"{yr}-05-01", f"{yr}-10-31")
    if len(items) < 3:
        items = _search(f"{yr-1}-05-01", f"{yr+1}-10-31")  # widen if sparse
    ds = odc.stac.load(items, bands=["nir08", "swir22", "qa_pixel"], bbox=[W, S, E, N],
                       crs="EPSG:4326", resolution=RES, groupby="solar_day",
                       chunks={"x": 1024, "y": 1024})
    clear = qa_clear(ds.qa_pixel.astype("uint16"))
    nir = (ds.nir08.where(ds.nir08 > 0) * 0.0000275 - 0.2).where(clear)
    swir = (ds.swir22.where(ds.swir22 > 0) * 0.0000275 - 0.2).where(clear)
    nbr = (nir - swir) / (nir + swir)
    comp = nbr.median(dim="time", skipna=True).compute().rio.write_crs("EPSG:4326")
    comp.rio.to_raster(OUT / f"nbr_{yr}.tif", compress="deflate")
    valid = float(np.isfinite(comp.values).mean()) * 100
    print(f"  nbr_{yr}.tif: {len(items)} scenes -> {tuple(comp.shape)}  {valid:.0f}% cloud-free valid", flush=True)


if __name__ == "__main__":
    years = [int(sys.argv[1])] if len(sys.argv) > 1 else EPOCHS
    for y in years:
        epoch_nbr(y)
    print("NBR_DONE", flush=True)
