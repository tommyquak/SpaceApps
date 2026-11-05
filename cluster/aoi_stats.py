import glob, os
import numpy as np
from osgeo import gdal
gdal.UseExceptions()
rows = []
for cor in sorted(glob.glob(os.path.expanduser("~/pilot/pairs/*/") + "*_5rlks_28alks.cor.geo.aoi.tif")):
    pair = cor.split("/pairs/")[1].split("/")[0]
    d0, d1 = pair.split("_")
    from datetime import date
    span = (date.fromisoformat(d1) - date.fromisoformat(d0)).days
    ds = gdal.Open(cor)
    band = ds.GetRasterBand(ds.RasterCount)
    a = band.ReadAsArray()
    gt = ds.GetGeoTransform()
    x0 = int((114.21 - gt[0]) / gt[1]); x1 = int((114.83 - gt[0]) / gt[1])
    y0 = int((-2.42 - gt[3]) / gt[5]); y1 = int((-3.0 - gt[3]) / gt[5])
    aoi = a[min(y0,y1):max(y0,y1), x0:x1]
    v = aoi[aoi > 0]
    if v.size == 0:
        print(f"{pair} span={span}d  NO VALID AOI DATA"); continue
    print(f"{pair} span={span}d  aoi_median={np.median(v):.3f}  frac>=0.3={float((v>=0.3).mean()):.1%}  frac>=0.5={float((v>=0.5).mean()):.1%}")
