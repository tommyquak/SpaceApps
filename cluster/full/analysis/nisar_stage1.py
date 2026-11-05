"""NISAR GUNW stage 1: what is in the files, and is the peat coherent at 12 days?

Walks each provisional GUNW HDF5, locates the coherence, unwrapped-phase and
connected-component layers by name (the provisional layout is not assumed),
reports grid CRS and extent, and measures coherence over the whole frame and
inside the Tay strip, for direct comparison with ALOS-2 at 14 days (0.22).
"""
import glob
import os
import sys
import numpy as np
import h5py
from osgeo import gdal, osr

D = os.path.expanduser(os.environ.get("GUNW_DIR", "/mnt/scratch/t/tommy/nisar/gunw"))
TAY = os.path.expanduser("~/full/tay/InSAR_LongTermRate.tif")
gdal.UseExceptions()


def find(h, names):
    hits = {}
    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            base = name.split("/")[-1]
            for n in names:
                if base == n and n not in hits:
                    hits[n] = name
    h.visititems(visit)
    return hits


def tree(h, maxdepth=6):
    out = []
    def visit(name, obj):
        if name.count("/") <= maxdepth:
            if isinstance(obj, h5py.Dataset):
                out.append(f"  {name}  {obj.shape} {obj.dtype}")
            elif name.count("/") <= 4:
                out.append(f"  {name}/")
    h.visititems(visit)
    return out


files = sorted(glob.glob(os.path.join(D, "*.h5")))
print(f"{len(files)} GUNW files under {D}")
if not files:
    sys.exit(1)

WANT = ["coherenceMagnitude", "unwrappedPhase", "connectedComponents", "ionospherePhaseScreen",
        "xCoordinates", "yCoordinates", "projection", "epsg", "xCoordinateSpacing", "yCoordinateSpacing",
        "troposphericPhaseScreen", "waterMask", "layoverShadowMask"]
first = True
for f in files:
    with h5py.File(f, "r") as h:
        if first:
            print("\n--- HDF5 layout (first file, depth<=6) ---")
            for line in tree(h)[:120]:
                print(line)
            first = False
        hits = find(h, WANT)
        print(f"\n=== {os.path.basename(f)[:60]}")
        for k, v in hits.items():
            print(f"   {k:24s} -> {v}")
        if "coherenceMagnitude" not in hits:
            print("   no coherence layer found by name; skipping stats")
            continue
        coh = h[hits["coherenceMagnitude"]][:].astype(np.float32)
        ok = np.isfinite(coh) & (coh > 0)
        print(f"   coherence grid {coh.shape}; valid {ok.mean()*100:.1f}%; mean {coh[ok].mean():.3f} median {np.median(coh[ok]):.3f}"
              f"; >0.3: {100*(coh[ok] > 0.3).mean():.1f}%  >0.5: {100*(coh[ok] > 0.5).mean():.1f}%")
        # georeference
        try:
            x = h[hits["xCoordinates"]][:]; y = h[hits["yCoordinates"]][:]
            epsg = None
            for key in ("epsg", "projection"):
                if key in hits:
                    v = h[hits[key]][()]
                    try:
                        epsg = int(v)
                    except Exception:
                        s = v.decode() if isinstance(v, bytes) else str(v)
                        epsg = int(s.split(":")[-1]) if ":" in s else None
            print(f"   grid x {x.min():.1f}..{x.max():.1f} ({x.size}), y {y.min():.1f}..{y.max():.1f} ({y.size}), spacing {abs(x[1]-x[0]):.1f} m, epsg {epsg}")
            if epsg:
                # Tay raster onto this grid, then coherence inside the strip
                ny, nx = coh.shape
                mem = gdal.GetDriverByName("MEM").Create("", nx, ny, 1, gdal.GDT_Float32)
                mem.SetGeoTransform((float(x.min()), float(x[1]-x[0]), 0.0, float(y.max()), 0.0, float(y[1]-y[0]) if y[1] < y[0] else -abs(float(y[1]-y[0]))))
                srs = osr.SpatialReference(); srs.ImportFromEPSG(epsg); mem.SetProjection(srs.ExportToWkt())
                mem.GetRasterBand(1).Fill(float("nan"))
                src = gdal.Open(TAY)
                gdal.ReprojectImage(src, mem, src.GetProjection(), srs.ExportToWkt(), gdal.GRA_Average)
                tay = mem.GetRasterBand(1).ReadAsArray().astype(np.float32); tay[tay == 0] = np.nan
                strip = np.isfinite(tay) & ok
                if strip.sum() > 100:
                    print(f"   Tay strip: {strip.sum()} px; coherence mean {coh[strip].mean():.3f} median {np.median(coh[strip]):.3f}"
                          f"; >0.3: {100*(coh[strip] > 0.3).mean():.1f}%  >0.5: {100*(coh[strip] > 0.5).mean():.1f}%")
                else:
                    print(f"   Tay strip not inside this frame ({strip.sum()} px)")
        except Exception as e:
            print("   georeference step failed:", repr(e))
print("\nSTAGE1_DONE")
