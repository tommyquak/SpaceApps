"""Apply Tay et al. 2025's own quality metric to our stack.

Tay excludes pixels whose MEAN phase standard deviation exceeds 0.7 cm, with
sigma_phi^2 estimated from coherence gamma per pixel per interferogram:
    sigma_phi = (lambda / 4 pi) * sqrt((1 - gamma^2) / (2 N gamma^2))
[E] Cramer-Rao bound for multilooked phase (Rodriguez & Martin 1992; Hanssen 2001).
This is NOT MintPy temporalCoherence. We have been gating on the wrong metric.

Reports the retained fraction under Tay's rule with nominal looks (what Tay
almost certainly used) and with our measured effective looks, plus the overlap
of the Tay-passing set with Tay's published data and with our tempCoh gates.
"""
import os
import numpy as np
import h5py
from osgeo import gdal, osr

W = os.path.expanduser("~/full/mintpy_tay")
TAY = os.path.expanduser("~/full/tay/InSAR_LongTermRate.tif")
OUT = os.path.expanduser("~/full/sigphi")
os.makedirs(OUT, exist_ok=True)
N_NOMINAL = 17 * 56
N_EFFECTIVE = float(os.environ.get("N_EFF", "41"))
THRESH_CM = 0.7
ROWS = 64

stack = os.path.join(W, "inputs", "ifgramStack.h5")
with h5py.File(os.path.join(W, "velocity.h5"), "r") as f:
    att = dict(f.attrs)
wl = float(att.get("WAVELENGTH", 0.2384))
print(f"wavelength {wl:.5f} m  (from attrs: {'WAVELENGTH' in att})")
k = wl / (4.0 * np.pi) * 100.0  # metres -> cm

with h5py.File(stack, "r") as f:
    coh = f["coherence"]
    npair, ny, nx = coh.shape
    print(f"pairs={npair} grid={ny}x{nx}")
    mean_sig_nom = np.full((ny, nx), np.nan, dtype=np.float64)
    mean_sig_eff = np.full((ny, nx), np.nan, dtype=np.float64)
    mean_gamma = np.full((ny, nx), np.nan, dtype=np.float64)
    for r0 in range(0, ny, ROWS):
        r1 = min(ny, r0 + ROWS)
        g = coh[:, r0:r1, :].astype(np.float64)
        g = np.clip(g, 1e-3, 0.9999)
        ratio = (1.0 - g * g) / (g * g)
        mean_sig_nom[r0:r1] = k * np.mean(np.sqrt(ratio / (2.0 * N_NOMINAL)), axis=0)
        mean_sig_eff[r0:r1] = k * np.mean(np.sqrt(ratio / (2.0 * N_EFFECTIVE)), axis=0)
        mean_gamma[r0:r1] = np.mean(g, axis=0)
        if (r0 // ROWS) % 5 == 0:
            print(f"  rows {r0}-{r1} done", flush=True)

np.save(os.path.join(OUT, "mean_sigphi_nominal_cm.npy"), mean_sig_nom.astype(np.float32))
np.save(os.path.join(OUT, "mean_sigphi_effective_cm.npy"), mean_sig_eff.astype(np.float32))
np.save(os.path.join(OUT, "mean_gamma.npy"), mean_gamma.astype(np.float32))

with h5py.File(os.path.join(W, "temporalCoherence.h5"), "r") as f:
    tc = f["temporalCoherence"][:].astype(np.float64)

x0 = float(att["X_FIRST"]); dx = float(att["X_STEP"])
y0 = float(att["Y_FIRST"]); dy = float(att["Y_STEP"])
mem = gdal.GetDriverByName("MEM").Create("", nx, ny, 1, gdal.GDT_Float32)
mem.SetGeoTransform((x0, dx, 0.0, y0, 0.0, dy))
srs = osr.SpatialReference(); srs.ImportFromEPSG(4326)
mem.SetProjection(srs.ExportToWkt())
mem.GetRasterBand(1).Fill(float("nan"))
src = gdal.Open(TAY)
gdal.ReprojectImage(src, mem, src.GetProjection(), srs.ExportToWkt(), gdal.GRA_Average)
tay = mem.GetRasterBand(1).ReadAsArray().astype(np.float64)
tay[tay == 0] = np.nan
tay_ok = np.isfinite(tay)

total = ny * nx
print(f"\nmean gamma over frame: {np.nanmean(mean_gamma):.4f}  median {np.nanmedian(mean_gamma):.4f}")
print(f"mean gamma in Tay-data strip: {np.nanmean(mean_gamma[tay_ok]):.4f}")
print(f"mean gamma outside strip:     {np.nanmean(mean_gamma[~tay_ok]):.4f}")

for label, sig in (("NOMINAL N=952", mean_sig_nom), (f"EFFECTIVE N={N_EFFECTIVE:g}", mean_sig_eff)):
    keep = sig <= THRESH_CM
    print(f"\n[{label}] Tay rule sigma_phi <= {THRESH_CM} cm")
    print(f"   retained {keep.sum():8d} px = {100*keep.sum()/total:6.2f}% of frame")
    print(f"   median sigma_phi {np.nanmedian(sig):.3f} cm;  p10 {np.nanpercentile(sig,10):.3f}  p90 {np.nanpercentile(sig,90):.3f}")
    both = keep & tay_ok
    print(f"   retained AND Tay has data: {both.sum():7d} px "
          f"({100*both.sum()/max(tay_ok.sum(),1):.1f}% of Tay strip passes our Tay-rule)")
    for t in (0.3, 0.5, 0.7):
        print(f"   of Tay-rule-retained px, fraction with tempCoh>={t}: "
              f"{100*np.mean(tc[keep] >= t):5.1f}%")
    np.save(os.path.join(OUT, f"mask_tayrule_{label.split()[0].lower()}.npy"), keep)

# Correlation with Tay under Tay's own mask (nominal), on the strip.
with h5py.File(os.path.join(W, "velocity.h5"), "r") as f:
    v = f["velocity"][:].astype(np.float64) * 1000.0
with h5py.File(os.path.join(W, "velocityClosureCor.h5"), "r") as f:
    vz = f["velocity"][:].astype(np.float64) * 1000.0
keep = (mean_sig_nom <= THRESH_CM) & tay_ok & np.isfinite(v) & np.isfinite(vz)
print(f"\nTay-rule (nominal) AND Tay data: n={keep.sum()}")
if keep.sum() > 100:
    print(f"   r(pre-Zheng, Tay)  = {np.corrcoef(v[keep], tay[keep])[0,1]:+.3f}")
    print(f"   r(post-Zheng, Tay) = {np.corrcoef(vz[keep], tay[keep])[0,1]:+.3f}")
    print(f"   our std pre {v[keep].std():.2f} post {vz[keep].std():.2f}; Tay std {1000*tay[keep].std():.2f} mm/yr")
print("SIGPHI_DONE")
