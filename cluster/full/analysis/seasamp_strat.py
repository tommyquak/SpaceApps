"""Is there a clean subset? Amplitude and month of minimum on Tay's footprint,
stratified by the fit-residual RMS quartile, against the noise expectation for
a sinusoid fitted to N epochs of white noise, amp_noise ~ RMS * sqrt(2 / N)."""
import os
import numpy as np
from osgeo import gdal, osr
gdal.UseExceptions()
O = os.path.expanduser("~/full/analysis_out"); TAY = os.path.expanduser("~/full/tay/InSAR_LongTermRate.tif")
amp = gdal.Open(f"{O}/bw60_seasonal_amp_mm.tif"); A = amp.GetRasterBand(1).ReadAsArray()
dsM = gdal.Open(f"{O}/bw60_seasonal_min_month.tif"); M = dsM.GetRasterBand(1).ReadAsArray()
dsR = gdal.Open(f"{O}/bw60_fit_resid_rms_mm.tif"); R = dsR.GetRasterBand(1).ReadAsArray()
mem = gdal.GetDriverByName("MEM").Create("", amp.RasterXSize, amp.RasterYSize, 1, gdal.GDT_Float32); mem.SetGeoTransform(amp.GetGeoTransform()); mem.SetProjection(amp.GetProjection()); mem.GetRasterBand(1).Fill(float("nan"))
src = gdal.Open(TAY); gdal.ReprojectImage(src, mem, src.GetProjection(), amp.GetProjection(), gdal.GRA_Average)
t = mem.GetRasterBand(1).ReadAsArray(); t[t == 0] = np.nan; strip = np.isfinite(t) & np.isfinite(A) & np.isfinite(R)
N = 93
print(f"footprint px {strip.sum()}; residual RMS quartiles: {np.percentile(R[strip], [25, 50, 75]).round(0)} mm")
q = np.percentile(R[strip], [0, 25, 50, 75, 100])
for k in range(4):
    m = strip & (R >= q[k]) & (R < q[k + 1] + (1e-6 if k == 3 else 0))
    a = A[m]; r = R[m]; mm = M[m]
    noise = np.median(r) * np.sqrt(2.0 / N)
    hist = np.histogram(mm, bins=np.arange(1, 14))[0]; peak = int(np.argmax(hist)) + 1; frac_peak = hist.max() / max(hist.sum(), 1)
    print(f"  Q{k+1} RMS {q[k]:.0f}-{q[k+1]:.0f} mm: n={m.sum():6d} amplitude median {np.median(a):.1f} mm, noise expectation {noise:.1f} mm, ratio {np.median(a)/noise:.2f}; month of minimum mode {peak} ({100*frac_peak:.0f}% of px), median {np.median(mm):.1f}")
print("STRAT_DONE")
