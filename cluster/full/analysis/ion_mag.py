import glob, os, numpy as np, re
# Does the ion correction that MAIN is applying right now exceed snaphu's
# recoverable gradient? Any |grad| above pi rad/pixel cannot be unwrapped.
for d, tag in ((\
    "/mnt/scratch/t/tommy/full/dates_ion_ml15", "ml15 (15x84, RUNNING NOW)"),
    ("/mnt/scratch/t/tommy/full/dates_ion", "ml2 (5x28, previous run)")):
    fs = sorted(glob.glob(d + "/filt_ion_*.ion"))
    if not fs:
        print(tag, "-> no grids"); continue
    # dims from the .vrt
    v = fs[0] + ".vrt"
    txt = open(v).read() if os.path.exists(v) else ""
    m = re.search(r'rasterXSize="(\d+)"\s+rasterYSize="(\d+)"', txt)
    if not m:
        m = re.search(r'rasterYSize="(\d+)"\s+rasterXSize="(\d+)"', txt)
        W, H = int(m.group(2)), int(m.group(1))
    else:
        W, H = int(m.group(1)), int(m.group(2))
    means, gr, ga, ptps = [], [], [], []
    for f in fs:
        a = np.fromfile(f, dtype=np.float32)
        if a.size != W * H: continue
        a = a.reshape(H, W)
        means.append(float(np.nanmean(a)))
        ptps.append(float(np.nanmax(a) - np.nanmin(a)))
        gr.append(float(np.nanmedian(np.abs(np.diff(a, axis=1)))))
        ga.append(float(np.nanmedian(np.abs(np.diff(a, axis=0)))))
    means, gr, ga, ptps = map(np.array, (means, gr, ga, ptps))
    print("\n===", tag, " grid %dx%d, %d dates" % (W, H, len(means)))
    print("  per-date mean:   median %12.1f  min %12.1f  max %12.1f rad" % (
        np.median(means), means.min(), means.max()))
    print("  per-date ptp:    median %12.1f  max %12.1f rad" % (np.median(ptps), ptps.max()))
    print("  |grad| range:    median %8.3f  max %8.3f rad/px   (pi = 3.142)" % (
        np.median(gr), gr.max()))
    print("  |grad| azimuth:  median %8.3f  max %8.3f rad/px" % (np.median(ga), ga.max()))
    over = int(((gr > np.pi) | (ga > np.pi)).sum())
    print("  dates with median gradient over pi/px: %d of %d" % (over, len(means)))
    # physical scale: L-band ~13.4 rad per TECU
    print("  implied TEC from per-date mean: median %.0f TECU (real ionosphere peaks ~200)" % (
        abs(np.median(means)) / 13.4))
