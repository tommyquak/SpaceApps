"""Round-11 self-test for crop_pack.py on realistic frame-3700-like numbers:
swath origins not on the 5-quantum, microsecond-quantised sensingStart at a
non-round PRF, 2-band BIL crop, pack identity checks, dangling-link refusal.
Round 12 added the overflow side file and the heavy-tailed pair; round 13
split the side-file refusals into missing, wrong-sized and foreign-but-right-
sized, so the int16-limit ownership check is exercised. Fits in the
login node's 1 GB virtual-memory limit (single-threaded BLAS, fixtures freed
between phases); set TMPDIR to a directory with quota for the ~300 MB of
fixtures."""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import datetime, json, sys, subprocess, tempfile, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import crop_pack as cp
# fixtures (~300 MB) go to a scratch directory, never next to the sources
T = tempfile.mkdtemp(prefix="crop_pack_selftest_")
os.chdir(T)
print("fixtures in", T)

# --- geometry fixture (resampled grid = swath 2 grid: 16.104 MHz, PRF 3314.512)
rps = 299792458.0 / (2 * 16.104e6)          # 9.3079 m
ali = 1.0 / 3314.512
t_min = datetime.datetime(2018, 11, 23, 3, 0, 0)
r_min = 815141.0
# synthetic swath origins (samples from the frame origin) and widths: NOT the
# real frame-3700 geometry (s3 at 4100, s4 at 9412; tests/golden)
off0 = {"2": 0, "3": 7623, "4": 15206}
N = {"2": 8016, "3": 7952, "4": 7900}
k_s = {"2": 42, "3": 14, "4": 0}
NL = 1400                                     # lines per swath (test size)
sw = {}
for sn in ("2", "3", "4"):
    ts = t_min + datetime.timedelta(seconds=k_s[sn] * ali)   # upstream quantisation (1 us)
    sw[sn] = {"orig": {"startingRange": r_min + off0[sn] * rps + 0.37 * rps,
                       "sensingStart": (ts - datetime.timedelta(microseconds=173)).strftime(cp.TFMT),
                       "rangePixelSize": rps * 1.0003, "azimuthLineInterval": ali * 0.9997,
                       "numberOfSamples": N[sn], "numberOfLines": NL,
                       "rangeSamplingRate": 16.104e6, "rangeBandwidth": 14e6},
              "res": {"startingRange": r_min + off0[sn] * rps, "sensingStart": ts.strftime(cp.TFMT),
                      "rangePixelSize": rps, "azimuthLineInterval": ali,
                      "numberOfSamples": N[sn], "numberOfLines": NL,
                      "rangeSamplingRate": 16.104e6, "rangeBandwidth": 14e6},
              "samples_from_origin": off0[sn], "lines_from_origin": k_s[sn]}
frame_w = max(off0[s] + N[s] for s in sw)                     # 23106
t_end = max(cp._t(sw[s]["res"]["sensingStart"]) + datetime.timedelta(seconds=NL * ali) for s in sw)
frame_l = int(round((t_end - t_min).total_seconds() / ali) / 14)
geom = {"ref_date": "181123", "frame": "3700", "swaths": sw,
        "frame_origin": {"startingRange": r_min, "sensingStart": t_min.strftime(cp.TFMT),
                         "rangePixelSize": rps, "azimuthLineInterval": ali},
        "frame_ml1": {"width": frame_w, "length": frame_l}, "subband": {}}
json.dump(geom, open("geom.json", "w"), indent=1)

# --- k_s recovery on the microsecond-quantised timestamps (old tol 1e-3 failed at k=42)
for sn in sw:
    v = (cp._t(sw[sn]["res"]["sensingStart"]) - t_min).total_seconds() / ali
    k = cp._lines_from_dt((cp._t(sw[sn]["res"]["sensingStart"]) - t_min).total_seconds(), ali)
    assert k == k_s[sn], (sn, k, v)
    print(f"s{sn}: k={k} residual {abs(v - k):.2e} lines")
try:
    cp._lines_from_dt(43.4 * ali, ali); raise AssertionError("non-integer accepted")
except SystemExit as e:
    print("non-integer line offset rejected:", str(e)[:60])

# --- lat/lon rasters: AOI in the middle, ml1 grid
W, L = frame_w, frame_l
lat = np.memmap("lat", dtype="float64", mode="w+", shape=(L, W))
lon = np.memmap("lon", dtype="float64", mode="w+", shape=(L, W))
rows = np.arange(L)[:, None]; cols = np.arange(W)[None, :]
lat[:] = -2.0 - 2.4 * rows / L; lon[:] = 113.0 + 3.0 * cols / W    # AOI lon 113.8-115.2 -> cols ~6160..16945
lat.flush(); lon.flush(); del lat, lon
win = cp.derive_window("lat", "lon", W, L, geom_path="geom.json", out="win.json")
lo_f = min(w["frame_sampOff"] for w in win["swaths"].values())
hi_f = max(w["frame_sampOff"] + w["width"] for w in win["swaths"].values())
assert (lo_f, hi_f) == (win["sampOff_frame"], win["sampOff_frame"] + win["cropWidth"])
for sn, w in win["swaths"].items():
    assert w["lineOff"] == win["lineOff"] - k_s[sn] and w["lineOff"] % 14 == 0
    assert w["width"] > 0 and w["frame_sampOff"] == off0[sn] + w["sampOff"]
print("per-swath widths mod 5:", {sn: w["width"] % 5 for sn, w in win["swaths"].items()},
      "sampOff mod 5:", {sn: w["sampOff"] % 5 for sn, w in win["swaths"].items()})
# mosaic_parameter re-run equivalent: frame width from per-swath edges
ends = [(off0[s] + w["sampOff"]) for s, w in win["swaths"].items()]
ends2 = [(off0[s] + w["sampOff"] + w["width"]) for s, w in win["swaths"].items()]
assert int(round((max(ends2) - min(ends)) * rps / rps)) == win["cropWidth"]

# --- pack with path identity checks: s3 of secondary 181012
w3, d3 = win["swaths"]["3"], sw["3"]
os.makedirs("dates_resampled/181012/f1_3700/s3", exist_ok=True)
os.makedirs("dates/181012/f1_3700/s3", exist_ok=True)
os.makedirs("dates/181012/f1_3700/s2", exist_ok=True)
rng = np.random.default_rng(1)
slc = (rng.standard_normal((NL, N["3"])) + 1j * rng.standard_normal((NL, N["3"]))).astype(np.complex64)
slc.tofile("dates_resampled/181012/f1_3700/s3/181012.slc")
# synthetic cull.off (writeOffset layout) on the original grid, total offset ~ smooth
pts = []
for y in range(20, NL - 20, 100):
    for x in range(50, N["3"] - 50, 400):
        dx = 3.2 + 1e-4 * x - 2e-4 * y; dy = -1.1 + 5e-5 * x + 1e-4 * y
        pts.append(f"{x:8d} {dx:10.4f} {y:8d} {dy:10.4f} {25.0:10.5f} {0.001:10.6f} {0.002:10.6f} {0.0:10.6f}")
open("dates/181012/f1_3700/s3/cull.off", "w").write("\n".join(pts) + "\n")
open("dates/181012/f1_3700/s2/cull.off", "w").write("\n".join(pts) + "\n")
json.dump(geom, open("geom.json", "w"), indent=1)
cp.HOME_GEOM, cp.HOME_WINDOW = os.path.abspath("geom.json"), os.path.abspath("win.json")
base = ["crop_pack.py", "pack", "dates_resampled/181012/f1_3700/s3/181012.slc", "s3.pack", "s3.json", "3"]
def run(args):
    try:
        cp.main(args); return "ok"
    except SystemExit as e:
        return f"exit: {str(e)[:90]}"
print("wrong-swath cull  ->", run(base + ["--cull", "dates/181012/f1_3700/s2/cull.off"]))
print("wrong-date cull   ->", run(base + ["--cull", "dates/181012/f1_3700/s3/cull.off", "--date", "181123"]))
print("secondary --ref   ->", run(base + ["--reference"]))
os.makedirs("dates/181123/f1_3700/s3", exist_ok=True)
open("dates/181123/f1_3700/s3/cull.off", "w").write("\n".join(pts) + "\n")
os.makedirs("work/181012", exist_ok=True)
os.chdir("work/181012")                                     # per-date job cwd
up = ["crop_pack.py", "pack", "../../dates_resampled/181012/f1_3700/s3/181012.slc", "s3.pack", "s3.json", "3"]
r = run(up + ["--cull", "../../dates/181123/f1_3700/s3/cull.off"])
print("cwd-date cull     ->", r); assert r.startswith("exit")
os.chdir(T)
os.makedirs("work/181012/dates/181123/f1_3700/s3", exist_ok=True)   # <date> lane above dates/
open("work/181012/dates/181123/f1_3700/s3/cull.off", "w").write("\n".join(pts) + "\n")
r = run(base + ["--cull", "work/181012/dates/181123/f1_3700/s3/cull.off"])
print("lane-date cull    ->", r); assert r.startswith("exit")
r = run(base + ["--cull", "dates/181012/f1_3700/s3/offset.off"])
print("not cull.off      ->", r); assert r.startswith("exit")
print("right cull        ->", run(base + ["--cull", "dates/181012/f1_3700/s3/cull.off"]))
sc = json.load(open("s3.json"))
assert sc["date"] == "181012" and sc["swath"] == 3 and sc["offset_polys"]["fit"] in ("isce", "numpy")
assert os.path.getsize("s3.pack") == sc["width"] * sc["length"] * 4
cp.unpack_slc("s3.pack", "s3.json", "s3.unp")
unp = np.fromfile("s3.unp", dtype=np.complex64).reshape(sc["length"], sc["width"])
ref = slc[w3["lineOff"]:w3["lineOff"] + w3["length"], w3["sampOff"]:w3["sampOff"] + w3["width"]]
dphi = np.angle(unp * np.conj(ref)); print("pack/unpack phase rms %.2e rad" % np.sqrt(np.mean(dphi ** 2)))
over = (np.abs(np.round(ref.real / sc["scale"])) > 32767) | (np.abs(np.round(ref.imag / sc["scale"])) > 32767)
n_over = int(over.sum())
assert sc["overflow_count"] == n_over, (sc["overflow_count"], n_over)
assert os.path.getsize("s3.pack.ovf") == 16 * n_over
print("gaussian crop: %d overflow pixels restored exactly" % n_over)
assert np.array_equal(unp[over], ref[over])
bright = np.abs(ref) > 0.01 * sc["peak"]                     # quantisation floor: |dphi| ~ scale/(sqrt(2)|z|)
print("phase err above 1 percent of peak: max %.1e rad" % np.abs(dphi[bright]).max())
assert np.abs(dphi[bright]).max() < 5e-3 and np.sqrt(np.mean(dphi ** 2)) < 1e-3

# --- heavy-tailed pair: speckle plus point targets 5-50x the int16 limit with a
# different phase per date; the packed pair must reproduce the float pair per
# pixel on the targets AND in every 5x28 multilooked interferogram cell
del unp, ref, dphi, over, bright
Lc, Wc = 28 * 20, 5 * 1400                                  # 3.9 M px
n_t = 60                                                     # 1.5e-5 << 1 - pctl, peak stays in the speckle tail
ty = rng.integers(0, Lc, n_t); tx = rng.integers(0, Wc, n_t)
gain = rng.uniform(5.0, 50.0, n_t) * 32767 * sc["scale"]
pair = {}
n_ovf = {}
for name, seed in (("ref", 11), ("sec", 12)):
    g = np.random.default_rng(seed)
    z = np.empty((Lc, Wc), dtype=np.complex64)
    z.real = g.standard_normal((Lc, Wc), dtype=np.float32)
    z.imag = g.standard_normal((Lc, Wc), dtype=np.float32)
    z[ty, tx] = (gain * np.exp(1j * g.uniform(-np.pi, np.pi, n_t))).astype(np.complex64)
    z[:3, :] = 0                                             # zero-padded lines stay zero
    z.tofile(f"{name}.slc")
    s = cp.pack_slc(f"{name}.slc", f"{name}.pack", f"{name}.json", 0, 0, Lc, Wc, Wc, Lc, sc["offset_polys"], swath=3,
                    date="181012" if name == "sec" else "181123")
    scn = cp.unpack_slc(f"{name}.pack", f"{name}.json", f"{name}.unp")
    u = np.fromfile(f"{name}.unp", dtype=np.complex64).reshape(Lc, Wc)
    assert scn["overflow_count"] >= n_t and scn["max_amp"] >= 5 * 32767 * scn["scale"]
    assert np.array_equal(u[ty, tx], z[ty, tx]) and np.array_equal(u[:3], z[:3])
    e = np.abs(np.angle(u * np.conj(z)))
    print(f"{name}: {scn['overflow_count']} overflow px, target phase err {e[ty, tx].max():.1e} rad, "
          f"amp err {np.abs(np.abs(u[ty, tx]) - np.abs(z[ty, tx])).max():.1e}")
    assert e[ty, tx].max() < 1e-2
    del e
    pair[name] = (z, u)
    n_ovf[name] = int(scn["overflow_count"])
def look(a, la=28, lr=5):
    a = a[:a.shape[0] // la * la, :a.shape[1] // lr * lr]
    return a.reshape(a.shape[0] // la, la, a.shape[1] // lr, lr).sum(axis=(1, 3))
ig_f = look(pair["ref"][0] * np.conj(pair["sec"][0]))
ig_p = look(pair["ref"][1] * np.conj(pair["sec"][1]))
cell_err = np.abs(np.angle(ig_p * np.conj(ig_f)))
hot = look(np.zeros((Lc, Wc)) + np.isin(np.arange(Lc * Wc).reshape(Lc, Wc), ty * Wc + tx)) > 0
print("5x28 igram cell phase err: max %.1e rad (target cells %.1e), amp ratio range %.6f..%.6f" %
      (cell_err.max(), cell_err[hot].max(), (np.abs(ig_p) / np.abs(ig_f)).min(), (np.abs(ig_p) / np.abs(ig_f)).max()))
assert cell_err.max() < 1e-2 and cell_err[hot].max() < 1e-3
del pair, ig_f, ig_p, cell_err, hot
# unpack refuses a missing, wrong-sized or foreign side file, and each case must
# be caught by its own check: the size test cannot stand in for the ownership test
os.replace("sec.pack.ovf", "sec.pack.ovf.bak")
try:
    cp.unpack_slc("sec.pack", "sec.json", "x"); raise AssertionError("missing ovf accepted")
except SystemExit as e:
    assert "missing or not" in str(e), e
    print("missing ovf refused:", e)
import shutil; shutil.copy("ref.pack.ovf", "sec.pack.ovf")
assert n_ovf["ref"] != n_ovf["sec"], n_ovf
try:
    cp.unpack_slc("sec.pack", "sec.json", "x"); raise AssertionError("wrong-sized ovf accepted")
except SystemExit as e:
    assert "missing or not" in str(e), e
    print("wrong-sized ovf refused:", e)
# right byte count, indices sorted and in range, but pointing at pixels that are
# not at the int16 limit: only the ownership check stands between it and a wrong
# unpack, so build it from pixels known to be below the limit
q = np.fromfile("sec.pack", dtype=np.int16).reshape(-1, 2)
free = np.flatnonzero(np.abs(q).max(axis=1) != 32767)[:n_ovf["sec"]].astype(np.int64)
assert free.size == n_ovf["sec"] and free.size > 0
with open("sec.pack.ovf", "wb") as f:
    free.tofile(f); np.zeros(free.size, dtype=np.complex64).tofile(f)
del q, free
assert os.path.getsize("sec.pack.ovf") == 16 * n_ovf["sec"]
try:
    cp.unpack_slc("sec.pack", "sec.json", "x"); raise AssertionError("foreign ovf accepted")
except SystemExit as e:
    assert "does not belong to" in str(e) and "int16 limit" in str(e), e
    print("foreign ovf refused:", e)
os.replace("sec.pack.ovf.bak", "sec.pack.ovf")
# polynomial recovery at the sidecar's evaluation coordinates: TOTAL offset on the resampled grid
p = sc["offset_polys"]
i, k = np.meshgrid(np.arange(0, w3["length"], 97), np.arange(0, w3["width"], 211), indexing="ij")
az = w3["lineOff"] + i + 1; rg = w3["sampOff"] + k + 1
o, r = d3["orig"], d3["res"]
# the applied polynomial is a function of the resampled coordinate x1 (= rg, as
# resamp_slc evaluates it): invert upstream's transform to the cull-grid (x, y),
# apply the synthetic offset field, and take the TOTAL offset (x + dx) - x1
x_c = (rg - (o["startingRange"] - r["startingRange"]) / r["rangePixelSize"]) * r["rangePixelSize"] / o["rangePixelSize"]
y_c = (az - (cp._t(o["sensingStart"]) - cp._t(r["sensingStart"])).total_seconds() / r["azimuthLineInterval"]) \
      * r["azimuthLineInterval"] / o["azimuthLineInterval"]
dx = 3.2 + 1e-4 * x_c - 2e-4 * y_c
dy = -1.1 + 5e-5 * x_c + 1e-4 * y_c
true_rg = (x_c + dx) - rg
true_az = (y_c + dy) - az
err_az = np.abs(cp.eval_poly(p["az"], az, rg) - true_az)
print("az poly max |err| %.2e px" % err_az.max())
assert err_az.max() < 1e-6
err = np.abs(cp.eval_poly(p["rg"], az, rg) - true_rg)
print("rg poly max |err| %.2e px" % err.max())
assert err.max() < 1e-6

# --- 2-band BIL crop (los), double and byte, and a dangling link
los = rng.standard_normal((L, 2, W)).astype(np.float32); los.tofile("los.rdr")
cp.crop_bil_raster("los.rdr", "los_crop.rdr", win["ml1"]["rowOff"], win["ml1"]["rows"],
                   win["ml1"]["colOff"], win["ml1"]["cols"], W, L, bands=2, dtype="float32")
got = np.fromfile("los_crop.rdr", dtype=np.float32).reshape(win["ml1"]["rows"], 2, win["ml1"]["cols"])
exp = los[win["ml1"]["rowOff"]:win["ml1"]["rowOff"] + win["ml1"]["rows"], :,
          win["ml1"]["colOff"]:win["ml1"]["colOff"] + win["ml1"]["cols"]]
assert np.array_equal(got, exp); print("2-band BIL crop exact")
hgt = rng.standard_normal((L, W)); hgt.tofile("hgt.rdr")
cp.crop_bil_raster("hgt.rdr", "hgt_crop.rdr", win["ml1"]["rowOff"], win["ml1"]["rows"],
                   win["ml1"]["colOff"], win["ml1"]["cols"], W, L, bands=1, dtype="double")
assert np.array_equal(np.fromfile("hgt_crop.rdr", dtype="float64").reshape(win["ml1"]["rows"], -1),
                      hgt[win["ml1"]["rowOff"]:win["ml1"]["rowOff"] + win["ml1"]["rows"],
                          win["ml1"]["colOff"]:win["ml1"]["colOff"] + win["ml1"]["cols"]]); print("double crop exact")
wbd = rng.integers(-1, 1, (L, W)).astype(np.int8); wbd.tofile("wbd.rdr")
cp.crop_bil_raster("wbd.rdr", "wbd_crop.rdr", 0, 14, 5, 25, W, L, bands=1, dtype="byte")
assert np.array_equal(np.fromfile("wbd_crop.rdr", dtype=np.int8).reshape(14, 25), wbd[0:14, 5:30]); print("byte crop exact")
try:
    cp.crop_bil_raster("hgt.rdr", "x", 0, 14, 5, 25, W, L, bands=1, dtype="float32"); raise AssertionError
except SystemExit as e:
    print("dtype/size mismatch refused:", str(e)[:50])
if os.path.lexists("dangling.off"): os.remove("dangling.off")
os.symlink("nonexistent_rg.off", "dangling.off")
try:
    cp.crop_bil_raster("dangling.off", "x", 0, 14, 5, 25, W, L); raise AssertionError
except SystemExit as e:
    print("dangling link refused:", str(e)[:70])
# --- render: cropoff names map to ISCE types, bad names and sizes refused before any xml
for bad in (("los_crop.rdr", win["ml1"]["cols"], win["ml1"]["rows"], "int16", 2, "BIL"),
            ("los_crop.rdr", win["ml1"]["cols"], win["ml1"]["rows"], "float32", 2, "BIQ"),
            ("los_crop.rdr", win["ml1"]["cols"], win["ml1"]["rows"], "float32", 1, "BIL")):
    try:
        cp.render_isce_xml(*bad); raise AssertionError(bad)
    except SystemExit as e:
        print("render refused:", str(e)[:60])
assert not os.path.exists("los_crop.rdr.xml")
try:
    import isce  # noqa: F401
    have_isce = True
except ImportError:
    have_isce = False
if have_isce:
    for (fn, dt, bands, off) in (("los_crop.rdr", "float32", 2, 4), ("hgt_crop.rdr", "double", 1, 8),
                                 ("s3.unp", "CFLOAT", 1, 8)):
        w_, l_ = (win["ml1"]["cols"], win["ml1"]["rows"]) if fn != "s3.unp" else (sc["width"], sc["length"])
        got = cp.render_isce_xml(fn, w_, l_, dt, bands, "BIL")
        vrt = open(fn + ".vrt").read(); xml = open(fn + ".xml").read()
        assert f"<PixelOffset>{off}</PixelOffset>" in vrt and f"<LineOffset>{w_ * off * bands}</LineOffset>" in vrt, vrt
        assert f"<value>{got}</value>" in xml and vrt.count("<VRTRasterBand") == bands
        print(f"render {fn} {dt} -> {got} x{bands}: xml+vrt ok")
    # rendered from elsewhere with a directory-qualified path: the xml must
    # still hold the bare basename (a stored path binds the data to the tree it
    # was rendered in; look_geom's Looks opens whatever the xml names) and the
    # xml/vrt must land beside the binary
    os.makedirs("insar", exist_ok=True)
    shutil.copy("los_crop.rdr", os.path.join("insar", "d_1rlks_14alks.los"))
    rel = os.path.join("insar", "d_1rlks_14alks.los")
    cp.render_isce_xml(rel, win["ml1"]["cols"], win["ml1"]["rows"], "float32", 2, "BIL")
    assert os.getcwd() == T, os.getcwd()
    xml = open(rel + ".xml").read()
    assert "<value>d_1rlks_14alks.los</value>" in xml, xml
    assert os.path.isfile(rel + ".vrt") and "<SourceFilename relativeToVRT=\"1\">d_1rlks_14alks.los" in open(rel + ".vrt").read()
    print("render from another directory: xml holds the basename, xml+vrt beside the binary")
else:
    print("render isce test skipped (no isce)")
print("SELFTEST_OK")
