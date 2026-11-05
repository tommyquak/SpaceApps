"""Cluster-only end-to-end check of the SEED metadata path on real ISCE products:
crop_pack geometry -> window -> rewrite -> mosaic_parameter re-run -> verify.

The smoke run persisted only the RAW 181012 tree, so the pre-crop resampled
reference tree is rebuilt here from it with the frame-3700 numbers the smoke
run printed (origins 0/4100/9412 samples, swath starts 42/14/0 lines, dims
172356x4440, 172355x5670, 172355x5794 on the s2 grid) and saved with the
authentic saveTrack, i.e. genuine ISCE xmls. Needs the smoke gate tree
(~/full/stack_smoke_artifacts/gate, override with SMOKE_GATE) and ISCE_STACK;
skips otherwise. Metadata only: no SLC binaries are touched.
"""
import datetime, json, os, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import crop_pack as cp

def _skip(msg):
    """Skip cleanly under BOTH pytest and direct execution.

    A module-level `raise SystemExit` aborts pytest's COLLECTION of the entire
    suite, not just this file, so the whole test run reported an error even
    though this cluster-only check was simply not applicable. pytest.skip with
    allow_module_level=True is the supported way to skip at import time; when
    pytest is absent (running this file as a plain script, which is how it is
    used on the cluster) fall back to the original SystemExit(0).
    """
    print(msg)
    # Detect whether we are UNDER pytest, not merely whether pytest is
    # installed. The previous version imported pytest and called pytest.skip
    # unconditionally, which raises Skipped outside a test runner and so broke
    # direct execution on any machine that happens to have pytest available.
    if "pytest" in sys.modules:
        import pytest
        pytest.skip(msg, allow_module_level=True)
    raise SystemExit(0)

G = os.environ.get("SMOKE_GATE", os.path.expanduser("~/full/stack_smoke_artifacts/gate"))
W = os.path.join(os.environ.get("TMPDIR", "/tmp"), "seed_metadata_check")
DATE = "181012"
if not os.path.isfile(os.path.join(G, "dates", DATE, DATE + ".track.xml")):
    _skip("seed metadata test skipped (no smoke gate tree at %s)" % G)
if "ISCE_STACK" not in os.environ:
    _skip("seed metadata test skipped (ISCE_STACK not set)")
loadTrack, saveTrack = cp._stackpulic("loadTrack", "saveTrack")

shutil.rmtree(W, ignore_errors=True)
os.makedirs(W)
shutil.copytree(os.path.join(G, "dates"), os.path.join(W, "dates"))

# swath_offset.txt in the upstream format (runSwathOffset.py:157-159); the
# matching column is what resample_common_grid and build_geometry read
mo = os.path.join(W, "dates", DATE, "f1_3700", "mosaic")
os.makedirs(mo, exist_ok=True)
RG_M = [0.0, -4100.224, -5312.159]
AZ_M = [0.0, 20.272, 17.607]
RG_G = [0.0, -4100.421, -5312.301]
AZ_G = [0.0, 19.466, 19.360]
with open(os.path.join(mo, "swath_offset.txt"), "w") as f:
    f.write("offset type       i     geometrical           match      difference\n")
    f.write("+" * 67 + "\n")
    for i, (a, b) in enumerate(zip(RG_G, RG_M)):
        f.write("range offset     {:2d}   {:13.3f}   {:13.3f}   {:13.3f}\n".format(i, a, b, a - b))
    for i, (a, b) in enumerate(zip(AZ_G, AZ_M)):
        f.write("azimuth offset   {:2d}   {:13.3f}   {:13.3f}   {:13.3f}\n".format(i, a, b, a - b))

# pre-crop resampled reference tree on the s2 grid
raw = loadTrack(os.path.join(W, "dates", DATE), DATE)
s2 = raw.frames[0].swaths[0]
RPS, ALI = float(s2.rangePixelSize), float(s2.azimuthLineInterval)
R0, T0 = float(s2.startingRange), s2.sensingStart
GRID = {"2": (0, 42, 4440, 172356), "3": (4100, 14, 5670, 172355), "4": (9412, 0, 5794, 172355)}
res = loadTrack(os.path.join(W, "dates", DATE), DATE)
for sw in res.frames[0].swaths:
    off, k, ns, nl = GRID[str(int(sw.swathNumber))]
    sw.startingRange = R0 + off * RPS
    sw.sensingStart = T0 + datetime.timedelta(seconds=k * ALI)
    sw.rangePixelSize, sw.azimuthLineInterval = RPS, ALI
    sw.numberOfSamples, sw.numberOfLines = ns, nl
rdir = os.path.join(W, "dates_resampled", DATE)
for sn in GRID:
    os.makedirs(os.path.join(rdir, "f1_3700", "s" + sn))   # stackDateStatistics globs f1_*/s*

cwd = os.getcwd()
os.chdir(rdir)
saveTrack(res, DATE)
os.chdir(cwd)

GEOM = os.path.join(W, "ref_geometry.json")
WIN = os.path.join(W, "crop_window.json")
g = cp.build_geometry(os.path.join(W, "dates"), os.path.join(W, "dates_resampled"), DATE, out=GEOM)
assert {sn: d["samples_from_origin"] for sn, d in g["swaths"].items()} == {"2": 0, "3": 4100, "4": 9412}
assert {sn: d["lines_from_origin"] for sn, d in g["swaths"].items()} == {"2": 42, "3": 14, "4": 0}
print("geometry ok:", {sn: (d["samples_from_origin"], d["lines_from_origin"]) for sn, d in g["swaths"].items()},
      "frame ml1", g["frame_ml1"])

# window: the real derive needs the ml1 lat/lon rasters (1.5 GB each, not
# persisted), so build the same structure with the same helpers
lineOff, cropLen = 56, 28 * 100
sampOff, cropWidth = 0, 15205
win = {"lineOff": lineOff, "cropLen": cropLen, "sampOff_frame": sampOff, "cropWidth": cropWidth,
       "swaths": {}}
for sn, d in g["swaths"].items():
    rw = cp.swath_window(win, d["res"]["startingRange"], g["frame_origin"]["startingRange"],
                         g["frame_origin"]["rangePixelSize"], d["res"]["numberOfSamples"], sn)
    lo = cp.swath_line_offset(win, d["res"]["sensingStart"], g["frame_origin"]["sensingStart"],
                              g["frame_origin"]["azimuthLineInterval"], d["res"]["numberOfLines"], sn)
    win["swaths"][sn] = {"lineOff": lo, "length": cropLen, "sampOff": rw["sampOff"],
                         "width": rw["width"], "frame_sampOff": d["samples_from_origin"] + rw["sampOff"],
                         "lines_from_origin": d["lines_from_origin"]}
lo_f = min(w["frame_sampOff"] for w in win["swaths"].values())
hi_f = max(w["frame_sampOff"] + w["width"] for w in win["swaths"].values())
assert (lo_f, hi_f) == (sampOff, sampOff + cropWidth), (lo_f, hi_f)
json.dump(win, open(WIN, "w"), indent=1)
print("window ok:", {sn: (w["lineOff"], w["sampOff"], w["width"]) for sn, w in win["swaths"].items()})

r = cp.rewrite_resampled(rdir, DATE, geom_path=GEOM, window_path=WIN)
# the rewrite must be refused a second time (the tree no longer matches the geometry)
try:
    cp.rewrite_resampled(rdir, DATE, geom_path=GEOM, window_path=WIN)
    raise AssertionError("second rewrite accepted")
except SystemExit as e:
    assert "already rewritten" in str(e), e
    print("second rewrite refused:", str(e)[:90])

mp = os.path.join(os.environ["ISCE_STACK"], "alosStack", "mosaic_parameter.py")
subprocess.run([sys.executable, mp, "-idir", "dates_resampled", "-ref_date", DATE,
                "-sec_date", DATE, "-nrlks1", "1", "-nalks1", "14"], cwd=W, check=True,
               stdout=subprocess.DEVNULL)
v = cp.verify_cropped(rdir, DATE, geom_path=GEOM, window_path=WIN)
assert v == {"width": cropWidth, "length": cropLen // 14}, v

# a window that disagrees with the tree must fail verify, not pass silently
bad = json.load(open(WIN)); bad["cropLen"] = cropLen + 28
json.dump(bad, open(WIN + ".bad", "w"))
try:
    cp.verify_cropped(rdir, DATE, geom_path=GEOM, window_path=WIN + ".bad")
    raise AssertionError("verify accepted a wrong window")
except SystemExit as e:
    print("verify catches a wrong window:", str(e)[:100])
shutil.rmtree(W, ignore_errors=True)
print("SEED_METADATA_OK")
