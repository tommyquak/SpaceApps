"""The validator has to catch a bad run, so check the checker: build a small
but structurally real persisted date, confirm check_run passes it, then break
it one way at a time and confirm each break is caught. Runs anywhere (numpy
only); the ISCE xml-load step inside check_run skips itself without ISCE."""
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import crop_pack as cp

D, REF, FRAME = "181012", "181123", "3700"
T = tempfile.mkdtemp(prefix="check_run_test_")
P = os.path.join(T, "stack")
LINE_OFF, CROP_LEN, SAMP_OFF, CROP_W = 56, 280, 0, 300
SW = {"2": {"from": 0, "lineOff": 56, "res_t0": "2018-11-23T03:00:00.000000"},
      "3": {"from": 100, "lineOff": 42, "res_t0": "2018-11-23T03:00:00.004224"},
      "4": {"from": 200, "lineOff": 28, "res_t0": "2018-11-23T03:00:00.008448"}}
IN_W, IN_L = 200, 400
ALI, RPS, R0 = 0.0003017035388618294, 8.582533435738496, 734356.0
polys = {"rg": cp._poly_dict([[0.4]]), "az": cp._poly_dict([[-0.2]]),
         "fit": "numpy", "n": 12, "source": "synthetic", "eval": cp.POLY_EVAL}


def build():
    shutil.rmtree(P, ignore_errors=True)
    os.makedirs(os.path.join(P, "dates", REF, "f1_" + FRAME, "mosaic"))
    open(os.path.join(P, "GATE_REF_PUBLISHED"), "w").close()
    for f in (os.path.join(P, "dates", REF, REF + ".track.xml"),
              os.path.join(P, "dates", REF, "f1_" + FRAME, REF + ".frame.xml"),
              os.path.join(P, "dates", REF, "f1_" + FRAME, "mosaic", "swath_offset.txt")):
        open(f, "w").write("<placeholder/>\n")
    os.makedirs(os.path.join(P, "dates", D, "f1_" + FRAME))
    open(os.path.join(P, "dates", D, D + ".track.xml"), "w").write("<placeholder/>\n")
    open(os.path.join(P, "dates", D, "f1_" + FRAME, D + ".frame.xml"), "w").write("<placeholder/>\n")
    geom = {"ref_date": REF, "frame": FRAME, "frame_ml1": {"width": CROP_W, "length": CROP_LEN // 14},
            "frame_origin": {"startingRange": R0, "sensingStart": "2018-11-23T03:00:00.000000",
                             "rangePixelSize": RPS, "azimuthLineInterval": ALI},
            "swaths": {}}
    win = {"lineOff": LINE_OFF, "cropLen": CROP_LEN, "sampOff_frame": SAMP_OFF, "cropWidth": CROP_W,
           "ml1": {"rowOff": LINE_OFF // 14, "rows": CROP_LEN // 14, "colOff": SAMP_OFF, "cols": CROP_W},
           "swaths": {}}
    for sn, d in SW.items():
        geom["swaths"][sn] = {
            "orig": {"startingRange": R0 + d["from"] * RPS, "sensingStart": d["res_t0"],
                     "rangePixelSize": RPS, "azimuthLineInterval": ALI,
                     "numberOfSamples": IN_W, "numberOfLines": IN_L},
            "res": {"startingRange": R0 + d["from"] * RPS, "sensingStart": d["res_t0"],
                    "rangePixelSize": RPS, "azimuthLineInterval": ALI,
                    "numberOfSamples": IN_W, "numberOfLines": IN_L},
            "samples_from_origin": d["from"], "lines_from_origin": LINE_OFF - d["lineOff"]}
        win["swaths"][sn] = {"lineOff": d["lineOff"], "length": CROP_LEN, "sampOff": 0,
                             "width": 100, "frame_sampOff": d["from"],
                             "lines_from_origin": LINE_OFF - d["lineOff"]}
    json.dump(geom, open(os.path.join(P, "ref_geometry.json"), "w"), indent=1)
    json.dump(win, open(os.path.join(P, "crop_window.json"), "w"), indent=1)

    rng = np.random.default_rng(7)
    od = os.path.join(P, "dates_resampled", D, "f1_" + FRAME)
    for sn, d in SW.items():
        sd = os.path.join(od, "s" + sn)
        os.makedirs(sd)
        src = os.path.join(T, f"{sn}.slc")
        z = np.empty((IN_L, IN_W), dtype=np.complex64)
        z.real = rng.standard_normal((IN_L, IN_W), dtype=np.float32)
        z.imag = rng.standard_normal((IN_L, IN_W), dtype=np.float32)
        # a realistic heavy tail: a sparse set of scatterers 5-15x the speckle,
        # so P99.995 (the packing scale) still sits in the speckle and only the
        # brightest targets overflow the int16 range
        ty = rng.integers(0, IN_L, 40); tx = rng.integers(0, IN_W, 40)
        z[ty, tx] *= rng.uniform(5.0, 15.0, 40)
        z.tofile(src)
        cp.pack_slc(src, os.path.join(sd, D + ".slc"), os.path.join(sd, D + ".slc.json"),
                    d["lineOff"], 0, CROP_LEN, 100, IN_W, IN_L, polys, swath=int(sn), date=D,
                    frame_lineOff=LINE_OFF, frame_sampOff=d["from"])
    ins = os.path.join(P, "dates_resampled", D, "insar")
    os.makedirs(ins)
    off = os.path.join(ins, f"{D}_1rlks_14alks_rg_rect.off")
    (rng.standard_normal((CROP_LEN // 14, CROP_W)).astype(np.float32) * 3).tofile(off)
    open(off + ".xml", "w").write("<placeholder/>\n")


def run():
    env = dict(os.environ, PERSIST=P, REF_DATE=REF, FRAME=FRAME)
    r = subprocess.run([sys.executable, os.path.join(ROOT, "check_run.py"), "date", D],
                       capture_output=True, text=True, env=env)
    return r.returncode, r.stdout + r.stderr


build()
rc, out = run()
assert rc == 0, "checker rejected a good date:\n" + out
print(out.strip().splitlines()[-1])

CASES = []


def case(name, fn, expect):
    build()
    fn()
    rc, out = run()
    assert rc == 1, f"{name}: checker PASSED a broken date:\n{out}"
    assert expect in out, f"{name}: expected {expect!r}, got:\n{out}"
    CASES.append(name)
    print(f"  caught {name}: {[l.strip() for l in out.splitlines() if l.strip().startswith('-')][0][:96]}")


def _rewrite_window():
    w = json.load(open(os.path.join(P, "crop_window.json")))
    w["swaths"]["2"]["width"] = 95          # the window moved after this crop was packed
    w["swaths"]["3"]["frame_sampOff"] = 95
    w["swaths"]["3"]["width"] = 105
    json.dump(w, open(os.path.join(P, "crop_window.json"), "w"), indent=1)


sd2 = os.path.join(P, "dates_resampled", D, "f1_" + FRAME, "s2")
case("truncated sidecar",
     lambda: open(os.path.join(sd2, D + ".slc.json"), "r+").truncate(120),
     "does not parse")
case("zeroed pack",
     lambda: np.zeros(os.path.getsize(os.path.join(sd2, D + ".slc")) // 2,
                      dtype=np.int16).tofile(os.path.join(sd2, D + ".slc")),
     "zero")
case("truncated pack",
     lambda: open(os.path.join(sd2, D + ".slc"), "r+").truncate(1024),
     "sidecar says")
case("foreign overflow file",
     lambda: shutil.copy(os.path.join(P, "dates_resampled", D, "f1_" + FRAME, "s3", D + ".slc.ovf"),
                         os.path.join(sd2, D + ".slc.ovf")),
     "overflow")
case("crop taken with a stale window",
     lambda: _rewrite_window(),
     "stale crop?")
case("missing overflow file",
     lambda: os.remove(os.path.join(sd2, D + ".slc.ovf")),
     "overflow pixels")
case("wrong-sized rg_rect.off",
     lambda: open(os.path.join(P, "dates_resampled", D, "insar",
                               f"{D}_1rlks_14alks_rg_rect.off"), "r+").truncate(64),
     "rg_rect.off is")
case("rg_rect.off without its xml",
     lambda: os.remove(os.path.join(P, "dates_resampled", D, "insar",
                                    f"{D}_1rlks_14alks_rg_rect.off.xml")),
     "no xml")


shutil.rmtree(T, ignore_errors=True)
print(f"CHECK_RUN_TEST_OK ({len(CASES)} breakages caught)")
