"""Wave-1 crop/pack tool for the alosStack pivot (spec: PIVOT_ALOSSTACK.md v5,
amended by audit round 10).

Implements, per the certified contract:
  sec 1: one absolute azimuth window for all swaths; per-swath range windows
         derived from absolute frame-level edges; divisibility lineOff%28,
         cropLen%28, sampOff%5, cropWidth%5 enforced at derivation time.
         Each resampled swath starts k_s*14 lines after the frame origin
         (resample_common_grid snapping), so the swath-local crop start is
         lineOff_s = lineOff - 14*k_s (multiple of 14, not necessarily 28).
         Range origins snap to nrlks1 = 1, so the swath-local sampOff_s and
         width_s are unquantised integers; only the frame-level sampOff and
         cropWidth carry the 5-sample quantum, and the union of the per-swath
         windows must equal the frame window (Amendment A, round 11).
         The window is clipped to the common azimuth extent of all swaths
         and widened in range so the s2/s3 and s3/s4 overlap strips stay
         whole; a swath falling entirely outside the window is an error.
  sec 3: each date's APPLIED range/azimuth offset polynomials, reproduced the
         way resample_common_grid.py builds them (secondary: cull.off mapped
         onto the resampled swath grid, TOTAL offset, OffsetField
         .getFitPolynomials; reference: the analytic grid-snap polynomial),
         serialised in the Poly2D convention so wave-3 evaluates them exactly
         as resamp_slc.f90 did.
  sec 4A: complex-int16 packing with a JSON sidecar as the single source of
         truth {swath, width, length, sampOff, lineOff, scale, offset_polys}.
  sec 4B/4F: ml1-window crops of BIL rasters (rg_rect.off float32; lat/lon/
         hgt double; wbd byte; los 2-band float32) + rewritten xml.
render_isce_xml renders .slc.xml/.vrt via isceobj; block jobs call it after
unpack, before pair_up (sec 5b).

SEED-time order (all BEFORE the sec-2 metadata rewrite, because the
resampled xmls must still describe the pre-crop grid):
  geometry -> derive -> pack --reference (per swath) -> sec-2 rewrite
Wave-1 date jobs: pack --cull <cull.off> per swath.
Dimensions always come from ref_geometry.json, never from file sizes; file
sizes are only ever checked against them.
"""
import datetime
import json
import os
import re
import sys

import numpy as np

FULL = os.path.expanduser("~/full")
HOME_WINDOW = os.path.join(FULL, "crop_window.json")
HOME_GEOM = os.path.join(FULL, "ref_geometry.json")
AOI = (113.8, -3.4, 115.2, -2.34)  # w s e n (campaign AOI, frame-clipped)
NALKS1, NRLKS1, NALKS2, NRLKS2 = 14, 1, 2, 5
LINE_Q = NALKS1 * NALKS2   # 28: frame-level lineOff and cropLen quantum
SAMP_Q = NRLKS1 * NRLKS2   # 5:  sampOff and cropWidth quantum
TFMT = "%Y-%m-%dT%H:%M:%S.%f"
POLY_EVAL = ("Poly2D: value = sum_ij coeffs[i][j] * ((az - meanAzimuth)/normAzimuth)**i"
             " * ((rg - meanRange)/normRange)**j. az, rg are the 1-based line and"
             " sample of the RESAMPLED SWATH grid, exactly as resamp_slc.f90 evaluates"
             " them (r_at=j, r_rt=i); for 0-based cropped (i, k) of this sidecar:"
             " az = lineOff + i + 1, rg = sampOff + k + 1. Ramp restoration for a"
             " subband regenerated from the CROPPED file (rg_filter demodulates by"
             " the cropped 0-based index k while the authentic subband was"
             " demodulated on the original grid at k + sampOff + off_rg): multiply"
             " by exp(-1j*2*pi*bc*(off_rg(az, rg) + sampOff)); the sampOff term is a"
             " per-swath constant that cancels in interferograms")


def _snap_down(v, q):
    return int(v) // q * q


def _snap_up(v, q):
    return -(-int(v) // q) * q


def _t(s):
    return datetime.datetime.strptime(s, TFMT)


def _near_int(v, tol=1e-3):
    """Integer sample offset from a float64 range difference (range sites only)."""
    r = int(round(v))
    if abs(v - r) > tol:
        raise SystemExit(f"expected an integer pixel offset, got {v}")
    return r


def _lines_from_dt(dt_seconds, ali, tol_seconds=1.5e-6):
    """Integer line count from a datetime difference. resample_common_grid
    stores sensingStart = t_min + datetime.timedelta(seconds=k*14*ali), which
    timedelta (and the ISCE xml %f timestamp) quantise to 1 us, so the
    residual must be tested in seconds: 0.5e-6*prf lines exceeds a 1e-3 line
    tolerance at ALOS-2 ScanSAR PRFs (2270-3314 Hz)."""
    k = int(round(dt_seconds / ali))
    if abs(dt_seconds - k * ali) > tol_seconds:
        raise SystemExit(f"expected an integer line offset, got {dt_seconds / ali} lines "
                         f"(residual {dt_seconds - k * ali:.3e} s > {tol_seconds:.1e} s)")
    return k


# ---------------------------------------------------------------- geometry

def _swath_dict(sw):
    return {"startingRange": float(sw.startingRange),
            "sensingStart": sw.sensingStart.strftime(TFMT),
            "rangePixelSize": float(sw.rangePixelSize),
            "azimuthLineInterval": float(sw.azimuthLineInterval),
            "numberOfSamples": int(sw.numberOfSamples),
            "numberOfLines": int(sw.numberOfLines),
            "rangeSamplingRate": float(sw.rangeSamplingRate),
            "rangeBandwidth": float(sw.rangeBandwidth)}


def _stackpulic(*names):
    """Import names from the alosStack StackPulic module. The scripts import it
    from their own directory, so the job environment only has
    $ISCE_STACK/alosStack on PATH (Amendment A.10)."""
    import isce  # noqa: F401
    try:
        mod = __import__("StackPulic", fromlist=list(names))
    except ImportError:
        import shutil
        cand = [os.path.join(os.environ.get("ISCE_STACK", ""), "alosStack")]
        w = shutil.which("resample_common_grid.py")
        if w:
            cand.append(os.path.dirname(os.path.realpath(w)))
        for c in cand:
            if os.path.isfile(os.path.join(c, "StackPulic.py")):
                sys.path.insert(0, c)
                break
        else:
            raise SystemExit("StackPulic.py not importable; set ISCE_STACK or "
                             "put $ISCE_STACK/alosStack on PYTHONPATH")
        mod = __import__("StackPulic", fromlist=list(names))
    return tuple(getattr(mod, n) for n in names)


def build_geometry(idir, odir, ref_date, out=None):
    """Per-swath reference geometry: the offset-ADJUSTED original swath
    parameters and the resampled (pre-crop) swath parameters, computed the way
    resample_common_grid.py computes them (single-frame stacks only), plus the
    frame origin as mosaic_parameter.py defines it and the subband parameters.
    Run at SEED time before the sec-2 rewrite. Requires the ISCE environment.
    """
    loadTrack, stackDateStatistics, subbandParameters = _stackpulic(
        "loadTrack", "stackDateStatistics", "subbandParameters")
    dateDirs, dates, frames, swaths, iref = stackDateStatistics(idir, ref_date)
    if len(frames) != 1:
        raise SystemExit("build_geometry: single-frame stacks only "
                         "(frame_offset.txt accumulation is not replicated)")
    frameDir = "f1_{}".format(frames[0])
    raw = loadTrack(dateDirs[iref], ref_date)
    res = loadTrack(os.path.join(os.path.abspath(odir), ref_date), ref_date)
    fr, fs = raw.frames[0], res.frames[0]
    nsw = len(fr.swaths)
    if nsw != len(swaths) or nsw != len(fs.swaths):
        raise SystemExit("build_geometry: swath count mismatch between raw, "
                         "resampled and directory listing")
    for j in range(nsw):
        if int(fr.swaths[j].swathNumber) != swaths[j] or int(fs.swaths[j].swathNumber) != swaths[j]:
            raise SystemExit("build_geometry: swath order mismatch")
    if nsw > 1:
        t0s = {x.sensingStart for x in fs.swaths}
        nls = {int(x.numberOfLines) for x in fs.swaths}
        if len(t0s) == 1 and len(nls) == 1:
            raise SystemExit("build_geometry: resampled swaths already share sensingStart "
                             "and numberOfLines; this looks post-rewrite. Geometry must be "
                             "built from the pre-crop resampled xmls")

    # matching swath offsets, read exactly as resample_common_grid.py reads them
    rg_m, az_m = [], []
    if nsw >= 2:
        with open(os.path.join(dateDirs[iref], frameDir, "mosaic", "swath_offset.txt")) as f:
            for ln in f:
                if "range offset" in ln:
                    rg_m.append(float(ln.split()[4]))
                if "azimuth offset" in ln:
                    az_m.append(float(ln.split()[4]))
    else:
        rg_m, az_m = [0.0], [0.0]
    if len(rg_m) != nsw or len(az_m) != nsw:
        raise SystemExit("build_geometry: swath_offset.txt row count != swath count")

    sw0 = fr.swaths[0]
    per = {}
    for j, sw in enumerate(fr.swaths):
        adj_r0 = sw0.startingRange - sum(rg_m[0:j + 1]) * sw0.rangePixelSize
        adj_t0 = sw0.sensingStart - datetime.timedelta(
            seconds=sum(az_m[0:j + 1]) * sw0.azimuthLineInterval)
        o = _swath_dict(sw)
        o["startingRange"] = float(adj_r0)
        o["sensingStart"] = adj_t0.strftime(TFMT)
        per[str(int(sw.swathNumber))] = {"orig": o, "res": _swath_dict(fs.swaths[j])}

    # frame origin and ml1 frame dims as mosaic_parameter.py defines them
    r0 = min(x.startingRange for x in fs.swaths)
    t0 = min(x.sensingStart for x in fs.swaths)
    r1 = max(x.startingRange + x.rangePixelSize * x.numberOfSamples for x in fs.swaths)
    t1 = max(x.sensingStart + datetime.timedelta(seconds=x.azimuthLineInterval * x.numberOfLines)
             for x in fs.swaths)
    rps = float(fs.swaths[0].rangePixelSize)
    ali = float(fs.swaths[0].azimuthLineInterval)
    for x in fs.swaths:
        if abs(x.rangePixelSize - rps) > 1e-9 or abs(x.azimuthLineInterval - ali) > 1e-12:
            raise SystemExit("build_geometry: resampled swaths do not share pixel sizes")
    frame_ml1 = {"width": int(round((r1 - r0) / rps) / NRLKS1),
                 "length": int(round((t1 - t0).total_seconds() / ali) / NALKS1)}
    for sn, d in per.items():
        d["samples_from_origin"] = _near_int((d["res"]["startingRange"] - r0) / rps)
        k = _lines_from_dt((_t(d["res"]["sensingStart"]) - t0).total_seconds(), ali)
        if k < 0 or k % NALKS1 != 0:
            raise SystemExit(f"build_geometry: swath {sn} starts {k} lines after the frame "
                             f"origin; expected a non-negative multiple of {NALKS1}")
        d["lines_from_origin"] = k

    wl, bw, fc, pre = subbandParameters(res)
    out = out or HOME_GEOM
    geom = {"ref_date": ref_date, "frame": frames[0], "swaths": per,
            "frame_origin": {"startingRange": float(r0), "sensingStart": t0.strftime(TFMT),
                             "rangePixelSize": rps, "azimuthLineInterval": ali},
            "frame_ml1": frame_ml1,
            "subband": {"radarWavelength": [float(v) for v in wl],
                        "bandWidth": [float(v) for v in bw],
                        "frequencyCenter": [float(v) for v in fc],
                        "prefix": list(pre)}}
    with open(out, "w") as f:
        json.dump(geom, f, indent=1)
    print("geometry: frame ml1 {}x{} (w x l); swath line offsets {}".format(
        frame_ml1["width"], frame_ml1["length"],
        {sn: d["lines_from_origin"] for sn, d in per.items()}))
    return geom


def load_geometry(path=None):
    return json.load(open(path or HOME_GEOM))


# ------------------------------------------------------------------ window

def derive_window(ref_ml1_lat, ref_ml1_lon, ml1_width, ml1_length,
                  geom_path=None, out=None):
    """Find the full-res frame-level crop window bounding the AOI, then
    validate and adjust it against every swath (spec sec 1):
      - clipped to the azimuth extent common to all resampled swaths;
      - widened so any s/s+1 overlap strip it touches is included whole;
      - per-swath windows derived from the absolute frame edges.
    ref_ml1_lat/lon: the reference 1rlks_14alks lat/lon rasters (float64),
    dims (ml1_length, ml1_width) which must equal the frame's ml1 dims.
    """
    g = load_geometry(geom_path)
    out = out or HOME_WINDOW
    if (ml1_width, ml1_length) != (g["frame_ml1"]["width"], g["frame_ml1"]["length"]):
        raise SystemExit("derive_window: ml1 dims {}x{} != frame ml1 dims {}x{} from geometry".format(
            ml1_width, ml1_length, g["frame_ml1"]["width"], g["frame_ml1"]["length"]))
    for p in (ref_ml1_lat, ref_ml1_lon):
        if os.path.getsize(p) != ml1_width * ml1_length * 8:
            raise SystemExit(f"derive_window: {p} size != {ml1_width}x{ml1_length} float64")
    lat = np.memmap(ref_ml1_lat, dtype="float64", mode="r", shape=(ml1_length, ml1_width))
    lon = np.memmap(ref_ml1_lon, dtype="float64", mode="r", shape=(ml1_length, ml1_width))
    inb = (lon >= AOI[0]) & (lon <= AOI[2]) & (lat >= AOI[1]) & (lat <= AOI[3])
    if not inb.any():
        raise SystemExit("AOI does not intersect the reference geometry")
    rows = np.where(inb.any(axis=1))[0]
    cols = np.where(inb.any(axis=0))[0]
    # full-res, quantized outward so the AOI is fully contained
    lineOff = _snap_down(rows[0] * NALKS1, LINE_Q)
    lineEnd = _snap_up((rows[-1] + 1) * NALKS1, LINE_Q)
    sampOff = _snap_down(cols[0], SAMP_Q)
    sampEnd = _snap_up(cols[-1] + 1, SAMP_Q)
    aoi_window = [int(lineOff), int(lineEnd), int(sampOff), int(sampEnd)]

    sw = {sn: d for sn, d in g["swaths"].items()}
    order = sorted(sw, key=lambda sn: sw[sn]["samples_from_origin"])

    # azimuth: clip to the extent every swath covers (frame-grid lines)
    az_lo = max(d["lines_from_origin"] for d in sw.values())
    az_hi = min(d["lines_from_origin"] + d["res"]["numberOfLines"] for d in sw.values())
    lineOff = max(lineOff, _snap_up(az_lo, LINE_Q))
    lineEnd = min(lineEnd, _snap_down(az_hi, LINE_Q))
    if lineEnd - lineOff < LINE_Q:
        raise SystemExit("derive_window: no azimuth window survives clipping to the common swath extent")
    az_clipped = [lineOff, lineEnd] != aoi_window[:2]

    # range: keep every overlap strip the window touches whole (frame-grid samples)
    ext = {sn: (sw[sn]["samples_from_origin"],
                sw[sn]["samples_from_origin"] + sw[sn]["res"]["numberOfSamples"]) for sn in sw}
    frame_w = max(e[1] for e in ext.values())
    strips = []
    for a, b in zip(order, order[1:]):
        lo, hi = ext[b][0], ext[a][1]
        if hi <= lo:
            raise SystemExit(f"derive_window: swaths {a}/{b} do not overlap in range")
        strips.append((a, b, lo, hi))
    for a, b, lo, hi in strips:
        if sampOff < hi and sampEnd > lo:  # touched: include whole
            sampOff = min(sampOff, _snap_down(lo, SAMP_Q))
            sampEnd = max(sampEnd, _snap_up(hi, SAMP_Q))
    sampOff = max(sampOff, 0)
    sampEnd = min(sampEnd, _snap_down(frame_w, SAMP_Q))
    rg_widened = [sampOff, sampEnd] != aoi_window[2:]

    win = {
        "lineOff": int(lineOff), "cropLen": int(lineEnd - lineOff),
        "sampOff_frame": int(sampOff), "cropWidth": int(sampEnd - sampOff),
        "ml1": {"rowOff": int(lineOff // NALKS1), "rows": int((lineEnd - lineOff) // NALKS1),
                "colOff": int(sampOff), "cols": int(sampEnd - sampOff)},
        "ml2": {"rowOff": int(lineOff // LINE_Q), "rows": int((lineEnd - lineOff) // LINE_Q),
                "colOff": int(sampOff // SAMP_Q), "cols": int(sampEnd - sampOff) // SAMP_Q},
        "aoi": AOI, "aoi_window": aoi_window,
        "az_clipped": bool(az_clipped), "rg_widened": bool(rg_widened),
        "swaths": {},
    }
    assert win["lineOff"] % LINE_Q == 0 and win["cropLen"] % LINE_Q == 0
    assert win["sampOff_frame"] % SAMP_Q == 0 and win["cropWidth"] % SAMP_Q == 0
    for sn in order:
        d = sw[sn]
        rw = swath_window(win, d["res"]["startingRange"], g["frame_origin"]["startingRange"],
                          g["frame_origin"]["rangePixelSize"], d["res"]["numberOfSamples"], sn)
        lo_s = swath_line_offset(win, d["res"]["sensingStart"], g["frame_origin"]["sensingStart"],
                                 g["frame_origin"]["azimuthLineInterval"], d["res"]["numberOfLines"], sn)
        win["swaths"][sn] = {"lineOff": lo_s, "length": win["cropLen"],
                             "sampOff": rw["sampOff"], "width": rw["width"],
                             "frame_sampOff": d["samples_from_origin"] + rw["sampOff"],
                             "lines_from_origin": d["lines_from_origin"]}
    lo_f = min(w["frame_sampOff"] for w in win["swaths"].values())
    hi_f = max(w["frame_sampOff"] + w["width"] for w in win["swaths"].values())
    if (lo_f, hi_f) != (win["sampOff_frame"], win["sampOff_frame"] + win["cropWidth"]):
        # mosaic_parameter's re-run and swathMosaic's outWidth both depend on
        # the per-swath edges reproducing the frame window exactly
        raise SystemExit(f"derive_window: per-swath windows span [{lo_f},{hi_f}) != frame window "
                         f"[{win['sampOff_frame']},{win['sampOff_frame'] + win['cropWidth']})")
    for a, b, lo, hi in strips:  # the overlap strips must survive whole
        wa, wb = win["swaths"][a], win["swaths"][b]
        if not (wa["frame_sampOff"] <= lo and wa["frame_sampOff"] + wa["width"] >= hi
                and wb["frame_sampOff"] <= lo and wb["frame_sampOff"] + wb["width"] >= hi):
            raise SystemExit(f"derive_window: overlap strip {a}/{b} not whole after cropping")
    with open(out, "w") as f:
        json.dump(win, f, indent=1)
    print(f"window: lines {win['lineOff']}+{win['cropLen']}, samples "
          f"{win['sampOff_frame']}+{win['cropWidth']} (full-res, frame grid); "
          f"az_clipped={az_clipped} rg_widened={rg_widened}; per swath "
          + ", ".join(f"s{sn}: lines {w['lineOff']}+{w['length']} samples {w['sampOff']}+{w['width']}"
                      for sn, w in win["swaths"].items()))
    return win


def load_window(path=None):
    return json.load(open(path or HOME_WINDOW))


# ------------------------------------------------- sec-2 metadata rewrite

def rewrite_resampled(res_ref_dir, ref_date, geom_path=None, window_path=None):
    """Spec sec 2: rewrite the ONE target, dates_resampled/<ref>/<ref>.track
    .xml and f1_*/<ref>.frame.xml, to the cropped geometry. Per swath only:
    numberOfSamples = width_s, numberOfLines = cropLen, startingRange +=
    sampOff_s*rangePixelSize, sensingStart += lineOff_s*azimuthLineInterval
    (which is T0 for every swath, Amendment A.1). Everything else, including
    every frame- and track-level field, is left to the mosaic_parameter.py
    re-run the caller must do next; nothing here is hand-written.

    Runs AFTER pack --reference: the packing needs the pre-crop dims."""
    g = load_geometry(geom_path)
    win = load_window(window_path)
    if ref_date != g["ref_date"]:
        raise SystemExit(f"rewrite: {ref_date} is not the geometry's reference date {g['ref_date']}")
    loadTrack, saveTrack = _stackpulic("loadTrack", "saveTrack")
    d = os.path.abspath(res_ref_dir)
    track = loadTrack(d, ref_date)
    if len(track.frames) != 1:
        raise SystemExit("rewrite: single-frame stacks only")
    swaths = track.frames[0].swaths
    if len(swaths) != len(g["swaths"]):
        raise SystemExit(f"rewrite: {len(swaths)} swaths in the track, "
                         f"{len(g['swaths'])} in the geometry")
    t0 = _t(g["frame_origin"]["sensingStart"]) + datetime.timedelta(
        seconds=win["lineOff"] * g["frame_origin"]["azimuthLineInterval"])
    plan = []
    for sw in swaths:
        sn = str(int(sw.swathNumber))
        if sn not in win["swaths"]:
            raise SystemExit(f"rewrite: swath {sn} is not in the window file")
        w, gs = win["swaths"][sn], g["swaths"][sn]
        # the tree must still be the pre-crop one the geometry was built from
        if int(sw.numberOfSamples) != int(gs["res"]["numberOfSamples"]) or \
                int(sw.numberOfLines) != int(gs["res"]["numberOfLines"]):
            raise SystemExit(
                f"rewrite: swath {sn} is {sw.numberOfSamples}x{sw.numberOfLines} but the geometry "
                f"recorded {gs['res']['numberOfSamples']}x{gs['res']['numberOfLines']} "
                "(already rewritten, or not the tree the geometry was built from)")
        if abs(sw.startingRange - gs["res"]["startingRange"]) > 1e-6 or \
                abs((sw.sensingStart - _t(gs["res"]["sensingStart"])).total_seconds()) > 1.5e-6:
            raise SystemExit(f"rewrite: swath {sn} origin does not match the geometry")
        if w["lineOff"] < 0 or w["lineOff"] + w["length"] > int(sw.numberOfLines) or \
                w["sampOff"] < 0 or w["sampOff"] + w["width"] > int(sw.numberOfSamples):
            raise SystemExit(f"rewrite: swath {sn} window is outside the resampled swath")
        r0 = sw.startingRange + w["sampOff"] * sw.rangePixelSize
        ts = sw.sensingStart + datetime.timedelta(seconds=w["lineOff"] * sw.azimuthLineInterval)
        # sec 2: sensingStart_s = T0 exactly for every swath
        if abs((ts - t0).total_seconds()) > 1.5e-6:
            raise SystemExit(f"rewrite: swath {sn} cropped sensingStart {ts} != window T0 {t0}")
        plan.append((sw, sn, r0, ts, w))
    for sw, sn, r0, ts, w in plan:
        sw.startingRange = float(r0)
        sw.sensingStart = ts
        sw.numberOfSamples = int(w["width"])
        sw.numberOfLines = int(w["length"])
    cwd = os.getcwd()
    try:
        os.chdir(d)
        saveTrack(track, ref_date)
    finally:
        os.chdir(cwd)
    back = loadTrack(d, ref_date)
    for sw in back.frames[0].swaths:
        sn = str(int(sw.swathNumber))
        w = win["swaths"][sn]
        if int(sw.numberOfSamples) != w["width"] or int(sw.numberOfLines) != w["length"] or \
                abs((sw.sensingStart - t0).total_seconds()) > 1.5e-6:
            raise SystemExit(f"rewrite: swath {sn} did not persist as written")
    print("rewrite: {} swaths cropped to {} lines from T0 {}; samples {}".format(
        len(plan), win["cropLen"], t0.strftime(TFMT),
        {sn: w["width"] for sn, w in win["swaths"].items()}))
    return {"T0": t0.strftime(TFMT), "swaths": {sn: {"startingRange": r0, "width": w["width"],
                                                     "length": w["length"]}
                                                for _, sn, r0, _, w in plan}}


def verify_cropped(res_ref_dir, ref_date, geom_path=None, window_path=None):
    """Spec sec 1 assertion, run AFTER the rewrite and the mosaic_parameter.py
    re-run: the cropped frame must be cropWidth x cropLen/NALKS1 at ml1, its
    origin the window origin, and the track must agree with the frame."""
    g = load_geometry(geom_path)
    win = load_window(window_path)
    loadTrack, = _stackpulic("loadTrack")
    track = loadTrack(os.path.abspath(res_ref_dir), ref_date)
    fr = track.frames[0]
    t0 = _t(g["frame_origin"]["sensingStart"]) + datetime.timedelta(
        seconds=win["lineOff"] * g["frame_origin"]["azimuthLineInterval"])
    r0 = g["frame_origin"]["startingRange"] + win["sampOff_frame"] * g["frame_origin"]["rangePixelSize"]
    want_w, want_l = win["cropWidth"] // NRLKS1, win["cropLen"] // NALKS1
    bad = []
    if int(fr.numberOfSamples) != want_w or int(fr.numberOfLines) != want_l:
        bad.append(f"frame ml1 {fr.numberOfSamples}x{fr.numberOfLines} != {want_w}x{want_l}")
    if abs(fr.startingRange - r0) > 1e-6:
        bad.append(f"frame startingRange {fr.startingRange} != {r0}")
    if abs((fr.sensingStart - t0).total_seconds()) > 1.5e-6:
        bad.append(f"frame sensingStart {fr.sensingStart} != {t0}")
    if int(track.numberOfSamples) != int(fr.numberOfSamples) or \
            int(track.numberOfLines) != int(fr.numberOfLines):
        bad.append("track dims != frame dims (single-frame stack)")
    if abs(track.startingRange - fr.startingRange) > 1e-6 or \
            abs((track.sensingStart - fr.sensingStart).total_seconds()) > 1.5e-6:
        bad.append("track origin != frame origin (single-frame stack)")
    for sw in fr.swaths:
        sn = str(int(sw.swathNumber))
        w = win["swaths"][sn]
        if int(sw.numberOfSamples) != w["width"] or int(sw.numberOfLines) != w["length"]:
            bad.append(f"swath {sn} {sw.numberOfSamples}x{sw.numberOfLines} != "
                       f"{w['width']}x{w['length']}")
    if bad:
        raise SystemExit("verify: " + "; ".join(bad))
    print(f"verify: cropped frame ml1 {want_w}x{want_l}, origin {r0:.3f} m / {t0.strftime(TFMT)}, "
          f"{len(fr.swaths)} swaths match the window")
    return {"width": want_w, "length": want_l}


def swath_window(win, swath_startingRange, frame_startingRange, rangePixelSize, swath_samples, sn=""):
    """Per-swath range window from the absolute frame edges (spec sec 1)."""
    off0 = _near_int((swath_startingRange - frame_startingRange) / rangePixelSize)
    s0 = max(win["sampOff_frame"] - off0, 0)
    s1 = min(win["sampOff_frame"] + win["cropWidth"] - off0, int(swath_samples))
    if s1 <= s0:
        raise SystemExit(f"swath {sn} lies entirely outside the range window; the contract "
                         "persists and rewrites every swath, so the window must intersect all")
    # no per-swath quantum: range origins snap to nrlks1 = 1, so sampOff_s and
    # width_s are any integers; the 5-sample quantum is frame-level only
    return {"sampOff": int(s0), "width": int(s1 - s0)}


def swath_line_offset(win, swath_sensingStart, frame_sensingStart, azimuthLineInterval, swath_lines, sn=""):
    """Per-swath azimuth start (spec sec 1/2): lineOff_s = lineOff - 14*k_s."""
    k = _lines_from_dt((_t(swath_sensingStart) - _t(frame_sensingStart)).total_seconds(), azimuthLineInterval)
    if k < 0 or k % NALKS1 != 0:
        raise SystemExit(f"swath {sn} starts {k} lines after the frame origin; expected a "
                         f"non-negative multiple of {NALKS1}")
    lo = win["lineOff"] - k
    if lo < 0 or lo + win["cropLen"] > int(swath_lines):
        raise SystemExit(f"swath {sn}: azimuth window {lo}+{win['cropLen']} outside "
                         f"[0, {swath_lines})")
    if lo % NALKS1:
        raise SystemExit(f"swath {sn}: lineOff_s {lo} not a multiple of {NALKS1}")
    return int(lo)


# ------------------------------------------------------- offset polynomials

def _poly_dict(coeffs, meanAzimuth=0.0, normAzimuth=1.0, meanRange=0.0, normRange=1.0):
    return {"coeffs": [[float(c) for c in row] for row in coeffs],
            "meanAzimuth": float(meanAzimuth), "normAzimuth": float(normAzimuth),
            "meanRange": float(meanRange), "normRange": float(normRange)}


def _from_poly2d(p):
    return _poly_dict(p.getCoeffs(), p.getMeanAzimuth(), p.getNormAzimuth(),
                      p.getMeanRange(), p.getNormRange())


def eval_poly(p, az, rg):
    """Evaluate a serialised polynomial at (az, rg); see POLY_EVAL for the
    coordinate convention. Vectorised over numpy inputs."""
    y = (np.asarray(az, dtype=np.float64) - p["meanAzimuth"]) / p["normAzimuth"]
    x = (np.asarray(rg, dtype=np.float64) - p["meanRange"]) / p["normRange"]
    out = np.zeros(np.broadcast(y, x).shape, dtype=np.float64)
    for i, row in enumerate(p["coeffs"]):
        for j, c in enumerate(row):
            if c != 0.0:
                out = out + c * (y ** i) * (x ** j)
    return out


def _read_cull(path):
    """cull.off rows as written by Alos2ProcPublic.writeOffset:
    x dx y dy snr sigmax sigmay sigmaxy (x, y full-res original-grid pixels)."""
    rows = []
    for ln in open(path):
        if ln.strip() == "" or "*" in ln:
            continue
        p = ln.split()
        rows.append([float(v) for v in p[:8]])
    if len(rows) < 6:
        raise SystemExit(f"{path}: {len(rows)} offsets, need at least 6 for the order-2 fit")
    return np.array(rows, dtype=np.float64)


def fit_offset_polys(cull_off_path, orig, res):
    """Reproduce the offset polynomials resample_common_grid.py applied to a
    secondary date's swath: map each cull.off point onto the resampled swath
    grid, fit the TOTAL offset (grid snap + misregistration) with
    OffsetField.getFitPolynomials(2, 2, maxOrder=True), same weights, same
    azimuth mean shift. orig/res: the swath's offset-adjusted original and
    resampled parameter dicts from ref_geometry.json."""
    a = _read_cull(cull_off_path)
    x, dx, y, dy, snr = a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4]
    x1 = (x * orig["rangePixelSize"] / res["rangePixelSize"]
          + (orig["startingRange"] - res["startingRange"]) / res["rangePixelSize"])
    y1 = (y * orig["azimuthLineInterval"] / res["azimuthLineInterval"]
          + (_t(orig["sensingStart"]) - _t(res["sensingStart"])).total_seconds()
          / res["azimuthLineInterval"])
    dxt = (x + dx) - x1
    dyt = (y + dy) - y1
    try:
        import isce  # noqa: F401
        from isceobj.Location.Offset import OffsetField, Offset
    except ImportError:
        if os.environ.get("CROP_PACK_NUMPY_FIT") != "1":
            raise SystemExit("fit_offset_polys: ISCE not importable; set CROP_PACK_NUMPY_FIT=1 "
                             "to use the numpy replica of offsetpoly.f90 (recorded in the sidecar)")
        rg, az = _numpy_offsetpoly(x1, y1, dxt, dyt, snr)
        fit = "numpy"
    else:
        field = OffsetField()
        for i in range(len(x1)):
            o = Offset()
            o.setCoordinate(float(x1[i]), float(y1[i]))
            o.setOffset(float(dxt[i]), float(dyt[i]))
            o.setSignalToNoise(float(snr[i]))
            o.setCovariance(float(a[i, 5]), float(a[i, 6]), float(a[i, 7]))
            field.addOffset(o)
        azp, rgp = field.getFitPolynomials(rangeOrder=2, azimuthOrder=2, maxOrder=True, usenumpy=False)
        rg, az = _from_poly2d(rgp), _from_poly2d(azp)
        fit = "isce"
    return {"rg": rg, "az": az, "fit": fit, "n": int(len(x1)), "source": os.path.abspath(cull_off_path),
            "eval": POLY_EVAL}


def _numpy_offsetpoly(x, y, dx, dy, snr):
    """Replica of stdproc offsetpoly.f90 with 6 coefficients: weighted least
    squares, sigma = 1 + 1/snr, azimuth shifted by its minimum, coefficient
    order const, rg, az, rg*az, rg**2, az**2 unpacked as getFitPolynomials does."""
    azmin = float(np.min(y))
    ys = y - azmin
    sig = 1.0 + 1.0 / snr
    A = np.column_stack([np.ones_like(x), x, ys, x * ys, x * x, ys * ys]) / sig[:, None]
    out = []
    for off in (dx, dy):
        v, *_ = np.linalg.lstsq(A, off / sig, rcond=None)
        out.append(_poly_dict([[v[0], v[1], v[4]], [v[2], v[3], 0.0], [v[5], 0.0, 0.0]],
                              meanAzimuth=azmin))
    return out[0], out[1]


def reference_polys(orig, res):
    """The analytic polynomials resample_common_grid.py applies to the reference
    date (idate == dateIndexReference branch): grid-snap constant plus
    pixel-size stretch, rangeOrder 1 / azimuthOrder 0 and 0 / 1."""
    rg = _poly_dict([[(res["startingRange"] - orig["startingRange"]) / orig["rangePixelSize"],
                      res["rangePixelSize"] / orig["rangePixelSize"] - 1.0]])
    az = _poly_dict([[(_t(res["sensingStart"]) - _t(orig["sensingStart"])).total_seconds()
                      / orig["azimuthLineInterval"]],
                     [res["azimuthLineInterval"] / orig["azimuthLineInterval"] - 1.0]])
    return {"rg": rg, "az": az, "fit": "analytic", "n": 0, "source": "reference", "eval": POLY_EVAL}


# -------------------------------------------------------------- pack/unpack

def pack_slc(slc_path, out_path, sidecar_path, lineOff, sampOff, length, width,
             in_width, in_length, offset_polys, swath=None, frame_lineOff=None,
             frame_sampOff=None, pctl=99.995, date=None):
    """Crop a complex64 SLC (raw BIL binary, no headers, in_width samples per
    line, in_length lines) and pack it as int16 I,Q with a JSON sidecar."""
    if offset_polys is None:
        raise SystemExit("pack_slc: offset_polys is required (spec sec 3); no sidecar ships without them")
    if not (0 <= lineOff and 0 <= sampOff and lineOff + length <= in_length and sampOff + width <= in_width):
        raise SystemExit(f"pack_slc: window lines {lineOff}+{length} samples {sampOff}+{width} "
                         f"outside the {in_width}x{in_length} input")
    if length <= 0 or width <= 0:
        raise SystemExit("pack_slc: empty window")
    expect = in_width * in_length * 8
    if os.path.getsize(slc_path) != expect:
        raise SystemExit(f"pack_slc: {slc_path} is {os.path.getsize(slc_path)} bytes, "
                         f"expected {expect} for {in_width}x{in_length} complex64")
    n_line_bytes = in_width * 8
    out = np.empty((length, width), dtype=np.complex64)
    with open(slc_path, "rb") as f:
        for i in range(length):
            f.seek((lineOff + i) * n_line_bytes + sampOff * 8)
            buf = f.read(width * 8)
            if len(buf) != width * 8:
                raise SystemExit(f"pack_slc: short read at line {lineOff + i}")
            out[i] = np.frombuffer(buf, dtype=np.complex64)
    if not np.isfinite(out).all():
        raise SystemExit(f"pack_slc: {slc_path} window holds non-finite samples")
    amp = np.abs(out)
    nz = amp > 0
    peak = float(np.percentile(amp[nz], pctl)) if nz.any() else 1.0
    max_amp = float(amp.max()) if nz.any() else 0.0
    scale = peak / 32000.0 if peak > 0 else 1.0
    # int16 quantum from the pctl amplitude; the few pixels above the int16
    # range (the brightest scatterers, by construction ~1 - pctl of the swath)
    # are NOT clipped or saturated (either would change their phase or their
    # weight in a multilooked cell): they are stored exactly (complex64) in
    # the .ovf side file and the int16 pair holds the clipped value only as a
    # placeholder. unpack_slc restores them.
    packed = np.empty((length, width, 2), dtype=np.int16)
    q = np.round(out.real / scale)
    ovf = np.abs(q) > 32767
    packed[..., 0] = np.clip(q, -32767, 32767)
    q = np.round(out.imag / scale)
    ovf |= np.abs(q) > 32767
    packed[..., 1] = np.clip(q, -32767, 32767)
    del q
    idx = np.flatnonzero(ovf).astype(np.int64)
    vals = out.reshape(-1)[idx]
    del ovf
    ovf_path = out_path + ".ovf"
    tmp = ovf_path + ".tmp"
    with open(tmp, "wb") as f:
        idx.tofile(f)
        vals.tofile(f)
    os.replace(tmp, ovf_path)
    tmp = out_path + ".tmp"
    packed.tofile(tmp)
    os.replace(tmp, out_path)
    # the sidecar is the contract's single source of truth for this crop, so it
    # is published atomically like the binaries: a lane killed by the SLURM time
    # limit mid-write would otherwise leave a truncated sidecar beside two
    # complete files, and the crop cannot be read back without it
    tmp = sidecar_path + ".tmp"
    with open(tmp, "w") as f:
        json.dump({"swath": swath, "width": int(width), "length": int(length),
                   "sampOff": int(sampOff), "lineOff": int(lineOff),
                   "frame_lineOff": frame_lineOff, "frame_sampOff": frame_sampOff,
                   "in_width": int(in_width), "in_length": int(in_length),
                   "scale": scale, "dtype": "int16 I,Q", "pctl": pctl, "peak": peak,
                   "max_amp": max_amp, "nonzero_fraction": float(nz.mean()),
                   "overflow_count": int(idx.size), "overflow_file": os.path.basename(ovf_path),
                   "overflow_layout": "int64 row-major flat indices then complex64 values, "
                                      "exact samples above the int16 range",
                   "source": os.path.abspath(slc_path), "date": date,
                   "offset_polys": offset_polys}, f, indent=1)
    os.replace(tmp, sidecar_path)
    return scale


def unpack_slc(packed_path, sidecar_path, out_path):
    """int16 -> complex64 for block-job staging (spec 5b); the overflow
    pixels come back exactly from the .ovf side file."""
    sc = json.load(open(sidecar_path))
    length, width = sc["length"], sc["width"]
    expect = length * width * 4
    if os.path.getsize(packed_path) != expect:
        raise SystemExit(f"unpack_slc: {packed_path} is {os.path.getsize(packed_path)} bytes, "
                         f"expected {expect} for {width}x{length} int16 I,Q")
    if "overflow_count" not in sc or "overflow_file" not in sc:
        raise SystemExit(f"unpack_slc: {sidecar_path} has no overflow record (pre-round-12 pack?)")
    n = int(sc["overflow_count"])
    ovf_path = os.path.join(os.path.dirname(os.path.abspath(packed_path)), sc["overflow_file"])
    if not os.path.isfile(ovf_path) or os.path.getsize(ovf_path) != n * 16:
        raise SystemExit(f"unpack_slc: {ovf_path} missing or not {n * 16} bytes for "
                         f"{n} overflow pixels")
    raw = np.fromfile(packed_path, dtype=np.int16).reshape(length, width, 2)
    out = np.empty((length, width), dtype=np.complex64)
    out.real = raw[..., 0]
    out.imag = raw[..., 1]
    out *= np.float32(sc["scale"])
    if n:
        with open(ovf_path, "rb") as f:
            idx = np.fromfile(f, dtype=np.int64, count=n)
            vals = np.fromfile(f, dtype=np.complex64, count=n)
        if idx.size != n or vals.size != n or idx.min() < 0 or idx.max() >= length * width \
                or np.any(np.diff(idx) <= 0):
            raise SystemExit(f"unpack_slc: {ovf_path} indices are not a sorted subset of the crop")
        if not np.all(np.abs(raw.reshape(-1, 2)[idx]).max(axis=1) == 32767):
            raise SystemExit(f"unpack_slc: {ovf_path} does not belong to {packed_path} "
                             "(overflow pixels are not at the int16 limit)")
        out.reshape(-1)[idx] = vals
    del raw
    tmp = out_path + ".tmp"
    out.tofile(tmp)
    os.replace(tmp, out_path)
    return sc


ITEMSIZE = {"float32": 4, "float": 4, "double": 8, "float64": 8, "byte": 1, "int8": 1}
ISCE_DTYPE = {"float32": "FLOAT", "float": "FLOAT", "double": "DOUBLE", "float64": "DOUBLE",
              "byte": "BYTE", "int8": "BYTE"}


def crop_bil_raster(in_path, out_path, rowOff, rows, colOff, cols, in_width, in_length,
                    bands=1, dtype="float32"):
    """Crop a BIL raster by window, band by band. dtype: float32 (rg_rect.off,
    los 2-band), double (ml1 lat/lon/hgt), byte (wbd). Symlinks (rg_rect.off is
    a link to rg.off under an identity affine) are dereferenced."""
    if dtype not in ITEMSIZE:
        raise SystemExit(f"crop_bil_raster: unknown dtype {dtype}; one of {sorted(ITEMSIZE)}")
    isz = ITEMSIZE[dtype]
    src = os.path.realpath(in_path)
    if not os.path.isfile(src):
        raise SystemExit(f"crop_bil_raster: {in_path} -> {src} is not a regular file (dangling link?)")
    lb = in_width * isz * bands
    if os.path.getsize(src) != lb * in_length:
        raise SystemExit(f"crop_bil_raster: {in_path} size {os.path.getsize(src)} != "
                         f"{in_width}x{in_length}x{bands} {dtype}")
    if not (0 <= rowOff and 0 <= colOff and rowOff + rows <= in_length and colOff + cols <= in_width):
        raise SystemExit(f"crop_bil_raster: window rows {rowOff}+{rows} cols {colOff}+{cols} "
                         f"outside {in_width}x{in_length}")
    tmp = out_path + ".tmp"
    with open(src, "rb") as fi, open(tmp, "wb") as fo:
        for i in range(rows):
            base = (rowOff + i) * lb
            for b in range(bands):
                fi.seek(base + b * in_width * isz + colOff * isz)
                fo.write(fi.read(cols * isz))
    os.replace(tmp, out_path)


crop_float_raster = crop_bil_raster  # former name (float32, 1 band)


ISCE_ITEMSIZE = {"BYTE": 1, "SHORT": 2, "INT": 4, "FLOAT": 4, "DOUBLE": 8, "CFLOAT": 8, "CDOUBLE": 16}


def render_isce_xml(binary_path, width, length, dtype="CFLOAT", bands=1, scheme="BIL"):
    """Render .xml/.vrt for an unpacked/cropped ISCE image (spec 5b). dtype is
    a cropoff name (float32/double/byte, ...) or an ISCE name; validated, with
    the binary size, before anything is written (renderHdr writes the .xml
    first and only then renderVRT looks the type up).

    The xml records the BARE BASENAME with cwd = the binary's directory, as
    upstream create_xml does: consumers that open the data through the xml
    (look_geom hands the loaded los image to mroipac Looks, which calls
    createImage() without resetting the filename) would otherwise follow the
    render-time path and read the persist tree instead of the staged
    node-local copy (Amendment A.17)."""
    dt = ISCE_DTYPE.get(str(dtype).lower(), str(dtype).upper())
    if dt not in ISCE_ITEMSIZE:
        raise SystemExit(f"render: unknown dtype {dtype}; one of {sorted(ISCE_DTYPE)} or "
                         f"{sorted(ISCE_ITEMSIZE)}")
    sch = str(scheme).upper()
    if sch not in ("BIL", "BIP", "BSQ"):
        raise SystemExit(f"render: unknown scheme {scheme}; BIL, BIP or BSQ")
    if width <= 0 or length <= 0 or bands <= 0:
        raise SystemExit(f"render: bad dims {width}x{length}x{bands}")
    expect = width * length * bands * ISCE_ITEMSIZE[dt]
    if not os.path.isfile(binary_path) or os.path.getsize(binary_path) != expect:
        raise SystemExit(f"render: {binary_path} is not {expect} bytes "
                         f"({width}x{length}x{bands} {dt})")
    import isce  # noqa: F401
    import isceobj
    d, base = os.path.split(os.path.abspath(binary_path))
    img = isceobj.createImage()
    img.setFilename(base)
    img.setWidth(width)
    img.setLength(length)
    img.dataType = dt
    img.bands = bands
    img.scheme = sch
    img.setAccessMode("READ")
    cwd = os.getcwd()
    try:
        os.chdir(d)
        img.renderHdr()
    finally:
        os.chdir(cwd)
    return dt


# --------------------------------------------------------------------- CLI

USAGE = """crop_pack.py commands
  geometry <idir> <resampled_odir> <ref_date>
  derive   <ref_ml1_lat> <ref_ml1_lon> <ml1_width> <ml1_length>
  pack     <slc> <out> <sidecar> <swath> (--cull <cull.off> | --reference) [--date <yymmdd>]
           (date and swath are taken from the SLC path .../<date>/f1_*/s<n>/<date>.slc
           unless --date is given; the cull.off path must belong to the same date/swath)
  unpack   <packed> <sidecar> <out>
  cropoff  <in> <out> <rowOff> <rows> <colOff> <cols> <in_width> <in_length> [bands] [dtype]
           (dtype float32|double|byte, default float32; BIL)
  render   <binary> <width> <length> [dtype] [bands] [scheme]
           (dtype a cropoff name or ISCE FLOAT|DOUBLE|BYTE|CFLOAT, default CFLOAT;
           bands default 1; scheme BIL|BIP|BSQ, default BIL; the binary size
           must match)
  rewrite  <dates_resampled/<ref>> <ref_date>
           (sec-2 metadata rewrite to the cropped geometry, per swath only;
           re-run mosaic_parameter.py on the tree afterwards)
  verify   <dates_resampled/<ref>> <ref_date>
           (sec-1 assertion after the rewrite and the mosaic_parameter re-run)"""


def _slc_identity(slc_path, argv):
    """(date, swath) of an SLC from its path or --date; None where unknown."""
    date = argv[argv.index("--date") + 1] if "--date" in argv else None
    base = os.path.basename(slc_path)
    m = re.match(r"^(\d{6})(_lower|_upper)?\.slc$", base)
    if m and date is None:
        date = m.group(1)
    if m and date is not None and m.group(1) != date:
        raise SystemExit(f"pack: --date {date} does not match the SLC name {base}")
    parent = os.path.basename(os.path.dirname(os.path.abspath(slc_path)))
    m2 = re.match(r"^s(\d+)$", parent)
    return date, (str(int(m2.group(1))) if m2 else None)


def _cull_belongs_to(cull, date, sn):
    """True when cull is .../<date>/f<i>_<frame>/s<sn>/cull.off by position
    (the upstream layout, estimate_slc_offset.py); membership of the path
    components would also match a <date> or s<n> directory anywhere above,
    including the job cwd."""
    parts = os.path.normpath(os.path.abspath(cull)).split(os.sep)
    return (len(parts) >= 4 and parts[-1] == "cull.off" and parts[-2] == f"s{sn}"
            and re.fullmatch(r"f\d+_\d+", parts[-3]) is not None and parts[-4] == date)


def main(argv):
    if len(argv) < 2:
        raise SystemExit(USAGE)
    cmd = argv[1]
    if cmd == "geometry":
        build_geometry(argv[2], argv[3], argv[4])
    elif cmd == "derive":
        derive_window(argv[2], argv[3], int(argv[4]), int(argv[5]))
    elif cmd == "pack":
        slc, out, sidecar, sn = argv[2], argv[3], argv[4], str(int(argv[5]))
        win, g = load_window(), load_geometry()
        if sn not in win["swaths"] or sn not in g["swaths"]:
            raise SystemExit(f"pack: swath {sn} not in the window/geometry files")
        w, d = win["swaths"][sn], g["swaths"][sn]
        date, path_sn = _slc_identity(slc, argv)
        if date is None:
            raise SystemExit("pack: cannot tell the date from the SLC path; give --date <yymmdd>")
        if path_sn is not None and path_sn != sn:
            raise SystemExit(f"pack: SLC path is swath s{path_sn} but swath {sn} was given")
        if "--cull" in argv:
            if date == g["ref_date"]:
                raise SystemExit(f"pack: {date} is the reference date; use --reference")
            cull = argv[argv.index("--cull") + 1]
            if not _cull_belongs_to(cull, date, sn):
                raise SystemExit(f"pack: {cull} does not belong to date {date} swath s{sn} "
                                 "(expected .../<date>/f1_*/s<n>/cull.off)")
            polys = fit_offset_polys(cull, d["orig"], d["res"])
        elif "--reference" in argv:
            if date != g["ref_date"]:
                raise SystemExit(f"pack: --reference given for {date}, but the reference date "
                                 f"is {g['ref_date']}")
            polys = reference_polys(d["orig"], d["res"])
        else:
            raise SystemExit("pack: give --cull <cull.off> (secondary) or --reference")
        pack_slc(slc, out, sidecar, w["lineOff"], w["sampOff"], w["length"], w["width"],
                 d["res"]["numberOfSamples"], d["res"]["numberOfLines"], polys, swath=int(sn),
                 frame_lineOff=win["lineOff"], frame_sampOff=w["frame_sampOff"], date=date)
    elif cmd == "unpack":
        unpack_slc(argv[2], argv[3], argv[4])
    elif cmd == "cropoff":
        crop_bil_raster(argv[2], argv[3], int(argv[4]), int(argv[5]), int(argv[6]), int(argv[7]),
                        int(argv[8]), int(argv[9]), int(argv[10]) if len(argv) > 10 else 1,
                        argv[11] if len(argv) > 11 else "float32")
    elif cmd == "rewrite":
        rewrite_resampled(argv[2], argv[3])
    elif cmd == "verify":
        verify_cropped(argv[2], argv[3])
    elif cmd == "render":
        render_isce_xml(argv[2], int(argv[3]), int(argv[4]),
                        argv[5] if len(argv) > 5 else "CFLOAT",
                        int(argv[6]) if len(argv) > 6 else 1,
                        argv[7] if len(argv) > 7 else "BIL")
    else:
        raise SystemExit(f"unknown command {cmd}\n{USAGE}")


if __name__ == "__main__":
    main(sys.argv)
