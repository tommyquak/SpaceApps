"""Post-run validity check for the alosStack fleet (spec secs 1-5).

Every job runs this on its own output before it may mark itself DONE, and it
is safe to re-run afterwards on any persisted date. It checks CONTENT, not
presence: a directory that merely looks finished fails here.

  check_run.py seed                 the SEED bundle (manifest E, F, G + ref crops)
  check_run.py ion <block>          one persisted ion block (what ion_ls consumes)
  check_run.py dates_ion            the per-date ion grids ion_ls wrote
  check_run.py date <YYMMDD>        one wave-1 persisted date (manifest A, B, C)
  check_run.py bundle               the once-bundle only (what wave-1 lanes read)

Exit 0 and print <STAGE>_CHECK_OK, or exit 1 listing every failure. ISCE is
used when importable (xml loadability) and skipped otherwise.
"""
import json
import os
import sys

import numpy as np

FULL = os.path.expanduser("~/full")
PERSIST = os.environ.get("PERSIST", "/mnt/scratch/t/tommy/full/stack")
ION_OUT = os.environ.get("ION_OUT", "/mnt/scratch/t/tommy/full/pairs_ion")
REF = os.environ.get("REF_DATE", "181123")
FRAME = os.environ.get("FRAME", "3700")
NRLKS1, NALKS1 = 1, 14
ML1 = "_{}rlks_{}alks".format(NRLKS1, NALKS1)
NRLKS2, NALKS2 = 5, 2
NRLKS_ION, NALKS_ION = 80, 32
ML2 = "_{}rlks_{}alks".format(NRLKS1 * NRLKS2, NALKS1 * NALKS2)
# ion_cal products live on the ION grid, not ml2: ion_ls.py:171 builds the name
# from nrlks1*nrlks_ion and nalks1*nalks_ion
ML_ION = "_{}rlks_{}alks".format(NRLKS1 * NRLKS_ION, NALKS1 * NALKS_ION)


def _xml_dims(path):
    """width, length from an ISCE xml, so the checker never has to assume a grid."""
    import xml.etree.ElementTree as ET
    d = {}
    root = ET.parse(path + ".xml").getroot()
    for prop in root.iter("property"):
        n = (prop.get("name") or "").lower()
        if n in ("width", "length"):
            v = prop.find("value")
            if v is not None:
                d[n] = int(v.text.strip())
    return d.get("width"), d.get("length")
AOI = (113.8, -3.4, 115.2, -2.34)   # w s e n
SIDECAR_KEYS = ("swath", "width", "length", "sampOff", "lineOff", "frame_lineOff",
                "frame_sampOff", "in_width", "in_length", "scale", "pctl", "peak",
                "max_amp", "nonzero_fraction", "overflow_count", "overflow_file",
                "source", "date", "offset_polys")

fails = []


def bad(msg):
    fails.append(msg)


def need(path, what):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        bad(f"{what}: {path} missing or empty")
        return False
    return True


def check_sidecar(sd, d, sn, win, band=""):
    """One packed swath: sidecar completeness, sizes, provenance, polynomials,
    then the pixels themselves. band is "" for the full band, or _lower/_upper
    for the authentic subbands wave 1 persists for the ion pass."""
    tag = f"s{sn}{band}"
    j = os.path.join(sd, f"{d}{band}.slc.json")
    p = os.path.join(sd, f"{d}{band}.slc")
    if not need(j, f"{tag} sidecar") or not need(p, f"{tag} pack"):
        return
    try:
        sc = json.load(open(j))
    except Exception as e:                      # a truncated sidecar is the failure mode
        bad(f"{tag} sidecar does not parse: {e}")
        return
    for k in SIDECAR_KEYS:
        if k not in sc:
            bad(f"{tag} sidecar has no {k}")
    if fails:
        return
    if str(sc["date"]) != d or str(sc["swath"]) != str(sn):
        bad(f"{tag} sidecar says date {sc['date']} swath {sc['swath']}")
    want = sc["width"] * sc["length"] * 4
    got = os.path.getsize(p)
    if got != want:
        bad(f"{tag} pack is {got} B, sidecar says {want} ({sc['width']}x{sc['length']} int16 I,Q)")
        return
    if win is not None:
        w = win["swaths"].get(str(sn))
        if w is None:
            bad(f"{tag} is not in the crop window")
        else:
            for k in ("width", "length", "sampOff", "lineOff"):
                if int(sc[k]) != int(w[k]):
                    bad(f"{tag} sidecar {k}={sc[k]} but the window says {w[k]} (stale crop?)")
    pol = sc["offset_polys"]
    if pol.get("fit") not in ("isce", "analytic", "numpy"):
        bad(f"{tag} polynomials have fit={pol.get('fit')}")
    for nm in ("rg", "az"):
        c = np.array(pol[nm]["coeffs"], dtype=float)
        if c.size == 0 or not np.isfinite(c).all():
            bad(f"{tag} {nm} polynomial is empty or non-finite")
    if sc["scale"] <= 0 or not np.isfinite(sc["scale"]):
        bad(f"{tag} scale is {sc['scale']}")
    if not (0.0 < sc["nonzero_fraction"] <= 1.0):
        bad(f"{tag} nonzero_fraction is {sc['nonzero_fraction']}")

    # the overflow side file must belong to THIS pack: same rule unpack enforces
    ovf = os.path.join(sd, sc["overflow_file"])
    n = int(sc["overflow_count"])
    if not os.path.isfile(ovf) or os.path.getsize(ovf) != n * 16:
        bad(f"{tag} {os.path.basename(ovf)} missing or not {n * 16} B for {n} overflow pixels")
    elif n:
        raw = np.memmap(p, dtype=np.int16, mode="r", shape=(sc["length"] * sc["width"], 2))
        with open(ovf, "rb") as f:
            idx = np.fromfile(f, dtype=np.int64, count=n)
            vals = np.fromfile(f, dtype=np.complex64, count=n)
        if idx.size != n or vals.size != n:
            bad(f"{tag} overflow file is short")
        elif idx.min() < 0 or idx.max() >= sc["length"] * sc["width"] or np.any(np.diff(idx) <= 0):
            bad(f"{tag} overflow indices are not a sorted subset of the crop")
        elif not np.all(np.abs(raw[idx]).max(axis=1) == 32767):
            bad(f"{tag} overflow file does not belong to this pack")
        elif not np.isfinite(vals.view(np.float32)).all():
            bad(f"{tag} overflow values are not finite")
        del raw

    # content: real SAR, not a zeroed or truncated crop. Sample stripes rather
    # than the whole 0.9 GB file.
    a = np.memmap(p, dtype=np.int16, mode="r", shape=(sc["length"], sc["width"], 2))
    rows = np.linspace(0, sc["length"] - 1, min(64, sc["length"])).astype(int)
    s = np.asarray(a[rows], dtype=np.float32)
    del a
    amp = np.hypot(s[..., 0], s[..., 1])
    nz = float((amp > 0).mean())
    if nz < 0.5:
        bad(f"{tag} sampled stripes are {100 * (1 - nz):.1f} percent zero (dead crop?)")
    if amp.max() == 0:
        bad(f"{tag} sampled stripes are all zero")
        return
    # the packing scale puts P99.995 of the swath at 32000 counts, so a healthy
    # crop has its bright tail in the thousands, not at 1 count or clipped flat
    p999 = float(np.percentile(amp[amp > 0], 99.9))
    if p999 < 10:
        bad(f"{tag} bright tail is only {p999:.1f} counts (wrong scale or empty swath)")
    print(f"    {tag}: {sc['width']}x{sc['length']}, {100 * nz:.1f} percent nonzero, "
          f"P99.9 {p999:.0f} counts, {n} overflow px, fit {pol['fit']}")


def check_geometry_crops(win):
    """Manifest F: the cropped ml1 geometry the block jobs stage."""
    gd = os.path.join(PERSIST, "dates_resampled", REF, "insar")
    rows, cols = win["ml1"]["rows"], win["ml1"]["cols"]
    for ext, dt, nb, itemsize in (("lat", "double", 1, 8), ("lon", "double", 1, 8),
                                  ("hgt", "double", 1, 8), ("wbd", "byte", 1, 1),
                                  ("los", "float32", 2, 4)):
        f = os.path.join(gd, f"{REF}{ML1}.{ext}")
        if not need(f, f"geometry {ext}"):
            continue
        want = rows * cols * nb * itemsize
        if os.path.getsize(f) != want:
            bad(f"geometry {ext} is {os.path.getsize(f)} B, expected {want} "
                f"({cols}x{rows}x{nb} {dt})")
            continue
        if not os.path.isfile(f + ".xml") or not os.path.isfile(f + ".vrt"):
            bad(f"geometry {ext} has no xml/vrt")
        else:
            x = open(f + ".xml").read()
            if f"<value>{os.path.basename(f)}</value>" not in x:
                bad(f"geometry {ext} xml does not carry the bare basename (not relocatable)")
    lat = np.memmap(os.path.join(gd, f"{REF}{ML1}.lat"), dtype="float64", mode="r",
                    shape=(rows, cols))
    lon = np.memmap(os.path.join(gd, f"{REF}{ML1}.lon"), dtype="float64", mode="r",
                    shape=(rows, cols))
    hgt = np.memmap(os.path.join(gd, f"{REF}{ML1}.hgt"), dtype="float64", mode="r",
                    shape=(rows, cols))
    r = np.linspace(0, rows - 1, min(64, rows)).astype(int)
    la, lo, hg = np.asarray(lat[r]), np.asarray(lon[r]), np.asarray(hgt[r])
    del lat, lon, hgt
    if not (np.isfinite(la).all() and np.isfinite(lo).all() and np.isfinite(hg).all()):
        bad("geometry lat/lon/hgt hold non-finite values")
    # the crop must actually cover the AOI it was derived for
    if la.max() < AOI[1] or la.min() > AOI[3] or lo.max() < AOI[0] or lo.min() > AOI[2]:
        bad(f"cropped geometry does not intersect the AOI: lat {la.min():.3f}..{la.max():.3f} "
            f"lon {lo.min():.3f}..{lo.max():.3f}")
    if hg.min() < -200 or hg.max() > 3000:
        bad(f"cropped heights {hg.min():.1f}..{hg.max():.1f} m are implausible for this AOI")
    print(f"    geometry {cols}x{rows}: lat {la.min():.3f}..{la.max():.3f}, "
          f"lon {lo.min():.3f}..{lo.max():.3f}, hgt {hg.min():.1f}..{hg.max():.1f} m")
    aff = os.path.join(gd, "affine_transform.txt")
    if need(aff, "affine_transform.txt"):
        # upstream format (rect_range_offset.py:70-74): line 0 is
        # "<nrlksSim> <nalksSim>", line 1 a bracketed comma-separated list of
        # the six affine parameters
        try:
            ln = open(aff).readlines()
            looks = [int(x) for x in ln[0].split()[:2]]
            vals = [float(x) for x in ln[1].strip().strip("[]").split(",")]
        except Exception as e:
            bad(f"affine_transform.txt does not parse the way rect_range_offset reads it: {e}")
        else:
            if len(looks) != 2 or min(looks) < 1:
                bad(f"affine_transform.txt sim looks are {looks}")
            if len(vals) != 6 or not all(np.isfinite(vals)):
                bad(f"affine_transform.txt does not hold 6 finite parameters: {vals}")
            # radar_dem_offset returns a near-identity transform on a good fit;
            # a wild one means the DEM/amplitude match failed
            elif abs(vals[0] - 1) > 0.05 or abs(vals[3] - 1) > 0.05 or max(abs(v) for v in vals[4:]) > 500:
                bad(f"affine_transform.txt is not a plausible near-identity fit: {vals}")
            else:
                print(f"    affine: scale {vals[0]:.6f}/{vals[3]:.6f}, shift "
                      f"{vals[4]:.3f}/{vals[5]:.3f} px, sim looks {looks}")


def _isnum(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def check_window(win, geom):
    """The window itself: what every crop in the fleet is measured against."""
    if win["lineOff"] % 28 or win["cropLen"] % 28:
        bad(f"window azimuth {win['lineOff']}+{win['cropLen']} is not on the 28 quantum")
    if win["sampOff_frame"] % 5 or win["cropWidth"] % 5:
        bad(f"window range {win['sampOff_frame']}+{win['cropWidth']} is not on the 5 quantum")
    lo = min(w["frame_sampOff"] for w in win["swaths"].values())
    hi = max(w["frame_sampOff"] + w["width"] for w in win["swaths"].values())
    if (lo, hi) != (win["sampOff_frame"], win["sampOff_frame"] + win["cropWidth"]):
        bad(f"per-swath windows span [{lo},{hi}) but the frame window is "
            f"[{win['sampOff_frame']},{win['sampOff_frame'] + win['cropWidth']})")
    for sn, w in win["swaths"].items():
        g = geom["swaths"][sn]
        if w["lineOff"] % NALKS1:
            bad(f"s{sn} lineOff {w['lineOff']} is not a multiple of {NALKS1}")
        if w["lineOff"] + w["length"] > g["res"]["numberOfLines"]:
            bad(f"s{sn} window overruns the resampled swath in azimuth")
        if w["sampOff"] + w["width"] > g["res"]["numberOfSamples"]:
            bad(f"s{sn} window overruns the resampled swath in range")
    if win["ml1"]["rows"] * NALKS1 != win["cropLen"] or win["ml1"]["cols"] != win["cropWidth"]:
        bad("the ml1 window does not match the full-res window")
    print(f"    window: lines {win['lineOff']}+{win['cropLen']}, samples "
          f"{win['sampOff_frame']}+{win['cropWidth']}, ml1 {win['ml1']['cols']}x{win['ml1']['rows']}")


def check_xmls(paths):
    try:
        sys.path.insert(0, FULL)
        import crop_pack as cp
        loadTrack, = cp._stackpulic("loadTrack")
    except Exception as e:
        print(f"    (xml load check skipped: {e})")
        return
    for d, date in paths:
        try:
            t = loadTrack(d, date)
            if len(t.frames[0].swaths) != 3:
                bad(f"{d}: track has {len(t.frames[0].swaths)} swaths, expected 3")
        except Exception as e:
            bad(f"{d}: {date} track does not load: {e}")


def check_bundle():
    print("bundle (manifest E, F, G):")
    ok = True
    for f, what in ((os.path.join(PERSIST, "GATE_REF_PUBLISHED"), "gate marker"),
                    (os.path.join(PERSIST, "dates", REF, f"{REF}.track.xml"), "ref track"),
                    (os.path.join(PERSIST, "dates", REF, f"f1_{FRAME}", f"{REF}.frame.xml"), "ref frame"),
                    (os.path.join(PERSIST, "dates", REF, f"f1_{FRAME}", "mosaic", "swath_offset.txt"),
                     "swath_offset.txt"),
                    (os.path.join(PERSIST, "ref_geometry.json"), "ref_geometry.json"),
                    (os.path.join(PERSIST, "crop_window.json"), "crop_window.json")):
        if what == "gate marker":
            if not os.path.exists(f):
                bad("gate marker missing")
                ok = False
            continue
        ok &= need(f, what)
    if not ok:
        return None, None
    geom = json.load(open(os.path.join(PERSIST, "ref_geometry.json")))
    win = json.load(open(os.path.join(PERSIST, "crop_window.json")))
    if geom["ref_date"] != REF:
        bad(f"ref_geometry.json is for {geom['ref_date']}, not {REF}")
    check_window(win, geom)
    check_xmls([(os.path.join(PERSIST, "dates", REF), REF),
                (os.path.join(PERSIST, "dates_resampled", REF), REF)])
    return geom, win


def check_seed():
    geom, win = check_bundle()
    if win is None:
        return
    print("seed reference crops (manifest A) and geometry (F):")
    for sn in ("2", "3", "4"):
        check_sidecar(os.path.join(PERSIST, "dates_resampled", REF, f"f1_{FRAME}", f"s{sn}"),
                      REF, sn, win)
    check_geometry_crops(win)


def check_date(d):
    geom, win = check_bundle()
    if win is None:
        return
    print(f"date {d} (manifest A, B, C):")
    od = os.path.join(PERSIST, "dates_resampled", d, f"f1_{FRAME}")
    for sn in ("2", "3", "4"):
        sd = os.path.join(od, f"s{sn}")
        check_sidecar(sd, d, sn, win)
        # authentic subbands, when wave 1 persisted them (PERSIST_SUBBANDS=1)
        for band in ("_lower", "_upper"):
            if os.path.exists(os.path.join(sd, f"{d}{band}.slc.json")):
                check_sidecar(sd, d, sn, win, band)
            elif os.environ.get("REQUIRE_SUBBANDS") == "1":
                bad(f"s{sn}{band} subband crop is missing (the ion pass needs it)")
    if d != REF:
        f = os.path.join(PERSIST, "dates_resampled", d, "insar", f"{d}{ML1}_rg_rect.off")
        rows, cols = win["ml1"]["rows"], win["ml1"]["cols"]
        if need(f, "rg_rect.off"):
            if os.path.getsize(f) != rows * cols * 4:
                bad(f"rg_rect.off is {os.path.getsize(f)} B, expected {rows * cols * 4}")
            else:
                a = np.memmap(f, dtype="float32", mode="r", shape=(rows, cols))
                r = np.linspace(0, rows - 1, min(64, rows)).astype(int)
                v = np.asarray(a[r])
                del a
                if not np.isfinite(v).all():
                    bad("rg_rect.off holds non-finite values")
                elif np.abs(v).max() > 10000:
                    bad(f"rg_rect.off max |offset| {np.abs(v).max():.1f} samples is implausible")
                else:
                    print(f"    rg_rect.off {cols}x{rows}: {v.min():.2f}..{v.max():.2f} samples")
            if not os.path.isfile(f + ".xml"):
                bad("rg_rect.off has no xml (diff_interferogram loads it through the xml)")
    check_xmls([(os.path.join(PERSIST, "dates", d), d)])


def check_ion(block):
    """One ion block: exactly the products ion_ls reads back (ion_ls.py:171-188,
    287), sized against their own xml and checked for finite, non-flat values."""
    pf = os.path.expanduser(f"~/full/ion_blocks/{block}.pairs")
    if not need(pf, f"ion block {block} pair list"):
        return
    pairs = [x.strip() for x in open(pf) if x.strip()]
    print(f"ion block {block}: {len(pairs)} pairs, products on the {ML_ION[1:]} grid")
    ok = 0
    stds = []
    dims = None
    for p in pairs:
        d = os.path.join(ION_OUT, p, "ion", "ion_cal")
        miss = [e for e in (".ion", ".std", ".win")
                if not os.path.isfile(os.path.join(d, f"filt_ion{ML_ION}{e}"))]
        if miss:
            bad(f"{p}: ion_cal has no filt_ion{ML_ION}{','.join(miss)}")
            continue
        base = os.path.join(d, f"filt_ion{ML_ION}.ion")
        if os.path.isfile(base + ".xml"):
            w, l = _xml_dims(base)
            if w and l:
                if dims is None:
                    dims = (w, l)
                elif dims != (w, l):
                    bad(f"{p}: grid {w}x{l} != {dims[0]}x{dims[1]} of the other pairs")
                for e in (".ion", ".std", ".win"):
                    f = os.path.join(d, f"filt_ion{ML_ION}{e}")
                    if os.path.getsize(f) != w * l * 4:
                        bad(f"{p}: filt_ion{ML_ION}{e} is {os.path.getsize(f)} B, "
                            f"expected {w * l * 4} for {w}x{l}")
        a = np.fromfile(base, dtype="float32")
        if a.size and not np.isfinite(a).all():
            bad(f"{p}: filt_ion{ML_ION}.ion holds non-finite values")
        elif a.size and float(np.abs(a).max()) == 0.0:
            bad(f"{p}: filt_ion{ML_ION}.ion is all zero")
        st = np.fromfile(os.path.join(d, f"filt_ion{ML_ION}.std"), dtype="float32")
        if st.size and np.isfinite(st).any():
            stds.append(float(np.nanmean(st[np.isfinite(st)])))
        ok += 1
    if pairs and not os.path.isfile(os.path.join(ION_OUT, pairs[0], f"{REF}.track.xml")):
        bad(f"{pairs[0]}: no {REF}.track.xml (ion_ls loads it from the first pair dir)")
    if dims:
        print(f"    ion grid {dims[0]}x{dims[1]}")
    if stds:
        v = sorted(stds)
        print(f"    filt_ion std: min {v[0]:.4f}, median {v[len(v) // 2]:.4f}, max {v[-1]:.4f} rad "
              f"({len(v)} pairs; the barrier excludes the tail)")
    print(f"    {ok}/{len(pairs)} pairs complete")


def check_dates_ion():
    """The barrier's product: one ion grid per date on the ml2 grid, which the
    MAIN pass subtracts. A silently empty or wrong-sized grid would corrupt
    every pair of that date."""
    ion_dir = os.environ.get("ION_DIR", "/mnt/scratch/t/tommy/full/dates_ion")
    win = json.load(open(os.path.join(PERSIST, "crop_window.json")))
    w2 = win["cropWidth"] // NRLKS2
    l2 = win["cropLen"] // (NALKS1 * NALKS2)
    if not os.path.isdir(ion_dir):
        bad(f"no dates_ion at {ion_dir}")
        return
    grids = sorted(f for f in os.listdir(ion_dir) if f.endswith(".ion"))
    if not grids:
        bad(f"{ion_dir} holds no .ion grids")
        return
    print(f"dates_ion: {len(grids)} grids, ml2 {w2}x{l2}")
    flat = 0
    for g in grids:
        f = os.path.join(ion_dir, g)
        if os.path.getsize(f) != w2 * l2 * 4:
            bad(f"{g} is {os.path.getsize(f)} B, expected {w2 * l2 * 4}")
            continue
        a = np.fromfile(f, dtype="float32")
        if not np.isfinite(a).all():
            bad(f"{g} holds non-finite values")
        elif float(a.std()) == 0.0:
            flat += 1
    if flat:
        # the zero-reference date is legitimately flat; more than one is not
        (bad if flat > 1 else print)(f"{flat} ion grids are flat (expect at most the zero date)")


def main(argv):
    if len(argv) < 2:
        raise SystemExit(__doc__)
    what = argv[1]
    if what == "seed":
        check_seed()
        tag = "SEED"
    elif what == "bundle":
        check_bundle()
        tag = "BUNDLE"
    elif what == "ion":
        if len(argv) < 3:
            raise SystemExit("check_run.py ion <block>")
        check_ion(argv[2])
        tag = f"ION_{argv[2]}"
    elif what == "dates_ion":
        check_dates_ion()
        tag = "DATES_ION"
    elif what == "date":
        if len(argv) < 3:
            raise SystemExit("check_run.py date <YYMMDD>")
        check_date(argv[2])
        tag = f"DATE_{argv[2]}"
    else:
        raise SystemExit(__doc__)
    if fails:
        print(f"{tag}_CHECK_FAILED ({len(fails)}):")
        for f in fails:
            print("  -", f)
        sys.exit(1)
    print(f"{tag}_CHECK_OK")


if __name__ == "__main__":
    main(sys.argv)
