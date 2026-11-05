"""Golden replay of crop_pack against the authentic alosStack smoke run.

Rebuilds the 181123 reference geometry the way resample_common_grid.py does
(offset-adjusted swath origins from swath_offset.txt, minima, 14-line
quantised resampled origins), then checks that

  * crop_pack.reference_polys reproduces the six analytic polynomials
    resamp_slc printed for the reference date (grid snap and stretch), and
  * crop_pack.fit_offset_polys on the three authentic 181012 cull.off files
    reproduces the applied secondary polynomials to within 1e-4 pixel over the
    whole resampled swath grid, against both the 6-digit printPoly2d output
    and the full-precision offsetpoly.f90 prints (all six fitted terms), with
    the azimuth mean shift pinned to the value Offset.getFitPolynomials would
    derive from the cull file rather than the one crop_pack stored.

Runs with the numpy replica (CROP_PACK_NUMPY_FIT=1) or with ISCE on the
cluster, where it also validates the isce branch.
"""
import datetime, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
os.environ.setdefault("CROP_PACK_NUMPY_FIT", "1")
import crop_pack as cp

REF = json.load(open(os.path.join(HERE, "golden", "smoke_181123_181012_ref.json")))
RPS = REF["rangePixelSize"]
ALI = REF["azimuthLineInterval"]
RG_M = REF["swath_offset_matching"]["range"]
AZ_M = REF["swath_offset_matching"]["azimuth"]
R0_S2 = REF["raw_s2_startingRange"]
T_MIN = cp._t(REF["printed_sensingStartMinimum"])
NALKS1 = REF["nalks1"]
# raw s2 sensingStart: the printed minimum is the s4 adjusted start
T_S2 = T_MIN + datetime.timedelta(seconds=sum(AZ_M) * ALI["2"])
RAW_DIMS = REF["reference_raw_dims"]            # 181123 raw swath [lines, samples]
DIMS_PRINTED = REF["resampled_dims_printed"]    # resamp_slc Output Image Dimensions


def geometry():
    """Reference geometry as build_geometry derives it, from the printed numbers."""
    per = {}
    for j, sn in enumerate(["2", "3", "4"]):
        adj_r0 = R0_S2 - sum(RG_M[0:j + 1]) * RPS
        adj_t0 = T_S2 - datetime.timedelta(seconds=sum(AZ_M[0:j + 1]) * ALI["2"])
        orig = {"startingRange": adj_r0, "sensingStart": adj_t0.strftime(cp.TFMT),
                "rangePixelSize": RPS, "azimuthLineInterval": ALI[sn]}
        offR = (adj_r0 - R0_S2) / RPS
        offA = (adj_t0 - T_MIN).total_seconds() / (ALI["2"] * NALKS1)
        res = {"startingRange": R0_S2 + round(offR) * RPS,
               "sensingStart": (T_MIN + datetime.timedelta(seconds=round(offA) * NALKS1 * ALI["2"])).strftime(cp.TFMT),
               "rangePixelSize": RPS, "azimuthLineInterval": ALI["2"]}
        # resampled extents as resample_common_grid.py lines 250-251 size them
        # from the reference date's raw swath; must match the resamp_slc prints
        nl, ns = RAW_DIMS[sn]
        dims = (round(nl * ALI[sn] / ALI["2"]), round(ns * RPS / RPS))
        assert list(dims) == DIMS_PRINTED[sn], (sn, dims, DIMS_PRINTED[sn])
        per[sn] = {"orig": orig, "res": res, "dims": dims}
    # the reconstruction must reproduce the printed minima, the k_s used by the
    # spec and the swath origins (samples from the frame origin 0 / 4100 / 9412)
    assert min(cp._t(d["res"]["sensingStart"]) for d in per.values()) == T_MIN
    assert min(d["res"]["startingRange"] for d in per.values()) == R0_S2
    ks = {sn: cp._lines_from_dt((cp._t(d["res"]["sensingStart"]) - T_MIN).total_seconds(), ALI["2"])
          for sn, d in per.items()}
    assert ks == {"2": 42, "3": 14, "4": 0}, ks
    so = {sn: cp._near_int((d["res"]["startingRange"] - R0_S2) / RPS) for sn, d in per.items()}
    assert so == {"2": 0, "3": 4100, "4": 9412}, so
    assert max(so[sn] + per[sn]["dims"][1] for sn in per) == 15206
    return per


def max_grid_diff(pa, pb, lines, samples):
    az = np.linspace(1, lines, 25)[:, None]
    rg = np.linspace(1, samples, 25)[None, :]
    return float(np.max(np.abs(cp.eval_poly(pa, az, rg) - cp.eval_poly(pb, az, rg))))


def main():
    per = geometry()
    fails = []
    for sn in ["2", "3", "4"]:
        lines, samples = per[sn]["dims"]
        o, r = per[sn]["orig"], per[sn]["res"]
        rp = cp.reference_polys(o, r)
        for name in ("rg", "az"):
            exp = np.array(REF["reference_polys_printed"][sn][name], dtype=float)
            got = np.array(rp[name]["coeffs"], dtype=float)
            # printed with 6 significant digits
            bad = np.abs(got - exp) > np.maximum(np.abs(exp) * 5e-6, 1e-9)
            print(f"s{sn} reference {name}: got {got.tolist()} printed {exp.tolist()}")
            if bad.any():
                fails.append(f"s{sn} reference {name} poly")
        cull = os.path.join(HERE, "golden", f"181012_s{sn}_cull.off")
        fp = cp.fit_offset_polys(cull, o, r)
        # the azimuth mean shift is not printed, but upstream derives it from
        # the cull points themselves (Offset.getFitPolynomials: azmin = min y1),
        # so pin it here instead of borrowing the value under test
        yc = cp._read_cull(cull)[:, 2]
        azmin = float(np.min(yc * o["azimuthLineInterval"] / r["azimuthLineInterval"]
                             + (cp._t(o["sensingStart"]) - cp._t(r["sensingStart"])).total_seconds()
                             / r["azimuthLineInterval"]))
        for name in ("rg", "az"):
            printed = np.array(REF["applied_polys_printed"][sn][name], dtype=float)
            got = np.array(fp[name]["coeffs"], dtype=float)
            rel = np.abs(got - printed) / np.maximum(np.abs(printed), 1e-300)
            rel[printed == 0] = np.abs(got[printed == 0])
            # upstream polynomial: all six offsetpoly.f90 terms at full precision
            fx = REF["offsetpoly_full_precision"][sn][name]
            c = printed.tolist()
            c[0][0], c[0][1], c[1][0], c[1][1] = fx["const"], fx["rg"], fx["az"], fx["cross"]
            c[0][2], c[2][0] = fx["rg2"], fx["az2"]
            got_ma = fp[name]["meanAzimuth"]
            if abs(got_ma - azmin) > 1e-9 * max(1.0, abs(azmin)):
                print(f"s{sn} fit {name}: meanAzimuth {got_ma!r} != cull azmin {azmin!r}")
                fails.append(f"s{sn} fit {name} meanAzimuth")
            upstream = cp._poly_dict(c, meanAzimuth=azmin)
            d = max_grid_diff(fp[name], upstream, lines, samples)
            allow = 1e-4
            print(f"s{sn} fit {name} ({fp['fit']}, n={fp['n']}): max |diff| over the swath grid "
                  f"vs upstream {d:.2e} px (allowed {allow:.2e}); worst coefficient rel err "
                  f"vs 6-digit print {rel.max():.1e}; meanAzimuth {got_ma:.6f}")
            if d > allow or rel.max() > 1e-4:
                fails.append(f"s{sn} fit {name} poly")
    if fails:
        print("GOLDEN_REPLAY_FAILED:", fails)
        sys.exit(1)
    print("GOLDEN_REPLAY_OK")


if __name__ == "__main__":
    main()
