"""Unpack an ALOS-1 PALSAR L1.1 (CEOS SLC) granule for ISCE2 stripmapStack.

Usage: p1_unpack_slc.py -i <granule dir with IMG-HH*/LED*> -o <SLC/<date>>

ISCE's own unpackFrame_ALOS.py stores only the frame in the shelve; every
downstream stripmapStack step (baselineGrid, geo2rdr, resampleSlc, topo,
geocode) also reads db['doppler'] and db['fmrate'], which is why the L1.1
route "broke geo2rdr" in May 2026. This mirrors unpackFrame_ALOS2.py: the
Doppler centroid and azimuth FM-rate polynomials come from the CEOS leader
(cross-track terms vs range pixel) that ALOS_SLC already parses.
"""
import argparse
import glob
import os
import shelve

import isce  # noqa: F401
from isceobj.Sensor import createSensor
from isceobj.Util import Poly1D


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", dest="indir", required=True)
    ap.add_argument("-o", dest="slcdir", required=True)
    a = ap.parse_args()
    indir, slcdir = a.indir.rstrip("/"), a.slcdir.rstrip("/")
    imgs = sorted(glob.glob(os.path.join(indir, "IMG-HH*")) + glob.glob(os.path.join(indir, "*", "IMG-HH*")))
    leds = sorted(glob.glob(os.path.join(indir, "LED*")) + glob.glob(os.path.join(indir, "*", "LED*")))
    if not imgs or not leds:
        raise SystemExit(f"no IMG-HH/LED under {indir}")
    os.makedirs(slcdir, exist_ok=True)
    date = os.path.basename(slcdir)
    obj = createSensor("ALOS_SLC")
    obj.configure()
    obj._leaderFile = leds[0]
    obj._imageFile = imgs[0]
    obj.output = os.path.join(slcdir, date + ".slc")
    obj.extractImage()
    obj.frame.getImage().renderHdr()
    obj.extractDoppler()
    poly = Poly1D.Poly1D(); poly.initPoly(order=len(obj.doppler_coeff) - 1); poly.setCoeffs(obj.doppler_coeff)
    fpoly = Poly1D.Poly1D(); fpoly.initPoly(order=len(obj.azfmrate_coeff) - 1); fpoly.setCoeffs(obj.azfmrate_coeff)
    with shelve.open(os.path.join(slcdir, "data")) as db:
        db["frame"] = obj.frame
        db["doppler"] = poly
        db["fmrate"] = fpoly
        # ALOS-1 leaders carry no facility record (ALOS-2's unpacker stores one); nothing downstream needs it
    print(f"unpacked {date}: {obj.frame.getNumberOfSamples()} x {obj.frame.getNumberOfLines()}, "
          f"doppler coeffs {obj.doppler_coeff}, fmrate {obj.azfmrate_coeff}")


if __name__ == "__main__":
    main()
