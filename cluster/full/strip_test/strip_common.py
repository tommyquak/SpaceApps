"""Shared pieces of the Tay-overlap strip experiment (cluster-side only).

Question: on identical inputs (the surviving 17x56 filtered interferograms cropped to
the strip where Tay's raster exists), does inverting in radar coordinates without a
deramp (Tay's order) correlate with Tay where our geocode-first, deramped chain did not?
Pre-registered: B reaches r >= +0.30 with a block-bootstrap CI clear of zero at
temporal coherence 0.5 -> the processing order and deramp were material defects;
both A' and B below +0.20 -> the limit is the data, not the order.
"""
import glob
import json
import os
import re

import numpy as np

HOME = os.path.expanduser("~")
PAIRS = f"{HOME}/full/pairs_tay"                       # 1650 pairs, filt_<p>_17rlks_56alks.int (770 x 960 CFLOAT)
ION = "/mnt/scratch/t/tommy/full/pairs_ion/190621-190830/ion/ion_cal"   # lat/lon at 80x448 (163 x 120)
REFSTACK = f"{HOME}/full/mintpy_tay/inputs/ifgramStack.h5"              # dates, bperp, radar attributes
GEOGEO = f"{HOME}/full/mintpy_tay/inputs/geometryGeo.h5"                # incidence/azimuth/height/slant on the 90 m grid
ROOT = "/mnt/scratch/t/tommy/strip_test"
OUT = f"{HOME}/full/strip_test"
W_FULL, L_FULL = 770, 960
STRIP = dict(lon_min=113.80, lon_max=114.50, lat_min=-2.90, lat_max=-2.20)   # Tay raster coverage
ML = "_17rlks_56alks"


def pair_list():
    cals = sorted(os.path.basename(d.rstrip("/")) for d in glob.glob(f"{PAIRS}/*/") if os.path.exists(f"{d}/DONE"))
    return cals


def cal_to_p(cal):
    r, s = cal[:10], cal[11:]
    return f"{r[2:4]}{r[5:7]}{r[8:10]}-{s[2:4]}{s[5:7]}{s[8:10]}"


def read_int(cal):
    p = cal_to_p(cal)
    fn = f"{PAIRS}/{cal}/filt_{p}{ML}.int"
    xml = open(fn + ".xml").read()
    w = int(re.search(r'name="width".*?<value>(\d+)</value>', xml, re.S).group(1))
    l = int(re.search(r'name="length".*?<value>(\d+)</value>', xml, re.S).group(1))
    return np.fromfile(fn, dtype=np.complex64).reshape(l, w)


def window():
    return json.load(open(f"{OUT}/window.json"))


def phase_consistency(z, half=2):
    """|mean z| / mean |z| in a (2*half+1)^2 window of the FILTERED interferogram: a
    biased-high coherence proxy used only for snaphu costs and inverse-variance weights."""
    from scipy.ndimage import uniform_filter
    k = 2 * half + 1
    num = np.abs(uniform_filter(z.real, k) + 1j * uniform_filter(z.imag, k))
    den = uniform_filter(np.abs(z), k)
    return np.where(den > 0, num / np.maximum(den, 1e-12), 0).astype(np.float32)
