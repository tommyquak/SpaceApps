"""Input and fitting checks for the provisional NISAR contrast experiment."""

from datetime import datetime

import h5py
import numpy as np


def read_pair_identity(path):
    """Read grouping metadata instead of inferring acquisition identity from a name."""
    with h5py.File(path, "r") as source:
        identification = source["science/LSAR/identification"]

        def scalar(name):
            value = identification[name][()]
            return value.decode() if isinstance(value, bytes) else value

        direction = str(scalar("orbitPassDirection")).upper()
        if direction not in ("ASCENDING", "DESCENDING"):
            raise ValueError(f"Unsupported orbit direction in {path}: {direction}")
        dates = []
        for name in ("referenceZeroDopplerStartTime", "secondaryZeroDopplerStartTime"):
            date = str(scalar(name))[:10]
            datetime.strptime(date, "%Y-%m-%d")
            dates.append(date.replace("-", ""))
        if dates[0] >= dates[1]:
            raise ValueError(f"Non-forward interferogram in {path}: {dates}")
        crid = str(scalar("compositeReleaseId"))
        if not crid.startswith("P"):
            raise ValueError(f"Expected provisional CRID in {path}, got {crid}")
        return dict(path=path, track=int(scalar("trackNumber")),
                    frame=int(scalar("frameNumber")), ad=direction[0],
                    crid=crid, d1=dates[0], d2=dates[1])


def select_chain(records, direction, track=None, frame=None, crid=None,
                 start=None, end=None):
    """Require one explicit product group and an uninterrupted, nonbranching chain."""
    candidates = [r for r in records if r["ad"] == direction
                  and (track is None or r["track"] == int(track))
                  and (frame is None or r["frame"] == int(frame))
                  and (crid is None or r["crid"] == crid)
                  and (start is None or r["d1"] >= start.replace("-", ""))
                  and (end is None or r["d2"] <= end.replace("-", ""))]
    groups = {(r["track"], r["frame"], r["crid"]) for r in candidates}
    if len(groups) != 1:
        raise ValueError(f"[{direction}] expected one track/frame/CRID group, got {sorted(groups)}; "
                         f"set NISAR_{direction}_TRACK, NISAR_{direction}_FRAME and NISAR_CRID")
    chain = sorted(candidates, key=lambda r: (r["d1"], r["d2"]))
    if len(chain) < 2:
        raise ValueError(f"[{direction}] fewer than two pairs in selected group")
    pairs = [(r["d1"], r["d2"]) for r in chain]
    if len(set(pairs)) != len(pairs):
        raise ValueError(f"[{direction}] duplicate/reprocessed date pairs; select one product per pair")
    for left, right in zip(chain, chain[1:]):
        if left["d2"] != right["d1"]:
            raise ValueError(f"[{direction}] gap or branching network at {left['d2']} -> {right['d1']}; "
                             f"select a contiguous interval with NISAR_{direction}_START/END")
    return chain


def validate_grids(grids):
    """Reject incompatible coordinate support before adding phase arrays."""
    first = grids[0]
    for index, grid in enumerate(grids):
        shape = grid["coh"].shape
        if shape != (len(grid["y"]), len(grid["x"])):
            raise ValueError(f"Pair {index}: raster shape does not match coordinate vectors")
        for name in ("unw", "cc", "ion", "ionu"):
            if grid[name] is None or grid[name].shape != shape:
                raise ValueError(f"Pair {index}: missing or mismatched {name} layer")
        for name in ("x", "y"):
            coordinates = grid[name]
            if len(coordinates) < 2 or not np.all(np.isfinite(coordinates)):
                raise ValueError(f"Pair {index}: invalid {name} coordinates")
            delta = np.diff(coordinates)
            if delta[0] == 0 or not np.allclose(delta, delta[0], rtol=0, atol=1e-6):
                raise ValueError(f"Pair {index}: irregular {name} coordinates")
            if coordinates.shape != first[name].shape or not np.array_equal(coordinates, first[name]):
                raise ValueError(f"Pair {index}: {name} grid differs from reference pair")
        if grid["epsg"] != first["epsg"]:
            raise ValueError(f"Pair {index}: CRS differs from reference pair")
        if not np.isfinite(grid["wl"]) or grid["wl"] <= 0 or not np.isclose(
                grid["wl"], first["wl"], rtol=1e-8, atol=0):
            raise ValueError(f"Pair {index}: invalid or inconsistent metadata wavelength")


def deramp(a, mask, minimum_stable_pixels=500):
    """Fit only the supplied stable support; never fall back to the target field."""
    if a.shape != mask.shape:
        raise ValueError("Stable-ground mask and phase shapes differ")
    ny, nx = a.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    yy = (yy - ny / 2) / ny
    xx = (xx - nx / 2) / nx
    valid = mask & np.isfinite(a)
    if valid.sum() < minimum_stable_pixels:
        raise ValueError(f"Insufficient stable-ground pixels: {valid.sum()} < {minimum_stable_pixels}")
    design = np.stack([np.ones(valid.sum()), xx[valid], yy[valid], xx[valid] ** 2,
                       yy[valid] ** 2, xx[valid] * yy[valid]], axis=1)
    coef, _, rank, _ = np.linalg.lstsq(design, a[valid], rcond=None)
    if rank != 6:
        raise ValueError(f"Stable-ground quadratic fit is rank deficient: {rank}/6")
    surface = (coef[0] + coef[1] * xx + coef[2] * yy + coef[3] * xx ** 2
               + coef[4] * yy ** 2 + coef[5] * xx * yy)
    return a - surface
