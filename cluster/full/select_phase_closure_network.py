#!/usr/bin/env python3
"""Find candidate MintPy phase-closure subnetworks for an isolated experiment.

[E] MintPy's implementation of Yunjun, Fattahi, and Amelung (2019),
doi:10.1016/j.cageo.2019.104331, solves integer ambiguities only inside a
connected region shared by every retained interferogram before transferring
the corrections to intersecting local components. A large sequential network
can therefore fail even when smaller sequential subnetworks have an anchor.

This program performs a read-only, streaming preflight. It tests sequential
bandwidths from small to large, reports the common-region area and number of
closure triangles, and can write the pair list for the largest passing bandwidth.
Passing this diagnostic does not establish peat coverage, correct unwrapping,
or an accurate deformation estimate. A common region may be an anchor outside
peat; component support at the target sites requires a separate assessment.
The scan respects ``dropIfgram`` unless explicitly told to include dropped pairs.
It does not alter any dataset in the input HDF5 file. Do not run it on a stack
that another process is writing.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import itertools
import json
import math
import os
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
from scipy import ndimage


def _fmt(value, spec=".6f"):
    """Format a statistic that may legitimately be unavailable.

    valid_count == 0 yields None for every closure statistic, and formatting
    None with a float spec raises
    TypeError: unsupported format string passed to NoneType.__format__.
    """
    if value is None:
        return "n/a (closure sampling incomplete or unavailable)"
    return format(value, spec)


@dataclass(frozen=True)
class Pair:
    stack_index: int
    date1: str
    date2: str
    gap: int

    @property
    def date12(self) -> str:
        return f"{self.date1}_{self.date2}"


@dataclass(frozen=True)
class BandwidthResult:
    """Structural candidate metrics; ``viable`` does not certify the phase."""

    bandwidth: int
    interferograms: int
    triangles: int
    common_pixels: int
    largest_region_pixels: int
    common_regions: int
    complete_sequential_levels: bool
    connected_date_graph: bool
    aoi_region_pixels: int
    viable: bool


def _decode(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def read_pairs(stack: h5py.File) -> tuple[list[str], list[Pair]]:
    """Return sorted acquisition dates and validated pair metadata."""
    if "date" not in stack:
        raise ValueError("ifgram stack has no date dataset")
    if stack["date"].ndim != 2 or stack["date"].shape[1] != 2:
        raise ValueError("date dataset must have shape (interferograms, 2)")

    raw_dates = [(_decode(row[0]), _decode(row[1])) for row in stack["date"][:]]
    dates = sorted({date for pair in raw_dates for date in pair})
    date_index = {date: index for index, date in enumerate(dates)}

    pairs: list[Pair] = []
    seen: set[tuple[str, str]] = set()
    for stack_index, (date1, date2) in enumerate(raw_dates):
        key = (date1, date2)
        if key in seen:
            raise ValueError(f"duplicate interferogram: {date1}_{date2}")
        seen.add(key)
        gap = date_index[date2] - date_index[date1]
        if gap <= 0:
            raise ValueError(f"dates are not chronological: {date1}_{date2}")
        pairs.append(Pair(stack_index, date1, date2, gap))
    return dates, pairs


def _retained_flags(stack: h5py.File, pair_count: int) -> np.ndarray:
    if "dropIfgram" not in stack:
        return np.ones(pair_count, dtype=bool)
    flags = stack["dropIfgram"][:]
    if flags.shape != (pair_count,) or not np.isin(flags, [0, 1]).all():
        raise ValueError("dropIfgram must contain one boolean flag per pair")
    # MintPy uses True for retained interferograms despite the dataset name.
    return flags.astype(bool)


def _reference_yx(stack: h5py.File, shape: tuple[int, int]) -> tuple[int, int]:
    try:
        values = [float(stack.attrs[key]) for key in ("REF_Y", "REF_X")]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("ifgram stack needs valid REF_Y/REF_X attributes") from exc
    if any(not math.isfinite(value) or not value.is_integer() for value in values):
        raise ValueError("REF_Y/REF_X must be finite integer pixel indices")
    ref_y, ref_x = map(int, values)
    if not (0 <= ref_y < shape[0] and 0 <= ref_x < shape[1]):
        raise ValueError("REF_Y/REF_X is outside the stack grid")
    return ref_y, ref_x


def _aoi_mask_from_bbox(shape, attrs, bbox):
    """Optional geographic filter for anchor regions, not a peat-validity test."""
    if bbox is None:
        return None
    west, south, east, north = bbox
    if not all(math.isfinite(value) for value in bbox):
        raise ValueError("--aoi-bbox coordinates must be finite")
    if not (-180 <= west < east <= 180 and -90 <= south < north <= 90):
        raise ValueError("--aoi-bbox requires ordered W S E N longitude/latitude bounds")
    if "UTM_ZONE" in attrs or any(
        _decode(attrs.get(key, "degrees")).lower() not in ("degree", "degrees")
        for key in ("X_UNIT", "Y_UNIT")
    ):
        raise ValueError("--aoi-bbox requires a longitude/latitude grid")
    try:
        y0 = float(attrs["Y_FIRST"])
        ys = float(attrs["Y_STEP"])
        x0 = float(attrs["X_FIRST"])
        xs = float(attrs["X_STEP"])
    except KeyError:
        raise ValueError(
            "--aoi-bbox needs X_FIRST/Y_FIRST/X_STEP/Y_STEP on the stack; "
            "this stack is not geocoded"
        )
    if not all(math.isfinite(value) for value in (y0, ys, x0, xs)) or not ys or not xs:
        raise ValueError("stack grid origin and steps must be finite, with nonzero steps")
    # build_stack2 records GDAL's upper-left pixel corner in X_FIRST/Y_FIRST.
    rows = y0 + ys * (np.arange(shape[0]) + 0.5)
    cols = x0 + xs * (np.arange(shape[1]) + 0.5)
    in_row = (rows <= max(north, south)) & (rows >= min(north, south))
    in_col = (cols >= min(west, east)) & (cols <= max(west, east))
    mask = np.outer(in_row, in_col)
    if not mask.any():
        raise ValueError("--aoi-bbox does not intersect the stack grid")
    return mask


def _largest_regions(
    mask: np.ndarray,
    min_region_pixels: int,
    erosion_size: int = 5,
    aoi_mask: np.ndarray | None = None,
    min_aoi_pixels: int = 0,
) -> tuple[np.ndarray, int, int]:
    """Return the largest retained region mask, its area, and region count.

    Four-neighbour connectivity and a 5 by 5 erosion survival check mirror the
    important parts of MintPy ``label_conn_comp``. Unlike MintPy, the explicit
    minimum area is also enforced when only one region exists, avoiding a small
    single-region edge case that later fails when 100 samples are requested.
    """
    structure = ndimage.generate_binary_structure(2, 1)
    labels, count = ndimage.label(mask, structure=structure)
    if count == 0:
        return np.zeros_like(mask, dtype=bool), 0, 0

    areas = np.bincount(labels.ravel())
    kept: list[tuple[int, int]] = []
    erosion_structure = np.ones((erosion_size, erosion_size), dtype=bool)
    for label in range(1, count + 1):
        area = int(areas[label])
        if area < min_region_pixels:
            continue
        region = labels == label
        if not ndimage.binary_erosion(region, structure=erosion_structure).any():
            continue
        kept.append((label, area))

    if not kept:
        return np.zeros_like(mask, dtype=bool), 0, 0
    kept.sort(key=lambda item: item[1], reverse=True)
    if aoi_mask is not None:
        # A smaller retained region may pass the requested anchor filter even
        # when the largest region is outside it.
        for label, area in kept:
            region = labels == label
            if int((region & aoi_mask).sum()) >= min_aoi_pixels:
                return region, area, len(kept)
        return np.zeros_like(mask, dtype=bool), 0, len(kept)
    largest_label, largest_area = kept[0]
    return labels == largest_label, largest_area, len(kept)


def _triangles(pairs: Iterable[Pair]) -> list[tuple[int, int, int]]:
    """Return local pair-column indices for all available date triangles."""
    pair_list = list(pairs)
    by_dates = {
        (pair.date1, pair.date2): local_index
        for local_index, pair in enumerate(pair_list)
    }
    dates = sorted({date for pair in pair_list for date in (pair.date1, pair.date2)})
    triangles: list[tuple[int, int, int]] = []
    for date1, date2, date3 in itertools.combinations(dates, 3):
        try:
            pair12 = by_dates[(date1, date2)]
            pair23 = by_dates[(date2, date3)]
            pair13 = by_dates[(date1, date3)]
        except KeyError:
            continue
        triangles.append((pair12, pair23, pair13))
    return triangles


def _date_graph_connected(dates: list[str], pairs: Iterable[Pair]) -> bool:
    if not dates:
        return False
    neighbours = {date: set() for date in dates}
    for pair in pairs:
        neighbours[pair.date1].add(pair.date2)
        neighbours[pair.date2].add(pair.date1)
    visited = {dates[0]}
    frontier = [dates[0]]
    while frontier:
        current = frontier.pop()
        for neighbour in neighbours[current] - visited:
            visited.add(neighbour)
            frontier.append(neighbour)
    return len(visited) == len(dates)


def scan_sequential_networks(
    stack_path: str | Path,
    min_bandwidth: int = 2,
    max_bandwidth: int = 20,
    min_region_pixels: int = 2500,
    aoi_bbox: tuple[float, float, float, float] | None = None,
    min_aoi_pixels: int = 500,
    include_dropped: bool = False,
) -> tuple[list[BandwidthResult], BandwidthResult | None, np.ndarray | None, list[Pair]]:
    """Scan nested sequential networks without loading the 3-D stack at once."""
    if min_bandwidth < 2:
        raise ValueError("min_bandwidth must be at least 2 to form triangles")
    if max_bandwidth < min_bandwidth:
        raise ValueError("max_bandwidth must be at least min_bandwidth")
    if min_region_pixels < 1:
        raise ValueError("min_region_pixels must be positive")
    if min_aoi_pixels < 1:
        raise ValueError("min_aoi_pixels must be positive")
    aoi_mask = None

    results: list[BandwidthResult] = []
    best: BandwidthResult | None = None
    best_mask: np.ndarray | None = None
    best_pairs: list[Pair] = []

    with h5py.File(stack_path, "r") as stack:
        if "connectComponent" not in stack:
            raise ValueError("ifgram stack has no connectComponent dataset")
        dates, pairs = read_pairs(stack)
        if len(dates) < 3:
            raise ValueError("at least three acquisition dates are required")
        conncomp = stack["connectComponent"]
        if conncomp.ndim != 3 or conncomp.shape[0] != len(pairs):
            raise ValueError("connectComponent and date dataset shapes disagree")
        flags = _retained_flags(stack, len(pairs))
        if not include_dropped:
            pairs = [pair for pair in pairs if flags[pair.stack_index]]

        common = np.ones(conncomp.shape[1:], dtype=bool)
        if aoi_bbox is not None:
            aoi_mask = _aoi_mask_from_bbox(
                conncomp.shape[1:], dict(stack.attrs), aoi_bbox
            )
        levels_complete = True
        pair_keys = {(pair.date1, pair.date2) for pair in pairs}

        for bandwidth in range(1, min(max_bandwidth, len(dates) - 1) + 1):
            expected_level = {
                (dates[index], dates[index + bandwidth])
                for index in range(len(dates) - bandwidth)
            }
            levels_complete = levels_complete and expected_level.issubset(pair_keys)

            for pair in pairs:
                if pair.gap != bandwidth:
                    continue
                values = conncomp[pair.stack_index, :, :]
                common &= np.isfinite(values) & (values > 0)

            kept_pairs = [pair for pair in pairs if pair.gap <= bandwidth]
            triangle_count = len(_triangles(kept_pairs))
            connected = _date_graph_connected(dates, kept_pairs)
            largest_mask, largest_area, region_count = _largest_regions(
                common,
                min_region_pixels=min_region_pixels,
                aoi_mask=aoi_mask,
                min_aoi_pixels=min_aoi_pixels,
            )
            if aoi_mask is None:
                aoi_area = -1
                aoi_ok = True
            else:
                aoi_area = int((largest_mask & aoi_mask).sum())
                aoi_ok = aoi_area >= min_aoi_pixels
            viable = (
                bandwidth >= min_bandwidth
                and levels_complete
                and connected
                and triangle_count > 0
                and largest_area >= min_region_pixels
                and aoi_ok
            )
            result = BandwidthResult(
                bandwidth=bandwidth,
                interferograms=len(kept_pairs),
                triangles=triangle_count,
                common_pixels=int(common.sum()),
                largest_region_pixels=largest_area,
                common_regions=region_count,
                complete_sequential_levels=levels_complete,
                connected_date_graph=connected,
                aoi_region_pixels=aoi_area,
                viable=viable,
            )
            results.append(result)
            if viable:
                best = result
                best_mask = largest_mask.copy()
                best_pairs = kept_pairs

            # Intersections are nested. Once empty, a wider network cannot
            # regain a common pixel, so avoid unnecessary HDF5 reads.
            if not common.any():
                break

    return results, best, best_mask, best_pairs


def sample_closure_stats(
    stack_path: str | Path,
    pairs: list[Pair],
    region_mask: np.ndarray,
    sample_count: int = 100,
    sample_seed: int = 0,
) -> dict[str, object]:
    """Measure closure only with a valid reference and complete finite samples.

    Missing data invalidate the diagnostic instead of silently removing the
    affected triangles from the denominator. Zero phase is a legitimate value
    when connectComponent is positive, including at an already referenced pixel.
    """
    if sample_count < 1:
        raise ValueError("sample_count must be positive")
    coords = np.argwhere(region_mask)
    if coords.size == 0:
        raise ValueError("cannot sample an empty common region")
    chosen = np.random.default_rng(sample_seed).choice(
        len(coords), min(sample_count, len(coords)), replace=False
    )
    chosen.sort()
    sample_coords = coords[chosen]
    triangles = _triangles(pairs)
    if not triangles:
        raise ValueError("selected network has no closure triangles")

    phase = np.empty((len(pairs), len(sample_coords)), dtype=np.float64)
    with h5py.File(stack_path, "r") as stack:
        if "unwrapPhase" not in stack or "connectComponent" not in stack:
            raise ValueError("closure sampling needs unwrapPhase and connectComponent")
        unwrap = stack["unwrapPhase"]
        conncomp = stack["connectComponent"]
        if unwrap.ndim != 3 or unwrap.shape != conncomp.shape:
            raise ValueError("unwrapPhase and connectComponent dataset shapes disagree")
        if region_mask.shape != unwrap.shape[1:]:
            raise ValueError("region mask does not match the stack grid")
        _, source_pairs = read_pairs(stack)
        if unwrap.shape[0] != len(source_pairs):
            raise ValueError("unwrapPhase and date dataset shapes disagree")
        if len({pair.stack_index for pair in pairs}) != len(pairs) or any(
            pair.stack_index < 0
            or pair.stack_index >= len(source_pairs)
            or source_pairs[pair.stack_index] != pair
            for pair in pairs
        ):
            raise ValueError("selected pair metadata does not match the source stack")
        ref_y, ref_x = _reference_yx(stack, unwrap.shape[1:])
        bad_ref = []
        bad_samples = []
        for local_index, pair in enumerate(pairs):
            ref_phase = float(unwrap[pair.stack_index, ref_y, ref_x])
            ref_component = float(conncomp[pair.stack_index, ref_y, ref_x])
            if (
                not math.isfinite(ref_phase)
                or not math.isfinite(ref_component)
                or ref_component <= 0
            ):
                bad_ref.append(pair.date12)
                phase[local_index, :] = np.nan
                continue
            invalid_count = 0
            for sample_index, (y, x) in enumerate(sample_coords):
                value = float(unwrap[pair.stack_index, int(y), int(x)])
                component = float(conncomp[pair.stack_index, int(y), int(x)])
                valid = (
                    math.isfinite(value)
                    and math.isfinite(component)
                    and component > 0
                    and math.isfinite(value - ref_phase)
                )
                phase[local_index, sample_index] = value - ref_phase if valid else np.nan
                invalid_count += not valid
            if invalid_count:
                bad_samples.append({"pair": pair.date12, "invalid_samples": invalid_count})

    pair12 = np.fromiter((triangle[0] for triangle in triangles), dtype=np.int64)
    pair23 = np.fromiter((triangle[1] for triangle in triangles), dtype=np.int64)
    pair13 = np.fromiter((triangle[2] for triangle in triangles), dtype=np.int64)
    with np.errstate(invalid="ignore", over="ignore"):
        closure = phase[pair12] + phase[pair23] - phase[pair13]
        wrapped = np.angle(np.exp(1j * closure))
        ambiguity = np.rint((closure - wrapped) / (2.0 * math.pi))
    finite = np.isfinite(ambiguity)
    bad = finite & (ambiguity != 0)
    valid_count = int(finite.sum())
    bad_per_sample = np.divide(
        bad.sum(axis=0),
        finite.sum(axis=0),
        out=np.full(len(sample_coords), np.nan),
        where=finite.sum(axis=0) > 0,
    )
    total_cells = int(ambiguity.size)
    complete = valid_count == total_cells and not bad_ref and not bad_samples
    return {
        "diagnostic_valid": complete,
        "samples": int(len(sample_coords)),
        "sample_seed": sample_seed,
        "sampling_method": "seeded random sampling without replacement",
        "reference_yx": [ref_y, ref_x],
        "triangles": int(len(triangles)),
        "valid_triangle_samples": valid_count,
        # Without these two the fraction silently describes only the survivors.
        "total_triangle_samples": total_cells,
        "finite_coverage_fraction": (
            float(valid_count / total_cells) if total_cells else 0.0
        ),
        "interferograms_with_invalid_reference": sorted(bad_ref),
        "interferograms_with_invalid_samples": bad_samples,
        "nonzero_integer_ambiguity_fraction": (
            float(bad.sum() / valid_count) if complete else None
        ),
        "median_bad_triplet_fraction_per_sample": (
            float(np.median(bad_per_sample)) if complete else None
        ),
        "max_abs_integer_cycles": (
            float(np.max(np.abs(ambiguity))) if complete else None
        ),
        "sample_yx": sample_coords.tolist(),
        "interpretation": (
            "Diagnostic only. A bad-triplet fraction does not prove MintPy's "
            "majority-rightness assumption, target-site support, or deformation "
            "accuracy. Incomplete finite coverage invalidates this diagnostic. "
            "Validate closure again after correction."
        ),
    }


def _print_results(results: list[BandwidthResult]) -> None:
    print(
        "bandwidth  ifgrams  triangles  common_px  largest_region  "
        "aoi_px  regions  complete  connected  viable"
    )
    for result in results:
        print(
            f"{result.bandwidth:9d}  {result.interferograms:7d}  "
            f"{result.triangles:9d}  {result.common_pixels:9d}  "
            f"{result.largest_region_pixels:14d}  "
            f"{('n/a' if result.aoi_region_pixels < 0 else str(result.aoi_region_pixels)):>6s}  "
            f"{result.common_regions:7d}  "
            f"{str(result.complete_sequential_levels):8s}  "
            f"{str(result.connected_date_graph):9s}  {result.viable}"
        )


def source_provenance(stack_path: str | Path) -> dict[str, object]:
    """Record the network and source identity without hashing the phase cube.

    File stat information also detects ordinary writes during the preflight.
    This is a metadata snapshot, not a cryptographic checksum of phase contents
    and not permission to read concurrently with a writer.
    """
    path = Path(stack_path).resolve()
    before = path.stat()
    with h5py.File(path, "r") as stack:
        dates, pairs = read_pairs(stack)
        flags = _retained_flags(stack, len(pairs))
        metadata = {
            "date_pairs": [pair.date12 for pair in pairs],
            "retained_flags": flags.tolist(),
            "attributes": {key: _decode(value) for key, value in stack.attrs.items()},
            "datasets": {
                name: {
                    "shape": list(stack[name].shape),
                    "dtype": str(stack[name].dtype),
                    "modification_time": _decode(
                        stack[name].attrs.get("MODIFICATION_TIME", "")
                    ),
                }
                for name in ("unwrapPhase", "connectComponent")
                if name in stack
            },
        }
    after = path.stat()
    if (before.st_ino, before.st_size, before.st_mtime_ns) != (
        after.st_ino, after.st_size, after.st_mtime_ns
    ):
        raise ValueError("source stack changed while its metadata was being read")
    snapshot = {
        "path": str(path),
        "file_size": after.st_size,
        "file_mtime_ns": after.st_mtime_ns,
        "file_inode": after.st_ino,
        "metadata": metadata,
    }
    return {
        "snapshot_sha256": hashlib.sha256(
            json.dumps(snapshot, sort_keys=True).encode("utf-8")
        ).hexdigest(),
        "path": str(path),
        "file_size": after.st_size,
        "file_mtime_ns": after.st_mtime_ns,
        "acquisitions": len(dates),
        "interferograms": len(pairs),
        "retained_interferograms": int(flags.sum()),
        "attributes": metadata["attributes"],
        "datasets": metadata["datasets"],
        "fingerprint_scope": "network, attributes, dataset metadata and file stat; phase contents are not hashed",
    }


def _check_output_paths(stack_path: str, outputs: list[str | None]) -> None:
    seen = {Path(stack_path).resolve()}
    for output in outputs:
        if output is None:
            continue
        path = Path(output)
        resolved = path.resolve()
        if resolved in seen:
            raise ValueError(f"output path collides with the input or another output: {path}")
        seen.add(resolved)
        if path.exists() or path.is_symlink():
            raise ValueError(f"output already exists; choose a new path to avoid stale artifacts: {path}")
        if not path.parent.is_dir():
            raise ValueError(f"output parent directory does not exist: {path.parent}")


def _write_new_atomic(path: str, data: bytes) -> None:
    """Publish a complete file and refuse to replace an existing destination."""
    output = Path(path)
    with tempfile.NamedTemporaryFile(dir=output.parent, prefix=f".{output.name}.", delete=False) as temp:
        temporary_path = Path(temp.name)
        try:
            temp.write(data)
            temp.flush()
            os.fsync(temp.fileno())
            os.link(temporary_path, output)
        finally:
            temporary_path.unlink()


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stack", help="MintPy inputs/ifgramStack.h5")
    parser.add_argument("--min-bandwidth", type=int, default=2)
    parser.add_argument("--max-bandwidth", type=int, default=20)
    parser.add_argument("--min-region-pixels", type=int, default=2500)
    parser.add_argument(
        "--include-dropped", action="store_true",
        help="explicitly explore pairs currently marked dropped; default respects dropIfgram",
    )
    parser.add_argument(
        "--aoi-bbox", type=float, nargs=4, metavar=("W", "S", "E", "N"),
        help="optional anchor-region filter in lon/lat; requires --min-aoi-pixels "
             "inside the box. Mineral anchors may be useful. This does not test "
             "connected-component support at peat or dam sites.",
    )
    parser.add_argument(
        "--min-aoi-pixels", type=int, default=500,
        help="minimum common-region pixels inside --aoi-bbox (default 500)",
    )
    parser.add_argument(
        "--write-region-mask",
        help="write the selected common-region mask as .npy for inspection",
    )
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--sample-seed", type=int, default=0)
    parser.add_argument(
        "--write-pairs",
        help="write experimental YYYYMMDD_YYYYMMDD pairs; requires --json and a new output path",
    )
    parser.add_argument("--json", dest="json_path", help="write a machine-readable report")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = create_parser()
    args = parser.parse_args(argv)
    if args.write_pairs and not args.json_path:
        parser.error("--write-pairs requires --json to record the source and selected-pair fingerprint")
    try:
        return _run(args)
    except (OSError, ValueError, KeyError) as exc:
        print(f"PHASE_CLOSURE_PREFLIGHT_FAILED: {exc}")
        return 3


def _run(args: argparse.Namespace) -> int:
    _check_output_paths(args.stack, [args.write_pairs, args.json_path, args.write_region_mask])
    source = source_provenance(args.stack)
    results, best, best_mask, best_pairs = scan_sequential_networks(
        args.stack,
        min_bandwidth=args.min_bandwidth,
        max_bandwidth=args.max_bandwidth,
        min_region_pixels=args.min_region_pixels,
        aoi_bbox=tuple(args.aoi_bbox) if args.aoi_bbox else None,
        min_aoi_pixels=args.min_aoi_pixels,
        include_dropped=args.include_dropped,
    )
    _print_results(results)

    report: dict[str, object] = {
        "schema_version": 2,
        "stack": str(Path(args.stack).resolve()),
        "source": source,
        "selector_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "parameters": vars(args),
        "min_region_pixels": args.min_region_pixels,
        "bandwidths": [asdict(result) for result in results],
        "selected": asdict(best) if best else None,
        "closure_samples": None,
        "status": "no_complete_sequential_candidate",
        "selection_scope": (
            "Largest passing complete sequential bandwidth, not an optimized network. "
            "Only a diagnostic candidate for an isolated phase-closure experiment. "
            "No conclusion about arbitrary subnetworks, peat support, or deformation accuracy."
        ),
    }
    exit_code = 2
    if best is None or best_mask is None:
        print("NO_VIABLE_COMPLETE_SEQUENTIAL_PHASE_CLOSURE_NETWORK")
        print("Only the tested complete sequential networks were ruled out.")
        if args.write_region_mask:
            print("NO_REGION_MASK_WRITTEN (no viable bandwidth)")
    else:
        report["status"] = "candidate_diagnostic_invalid"
        report["selected_pair_list_sha256"] = hashlib.sha256(
            ("\n".join(pair.date12 for pair in best_pairs) + "\n").encode("utf-8")
        ).hexdigest()
        exit_code = 3
        try:
            closure_stats = sample_closure_stats(
                args.stack, best_pairs, best_mask,
                sample_count=args.samples, sample_seed=args.sample_seed,
            )
            report["closure_samples"] = closure_stats
            print(f"CANDIDATE_BANDWIDTH={best.bandwidth}")
            print("SAMPLED_BAD_TRIPLET_FRACTION=" + _fmt(closure_stats["nonzero_integer_ambiguity_fraction"]))
            print("SAMPLED_FINITE_COVERAGE=" + _fmt(closure_stats["finite_coverage_fraction"], ".4f"))
            if closure_stats["diagnostic_valid"]:
                report["status"] = "candidate_diagnostic_complete"
                exit_code = 0
            else:
                print("CANDIDATE_DIAGNOSTIC_INVALID: invalid reference or incomplete phase samples; no pair list will be written")
        except ValueError as exc:
            report["diagnostic_error"] = str(exc)
            print(f"CANDIDATE_DIAGNOSTIC_INVALID: {exc}")
        print(
            "NOTE: this is an experimental network candidate. Complete sampling "
            "does not establish correct unwrapping or useful peat coverage."
        )

    if source_provenance(args.stack) != source:
        raise ValueError("source stack changed during preflight; no artifacts will be published")
    if args.json_path:
        _write_new_atomic(args.json_path, (json.dumps(report, indent=2, allow_nan=False) + "\n").encode("utf-8"))
        print(f"WROTE_REPORT={args.json_path}")
    if args.write_region_mask and best_mask is not None:
        buffer = io.BytesIO()
        np.save(buffer, best_mask, allow_pickle=False)
        _write_new_atomic(args.write_region_mask, buffer.getvalue())
        print(f"WROTE_REGION_MASK={args.write_region_mask} ({int(best_mask.sum())} px; diagnostic anchor mask)")
    if args.write_pairs:
        if exit_code == 0:
            _write_new_atomic(args.write_pairs, ("\n".join(pair.date12 for pair in best_pairs) + "\n").encode("utf-8"))
            print(f"WROTE_PAIR_LIST={args.write_pairs}")
        else:
            print("NO_PAIR_LIST_WRITTEN (no candidate with a valid closure diagnostic)")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
