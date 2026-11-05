"""Conservative cache and product provenance for the full MintPy workflow.

Records identify the retained network, reference, configuration and exact file
versions used by a stage. Large phase arrays are not content-hashed. File size,
nanosecond modification time and HDF5 metadata detect ordinary rewrites; this is
reproducibility bookkeeping, not a cryptographic integrity guarantee or a test
of scientific accuracy. Unknown legacy products are never silently certified.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import tempfile
from datetime import datetime, timezone

import h5py
import numpy as np


class ProvenanceError(ValueError):
    pass


def _text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def read_config(path):
    values = {}
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def file_version(path):
    path = Path(path)
    st = path.stat()
    if not path.is_file() or st.st_size == 0:
        raise ProvenanceError(f"missing or empty input: {path}")
    version = {"size": st.st_size, "mtime_ns": st.st_mtime_ns}
    if path.suffix in {".json", ".cfg", ".txt"} or path.name.startswith("."):
        version["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    return version


def stack_state(work):
    path = Path(work) / "inputs/ifgramStack.h5"
    version = file_version(path)
    with h5py.File(path, "r") as f:
        dates = _date_pairs(f)
        flags = f["dropIfgram"][:]
        if flags.ndim != 1 or not np.all(np.isin(flags, [0, 1])):
            raise ProvenanceError("dropIfgram must be a one-dimensional array of boolean/0/1 flags")
        retained = np.asarray(flags, dtype=bool)
        if len(dates) != retained.size or not retained.any():
            raise ProvenanceError("empty or inconsistent retained network")
        datasets = {}
        for name in ("unwrapPhase", "connectComponent", "unwrapPhase_phaseClosure"):
            if name in f:
                ds = f[name]
                datasets[name] = {"shape": list(ds.shape), "dtype": str(ds.dtype),
                                  "attrs": {k: _text(v) for k, v in ds.attrs.items()}}
        if "unwrapPhase" not in datasets:
            raise ProvenanceError("missing unwrapPhase")
        shape = datasets["unwrapPhase"]["shape"]
        if len(shape) != 3 or shape[0] != len(dates) or min(shape) < 1:
            raise ProvenanceError("unwrapPhase dimensions do not match the date-pair table")
        if "connectComponent" not in datasets or datasets["connectComponent"]["shape"] != shape:
            raise ProvenanceError("connectComponent dimensions do not match unwrapPhase")
        return {"file": version, "date": dates, "retained": retained.tolist(),
                "attrs": {k: _text(v) for k, v in f.attrs.items()}, "datasets": datasets}


def _date_pairs(handle):
    raw = handle["date"][:]
    if raw.ndim != 2 or raw.shape[1] != 2 or not raw.shape[0]:
        raise ProvenanceError("date must be a nonempty N x 2 pair table")
    dates = [[_text(v) for v in pair] for pair in raw]
    if len({tuple(pair) for pair in dates}) != len(dates):
        raise ProvenanceError("duplicate date pairs in stack")
    for a, b in dates:
        for value in (a, b):
            if not re.fullmatch(r"\d{8}", value):
                raise ProvenanceError("stack dates must use YYYYMMDD")
            datetime.strptime(value, "%Y%m%d")
        if a >= b:
            raise ProvenanceError("stack contains a reversed or self pair")
    return dates


def selected_pairs(path):
    pairs = []
    for lineno, raw in enumerate(Path(path).read_text().splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        match = re.fullmatch(r"(\d{8})[_-](\d{8})", line)
        if not match:
            raise ProvenanceError(f"invalid pair at {path}:{lineno}: {line!r}")
        pair = tuple(match.groups())
        for value in pair:
            datetime.strptime(value, "%Y%m%d")
        if pair[0] >= pair[1] or pair in pairs:
            raise ProvenanceError(f"reversed, self or duplicate pair at {path}:{lineno}")
        pairs.append(pair)
    if not pairs:
        raise ProvenanceError("empty pair list")
    return set(pairs)


def network_selection(work, pair_file=None):
    with h5py.File(Path(work) / "inputs/ifgramStack.h5", "r") as f:
        dates = [tuple(pair) for pair in _date_pairs(f)]
    keep = selected_pairs(pair_file) if pair_file else set(dates)
    unknown = keep - set(dates)
    if unknown:
        raise ProvenanceError(f"selected pairs absent from stack: {sorted(unknown)[:8]}")
    flags = [pair in keep for pair in dates]
    return flags


def _read_record(work, name):
    try:
        record = json.loads((Path(work) / name).read_text())
        if record["schema"] != 1:
            raise ValueError("unsupported schema")
        return record
    except (OSError, ValueError, KeyError) as exc:
        raise ProvenanceError(f"missing or invalid {name}; legacy evidence is unverified") from exc


def _write_record(work, name, record):
    record = {"schema": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
              "job_id": os.environ.get("SLURM_JOB_ID", "manual"), **record}
    fd, tmp = tempfile.mkstemp(prefix=f"{name}.", dir=work)
    try:
        with os.fdopen(fd, "w") as stream:
            json.dump(record, stream, indent=2, sort_keys=True)
            stream.write("\n")
        os.replace(tmp, Path(work) / name)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)
    return record


def _versions(work, names):
    return {name: file_version(Path(work) / name) for name in names}


def _phase_parameters(work):
    cfg = read_config(Path(work) / "fullApp.cfg")
    return {k: v for k, v in cfg.items()
            if k.startswith(("mintpy.unwrapError.", "mintpy.network.", "mintpy.reference."))}


def prepare_phase(work, pair_file=None):
    """Check before reuse and return fresh or the verified phase-cache status.

    A current cache need not re-enter MintPy's phase prerequisites. In
    particular, a clean no-correction resume must not toggle its effective
    unwrap method from no to phase_closure and back before inversion.
    """
    work = Path(work)
    state = stack_state(work)
    existing = ("unwrapPhase_phaseClosure" in state["datasets"] or
                any((work / name).exists() for name in
                    ("numTriNonzeroIntAmbiguity.h5", "maskConnComp.h5", "timeseries.h5")))
    if not existing:
        return "fresh"
    record = require_phase(work)
    if state["retained"] != network_selection(work, pair_file):
        raise ProvenanceError("requested network differs from cached network; use an isolated fresh stack")
    if record["parameters"] != _phase_parameters(work):
        raise ProvenanceError("phase/reference/network settings changed; use an isolated fresh stack")
    if record["status"] not in {"correction_available", "no_integer_ambiguities", "unrepaired"}:
        raise ProvenanceError("phase cache has an unsupported status")
    return record["status"]


def ambiguity_state(work):
    """Missing or invalid measurements remain unrepaired, never 'zero'."""
    path = Path(work) / "numTriNonzeroIntAmbiguity.h5"
    with h5py.File(path, "r") as f:
        keys = [key for key in f if isinstance(f[key], h5py.Dataset) and f[key].ndim == 2]
        if len(keys) != 1:
            raise ProvenanceError("ambiguity map must contain exactly one 2-D dataset")
        counts = f[keys[0]][:]
    finite = counts[np.isfinite(counts)]
    if not finite.size or np.any(finite < 0) or np.any(finite != np.floor(finite)):
        raise ProvenanceError("ambiguity map has no valid nonnegative integer measurements")
    return int(np.count_nonzero(finite))


def record_phase(work, pair_file=None):
    state = stack_state(work)
    if state["retained"] != network_selection(work, pair_file):
        raise ProvenanceError("MintPy retained network does not match requested pair list")
    if "unwrapPhase_phaseClosure" in state["datasets"]:
        if (state["datasets"]["unwrapPhase_phaseClosure"]["shape"] !=
                state["datasets"]["unwrapPhase"]["shape"]):
            raise ProvenanceError("phase-corrected dataset has incompatible dimensions")
        status, reason = "correction_available", "MintPy produced a correction; spatial validity still needs evaluation"
    else:
        try:
            count = ambiguity_state(work)
            status = "no_integer_ambiguities" if count == 0 else "unrepaired"
            reason = f"{count} measured pixels with nonzero integer ambiguities"
        except (OSError, ValueError, KeyError) as exc:
            status, reason = "unrepaired", f"cannot verify integer ambiguities: {exc}"
    names = [name for name in ("numTriNonzeroIntAmbiguity.h5", "maskConnComp.h5")
             if (Path(work) / name).exists()]
    return _write_record(work, ".phase_provenance.json", {
        "stack": state, "parameters": _phase_parameters(work), "status": status,
        "reason": reason, "artifacts": _versions(work, names)})


def require_phase(work, allow_unrepaired=True):
    record = _read_record(work, ".phase_provenance.json")
    if record["stack"] != stack_state(work) or record["artifacts"] != _versions(work, record["artifacts"]):
        raise ProvenanceError("phase products no longer match their network/reference/input version")
    if not allow_unrepaired and record["status"] == "unrepaired":
        raise ProvenanceError(record["reason"])
    return record


def zheng_network(work, bw, nl):
    """Mirror MintPy's active inversion network and all-row loop observations.

    Pinned MintPy get_closure_phase_index uses all stored pairs for the loops,
    but estimate_bias_timeseries uses retained pairs for its bandwidth network.
    """
    state = stack_state(work)
    all_pairs = {tuple(pair) for pair in state["date"]}
    kept = {tuple(pair) for pair, flag in zip(state["date"], state["retained"]) if flag}
    dates = sorted({date for pair in kept for date in pair})
    if not 1 <= bw < len(dates) or not 2 <= nl < len(dates):
        raise ProvenanceError("Zheng bandwidth/loop length exceed the retained acquisition span")
    expected = {(dates[i], dates[i + conn]) for conn in range(1, bw + 1)
                for i in range(len(dates) - conn)}
    if kept != expected:
        raise ProvenanceError(f"retained network is not the exact complete bandwidth-{bw} network")
    for conn in set(range(1, max(2, bw) + 1)) | {nl}:
        missing = [(dates[i], dates[i + conn]) for i in range(len(dates) - conn)
                   if (dates[i], dates[i + conn]) not in all_pairs]
        if missing:
            raise ProvenanceError(f"missing source closure observations for connection {conn}: {missing[:5]}")


def stamp_stage(work, stage, names, status="usable"):
    phase = require_phase(work)
    return _write_record(work, f".{stage}_provenance.json", {
        "stack": stack_state(work), "phase_status": phase["status"], "status": status,
        "artifacts": _versions(work, names)})


def require_stage(work, stage, allow_unrepaired=False):
    record = _read_record(work, f".{stage}_provenance.json")
    require_phase(work, allow_unrepaired=allow_unrepaired)
    if record["stack"] != stack_state(work) or record["artifacts"] != _versions(work, record["artifacts"]):
        raise ProvenanceError(f"{stage} products or dependencies changed after they were recorded")
    if not allow_unrepaired and record["phase_status"] == "unrepaired":
        raise ProvenanceError(f"{stage} was generated from unrepaired phase")
    if stage == "closureCor":
        require_stage(work, "inversion", allow_unrepaired)
        require_stage(work, "bias", allow_unrepaired)
    return record


def apply_inputs(work, allow_unrepaired=False):
    work = Path(work)
    require_stage(work, "inversion", allow_unrepaired)
    require_stage(work, "bias", allow_unrepaired)
    cfg = read_config(work / "smallbaselineApp.cfg")
    # This workflow currently supports SET followed by deramp. Guessing the
    # newest existing filename can select an obsolete or incompatible branch.
    for key in ("mintpy.ionosphericDelay.method", "mintpy.troposphericDelay.method", "mintpy.topographicResidual"):
        if cfg.get(key, "no").lower() not in {"no", "none", "false", "auto"}:
            raise ProvenanceError(f"unsupported correction order: {key}={cfg[key]}")
    setting = cfg.get("mintpy.solidEarthTides", "").lower()
    if setting not in {"yes", "true", "no", "false"}:
        raise ProvenanceError("solidEarthTides must be resolved explicitly")
    base = "timeseries_SET.h5" if setting in {"yes", "true"} else "timeseries.h5"
    file_version(work / base)
    ramp = cfg.get("mintpy.deramp", "")
    if ramp not in {"no", "linear", "quadratic", "linear_range", "linear_azimuth", "quadratic_range", "quadratic_azimuth"}:
        raise ProvenanceError(f"unsupported or unresolved deramp setting: {ramp!r}")
    mask = cfg.get("mintpy.deramp.maskFile", "auto")
    if mask == "auto":
        mask = "maskTempCoh.h5"
    if ramp != "no" and mask != "no":
        file_version(work / mask)
    with h5py.File(work / base, "r") as a, h5py.File(work / "timeseriesBias.h5", "r") as b:
        if a["timeseries"].shape != b["timeseries"].shape or not np.array_equal(a["date"][:], b["date"][:]):
            raise ProvenanceError("bias and base time series must have identical grids and acquisition dates")
        for key in ("X_FIRST", "Y_FIRST", "X_STEP", "Y_STEP", "UNIT"):
            if key not in a.attrs or key not in b.attrs or _text(a.attrs[key]) != _text(b.attrs[key]):
                raise ProvenanceError(f"bias and base metadata disagree or are missing: {key}")
    return base, ramp, mask


def select_velocity(work, allow_unverified=False):
    """Require explicit opt-in before comparing legacy or incomplete products."""
    work = Path(work)
    if (work / ".workflow_incomplete").exists():
        if allow_unverified and (work / "velocity.h5").exists():
            return work / "velocity.h5", {"status": "UNVERIFIED_DIAGNOSTIC",
                                          "reason": "latest workflow attempt is incomplete or still running"}
        raise ProvenanceError("latest workflow attempt is incomplete or still running")
    failures = []
    for name, stage in (("velocityClosureCor.h5", "closureCor"), ("velocity.h5", "inversion")):
        if not (work / name).exists():
            continue
        try:
            record = require_stage(work, stage)
            return work / name, record
        except (OSError, ValueError, KeyError) as exc:
            failures.append(f"{name}: {exc}")
    if allow_unverified and (work / "velocity.h5").exists():
        return work / "velocity.h5", {"status": "UNVERIFIED_DIAGNOSTIC", "reason": "; ".join(failures)}
    raise ProvenanceError("no verified velocity product: " + "; ".join(failures))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("selection", "prepare-phase", "record-phase", "require-phase", "zheng-network", "stamp", "require", "apply-inputs"))
    parser.add_argument("work", type=Path)
    parser.add_argument("--pairs")
    parser.add_argument("--stage")
    parser.add_argument("--artifact", action="append", default=[])
    parser.add_argument("--allow-unrepaired", action="store_true")
    parser.add_argument("--bw", type=int, default=20)
    parser.add_argument("--nl", type=int, default=20)
    args = parser.parse_args(argv)
    try:
        if args.action == "selection":
            flags = network_selection(args.work, args.pairs)
            print(",".join(str(i) for i, keep in enumerate(flags) if not keep) or "no")
        elif args.action == "prepare-phase":
            print(prepare_phase(args.work, args.pairs))
        elif args.action == "record-phase":
            result = record_phase(args.work, args.pairs)
            print(result["status"])
            print(result["reason"], file=sys.stderr)
        elif args.action == "require-phase":
            require_phase(args.work, args.allow_unrepaired)
        elif args.action == "zheng-network":
            zheng_network(args.work, args.bw, args.nl)
        elif args.action == "stamp":
            if not args.stage or not args.artifact:
                raise ProvenanceError("stamp requires a stage and artifacts")
            stamp_stage(args.work, args.stage, args.artifact)
        elif args.action == "require":
            require_stage(args.work, args.stage, args.allow_unrepaired)
        elif args.action == "apply-inputs":
            print("\n".join(apply_inputs(args.work, args.allow_unrepaired)))
    except (OSError, ValueError, KeyError) as exc:
        print(f"PROVENANCE_FAILED: {exc}", file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
