"""Tests for the read-only phase-closure network selector."""

import importlib.util
import hashlib
import json
import math
import os
import sys

import h5py
import numpy as np
import pytest


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MODULE_PATH = os.path.join(ROOT, "select_phase_closure_network.py")
SPEC = importlib.util.spec_from_file_location("select_phase_closure_network", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def _write_stack(path):
    dates = ["20200101", "20200113", "20200125", "20200206", "20200218"]
    pairs = []
    for gap in range(1, 4):
        for index in range(len(dates) - gap):
            pairs.append((dates[index], dates[index + gap], gap))

    shape = (len(pairs), 8, 8)
    phase = np.zeros(shape, dtype=np.float32)
    conncomp = np.ones(shape, dtype=np.uint8)
    for pair_index, (_, _, gap) in enumerate(pairs):
        if gap == 3:
            conncomp[pair_index] = 0

    # Add one regional two-cycle error. Referencing at (0, 0) preserves it in
    # the lower half of the synthetic image and creates non-zero closures.
    phase[0, 4:, :] = 4.0 * math.pi
    with h5py.File(path, "w") as stack:
        stack.create_dataset(
            "date",
            data=np.array([(date1, date2) for date1, date2, _ in pairs], dtype="S8"),
        )
        stack.create_dataset("unwrapPhase", data=phase)
        stack.create_dataset("connectComponent", data=conncomp)
        stack.attrs["REF_Y"] = "0"
        stack.attrs["REF_X"] = "0"


def test_selects_largest_viable_sequential_bandwidth(tmp_path):
    stack_path = tmp_path / "ifgramStack.h5"
    _write_stack(stack_path)

    results, selected, mask, pairs = MODULE.scan_sequential_networks(
        stack_path,
        min_bandwidth=2,
        max_bandwidth=3,
        min_region_pixels=16,
    )

    assert [result.bandwidth for result in results] == [1, 2, 3]
    assert selected is not None
    assert selected.bandwidth == 2
    assert selected.interferograms == 7
    assert selected.triangles == 3
    assert selected.connected_date_graph
    assert selected.complete_sequential_levels
    assert mask is not None and int(mask.sum()) == 64
    assert len(pairs) == 7
    assert results[-1].common_pixels == 0


def test_sampled_closure_detects_integer_ambiguity(tmp_path):
    stack_path = tmp_path / "ifgramStack.h5"
    _write_stack(stack_path)
    _, selected, mask, pairs = MODULE.scan_sequential_networks(
        stack_path,
        min_bandwidth=2,
        max_bandwidth=3,
        min_region_pixels=16,
    )

    assert selected is not None and mask is not None
    stats = MODULE.sample_closure_stats(stack_path, pairs, mask, sample_count=8)

    assert stats["triangles"] == 3
    assert stats["nonzero_integer_ambiguity_fraction"] > 0
    assert stats["max_abs_integer_cycles"] == 2.0
    assert stats["diagnostic_valid"] is True
    assert stats["finite_coverage_fraction"] == 1.0
    assert stats["interferograms_with_invalid_reference"] == []


def test_no_triangle_network_is_not_viable(tmp_path):
    stack_path = tmp_path / "ifgramStack.h5"
    _write_stack(stack_path)

    try:
        MODULE.scan_sequential_networks(
            stack_path,
            min_bandwidth=1,
            max_bandwidth=1,
            min_region_pixels=16,
        )
    except ValueError as exc:
        assert "at least 2" in str(exc)
    else:
        raise AssertionError("bandwidth 1 must be rejected")


def _write_unviable_stack(path):
    """A stack whose every interferogram is fully unwrapped-invalid.

    No bandwidth can produce a non-empty common region, so main() must take the
    bail path.
    """
    dates = ["20200101", "20200113", "20200125", "20200206", "20200218"]
    pairs = []
    for gap in range(1, 4):
        for index in range(len(dates) - gap):
            pairs.append((dates[index], dates[index + gap]))

    shape = (len(pairs), 8, 8)
    with h5py.File(path, "w") as stack:
        stack.create_dataset("date", data=np.array(pairs, dtype="S8"))
        stack.create_dataset("unwrapPhase", data=np.zeros(shape, dtype=np.float32))
        stack.create_dataset("connectComponent", data=np.zeros(shape, dtype=np.uint8))
        stack.attrs["REF_Y"] = "0"
        stack.attrs["REF_X"] = "0"


def test_bail_path_announces_that_no_region_mask_was_written(tmp_path, capsys):
    """Defect 20: the mask writer sits after the bail path's return.

    A caller that passes --write-region-mask on a stack with no viable
    bandwidth must still be told no mask was produced, otherwise a missing
    file is indistinguishable from a crash.
    """
    stack_path = tmp_path / "ifgramStack.h5"
    _write_unviable_stack(stack_path)
    mask_path = tmp_path / "region.npy"

    code = MODULE.main([
        str(stack_path),
        "--min-bandwidth", "2",
        "--max-bandwidth", "3",
        "--min-region-pixels", "16",
        "--write-region-mask", str(mask_path),
    ])

    out = capsys.readouterr().out
    assert code == 2, f"expected bail exit 2, got {code}"
    assert "NO_VIABLE_COMPLETE_SEQUENTIAL_PHASE_CLOSURE_NETWORK" in out
    assert "NO_REGION_MASK_WRITTEN" in out, (
        "bail path must announce the absent mask; it is unreachable in the "
        "writer block below the return"
    )
    assert not mask_path.exists(), "no mask file may be written on the bail path"


def _scan(path):
    return MODULE.scan_sequential_networks(
        path, min_bandwidth=2, max_bandwidth=3, min_region_pixels=16,
    )


def _cli_args(stack_path, tmp_path):
    return [
        str(stack_path), "--max-bandwidth", "3", "--min-region-pixels", "16",
        "--write-pairs", str(tmp_path / "pairs.txt"),
        "--json", str(tmp_path / "report.json"),
        "--write-region-mask", str(tmp_path / "anchor.npy"),
    ]


@pytest.mark.parametrize("invalid_reference", ["masked", "nan", "missing", "negative", "outside"])
def test_invalid_reference_blocks_pair_output_but_preserves_candidate_mask(
    tmp_path, invalid_reference,
):
    stack_path = tmp_path / "ifgramStack.h5"
    _write_stack(stack_path)
    with h5py.File(stack_path, "r+") as stack:
        if invalid_reference == "masked":
            stack["connectComponent"][0, 0, 0] = 0
        elif invalid_reference == "nan":
            stack["unwrapPhase"][0, 0, 0] = np.nan
        elif invalid_reference == "missing":
            del stack.attrs["REF_Y"]
        else:
            stack.attrs["REF_Y"] = "-1" if invalid_reference == "negative" else "8"

    assert MODULE.main(_cli_args(stack_path, tmp_path)) == 3
    report = json.loads((tmp_path / "report.json").read_text())
    assert report["status"] == "candidate_diagnostic_invalid"
    assert report["selected"]["bandwidth"] == 2
    assert not (tmp_path / "pairs.txt").exists()
    assert np.load(tmp_path / "anchor.npy").any()
    if invalid_reference in ("masked", "nan"):
        stats = report["closure_samples"]
        assert stats["diagnostic_valid"] is False
        assert stats["interferograms_with_invalid_reference"] == ["20200101_20200113"]
        assert stats["nonzero_integer_ambiguity_fraction"] is None
    else:
        assert report["diagnostic_error"]


def test_partial_finite_samples_do_not_report_a_survivor_fraction(tmp_path):
    stack_path = tmp_path / "ifgramStack.h5"
    _write_stack(stack_path)
    with h5py.File(stack_path, "r+") as stack:
        stack["unwrapPhase"][0, 1:, :] = np.nan
    _, selected, mask, pairs = _scan(stack_path)
    assert selected is not None

    stats = MODULE.sample_closure_stats(stack_path, pairs, mask, sample_count=64)

    assert 0 < stats["valid_triangle_samples"] < stats["total_triangle_samples"]
    assert stats["diagnostic_valid"] is False
    assert stats["interferograms_with_invalid_samples"] == [
        {"pair": "20200101_20200113", "invalid_samples": 56}
    ]
    assert stats["nonzero_integer_ambiguity_fraction"] is None
    assert stats["median_bad_triplet_fraction_per_sample"] is None


def test_current_dropped_pairs_are_only_reenabled_by_explicit_opt_in(tmp_path):
    stack_path = tmp_path / "ifgramStack.h5"
    _write_stack(stack_path)
    with h5py.File(stack_path, "r+") as stack:
        flags = np.ones(len(stack["date"]), dtype=bool)
        flags[0] = False
        stack.create_dataset("dropIfgram", data=flags)

    assert _scan(stack_path)[1] is None
    _, selected, _, pairs = MODULE.scan_sequential_networks(
        stack_path, max_bandwidth=3, min_region_pixels=16, include_dropped=True,
    )
    assert selected is not None and selected.bandwidth == 2
    assert any(pair.stack_index == 0 for pair in pairs)
    with h5py.File(stack_path, "r") as stack:
        assert not stack["dropIfgram"][0]


def test_success_records_exact_pair_fingerprint_and_leaves_source_unchanged(tmp_path):
    stack_path = tmp_path / "ifgramStack.h5"
    _write_stack(stack_path)
    original = stack_path.read_bytes()

    assert MODULE.main(_cli_args(stack_path, tmp_path)) == 0

    report = json.loads((tmp_path / "report.json").read_text())
    assert report["status"] == "candidate_diagnostic_complete"
    assert report["source"] == MODULE.source_provenance(stack_path)
    assert report["selected_pair_list_sha256"] == hashlib.sha256(
        (tmp_path / "pairs.txt").read_bytes()
    ).hexdigest()
    assert "phase contents are not hashed" in report["source"]["fingerprint_scope"]
    assert "peat support" in report["selection_scope"]
    assert stack_path.read_bytes() == original


def test_no_candidate_writes_scoped_report_and_no_pair_list(tmp_path):
    stack_path = tmp_path / "ifgramStack.h5"
    _write_unviable_stack(stack_path)

    assert MODULE.main(_cli_args(stack_path, tmp_path)) == 2

    report = json.loads((tmp_path / "report.json").read_text())
    assert report["status"] == "no_complete_sequential_candidate"
    assert "No conclusion about arbitrary subnetworks" in report["selection_scope"]
    assert not (tmp_path / "pairs.txt").exists()
    assert not (tmp_path / "anchor.npy").exists()


def test_output_cannot_replace_input_or_existing_artifacts(tmp_path):
    stack_path = tmp_path / "ifgramStack.h5"
    _write_stack(stack_path)
    original = stack_path.read_bytes()
    assert MODULE.main([str(stack_path), "--json", str(stack_path)]) == 3
    assert stack_path.read_bytes() == original

    (tmp_path / "pairs.txt").write_text("old candidate\n")
    assert MODULE.main(_cli_args(stack_path, tmp_path)) == 3
    assert (tmp_path / "pairs.txt").read_text() == "old candidate\n"
    assert not (tmp_path / "report.json").exists()


def test_output_paths_must_be_distinct(tmp_path):
    stack_path = tmp_path / "ifgramStack.h5"
    _write_stack(stack_path)
    output = tmp_path / "artifact"

    assert MODULE.main([
        str(stack_path), "--json", str(output), "--write-pairs", str(output),
    ]) == 3
    assert not output.exists()


def test_changed_source_prevents_artifact_publication(tmp_path, monkeypatch):
    stack_path = tmp_path / "ifgramStack.h5"
    _write_stack(stack_path)
    original_sampler = MODULE.sample_closure_stats

    def sample_then_change_source(*args, **kwargs):
        result = original_sampler(*args, **kwargs)
        with h5py.File(stack_path, "r+") as stack:
            stack.attrs["REF_Y"] = "1"
        return result

    monkeypatch.setattr(MODULE, "sample_closure_stats", sample_then_change_source)
    assert MODULE.main(_cli_args(stack_path, tmp_path)) == 3
    assert not (tmp_path / "pairs.txt").exists()
    assert not (tmp_path / "report.json").exists()
    assert not (tmp_path / "anchor.npy").exists()


def test_anchor_filter_can_select_smaller_common_region():
    mask = np.zeros((12, 24), dtype=bool)
    mask[:8, :8] = True
    mask[:6, 12:18] = True
    aoi = np.zeros_like(mask)
    aoi[:, 12:] = True

    selected, area, count = MODULE._largest_regions(
        mask, min_region_pixels=16, aoi_mask=aoi, min_aoi_pixels=30,
    )

    assert area == 36 and count == 2
    assert not selected[:, :8].any()
    assert (selected & aoi).sum() == 36


def test_bbox_uses_pixel_centres_and_rejects_projected_grids():
    attrs = {"X_FIRST": "100", "Y_FIRST": "1", "X_STEP": "1", "Y_STEP": "-1"}
    mask = MODULE._aoi_mask_from_bbox((2, 2), attrs, (100.1, 0.1, 100.9, 0.9))
    np.testing.assert_array_equal(mask, [[True, False], [False, False]])
    with pytest.raises(ValueError, match="longitude/latitude grid"):
        MODULE._aoi_mask_from_bbox((2, 2), dict(attrs, X_UNIT="meters"), (100, -1, 102, 1))


def test_pair_output_requires_provenance_report(tmp_path):
    with pytest.raises(SystemExit) as exc:
        MODULE.main([str(tmp_path / "stack.h5"), "--write-pairs", str(tmp_path / "pairs.txt")])
    assert exc.value.code == 2
