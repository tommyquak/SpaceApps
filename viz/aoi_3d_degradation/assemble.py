#!/usr/bin/env python3
"""Assemble the v4 flyover MP4 from a frames_v4 directory written by render.py.

One list entry per global frame index g. The camera is a single eased function
of g (render.camera_t_at) and never stops: a year's "hold" is only its data
layer persisting while the orbit continues. For every year of the schedule the
list takes its n frames; the last xfade frames before a next year are replaced
by crossfade frames: year Y's frame at g blended with year Y+1's pre-roll frame
rendered at the same g (same camera, same flicker tick), so the transition is a
pure data dissolve under the moving camera. The azimuth sequence of the final
list is checked to be strictly monotonic with neighbouring steps within 10
percent (render.check_azimuth_steps) and the maximum step is printed.

Usage (from the repo root):
  .venv/bin/python3 deliverables/aoi_3d_degradation/assemble.py \
      --frames deliverables/aoi_3d_degradation/frames_v4 \
      --out deliverables/aoi_3d_degradation/PeatGuard_AOI_3D_fires_yearly_1990_2026_SLIDE_v4.mp4

The crossfade blends the raw scenes (render.py writes <frame>.raw.png for every
hold and pre-roll frame) and composes the slide-mode HUD once per frame, the
year switching at the midpoint with a four-frame text fade, so the numerals,
headline and legend never double-expose. Slide mode only.

Encoding: 12 fps, libx264, yuv420p, CRF 18, High profile level 4.0, faststart,
1920x1080 (PowerPoint on Mac and Windows).
"""

import argparse
import os
import shutil
import subprocess
import sys

from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import render as R  # noqa: E402

FFMPEG = os.environ.get("FFMPEG", shutil.which("ffmpeg") or "ffmpeg")


def crossfade(a_path, b_path, w, out_path, year_a, year_b, k, n_x):
    """Blend the RAW scenes (no HUD) of the two frames, then compose the HUD once:
    year_a's text before the midpoint, year_b's after, fading out and in over
    four frames so nothing double-exposes."""
    a = Image.open(R.raw_path(a_path)).convert("RGB")
    b = Image.open(R.raw_path(b_path)).convert("RGB")
    scene = Image.blend(a, b, w)
    epoch = year_a if k < n_x // 2 else year_b
    R.compose_minimal(scene, epoch, hud_alpha=R.hud_alpha_for(k, n_x)).save(out_path, compress_level=3)


def build_list(frames_dir, xfade_dir):
    """One entry per global frame index g; returns (paths, gs)."""
    sch = R.schedule()
    years = R.YEARS
    entries, gs = [], []
    os.makedirs(xfade_dir, exist_ok=True)
    for idx, y in enumerate(years):
        s = sch[y]
        nxt = years[idx + 1] if idx + 1 < len(years) else None
        n_x = R.xfade_frames_for(nxt) if nxt is not None else 0
        for i in range(s["n"]):
            p = os.path.join(frames_dir, "%d_%03d.png" % (y, i))
            k = i - (s["n"] - n_x)
            if nxt is not None and k >= 0:
                q = os.path.join(frames_dir, "%d_pre_%02d.png" % (nxt, k))
                w = (k + 1) / float(n_x + 1)
                out = os.path.join(xfade_dir, "xfade_%d_%d_%02d.png" % (y, nxt, k))
                if not os.path.exists(out):
                    crossfade(p, q, w, out, y, nxt, k, n_x)
                p = out
            if not os.path.exists(p):
                raise SystemExit("missing frame: %s" % p)
            entries.append(p)
            gs.append(s["tick0"] + i)
    assert gs == list(range(R.total_frames())), "list does not cover every global frame once"
    return entries, gs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frames", default=os.path.join(R.DATA_DIR, "frames_v4"))
    ap.add_argument("--out", default=os.path.join(R.DATA_DIR, "PeatGuard_AOI_3D_fires_yearly_1990_2026_SLIDE_v4.mp4"))
    ap.add_argument("--preset", default="slow")
    args = ap.parse_args(argv)

    frames_dir = os.path.abspath(args.frames)
    xfade_dir = os.path.join(frames_dir, "xfade")
    R.set_spec("v4")
    entries, gs = build_list(frames_dir, xfade_dir)
    max_step, max_ratio = R.check_azimuth_steps(gs)
    print("azimuth: %d frames, max step %.4f deg per frame, neighbouring steps within %.2f percent, monotonic"
          % (len(gs), max_step, 100 * (max_ratio - 1)), file=sys.stderr)
    list_path = os.path.join(frames_dir, "list.txt")
    dur = 1.0 / R.FPS
    with open(list_path, "w") as f:
        for p in entries:
            f.write("file '%s'\nduration %.6f\n" % (p, dur))
        f.write("file '%s'\n" % entries[-1])  # concat demuxer: last file once more, no duration
    print("%d frames, %.1f s" % (len(entries), len(entries) * dur), file=sys.stderr)

    cmd = [FFMPEG, "-y", "-f", "concat", "-safe", "0", "-i", list_path,
           "-vf", "fps=%d,scale=1920:1080:flags=lanczos,format=yuv420p" % R.FPS,
           "-c:v", "libx264", "-preset", args.preset, "-crf", "18",
           "-profile:v", "high", "-level", "4.0", "-pix_fmt", "yuv420p",
           "-movflags", "+faststart", os.path.abspath(args.out)]
    subprocess.run(cmd, check=True)
    print(os.path.abspath(args.out))


if __name__ == "__main__":
    main()
