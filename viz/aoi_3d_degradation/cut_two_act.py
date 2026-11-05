#!/usr/bin/env python3
"""Two-act cut: PeatGuard_AOI_3D_degradation_then_fire_SLIDE.mp4 (about 34 s).

Act 1, degradation, 1990 to 2005 (global frames 0 to 84): the green base with the
pre-existing degraded red growing frame by frame, no flames and no char even
where fires occurred; headline "Peat degraded: N percent of the box". The red is
interpolated per frame between the masks of consecutive key years (1990, 1995,
2000, 2001 to 2006): cells new in the later mask fade in with the fraction of the
segment elapsed, so the Landsat-epoch steps and the yearly steps both grow
smoothly. The last 4 act 1 frames are dissolved with 2006 (rendered at the same
global frames) to open act 2.

Act 2, fire, 2006 to 2026 (render.py --spec twoact, rendered directly at 12 fps,
no frame reuse): red as of the year, black char from the 2001 onward fires (so
2006 opens with the 2001 to 2005 burn history), selective 3D flames, the
fire-on-degraded headline and the full legend; 13 frames per year, 23 for 2015
and 2019 (1.2 s longer), 21 for 2024 (12-frame dissolve into the imagery 2026),
51 for 2026 (3 s hold); 4-frame data dissolves elsewhere.

The camera is one eased function of the global frame index over the whole cut
(render.camera_t_at); it never stops, including through holds and dissolves. The
final list is checked to be one entry per global frame with a strictly monotonic
azimuth whose consecutive steps stay within 10 percent of each other.

Usage (from the repo root; renders act 1 and act 2, then assembles):
  .venv/bin/python3 deliverables/aoi_3d_degradation/cut_two_act.py
  options: --skip-act1 / --skip-act2 to reuse frames already in --work
"""

import argparse
import os
import shutil
import subprocess
import sys
import time

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import render as R  # noqa: E402

FFMPEG = os.environ.get("FFMPEG", shutil.which("ffmpeg") or "ffmpeg")
ACT1_KEYS = [1990, 1995, 2000, 2001, 2002, 2003, 2004, 2005, 2006]
# frames per act 1 segment (label = the segment's year); the 2005 segment carries
# the 4 boundary frames that dissolve into 2006
ACT1_SEGMENT_FRAMES = {1990: 15, 1995: 15, 2000: 9, 2001: 8, 2002: 8, 2003: 8, 2004: 9, 2005: 13}
ACT1_RAMP_2005 = 9      # the 2005 red ramp completes after this many frames; the rest hold the data


def render_act1(outdir):
    """Act 1 frames act1_NNN.png at global frames 0 .. prefix - 1 (raw scene saved
    for the frames dissolved into 2006). Returns the per-frame degraded percent."""
    R.set_spec("twoact")
    prefix = R.SPECS["twoact"]["prefix"]
    assert sum(ACT1_SEGMENT_FRAMES.values()) == prefix, (sum(ACT1_SEGMENT_FRAMES.values()), prefix)
    n_x = R.xfade_frames_for(R.YEARS[0])
    os.makedirs(outdir, exist_ok=True)
    masks = {y: R.degraded_mask(y).astype(np.float32) for y in ACT1_KEYS}
    R.BG_OVERRIDE = "#000000"
    pl, _, _ = R.build_scene(1990, with_fire=False)
    z_full, _ = R.build_z_full(1990)
    grid = R.build_surface(z_full)
    pcts = []
    g = 0
    for k, year in enumerate(ACT1_KEYS[:-1]):
        nxt = ACT1_KEYS[k + 1]
        n = ACT1_SEGMENT_FRAMES[year]
        ramp = ACT1_RAMP_2005 if year == 2005 else n
        for j in range(n):
            w = min(j / float(ramp), 1.0)
            red = masks[year] + w * np.clip(masks[nxt] - masks[year], 0.0, 1.0)
            pct = 100.0 * float(red.mean())
            tex = R.build_texture(year, None, red_field=red, with_char=False)
            pl.add_mesh(grid, texture=tex, smooth_shading=True, show_edges=False,
                        ambient=0.50, diffuse=0.62, specular=0.03, specular_power=12, name="terrain")
            t_f = time.time()
            frame = R.render_global_frame(pl, None, g)
            path = os.path.join(outdir, "act1_%03d.png" % g)
            if g >= prefix - n_x:
                Image.fromarray(frame).save(R.raw_path(path), compress_level=3)
            R.compose_minimal(frame, year, degraded_pct=pct).save(path, compress_level=3)
            pcts.append(pct)
            print("%s  %d  g=%d az=%.3f red=%.1f%%  %.2f s" % (os.path.basename(path), year, g,
                  R.camera_azimuth_deg(g), pct, time.time() - t_f), file=sys.stderr)
            g += 1
    assert g == prefix, g
    pl.close()
    with open(os.path.join(outdir, "act1_pct.txt"), "w") as f:
        f.write("\n".join("%.3f" % p for p in pcts) + "\n")
    return pcts


def render_act2(outdir, jobs):
    cmd = [sys.executable, os.path.join(HERE, "render.py"), "--all", "--spec", "twoact", "--jobs", str(jobs),
           "--minimal", "--bg", "#000000", "--outdir", outdir]
    subprocess.run(cmd, check=True)


def blend(a, b, w, out, year_a, year_b, k, n_x, pct_a=None):
    """Blend the raw scenes of two frames rendered at the same global frame and
    compose the HUD once (year_a's text before the midpoint, year_b's after,
    four-frame fade), see assemble.py."""
    if not os.path.exists(out):
        scene = Image.blend(Image.open(R.raw_path(a)).convert("RGB"), Image.open(R.raw_path(b)).convert("RGB"), w)
        before = k < n_x // 2
        R.compose_minimal(scene, year_a if before else year_b, degraded_pct=pct_a if before else None,
                          hud_alpha=R.hud_alpha_for(k, n_x)).save(out, compress_level=3)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", default=os.path.join(R.DATA_DIR, "frames_two_act"))
    ap.add_argument("--out", default=os.path.join(R.DATA_DIR, "PeatGuard_AOI_3D_degradation_then_fire_SLIDE.mp4"))
    ap.add_argument("--jobs", type=int, default=7)
    ap.add_argument("--skip-act1", action="store_true", help="reuse act 1 frames already in --work")
    ap.add_argument("--skip-act2", action="store_true", help="reuse act 2 frames already in --work")
    args = ap.parse_args(argv)

    work = os.path.abspath(args.work)
    xdir = os.path.join(work, "xfade")
    os.makedirs(xdir, exist_ok=True)
    R.set_spec("twoact")
    sch = R.schedule()
    years = R.YEARS
    prefix = R.SPECS["twoact"]["prefix"]
    total = R.total_frames()

    if not args.skip_act1:
        render_act1(work)
    with open(os.path.join(work, "act1_pct.txt")) as f:
        pcts = [float(x) for x in f.read().split()]
    if not args.skip_act2:
        render_act2(work, args.jobs)

    entries, gs = [], []
    first = years[0]
    n_x0 = R.xfade_frames_for(first)
    for g in range(prefix):
        p = os.path.join(work, "act1_%03d.png" % g)
        k = g - (prefix - n_x0)
        if k >= 0:
            q = os.path.join(work, "%d_pre_%02d.png" % (first, k))
            p = blend(p, q, (k + 1) / float(n_x0 + 1), os.path.join(xdir, "x_act1_%d_%02d.png" % (first, k)),
                      2005, first, k, n_x0, pct_a=pcts[g])
        entries.append(p)
        gs.append(g)
    for idx, y in enumerate(years):
        s = sch[y]
        nxt = years[idx + 1] if idx + 1 < len(years) else None
        n_x = R.xfade_frames_for(nxt) if nxt is not None else 0
        for i in range(s["n"]):
            p = os.path.join(work, "%d_%03d.png" % (y, i))
            k = i - (s["n"] - n_x)
            if nxt is not None and k >= 0:
                q = os.path.join(work, "%d_pre_%02d.png" % (nxt, k))
                p = blend(p, q, (k + 1) / float(n_x + 1), os.path.join(xdir, "x_%d_%d_%02d.png" % (y, nxt, k)),
                          y, nxt, k, n_x)
            entries.append(p)
            gs.append(s["tick0"] + i)
    assert gs == list(range(total)), "list does not cover every global frame once"
    for p in entries:
        if not os.path.exists(p):
            raise SystemExit("missing frame: %s" % p)
    max_step, max_ratio = R.check_azimuth_steps(gs)
    print("azimuth: %d frames, max step %.4f deg per frame, neighbouring steps within %.2f percent, monotonic"
          % (len(gs), max_step, 100 * (max_ratio - 1)), file=sys.stderr)

    list_path = os.path.join(work, "list.txt")
    dur = 1.0 / R.FPS
    with open(list_path, "w") as f:
        for p in entries:
            f.write("file '%s'\nduration %.6f\n" % (p, dur))
        f.write("file '%s'\n" % entries[-1])
    print("%d frames, %.1f s" % (len(entries), len(entries) * dur), file=sys.stderr)

    cmd = [FFMPEG, "-y", "-loglevel", "error", "-f", "concat", "-safe", "0", "-i", list_path,
           "-vf", "fps=%d,scale=1920:1080:flags=lanczos,format=yuv420p" % R.FPS,
           "-c:v", "libx264", "-preset", "slow", "-crf", "18", "-profile:v", "high", "-level", "4.0",
           "-pix_fmt", "yuv420p", "-movflags", "+faststart", os.path.abspath(args.out)]
    subprocess.run(cmd, check=True)
    print(os.path.abspath(args.out))


if __name__ == "__main__":
    main()
