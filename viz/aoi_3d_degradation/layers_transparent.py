#!/usr/bin/env python3
"""Separable transparent layers of the two-act cut for PowerPoint.

Four aligned render passes on the two-act schedule (render.SPECS["twoact"],
411 frames, 34.25 s, the same global camera function, lighting and frame count),
each with a transparent background (RGBA), so they stack pixel for pixel:

  land  base terrain only: green hillshade, rivers, the NICFI imagery in the
        2026 segment; no red, no char, no fire
  peat  the same terrain mesh textured with an RGBA texture whose alpha is the
        red pre-existing degraded layer plus the black char (alpha 0 elsewhere),
        so only the red and char float in 3D at the right place
  fire  only the flame and smoke meshes
  all   the full composite (land + peat + fire baked as in the cut)

Nothing but geometry: no dam sites, labels, north mark, legend, text or year,
and no cyan outline, so all = land over peat over fire exactly. Year dissolves
are premultiplied RGBA blends of the two years rendered at the same global
frame. Each pass is encoded twice at 12 fps with alpha: ProRes 4444
(<pass>_prores4444.mov, PowerPoint on Mac) and HEVC with alpha
(<pass>_hevc_alpha.mov, hvc1, PowerPoint 365 on Mac and Windows 10 and later),
plus the last 2026 frame as a transparent PNG (<pass>_2026_last.png).

Usage (from the repo root):
  .venv/bin/python3 deliverables/aoi_3d_degradation/layers_transparent.py --test   # 2 s check
  .venv/bin/python3 deliverables/aoi_3d_degradation/layers_transparent.py          # full set
Internal: --job act1 | --job <year> renders one job into --outdir.
"""

import argparse
import json
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
import cut_two_act as C  # noqa: E402
import pyvista as pv  # noqa: E402

FFMPEG = os.environ.get("FFMPEG", shutil.which("ffmpeg") or "ffmpeg")
PASSES = ["land", "peat", "fire", "all"]
OUT_DEFAULT = os.path.join(R.DATA_DIR, "layers_transparent")


# ----------------------------------------------------------------------------
# textures and scenes per pass
# ----------------------------------------------------------------------------
def land_texture(year):
    """Base terrain: green hillshade (imagery for the 2026 segment) and rivers."""
    zero = np.zeros((R.NROW, R.NCOL), np.float32)
    return R.build_texture(year, None, red_field=zero, with_char=False)


def peat_rgba(year, red_field=None, with_char=True):
    """RGBA texture: red degraded layer and black char over nothing."""
    deg = R.degraded_mask(year).astype(np.float32) if red_field is None else np.clip(red_field, 0, 1).astype(np.float32)
    a_r = (R.IMAGERY_DEGRADED_ALPHA if year == R.IMAGERY_EPOCH else R.DEGRADED_ALPHA) * deg
    char = R.build_char(year) if with_char else None
    a_c = np.zeros_like(a_r) if char is None else char.astype(np.float32)
    # premultiplied: red at a_r, then black char over it
    prem = (1.0 - a_c)[..., None] * a_r[..., None] * R.PALETTE["degraded"]
    alpha = a_c + (1.0 - a_c) * a_r
    rgb = np.where(alpha[..., None] > 1e-4, prem / np.maximum(alpha, 1e-4)[..., None], 0.0)
    rgba = np.concatenate([rgb, alpha[..., None] * 255.0], axis=-1)
    return pv.Texture(np.ascontiguousarray(np.clip(rgba, 0, 255).astype(np.uint8)))


def all_texture(year, fire, red_field=None, with_char=True):
    return R.build_texture(year, fire, red_field=red_field, with_char=with_char)


TERRAIN_KW = dict(smooth_shading=True, show_edges=False, ambient=0.50, diffuse=0.62, specular=0.03,
                  specular_power=12, name="terrain")


def new_plotter():
    pl = pv.Plotter(off_screen=True, window_size=R.WINDOW, lighting="none")
    pl.set_background("#000000")
    pl.enable_anti_aliasing("fxaa")
    pl.enable_depth_peeling(number_of_peels=6, occlusion_ratio=0.0)
    pl.add_light(pv.Light(position=(-60, 80, 90), focal_point=(0, 0, 0), light_type="scene light",
                          intensity=0.95, color="white"))
    pl.add_light(pv.Light(light_type="headlight", intensity=0.25))
    pl.camera.view_angle = R.CAM_VIEW_ANGLE
    return pl


class PassScene:
    """One pass of one year: plotter, terrain grid, optional fire animator."""

    def __init__(self, year, pass_name, with_fire):
        self.year, self.pass_name = year, pass_name
        self.z_full, _ = R.build_z_full(year)
        self.grid = R.build_surface(self.z_full)
        self.fire = R.build_fire_glow(year) if with_fire else None
        self.pl = new_plotter()
        self.anim = None
        if pass_name in ("fire", "all") and self.fire is not None and float(self.fire[0].max()) > R.FLAME_MIN_ALPHA:
            X, Y = R.work_xy()
            self.anim = R.FireAnimator(self.z_full[::R.DECIMATE, ::R.DECIMATE], X, Y, self.fire[0], self.fire[1])
            self.pl.add_mesh(self.anim.flame, scalars="rgba", rgba=True, lighting=False, name="flame")
            self.pl.add_mesh(self.anim.smoke, scalars="rgba", rgba=True, lighting=False, name="smoke")
        if pass_name != "fire":
            self.set_texture()

    def set_texture(self, red_field=None, with_char=True):
        if self.pass_name == "land":
            tex = land_texture(self.year)
        elif self.pass_name == "peat":
            tex = peat_rgba(self.year, red_field, with_char)
        else:
            tex = all_texture(self.year, self.fire, red_field, with_char)
        self.pl.add_mesh(self.grid, texture=tex, **TERRAIN_KW)

    def frame(self, g):
        R.set_camera(self.pl, R.camera_t_at(g))
        if self.anim is not None:
            self.anim.update(g)
        self.pl.render()
        return self.pl.screenshot(transparent_background=True, return_img=True, window_size=R.WINDOW)

    def close(self):
        self.pl.close()


def save_rgba(img, path):
    Image.fromarray(np.ascontiguousarray(img), "RGBA").save(path, compress_level=3)


# ----------------------------------------------------------------------------
# jobs
# ----------------------------------------------------------------------------
def job_year(year, outdir, passes, limit=None):
    """Act 2 year: pre-roll frames (at the global frames of the previous year's
    last dissolve frames) and the year's own frames, for each pass."""
    sch = R.schedule()
    s = sch[year]
    years = R.YEARS
    n_x = R.xfade_frames_for(year)
    if s["prev"] is not None:
        p = sch[s["prev"]]
        pre_g0 = p["tick0"] + p["n"] - n_x
    else:
        pre_g0 = s["tick0"] - n_x          # dissolve from act 1
    for pass_name in passes:
        d = os.path.join(outdir, pass_name)
        os.makedirs(d, exist_ok=True)
        sc = PassScene(year, pass_name, with_fire=True)
        t0 = time.time()
        for k in range(n_x):
            save_rgba(sc.frame(pre_g0 + k), os.path.join(d, "%d_pre_%02d.png" % (year, k)))
        n = s["n"] if limit is None else min(limit, s["n"])
        for i in range(n):
            save_rgba(sc.frame(s["tick0"] + i), os.path.join(d, "%d_%03d.png" % (year, i)))
        sc.close()
        print("%d %s: %d frames in %.1f s" % (year, pass_name, n_x + n, time.time() - t0), file=sys.stderr)


def job_act1(outdir, passes):
    """Act 1 (global frames 0 .. prefix - 1): red interpolated per frame, no
    char, no fire; the fire pass is fully transparent."""
    prefix = R.SPECS["twoact"]["prefix"]
    masks = {y: R.degraded_mask(y).astype(np.float32) for y in C.ACT1_KEYS}
    blank = np.zeros((R.WINDOW[1], R.WINDOW[0], 4), np.uint8)
    for pass_name in passes:
        d = os.path.join(outdir, pass_name)
        os.makedirs(d, exist_ok=True)
        if pass_name == "fire":
            for g in range(prefix):
                save_rgba(blank, os.path.join(d, "act1_%03d.png" % g))
            continue
        sc = PassScene(1990, pass_name, with_fire=False)
        t0 = time.time()
        g = 0
        for k, year in enumerate(C.ACT1_KEYS[:-1]):
            nxt = C.ACT1_KEYS[k + 1]
            n = C.ACT1_SEGMENT_FRAMES[year]
            ramp = C.ACT1_RAMP_2005 if year == 2005 else n
            for j in range(n):
                if pass_name != "land":
                    w = min(j / float(ramp), 1.0)
                    red = masks[year] + w * np.clip(masks[nxt] - masks[year], 0.0, 1.0)
                    sc.set_texture(red_field=red, with_char=False)
                save_rgba(sc.frame(g), os.path.join(d, "act1_%03d.png" % g))
                g += 1
        assert g == prefix, g
        sc.close()
        print("act1 %s: %d frames in %.1f s" % (pass_name, prefix, time.time() - t0), file=sys.stderr)


def run_jobs(outdir, passes, jobs, years):
    """outdir is the layer output directory (the children append frames/ themselves)."""
    base = [sys.executable, os.path.abspath(__file__), "--outdir", outdir, "--passes", ",".join(passes)]
    pending = ["act1"] + [str(y) for y in years]
    running, failed = [], []
    while pending or running:
        while pending and len(running) < jobs:
            j = pending.pop(0)
            log = open(os.path.join(outdir, "frames", "job_%s.log" % j), "w")
            running.append((j, subprocess.Popen(base + ["--job", j], stdout=log, stderr=subprocess.STDOUT), log))
        time.sleep(0.5)
        still = []
        for j, p, log in running:
            if p.poll() is None:
                still.append((j, p, log))
            else:
                log.close()
                print("job %s finished with code %d" % (j, p.returncode), file=sys.stderr)
                if p.returncode != 0:
                    failed.append(j)
        running = still
    return failed


# ----------------------------------------------------------------------------
# assembly
# ----------------------------------------------------------------------------
def blend_rgba(a_path, b_path, w, out):
    """Premultiplied dissolve of two RGBA frames."""
    if os.path.exists(out):
        return out
    a = np.asarray(Image.open(a_path).convert("RGBA"), np.float32)
    b = np.asarray(Image.open(b_path).convert("RGBA"), np.float32)
    pa, pb = a[..., :3] * a[..., 3:] / 255.0, b[..., :3] * b[..., 3:] / 255.0
    prem = (1 - w) * pa + w * pb
    alpha = (1 - w) * a[..., 3] + w * b[..., 3]
    rgb = np.where(alpha[..., None] > 0.5, prem * 255.0 / np.maximum(alpha, 0.5)[..., None], 0.0)
    out_img = np.concatenate([rgb, alpha[..., None]], axis=-1)
    Image.fromarray(np.clip(out_img, 0, 255).astype(np.uint8), "RGBA").save(out, compress_level=3)
    return out


def build_list(outdir, pass_name):
    d = os.path.join(outdir, pass_name)
    xdir = os.path.join(d, "xfade")
    os.makedirs(xdir, exist_ok=True)
    sch = R.schedule()
    years = R.YEARS
    prefix = R.SPECS["twoact"]["prefix"]
    entries, gs = [], []
    first = years[0]
    n_x0 = R.xfade_frames_for(first)
    for g in range(prefix):
        p = os.path.join(d, "act1_%03d.png" % g)
        k = g - (prefix - n_x0)
        if k >= 0:
            p = blend_rgba(p, os.path.join(d, "%d_pre_%02d.png" % (first, k)), (k + 1) / float(n_x0 + 1),
                           os.path.join(xdir, "x_act1_%d_%02d.png" % (first, k)))
        entries.append(p)
        gs.append(g)
    for idx, y in enumerate(years):
        s = sch[y]
        nxt = years[idx + 1] if idx + 1 < len(years) else None
        n_x = R.xfade_frames_for(nxt) if nxt is not None else 0
        for i in range(s["n"]):
            p = os.path.join(d, "%d_%03d.png" % (y, i))
            k = i - (s["n"] - n_x)
            if nxt is not None and k >= 0:
                p = blend_rgba(p, os.path.join(d, "%d_pre_%02d.png" % (nxt, k)), (k + 1) / float(n_x + 1),
                               os.path.join(xdir, "x_%d_%d_%02d.png" % (y, nxt, k)))
            entries.append(p)
            gs.append(s["tick0"] + i)
    assert gs == list(range(R.total_frames()))
    for p in entries:
        if not os.path.exists(p):
            raise SystemExit("missing frame: %s" % p)
    return entries, gs


def encode(entries, out_base):
    list_path = out_base + "_list.txt"
    with open(list_path, "w") as f:
        for p in entries:
            f.write("file '%s'\nduration %.6f\n" % (os.path.abspath(p), 1.0 / R.FPS))
        f.write("file '%s'\n" % os.path.abspath(entries[-1]))
    common = [FFMPEG, "-y", "-loglevel", "error", "-f", "concat", "-safe", "0", "-i", list_path,
              "-vf", "fps=%d,format=rgba" % R.FPS, "-r", str(R.FPS)]
    prores = out_base + "_prores4444.mov"
    subprocess.run(common + ["-c:v", "prores_ks", "-profile:v", "4444", "-pix_fmt", "yuva444p10le", prores], check=True)
    hevc = out_base + "_hevc_alpha.mov"
    subprocess.run(common + ["-c:v", "hevc_videotoolbox", "-alpha_quality", "0.9", "-tag:v", "hvc1",
                             "-pix_fmt", "bgra", "-b:v", "16M", hevc], check=True)
    return prores, hevc


def alpha_check(path):
    """Report the alpha statistics of a frame and write a checkerboard composite
    beside it (<name>_check.png) for a visual confirmation."""
    a = np.asarray(Image.open(path).convert("RGBA"))
    al = a[..., 3]
    stats = {"transparent": float((al == 0).mean()), "partial": float(((al > 0) & (al < 255)).mean()),
             "opaque": float((al == 255).mean())}
    H, W = al.shape
    chk = ((np.indices((H, W)).sum(0) // 40) % 2) * 70 + 110
    bg = Image.fromarray(np.stack([chk] * 3, -1).astype(np.uint8)).convert("RGBA")
    Image.alpha_composite(bg, Image.fromarray(a, "RGBA")).convert("RGB").save(path[:-4] + "_check.png")
    return stats


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=OUT_DEFAULT)
    ap.add_argument("--passes", default=",".join(PASSES))
    ap.add_argument("--jobs", type=int, default=7)
    ap.add_argument("--job", help="internal: act1 or a year")
    ap.add_argument("--test", action="store_true", help="2 s of 2015 for the passes, both codecs, alpha check")
    ap.add_argument("--skip-render", action="store_true")
    args = ap.parse_args(argv)
    passes = args.passes.split(",")
    R.set_spec("twoact")
    R.BG_OVERRIDE = "#000000"
    args.outdir = os.path.abspath(args.outdir)
    frames_dir = os.path.join(args.outdir, "frames")

    if args.job:
        os.makedirs(frames_dir, exist_ok=True)
        if args.job == "act1":
            job_act1(frames_dir, passes)
        else:
            job_year(int(args.job), frames_dir, passes)
        return

    if args.test:
        tdir = os.path.join(args.outdir, "test")
        os.makedirs(tdir, exist_ok=True)
        report = {}
        for pass_name in passes:
            d = os.path.join(tdir, pass_name)
            os.makedirs(d, exist_ok=True)
            sc = PassScene(2015, pass_name, with_fire=True)
            s = R.schedule()[2015]
            entries = []
            for i in range(2 * R.FPS):
                p = os.path.join(d, "f_%03d.png" % i)
                save_rgba(sc.frame(s["tick0"] + i), p)
                entries.append(p)
            sc.close()
            prores, hevc = encode(entries, os.path.join(tdir, pass_name))
            report[pass_name] = {"alpha": alpha_check(entries[0]), "prores": os.path.getsize(prores),
                                 "hevc": os.path.getsize(hevc), "check_png": entries[0][:-4] + "_check.png"}
        print(json.dumps(report, indent=1))
        return

    os.makedirs(frames_dir, exist_ok=True)
    if not args.skip_render:
        failed = run_jobs(args.outdir, passes, args.jobs, R.YEARS)
        if failed:
            raise SystemExit("failed jobs: %s" % failed)
    report = {}
    for pass_name in passes:
        entries, gs = build_list(frames_dir, pass_name)
        max_step, max_ratio = R.check_azimuth_steps(gs)
        prores, hevc = encode(entries, os.path.join(args.outdir, pass_name))
        last = os.path.join(args.outdir, "%s_2026_last.png" % pass_name)
        shutil.copyfile(entries[-1], last)
        report[pass_name] = {"frames": len(entries), "seconds": len(entries) / R.FPS, "max_azimuth_step_deg": max_step,
                             "prores4444": prores, "prores_bytes": os.path.getsize(prores),
                             "hevc_alpha": hevc, "hevc_bytes": os.path.getsize(hevc), "last_png": last,
                             "last_alpha": alpha_check(last)}
    with open(os.path.join(args.outdir, "report.json"), "w") as f:
        json.dump(report, f, indent=1)
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
