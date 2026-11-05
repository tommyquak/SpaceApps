# Shared scene-completeness predicate, sourced by download.sbatch (skip check)
# and verify.sbatch (readiness gate) so both key on the same condition. A scene
# is usable only with the full CEOS member set alos2App needs for swaths 2-4:
# volume, leader and trailer files plus the HH images for F2, F3 and F4.
scene_complete() {
  local d=$1
  [ -d "$d" ] || return 1
  ls "$d"/VOL-* >/dev/null 2>&1 &&
  ls "$d"/LED-* >/dev/null 2>&1 &&
  ls "$d"/TRL-* >/dev/null 2>&1 &&
  ls "$d"/IMG-HH-*-F2 >/dev/null 2>&1 &&
  ls "$d"/IMG-HH-*-F3 >/dev/null 2>&1 &&
  ls "$d"/IMG-HH-*-F4 >/dev/null 2>&1
}
