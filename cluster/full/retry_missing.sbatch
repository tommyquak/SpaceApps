#!/bin/bash
#SBATCH --job-name=retrypair
#SBATCH --partition=normal
#SBATCH --cpus-per-task=24
#SBATCH --mem=60G
#SBATCH --tmp=400000
#SBATCH --time=8:00:00
#SBATCH --output=/home/t/tommy/full/logs/retry-%j.log
# Insurance for the Zheng et al. 2022 closure-bias correction.
#
# A complete bandwidth-20 network over 93 dates is exactly 1650 pairs, which is
# exactly what we plan: zero redundancy. mintpy2.sbatch's gate requires the
# sequential network complete through connection level 20 and deliberately
# refuses to lower --nl, so ONE missing pair costs the whole correction. That is
# precisely how it was lost last time (19 unwrap failures broke levels 1-5).
#
# This regenerates only the pairs that are missing from pairs_tay, reusing
# main.sbatch's own per-pair machinery via a one-block override.
set -uo pipefail
OUT=${OUT:-$HOME/full/pairs_tay}
python3 - "$OUT" > $HOME/full/retry.pairs <<'PY'
import os, sys, glob, re
root = sys.argv[1]
PAT = re.compile(r"^(\d{4}-\d{2}-\d{2})_(\d{4}-\d{2}-\d{2})$")
done = set()
for d in glob.glob(root + "/*/"):
    m = PAT.match(os.path.basename(d.rstrip("/")))
    if m and os.path.exists(d + "DONE"):
        done.add((m.group(1).replace("-", "")[2:], m.group(2).replace("-", "")[2:]))
planned = []
for f in sorted(glob.glob(os.path.expanduser("~/full/main_blocks/*.pairs"))):
    for ln in open(f):
        p = ln.strip()
        if "-" in p:
            planned.append(tuple(p.split("-")))
missing = [p for p in planned if p not in done]
for a, b in missing:
    print("%s-%s" % (a, b))
sys.stderr.write("planned %d, done %d, missing %d\n" % (len(planned), len(done), len(missing)))
PY
N=$(wc -l < $HOME/full/retry.pairs)
echo "RETRY missing pairs: $N"
[ "$N" -eq 0 ] && { echo "RETRY_NOTHING_TO_DO network already complete"; exit 0; }
[ "$N" -gt 200 ] && { echo "RETRY_ABORT $N missing is too many for a retry; investigate the run instead"; exit 1; }

mkdir -p $HOME/full/main_blocks_retry
cp $HOME/full/retry.pairs $HOME/full/main_blocks_retry/0.pairs
awk -F- '{print $1"\n"$2}' $HOME/full/retry.pairs | sort -u > $HOME/full/main_blocks_retry/0.dates
echo "retry block: $(wc -l < $HOME/full/main_blocks_retry/0.pairs) pairs over $(wc -l < $HOME/full/main_blocks_retry/0.dates) dates"
echo "now submit: sbatch --array=0-0 --export=ALL,NRLKS2=17,NALKS2=4,USE_ION=1,\\"
echo "  ION_DIR=/mnt/scratch/t/tommy/full/dates_ion_v3,OUT=$OUT,BLOCKS=\$HOME/full/main_blocks_retry main.sbatch"
