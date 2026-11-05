#!/bin/bash
#SBATCH --job-name=strip-unw
#SBATCH --partition=normal
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --output=/home/t/tommy/full/logs/strip_unw-%A_%a.log
set -uo pipefail
export PATH=$HOME/micromamba/envs/insar/bin:$HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications:$PATH
export PROJ_DATA=$HOME/micromamba/envs/insar/share/proj ISCE_STACK=$HOME/micromamba/envs/insar/share/isce2
export PATH=$ISCE_STACK/alosStack:$PATH PYTHONPATH=$HOME/full/strip_test:$ISCE_STACK/alosStack:${PYTHONPATH:-}
cd $HOME/full/strip_test
export OMP_NUM_THREADS=1
B=${SLURM_ARRAY_TASK_ID:-0}; NB=${NBLOCKS:-8}
python -c "import shutil,sys; sys.exit(shutil.which('snaphu') is None)" || { echo "STRIP_UNW_ABORT no snaphu"; exit 1; }
PAIRS=$(python -c "from strip_common import pair_list; print(' '.join(pair_list()))" | tr ' ' '\n' | awk -v b=$B -v nb=$NB 'NR % nb == b')
echo "block $B: $(echo $PAIRS | wc -w) pairs"
printf '%s\n' $PAIRS | xargs -P 12 -I{} python strip_unwrap.py {} | tee /tmp/strip_unw_$B.out
ok=$(grep -c "^UNW_OK\|^SKIP" /tmp/strip_unw_$B.out); bad=$(grep -c "SNAPHU_FAIL\|Traceback" /tmp/strip_unw_$B.out)
echo "block $B: ok=$ok failed=$bad"; [ "$bad" -eq 0 ] && echo "STRIP_UNW_BLOCK_OK $B" || { echo "STRIP_UNW_BLOCK_FAIL $B"; exit 1; }
