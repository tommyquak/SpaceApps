#!/bin/bash
#SBATCH --job-name=alos2-shepherd
#SBATCH --partition=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --output=/home/t/tommy/pilot/logs/shepherd.log
RANGE=${1:-}
NEXT=${2:-}
[ -n "$RANGE" ] || { echo "no range"; exit 1; }
a=$(sbatch --parsable --array=$RANGE%10 ~/pilot/fanout.sbatch)
[ -n "$a" ] || { echo "SELF-HEAL-FAILED: fanout slice $RANGE submission returned no job id"; exit 1; }
echo "submitted fanout slice $RANGE as $a"
if [ -n "$NEXT" ]; then
  n3=""
  [ "$NEXT" = "25-48" ] && n3="49-70"
  b=$(sbatch --parsable --dependency=afterany:$a ~/pilot/shepherd.sbatch "$NEXT" "$n3")
  [ -n "$b" ] || { echo "SELF-HEAL-FAILED: shepherd chain for $NEXT returned no job id"; exit 1; }
  echo "chained shepherd for $NEXT"
fi
