#!/bin/bash
#SBATCH --job-name=alos2-verify
#SBATCH --partition=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/home/t/tommy/pilot/logs/verify.log
cd ~/pilot
. ~/pilot/lib_checks.sh
total=$(wc -l < urls.txt)
ok=0
while read -r date_tag _; do
  if scene_complete "ceos/$date_tag"; then
    ok=$((ok+1))
  else
    echo "$date_tag INCOMPLETE"
  fi
done < urls.txt
echo "=== $ok/$total scenes verified ==="
du -sh ~/pilot/ceos ~/pilot/dem 2>/dev/null
if [ "$ok" -eq "$total" ]; then
  touch ~/pilot/READY_FOR_ALOSSTACK; echo PILOT_DATA_READY
else
  r=$(cat ~/pilot/.dl_retries 2>/dev/null || echo 0)
  if [ "$r" -lt 2 ]; then
    echo "SELF-HEAL: resubmitting download array (retry $((r+1)))"
    a=$(sbatch --parsable ~/pilot/download.sbatch)
    [ -n "$a" ] || { echo "SELF-HEAL-FAILED: download array submission returned no job id"; exit 1; }
    b=$(sbatch --parsable --dependency=afterany:$a ~/pilot/verify.sbatch)
    [ -n "$b" ] || { echo "SELF-HEAL-FAILED: verify resubmission returned no job id"; exit 1; }
    echo $((r+1)) > ~/pilot/.dl_retries
  else
    echo "GAVE UP after 2 retries, needs human"
    exit 1
  fi
fi
