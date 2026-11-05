#!/bin/bash
#SBATCH --job-name=full-ctl
#SBATCH --partition=normal
#SBATCH --dependency=singleton
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:50:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --output=/home/t/tommy/full/logs/controller.log
#SBATCH --open-mode=append
# Rolling-window campaign controller (v2, post ultra-audit).
# Filesystem is the only state:
#   ceos/<date>/.complete      scene fully extracted (atomic: dir appears by rename)
#   ceos/<date>.part/          in-flight extraction (counted in window; reaped if stale)
#   pairs/<R>_<S>/DONE|FAILED  terminal markers   ATTEMPTS  science-attempt counter
#   pairs/<R>_<S>/QUEUED       jobid of in-flight pack   MISSING1  one-tick vanish grace
#   CAMPAIGN_DONE              mintpy submitted; controllers retire
# Resubmitting this script is always safe (singleton).
set -uo pipefail
cd ~/full
PACK=4; MAXJOBS=26; MAXSCENES=23; MAXDL=4
[ -f CAMPAIGN_DONE ] && { echo "$(date '+%F %T') campaign done, retiring"; exit 0; }
# Successor-at-top: the pending singleton successor survives any death of this
# generation. On retirement the successor exits at the CAMPAIGN_DONE check above.
SUCC=0
for k in 1 2 3; do sbatch controller.sbatch >/dev/null 2>&1 && { SUCC=1; break; }; sleep 30; done
[ $SUCC -eq 1 ] || echo "$(date '+%F %T') WARN successor submit failed, retrying in loop"
mapfile -t SC < <(awk '{print $1}' scenes_full.txt)
N=${#SC[@]}
t_end=$(( $(date +%s) + 9600 ))
while [ $(date +%s) -lt $t_end ]; do
  [ -f CAMPAIGN_DONE ] && { echo "$(date '+%F %T') campaign done, retiring"; exit 0; }
  [ $SUCC -eq 1 ] || { sbatch controller.sbatch >/dev/null 2>&1 && SUCC=1 || true; }
  # ---- queue snapshot: fail closed; a valid snapshot always contains this job ----
  squeue -u $USER -h -o '%i %j %K' > /tmp/q.$$ 2>>logs/controller.log \
    && grep -q "^${SLURM_JOB_ID} " /tmp/q.$$ \
    || { echo "$(date '+%F %T') squeue snapshot invalid, skipping tick"; rm -f /tmp/q.$$; sleep 240; continue; }
  njobs=$(( $(grep -c 'full-dl' /tmp/q.$$ || true) + $(grep -c 'full-pairs' /tmp/q.$$ || true) ))
  ndl=$(grep -c 'full-dl' /tmp/q.$$ || true)

  # ---- pass 1: pair states (two-tick vanish confirmation before charging an attempt) ----
  declare -A SCENE_OPEN=()
  READY=(); NDONE=0; NFAIL=0
  ln=0
  while read -r REF SEC; do
    ln=$((ln+1)); P=pairs/${REF}_${SEC}
    [ -f $P/DONE ] && { NDONE=$((NDONE+1)); continue; }
    [ -f $P/FAILED ] && { NFAIL=$((NFAIL+1)); continue; }
    SCENE_OPEN[$REF]=1; SCENE_OPEN[$SEC]=1
    if [ -f $P/QUEUED ]; then
      jid=$(cat $P/QUEUED)
      if grep -q "^$jid " /tmp/q.$$; then rm -f $P/MISSING1; scontrol release "$jid" 2>/dev/null || true; continue; fi
      if [ ! -f $P/MISSING1 ]; then touch $P/MISSING1; continue; fi
      scancel "$jid" 2>/dev/null || true
      rm -f $P/MISSING1 $P/QUEUED
      if [ ! -f logs/pack-$jid.log ]; then
        # job vanished without ever starting (launch failure, black-hole node):
        # release without charging a science attempt; READY path resubmits.
        echo "$(date '+%F %T') LAUNCHFAIL job $jid for ${REF}_${SEC}, released without attempt"
      else
      att=$(( $(cat $P/ATTEMPTS 2>/dev/null || echo 0) + 1 ))
      echo $att > $P/ATTEMPTS
      [ $att -ge 2 ] && { touch $P/FAILED; echo "$(date '+%F %T') PERMFAIL ${REF}_${SEC} after $att attempts"; NFAIL=$((NFAIL+1)); continue; }
      fi
    fi
    [ -f ceos/$REF/.complete ] && [ -f ceos/$SEC/.complete ] && READY+=("$ln")
  done < pairs_full.txt
  NPAIRS=$ln

  # ---- pass 2: reap stale partials, release scenes with no open pairs ----
  for dp in ceos/*.part; do
    [ -d "$dp" ] || continue
    [ -z "$(find "$dp" -mmin -120 -print -quit 2>/dev/null)" ] && { echo "$(date '+%F %T') REAP stale partial $dp"; rm -rf "$dp"; }
  done
  for d in $(ls -d ceos/*/ 2>/dev/null | xargs -n1 basename); do
    case "$d" in *.part) continue;; esac
    [ -z "${SCENE_OPEN[$d]:-}" ] && { echo "$(date '+%F %T') RELEASE scene $d"; rm -rf ceos/$d; }
  done

  # ---- pass 3: frontier downloads; window = landed dirs (incl. partials) + in-flight ----
  ndirs=$(ls -d ceos/*/ 2>/dev/null | wc -l)
  nsc=$(( ndirs + ndl ))
  if [ $nsc -lt $MAXSCENES ] && [ $ndl -lt $MAXDL ]; then
    for i in $(seq 0 $((N-1))); do
      d=${SC[$i]}
      [ -z "${SCENE_OPEN[$d]:-}" ] && continue
      [ -f ceos/$d/.complete ] && continue
      grep -q "full-dl $((i+1))\$" /tmp/q.$$ && continue
      sbatch --array=$((i+1)) dl2.sbatch >/dev/null 2>&1 && { echo "$(date '+%F %T') DL scene $d (idx $((i+1)))"; ndl=$((ndl+1)); nsc=$((nsc+1)); }
      { [ $ndl -ge $MAXDL ] || [ $nsc -ge $MAXSCENES ]; } && break
    done
  fi

  # ---- pass 4: packed pair submission ----
  while [ ${#READY[@]} -ge 1 ] && [ $njobs -lt $MAXJOBS ]; do
    CHUNK=("${READY[@]:0:$PACK}"); READY=("${READY[@]:${#CHUNK[@]}}")
    jid=$(sbatch --parsable -H fanout2.sbatch "${CHUNK[@]}" 2>>logs/controller.log) || break
    for ln2 in "${CHUNK[@]}"; do
      line=$(sed -n "${ln2}p" pairs_full.txt); P=pairs/${line%% *}_${line##* }
      mkdir -p $P; echo $jid > $P/QUEUED; rm -f $P/MISSING1
    done
    scontrol release "$jid" 2>>logs/controller.log || echo "$(date +%F_%T) WARN release failed for $jid, pass 1 heals"
    echo "$(date '+%F %T') PACK job $jid pairs: ${CHUNK[*]}"
    njobs=$((njobs+1))
  done

  # ---- pass 4b: orphan-held reaper (audit r4): a -H pack whose jid never
  # reached a QUEUED marker would sit JobHeldUser forever, eating a job slot ----
  declare -A KNOWN=()
  for qf in pairs/*/QUEUED; do [ -f "$qf" ] && KNOWN[$(cat $qf)]=1; done
  mkdir -p orphan
  squeue -h -u $USER -n full-pairs -t PD -O 'JobID,Reason' 2>/dev/null > /tmp/pd.$$
  while read -r oj reason; do
    [ "$reason" = "JobHeldUser" ] || continue
    [ -n "${KNOWN[$oj]:-}" ] && { rm -f orphan/$oj; continue; }
    if [ -f orphan/$oj ]; then
      scancel $oj 2>/dev/null || true; rm -f orphan/$oj
      echo "$(date '+%F %T') ORPHAN_HELD reaped $oj"
    else
      touch orphan/$oj
    fi
  done < /tmp/pd.$$
  rm -f /tmp/pd.$$
  for om in orphan/*; do
    [ -f "$om" ] || continue
    grep -q "^$(basename $om) " /tmp/q.$$ || rm -f "$om"
  done

  # ---- pass 5: finale (retry-safe, logged) ----
  NTERM=$((NDONE+NFAIL))
  echo "$(date '+%F %T') STATUS done=$NDONE fail=$NFAIL open=$((NPAIRS-NTERM)) scenes=$nsc jobs=$njobs"
  if [ $NTERM -ge $NPAIRS ]; then
    if sbatch mintpy2.sbatch >>logs/controller.log 2>&1; then
      scancel --name=full-dl -u $USER 2>/dev/null || true
      rm -rf ceos; mkdir -p ceos
      touch CAMPAIGN_DONE; echo "$(date '+%F %T') ALL PAIRS TERMINAL, scene store cleared, mintpy2 submitted"; exit 0
    else
      echo "$(date '+%F %T') mintpy2 sbatch FAILED, retrying next tick"
    fi
  fi
  rm -f /tmp/q.$$
  sleep 240
done
if [ $SUCC -ne 1 ]; then
  for k in 1 2 3 4 5; do sbatch controller.sbatch >/dev/null 2>&1 && { SUCC=1; break; }; sleep 60; done
  [ $SUCC -eq 1 ] || echo "$(date '+%F %T') CHAIN_AT_RISK no successor could be queued"
fi
echo "$(date '+%F %T') generation end, successor pending"
