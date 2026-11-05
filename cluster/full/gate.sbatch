#!/bin/bash
#SBATCH --job-name=full-gate
#SBATCH --partition=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --exclude=xcng[0-1],xgpj0,xcni[0-19],xcnc15,xcnc18
#SBATCH --output=/home/t/tommy/full/logs/gate.log
cd ~/full
nd=$(ls pairs/*/DONE 2>/dev/null | wc -l)
echo "test pack DONE count: $nd/4"
if [ "$nd" -ge 3 ]; then
  sbatch controller.sbatch && echo "GATE_PASSED controller launched"
else
  echo "GATE_FAILED only $nd test pairs done, campaign held for review"
fi
