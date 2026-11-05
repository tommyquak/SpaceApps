#!/bin/bash
#SBATCH --job-name=zhengbw4
#SBATCH --output=/home/t/tommy/full/logs/zhengbw4_%j.out
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
# Zheng et al. 2022 closure-phase bias for a BANDWIDTH-4 network, estimated
# from the full 1650-pair stack (Tay's remedy for short-baseline rates), written
# to its own directory so the bandwidth-20 products are untouched.
export PATH=$HOME/micromamba/envs/insar/bin:$PATH PROJ_LIB=$HOME/micromamba/envs/insar/share/proj TMPDIR=$HOME/full/tmp PYTHONUNBUFFERED=1
W=$HOME/full/mintpy_tay; O=$W/zheng_bw4; mkdir -p $O; cd $W || exit 1
closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 20 --bw 4 -a mask -o $O; echo "ZBW4_MASK_RC=$?"
closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 20 --bw 4 -a quick_estimate -c local --num-worker 16 -o $O; echo "ZBW4_QUICK_RC=$?"
closure_phase_bias.py -i inputs/ifgramStack.h5 --nl 20 --bw 4 -a estimate -c local --num-worker 16 -o $O; rc=$?; echo "ZBW4_EST_RC=$rc"
ls -la $O/*.h5 2>/dev/null; echo ZBW4_DONE
