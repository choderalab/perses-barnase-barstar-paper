#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 2:00
#
# Set output file
#BSUB -o run_make_traj_per_replica.%I.out
#
# Set error file
#BSUB -eo run_make_traj_per_replica.%I.stderr
#
# Specify node group
#BSUB -m "lj-gpu ll-gpu ln-gpu ly-gpu lx-gpu lu-gpu ld-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -R "rusage[mem=60]" -sp 25
#BSUB -n 1 -gpu "num=1/task:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "47.26[1-36]"

source ~/.bashrc
conda activate perses-paper5

sub_dir=26
out_dir=/data/chodera/zhangi/perses_benchmark/repex/perses-bnbs-paper-fourth-attempt/47/26/replicate_1
phase=complex

python /data/chodera/zhangi/perses_benchmark/perses-barnase-barstar-paper/scripts/05_analyze/run_make_traj_per_replica.py $out_dir $sub_dir $phase "$((${LSB_JOBINDEX}-1))"
