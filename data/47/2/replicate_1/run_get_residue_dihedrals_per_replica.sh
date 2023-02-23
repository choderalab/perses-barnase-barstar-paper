#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 1:00
#
# Set output file
#BSUB -o run_get_residue_dihedrals_per_replica.%I.out
#
# Set error file
#BSUB -eo run_get_residue_dihedrals_per_replica.%I.stderr
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
#BSUB -J "47.2d"

source ~/.bashrc
conda activate perses-paper5

sub_dir=2
out_dir=/data/chodera/zhangi/perses_benchmark/repex/perses-bnbs-paper-fourth-attempt/47/2/replicate_1
phase=complex

python /data/chodera/zhangi/perses_benchmark/perses-barnase-barstar-paper/scripts/05_analyze/get_residue_dihedrals_per_replica.py $out_dir $sub_dir $phase
