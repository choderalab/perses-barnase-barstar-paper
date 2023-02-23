#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 1:00
#
# Set output file
#BSUB -o generate_heatmap_data_50ns.out
#
# Set error file
#BSUB -eo generate_heatmap_data_50ns.stderr
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
#BSUB -J "47.25"

source ~/.bashrc
conda activate perses-paper5

main_dir=47
sub_dir=25
phase=complex

python /data/chodera/zhangi/perses_benchmark/perses-barnase-barstar-paper/scripts/05_analyze/generate_heatmap_data_50ns.py $main_dir $sub_dir $phase
