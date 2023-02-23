#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 24:00 
#
# Set output file
#BSUB -o analyze_dg_timeseries_complex.out
#
# Set error file
#BSUB -eo analyze_dg_timeseries_complex.stderr
#
# Specify node group
#BSUB -m "lj-gpu ll-gpu ln-gpu ly-gpu lx-gpu lu-gpu lt-gpu"
#BSUB -q cpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -R "rusage[mem=3]"
#BSUB -n 24
#
# job name (default = name of script file)
#BSUB -J "52.25cr0a"

source ~/.bashrc
conda activate perses-paper5

main_dir=52
sub_dir=25
replicate=0
phase=complex
starting_iteration=1000
total_iterations=50000

python /data/chodera/zhangi/perses_benchmark/perses-barnase-barstar-paper/scripts/05_analyze/analyze_dg_timeseries.py $main_dir $sub_dir $replicate $phase $starting_iteration $total_iterations
