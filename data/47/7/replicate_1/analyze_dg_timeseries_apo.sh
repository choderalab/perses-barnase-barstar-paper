#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 1:00 
#
# Set output file
#BSUB -o analyze_dg_timeseries_apo.out
#
# Set error file
#BSUB -eo analyze_dg_timeseries_apo.stderr
#
# Specify node group
#BSUB -q cpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -R "rusage[mem=3]"
#BSUB -n 24
#
# job name (default = name of script file)
#BSUB -J "47.7ar1a"

source ~/.bashrc
conda activate perses-paper5

main_dir=47
sub_dir=7
replicate=1
phase=apo
starting_iteration=1000
total_iterations=10000

python /data/chodera/zhangi/perses_benchmark/perses-barnase-barstar-paper/scripts/05_analyze/analyze_dg_timeseries.py $main_dir $sub_dir $replicate $phase $starting_iteration $total_iterations
