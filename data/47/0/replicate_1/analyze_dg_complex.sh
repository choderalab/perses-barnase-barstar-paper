#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 1:00 
#
# Set output file
#BSUB -o analyze_dg_complex.out
#
# Set error file
#BSUB -eo analyze_dg_complex.stderr
#
# Specify node group
#BSUB -q cpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -R "rusage[mem=3]"
#BSUB -n 24
#
# job name (default = name of script file)
#BSUB -J "47.0a"

source ~/.bashrc
conda activate perses-paper5

main_dir=47
sub_dir=0
replicate=1
phase=complex
total_iterations=50000

python /data/chodera/zhangi/perses_benchmark/perses-barnase-barstar-paper/scripts/05_analyze/analyze_dg.py $main_dir $sub_dir $replicate $phase $total_iterations --bootstrap

