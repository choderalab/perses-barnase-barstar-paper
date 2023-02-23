#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 3:00 
#
# Set output file
#BSUB -o run_equil_apo.out
#
# Set error file
#BSUB -eo run_equil_apo.stderr
#
# Specify node group
#BSUB -m "ly-gpu lx-gpu lu-gpu lt-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=3]" -sp 25
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "equil.34.apo"

outdir=/data/chodera/zhangi/perses_benchmark/repex/43/34/
phase=apo

source ~/.bashrc
conda activate perses-paper3

cd /data/chodera/zhangi/perses_benchmark/perses-barnase-barstar-paper/scripts/02_run_equilibration/
python run_equilibration.py $outdir $phase --gentle