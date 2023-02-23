#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 1:00 
#
# Set output file
#BSUB -o  generate_solvated_inputs.out
#
# Set error file
#BSUB -eo generate_solvated_inputs.stderr
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
#BSUB -J "inputs.25"

protein_filename=/data/chodera/zhangi/perses_benchmark/perses-barnase-barstar-paper/input_files/1brs_barstar_F29.pdb
outdir=/data/chodera/zhangi/perses_benchmark/repex/43/25/

source ~/.bashrc
conda activate perses-paper3

cd /data/chodera/zhangi/perses_benchmark/perses-barnase-barstar-paper/scripts/01_generate_solvated_inputs/
python generate_solvated_inputs.py $protein_filename $outdir --ligand_input /data/chodera/zhangi/perses_benchmark/perses-barnase-barstar-paper/input_files/1brs_barnase_renumbered.pdb
