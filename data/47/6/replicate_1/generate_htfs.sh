#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 2:00 
#
# Set output file
#BSUB -o  generate_htfs.out
#
# Set error file
#BSUB -eo generate_htfs.stderr
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
#BSUB -J "htf.6"

protein_filename=/data/chodera/zhangi/perses_benchmark/repex/43/3/apo_equilibrated.pdb
ligand_input_filename=/data/chodera/zhangi/perses_benchmark/repex/43/3/complex_equilibrated.pdb
resid=80
mutant_aa=ALA
outdir=/data/chodera/zhangi/perses_benchmark/repex/47/6/replicate_1/

source ~/.bashrc
conda activate perses-paper3

cd /data/chodera/zhangi/perses_benchmark/perses-barnase-barstar-paper/scripts/03_generate_htfs/
python generate_htfs.py $protein_filename $resid $mutant_aa $outdir --ligand_input $ligand_input_filename
