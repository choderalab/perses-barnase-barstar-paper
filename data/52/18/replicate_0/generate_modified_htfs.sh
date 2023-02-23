#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 2:00 
#
# Set output file
#BSUB -o  generate_modified_htfs.out
#
# Set error file
#BSUB -eo generate_modified_htfs.stderr
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
#BSUB -J "htf.18"

in_dir=/data/chodera/zhangi/perses_benchmark/repex/perses-bnbs-paper-fourth-attempt/47/18/replicate_1/
out_dir=/data/chodera/zhangi/perses_benchmark/repex/perses-bnbs-paper-fourth-attempt/52/18/replicate_0/

source ~/.bashrc
conda activate perses-paper5

cd /data/chodera/zhangi/perses_benchmark/perses-barnase-barstar-paper/scripts/03_generate_htfs/
python generate_modified_htfs.py $in_dir $out_dir --rest_radius 0.5