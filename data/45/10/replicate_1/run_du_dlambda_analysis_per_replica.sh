#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 3:00
#
# Set output file
#BSUB -o run_du_dlambda_analysis_replica.%I.out
#
# Set error file
#BSUB -eo run_du_dlambda_analysis_replica.%I.stderr
#
# Specify node group
#BSUB -m "lj-gpu ll-gpu ln-gpu ly-gpu lx-gpu lu-gpu ld-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -R "rusage[mem=20]" -sp 25
#BSUB -n 1 -gpu "num=1/task:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "45.10.comr1du[1-36]"

source ~/.bashrc
conda activate perses-paper5

sub_dir=10
out_dir=/data/chodera/zhangi/perses_benchmark/repex/perses-bnbs-paper-fourth-attempt/45/$sub_dir/replicate_1/
phase=complex

python /home/zhangi/choderalab/perses_benchmark/perses_protein_mutations/code/39_perses_paper/run_du_dlambda_analysis_per_replica.py $out_dir $sub_dir $phase "$((${LSB_JOBINDEX}-1))"
