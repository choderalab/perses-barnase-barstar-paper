#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 24:00
#
# Set output file
#BSUB -o run_repex_complex.out
#
# Set error file
#BSUB -eo run_repex_complex.stderr
#
# Specify node group
#BSUB -m "lj-gpu ll-gpu ln-gpu ly-gpu lx-gpu lu-gpu ld-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -R "rusage[mem=3]" -sp 25
#BSUB -n 2 -gpu "num=1/task:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "45.10.comr2"

source ~/.bashrc
conda activate perses-paper5

outdir=/data/chodera/zhangi/perses_benchmark/repex/perses-bnbs-paper-fourth-attempt/45/10/replicate_1/
phase=complex
n_states=36
n_cycles=50000
t_max=300

build_mpirun_configfile --configfilepath configfile_${phase} --hostfilepath hostfile_${phase} "python run_repex.py $outdir $phase $n_states $n_cycles $t_max"
mpiexec.hydra -f hostfile_${phase} -configfile configfile_${phase}
