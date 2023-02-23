#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 2:00
#
# Set output file
#BSUB -o run_repex_apo.out
#
# Set error file
#BSUB -eo run_repex_apo.stderr
#
# Specify node group
#BSUB -m "lj-gpu ll-gpu ln-gpu ly-gpu lx-gpu lu-gpu ld-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -R "rusage[mem=2]"
#BSUB -n 2 -gpu "num=1/task:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "46.20.apo"

source ~/.bashrc
conda activate perses-paper3

outdir=/data/chodera/zhangi/perses_benchmark/repex/46/20/replicate_0/
phase=apo
n_states=12
n_cycles=5000
t_max=300

build_mpirun_configfile --configfilepath configfile_${phase} --hostfilepath hostfile_${phase} "python /data/chodera/zhangi/perses_benchmark/perses-barnase-barstar-paper/scripts/04_run_repex/run_repex.py $outdir $phase $n_states $n_cycles $t_max"
mpiexec.hydra -f hostfile_${phase} -configfile configfile_${phase}
