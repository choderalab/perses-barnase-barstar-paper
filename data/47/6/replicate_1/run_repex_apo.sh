#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 30:00
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
#BSUB -J "47.6.apor1"

source ~/.bashrc
conda activate perses-paper5

outdir=/data/chodera/zhangi/perses_benchmark/repex/perses-bnbs-paper-fourth-attempt/47/6/replicate_1/
phase=apo
n_states=36
n_cycles=10000
t_max=300
restraint=None
force_constant=None

build_mpirun_configfile --configfilepath configfile_${phase} --hostfilepath hostfile_${phase} "python /data/chodera/zhangi/perses_benchmark/perses-barnase-barstar-paper/scripts/04_run_repex/run_repex.py $outdir $phase $n_states $n_cycles $t_max"
mpiexec.hydra -f hostfile_${phase} -configfile configfile_${phase}
