#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 72:00
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
#BSUB -R "rusage[mem=3]"
#BSUB -n 2 -gpu "num=1/task:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "49.12.com"

source ~/.bashrc
conda activate perses-paper5

outdir=/data/chodera/zhangi/perses_benchmark/repex/perses-bnbs-paper-fourth-attempt/49/12/replicate_0/
phase=complex
n_states=36
n_cycles=10000
t_max=400
restraint=None
force_constant=None

build_mpirun_configfile --configfilepath configfile_${phase} --hostfilepath hostfile_${phase} "python /data/chodera/zhangi/perses_benchmark/perses-barnase-barstar-paper/scripts/04_run_repex/run_repex.py $outdir $phase $n_states $n_cycles $t_max"
mpiexec.hydra -f hostfile_${phase} -configfile configfile_${phase}
