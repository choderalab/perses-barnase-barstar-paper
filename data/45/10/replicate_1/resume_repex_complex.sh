#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 36:00
#
# Set output file
#BSUB -o resume_repex_complex.out
#
# Set error file
#BSUB -eo resume_repex_complex.stderr
#
# Specify node group
#BSUB -m "lj-gpu ll-gpu ln-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -R "rusage[mem=3]" -sp 25
#BSUB -n 2 -gpu "num=1/task:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "45.10.comr1"
#BSUB -w "ended(9369370)"

source ~/.bashrc
conda activate perses-paper5

outdir=/data/chodera/zhangi/perses_benchmark/repex/perses-bnbs-paper-fourth-attempt/45/10/replicate_1/
phase=complex
total_iterations=100000

build_mpirun_configfile --configfilepath configfile_${phase} --hostfilepath hostfile_${phase} "python /data/chodera/zhangi/perses_benchmark/perses-barnase-barstar-paper/scripts/04_run_repex/resume_repex.py $outdir $phase $total_iterations"
mpiexec.hydra -f hostfile_${phase} -configfile configfile_${phase}
