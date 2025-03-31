#!/bin/bash
#SBATCH --job-name=gmgpolar
#SBATCH --output=../../slurm_output/slurm-%A-Weak-Scaling-Benchmark.out
#SBATCH --error=../../slurm_output/slurm-%A-Weak-Scaling-Benchmark.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 56
#SBATCH --threads-per-core=1
#SBATCH -t 1400
# #SBATCH --nodelist="be-cpu02"
#SBATCH --exclusive

# Adjust parameters in src/weak_scaling.cpp

srun --cpus-per-task=56 ./../../../build/weak_scaling
