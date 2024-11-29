#!/bin/bash
#SBATCH --job-name=gmgpolar
#SBATCH --output=../../slurm_output/slurm-%A-Strong-Scaling-Benchmark.out
#SBATCH --error=../../slurm_output/slurm-%A-Strong-Scaling-Benchmark.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 56
#SBATCH --threads-per-core=1
#SBATCH -t 1400
# #SBATCH --nodelist="be-cpu02"
#SBATCH --exclusive

# Adjust parameters in src/strong_scaling.cpp

srun --cpus-per-task=56 ./../../../build/strong_scaling
