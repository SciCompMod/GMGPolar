#!/bin/bash
#SBATCH --job-name=RoofLine
#SBATCH --output=slurm-%A-Roofline-Model.out
#SBATCH --error=slurm-%A-Roofline-Model.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH --threads-per-core=1
#SBATCH --time=0:30:00
#SBATCH --exclusive
#SBATCH --partition=naples128
#SBATCH --account=2476029

srun likwid-topology

# CPU name: AMD EPYC 7601 32-Core Processor      
# CPU type: AMD K17 (Zen) architecture
# CPU stepping: 2
# Sockets: 2
# Cores per socket: 32
# Threads per core: 2
# Level 1: 32 kB
# Level 2: 512 kB
# Level 3: 8 MB
# NUMA domains: 8

# So we have 2 x 32 hardware threads on this machine. 
# We use a stream size of 2 x 32 x 32 kB = 2048 kB, 
# 32kB for each of the 64 hardware threads so that each vector chunk fits into the L1 cache of one core.

# Upper Peak GFLOPS Estimate: 2.2 GHz x 32 FLOPS/cycle x 64 cores = 4505.6 GFLOPS
srun likwid-bench -t peakflops_avx -W N:2048kB:64 # GFlops/s: 1256.6
# Upper Peak GFLOPS Estimate: 2.2 GHz x 32 FLOPS/cycle x 64 cores = 4505.6 GFLOPS
srun likwid-bench -t peakflops_avx512 -W N:2048kB:64 # GFlops/s: 1256.6
# Remark: peakflops_avx512 is not available

# Bandwidth according to cluster documentation: GByte/s: 158.95 GiB/s
srun likwid-bench -t load_avx512 -W N:2GB:64 # GByte/s: 355.6