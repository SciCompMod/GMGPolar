#!/bin/bash
#SBATCH --job-name=RoofLine
#SBATCH --output=slurm-%A-Roofline-Model.out
#SBATCH --error=slurm-%A-Roofline-Model.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 56
#SBATCH --threads-per-core=1
#SBATCH --nodelist="be-cpu03"
#SBATCH -t 1440
#SBATCH --exclusive

srun likwid-topology

# CPU name: Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
# CPU type: Intel Skylake SP processor
# CPU stepping: 4
# Sockets: 4
# Cores per socket: 14
# Threads per core: 1
# Level 1: 32 kB
# Level 2: 1 MB
# Level 3: 19 MB
# NUMA domains: 4

# So we have 4 x 14 hardware threads on this machine. 
# We use a stream size of 4 x 14 x 32 kB = 1792 kB, 
# 32kB for each of the 56 hardware threads so that each vector chunk fits into the L1 cache of one core.

# Upper Peak GFLOPS Estimate: 2.6 GHz x 32 FLOPS/cycle x 56 cores = 4659.2 GFLOPS
srun likwid-bench -t peakflops_avx -W N:1792kB:56 # GFlops/s: 1256.6
# Remark: peakflops_avx512 is not available

# Bandwidth according to cluster documentation: GByte/s: 281 GiB/s
srun likwid-bench -t stream_mem_avx_fma -W N:2GB:56 # GByte/s: 272