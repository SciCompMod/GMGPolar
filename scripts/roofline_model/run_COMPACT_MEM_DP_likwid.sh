#!/bin/bash
#SBATCH --job-name=MEM_DP
#SBATCH --output=slurm-%A-FLOPS_DP_CARA.out
#SBATCH --error=slurm-%A-FLOPS_DP_CARA.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH --threads-per-core=1
#SBATCH --time=0:30:00
#SBATCH --exclusive
#SBATCH --partition=naples128
#SBATCH --account=2476029

# potentially run benchmark on machine
# srun --cpus-per-task=64 likwid-bench -i 1000 -t stream_avx_fma -w S0:500MB:64

core_list=( 1 2 4 8 16 32 64 )
for m in ${core_list[@]}; do
    let mminus1=m-1
    output_file="data/COMPACT_MEM_DP_${m}.txt"
    # for testing that pin works correctly, potentially use likwid-pin beforehand
    # srun --cpus-per-task=64 likwid-pin -C E:N:$m ./../../build/gmgpolar --verbose 1 --paraview 0 --maxOpenMPThreads $m --threadReductionFactor 1.0 --stencilDistributionMethod 1 --cacheDensityProfileCoefficients 1 --cacheDomainGeometry 0 --R0 1e-8 --Rmax 1.3 --nr_exp 4 --ntheta_exp -1 --anisotropic_factor 3 --divideBy2 5 --write_grid_file 0 --load_grid_file 0 --file_grid_radii _radii.txt --file_grid_angles _angles.txt --DirBC_Interior 1 --geometry 2 --kappa_eps 0.3 --delta_e 1.4 --problem 1 --alpha_coeff 2 --alpha_jump 0.6288100000000001 --beta_coeff 1 --FMG 0 --FMG_iterations 3 --FMG_cycle 2 --extrapolation 1 --maxLevels 7 --preSmoothingSteps 1 --postSmoothingSteps 1 --multigridCycle 0 --maxIterations 150 --residualNormType 0 --absoluteTolerance 1e-10 --relativeTolerance 1e-10
    srun --cpus-per-task=64 likwid-perfctr -f -m -C E:N:$m -g MEM_DP -o $output_file ./../../build/gmgpolar --verbose 1 --paraview 0 --maxOpenMPThreads $m --threadReductionFactor 1.0 --stencilDistributionMethod 1 --cacheDensityProfileCoefficients 1 --cacheDomainGeometry 0 --R0 1e-8 --Rmax 1.3 --nr_exp 4 --ntheta_exp -1 --anisotropic_factor 3 --divideBy2 5 --write_grid_file 0 --load_grid_file 0 --file_grid_radii _radii.txt --file_grid_angles _angles.txt --DirBC_Interior 1 --geometry 2 --kappa_eps 0.3 --delta_e 1.4 --problem 1 --alpha_coeff 2 --alpha_jump 0.6288100000000001 --beta_coeff 1 --FMG 0 --FMG_iterations 3 --FMG_cycle 2 --extrapolation 1 --maxLevels 7 --preSmoothingSteps 1 --postSmoothingSteps 1 --multigridCycle 0 --maxIterations 150 --residualNormType 0 --absoluteTolerance 1e-10 --relativeTolerance 1e-10
done;
