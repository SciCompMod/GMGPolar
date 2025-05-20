#!/bin/bash
#SBATCH --job-name=gmgpolar1.0
#SBATCH --output=slurm-%A-convergenceorder_GMGPolar_CARA.out
#SBATCH --error=slurm-%A-convergenceorder_GMGPolar_CARA.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH --threads-per-core=1
#SBATCH --time=0:30:00
#SBATCH --exclusive
#SBATCH --partition=naples128
#SBATCH --account=2476029

R0=1e-8
Rmax=1.3
nr_exp=4
ntheta_exp=-1
anisotropic_factor=3
divideBy2=6

DirBC_Interior=0 # Across-origin(0), u_D_Interior(1)

# Test Cases #
geometry=2 # Circular (0), Shafranov(1), Czarny(2), Culham (3)
kappa_eps=0.3 # Shafranov: 0.3, Czarny = 0.3, else unused
delta_e=1.4 # Shafranov: 0.2, Czarny = 1.4, else unused
problem=6 # CartesianR2(5), CartesianR6(7), PolarR6(6), RefinedRadius(4)
alpha_coeff=0 # Poisson(3), Sonnendrucker(0), Zoni(1), Zoni-Shifted(2)
beta_coeff=1 # Zero(0), Gyro - Alpha Inverse(1)

extrapolation=1
maxLevels=1e6
preSmoothingSteps=1
postSmoothingSteps=1
multigridCycle=0

maxIterations=300
residualNormType=0 # convergence_criterium: scaled L2 (0), scaled inf (1), L2 (2), inf (3)
absoluteTolerance=1e-14
relativeTolerance=1e-8

# Loop over different core counts
for cores in 1 2 4 8 16 32 64; do
    # Set the number of OpenMP threads
    export OMP_NUM_THREADS=$cores
    
    echo "Running with $cores cores/threads"
    
    # Run the simulation
    ./../../../GMGPolar_1.0/build/gmgpolar_simulation \
        --verbose 2 \
        --optimized 1 \
        --matrix_free 1 \
        --debug 0 \
        --compute_rho 1 \
        --check_error 1 \
        --write_radii_angles 0 \
        --nr_exp $nr_exp \
        --ntheta_exp $ntheta_exp \
        --fac_ani $anisotropic_factor \
        --level $maxLevels \
        --maxiter $maxIterations \
        --v1 $preSmoothingSteps \
        --v2 $postSmoothingSteps \
        --cycle $multigridCycle \
        --mod_pk $geometry \
        --extrapolation $extrapolation \
        --DirBC_Interior $DirBC_Interior \
        --divideBy2 $divideBy2 \
        --prob $problem \
        --alpha_coeff $alpha_coeff \
        --beta_coeff $beta_coeff \
        --openmp $cores \
        --res_norm $residualNormType \
        --R0 $R0 \
        --R $Rmax \
        --kappa_eps $kappa_eps \
        --delta_e $delta_e \
        --tol_bound_check $absoluteTolerance \
        --rel_red_conv $relativeTolerance
    
    echo "Completed run with $cores cores/threads"
    echo "--------------------------------------"
done