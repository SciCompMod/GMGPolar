#!/bin/bash
#SBATCH --job-name=gmgpolar
#SBATCH --output=Output/slurm-%A-p6-r4-dbt7-mpk2-s3-e1--N1-R1-maxC128.out
#SBATCH --error=Output/slurm-%A-p6-r4-dbt7-mpk2-s3-e1--N1-R1-maxC128.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 56
#SBATCH --threads-per-core=1
#SBATCH -t 1400
# #SBATCH --nodelist="be-cpu01"
#SBATCH --exclusive

verbose=1
maxOpenMPThreads=32
threadReductionFactor=1.0

R0=1e-8
Rmax=1.0
nr_exp=12
ntheta_exp=-1
anisotropic_factor=0
divideBy2=0

paraview=0
write_grid_file=0
load_grid_file=0
file_grid_radii="_radii.txt"
file_grid_angles="_angles.txt"

# Interior boundary condition: 
# 0: Across-origin
# 1: u_D_Interior
DirBC_Interior=0

### Custom Test Cases ###
geometry=2 # Circular (0), Shafranov(1), Czarny(2), Culham (3)
problem=2 # CartesianR2(0), CartesianR6(1), PolarR6(2), RefinedRadius(3)
alpha_coeff=1 # Poisson(0), Sonnendrucker(1), Zoni(2), Zoni-Shifted(3)
beta_coeff=1 # Zero(0), Gyro - Alpha Inverse(1)
# Remark: For RefinedRadius choose alpha_coeff=3, beta_coeff=1
# Remark: For Culham Geometry choose geometry=3, problem=2,3, alpha_coeff=3, beta_coeff=1
# Remark: For Culham Geometry the provided exactSolution may be incorrect.

# Full Multigrid Method:
# 0: Initial approximation is set to zero
# 1: Initial approximation obtained by nested iteration (recommended)
FMG=1
FMG_iterations=3
FMG_multigridCycle=2
# Extrapolation Method:
# 0: No extrapolation
# 1: Implicit extrapolation (recommended)
# 2: Implicit extrapolation with full grid smoothing (residuals cannot be used as convergence criteria)
# 3: Combination of both implicit extrapolation methods (usefull if FMG=0)
extrapolation=1
# Maximum number of multigrid levels:
maxLevels=7
# Number of smoothing steps:
preSmoothingSteps=1
postSmoothingSteps=1
# Multigrid Cycle:
# 0: V-Cycle
# 1: W-Cycle
# 2: F-Cycle
multigridCycle=0

# Convergence criteria:
maxIterations=150
residualNormType=0 # L2-Norm(0) = 0, Weighted L2-Norm(1), Infinity-Norm(2)
absoluteTolerance=1e-15
relativeTolerance=1e-15

# Define additional geometry parameters
kappa_eps=0.0
delta_e=0.0
if [ "$geometry" -eq 1 ]; then
    kappa_eps=0.3
    delta_e=0.2
elif [ "$geometry" -eq 2 ]; then
    kappa_eps=0.3
    delta_e=1.4
fi

# Set alpha_jump based on alpha_coeff value
if [ "$alpha_coeff" -eq 0 ]; then
    alpha_jump=$(echo "0.5 * $Rmax" | bc)
elif [ "$alpha_coeff" -eq 1 ]; then
    alpha_jump=$(echo "0.66 * $Rmax" | bc)
elif [ "$alpha_coeff" -eq 2 ]; then
    alpha_jump=$(echo "0.4837 * $Rmax" | bc)
elif [ "$alpha_coeff" -eq 3 ]; then
    alpha_jump=$(echo "0.7081 * $Rmax" | bc)
else
    echo "Invalid value for alpha_coeff: $alpha_coeff"
    exit 1
fi

srun ./build/gmgpolar --verbose $verbose --maxOpenMPThreads $maxOpenMPThreads --threadReductionFactor $threadReductionFactor --R0 $R0 --Rmax $Rmax --nr_exp $nr_exp --ntheta_exp $ntheta_exp --anisotropic_factor $anisotropic_factor --divideBy2 $divideBy2 --write_grid_file $write_grid_file --load_grid_file $load_grid_file --file_grid_radii "$file_grid_radii" --file_grid_angles "$file_grid_angles" --DirBC_Interior $DirBC_Interior --geometry $geometry --kappa_eps $kappa_eps --delta_e $delta_e --problem $problem --alpha_coeff $alpha_coeff --alpha_jump $alpha_jump --beta_coeff $beta_coeff --FMG $FMG --extrapolation $extrapolation --maxLevels $maxLevels --preSmoothingSteps $preSmoothingSteps --postSmoothingSteps $postSmoothingSteps --multigridCycle $multigridCycle --maxIterations $maxIterations --residualNormType $residualNormType --absoluteTolerance $absoluteTolerance --relativeTolerance $relativeTolerance