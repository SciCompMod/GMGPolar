#!/bin/bash
#SBATCH --job-name=gmgpolar
#SBATCH --output=Output/slurm-%A-p6-r4-dbt7-mpk2-s3-e1--N1-R1-maxC128.out
#SBATCH --error=Output/slurm-%A-p6-r4-dbt7-mpk2-s3-e1--N1-R1-maxC128.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 55
#SBATCH --threads-per-core=1
#SBATCH -t 28800
#SBATCH --exclusive

maxOpenMPThreads=10
finestLevelThreads=10
threadReductionFactor=1.0

R0=1e-5
Rmax=1.3
nr_exp=4
ntheta_exp=-1
anisotropic_factor=0
divideBy2=6

write_grid_file=0
load_grid_file=0
file_grid_radii="_radii.txt"
file_grid_angles="_angles.txt"

DirBC_Interior=0 # Across-origin(0), u_D_Interior(1)

# Test Cases #
geometry=1 # Circular (0), Shafranov(1), Czarny(2), Culham (3)
problem=0 # CartesianR2(0), CartesianR6(1), PolarR6(2), RefinedRadius(3)
alpha_coeff=1 # Poisson(0), Sonnendrucker(1), Zoni(2), Zoni-Shifted(3)
beta_coeff=0 # Zero(0), Gyro - Alpha Inverse(1)
# Remark: For Culham Geometry choose geometry=3, problem=2,3, alpha_coeff=3, beta_coeff=1

extrapolation=1
maxLevels=5
preSmoothingSteps=1
postSmoothingSteps=1
multigridCycle=0

maxIterations=150
residualNormType=0
absoluteTolerance=1e-100
relativeTolerance=1e-8

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
    alpha_jump=0.5
elif [ "$alpha_coeff" -eq 1 ]; then
    alpha_jump=0.66
elif [ "$alpha_coeff" -eq 2 ]; then
    alpha_jump=0.4837
elif [ "$alpha_coeff" -eq 3 ]; then
    alpha_jump=0.7081
else
    echo "Invalid value for alpha_coeff: $alpha_coeff"
    exit 1
fi

srun ./build/gmgpolar --maxOpenMPThreads $maxOpenMPThreads --finestLevelThreads $finestLevelThreads --threadReductionFactor $threadReductionFactor --R0 $R0 --Rmax $Rmax --nr_exp $nr_exp --ntheta_exp $ntheta_exp --anisotropic_factor $anisotropic_factor --divideBy2 $divideBy2 --write_grid_file $write_grid_file --load_grid_file $load_grid_file --file_grid_radii "$file_grid_radii" --file_grid_angles "$file_grid_angles" --DirBC_Interior $DirBC_Interior --geometry $geometry --kappa_eps $kappa_eps --delta_e $delta_e --problem $problem --alpha_coeff $alpha_coeff --alpha_jump $alpha_jump --beta_coeff $beta_coeff --extrapolation $extrapolation --maxLevels $maxLevels --preSmoothingSteps $preSmoothingSteps --postSmoothingSteps $postSmoothingSteps --multigridCycle $multigridCycle --maxIterations $maxIterations --residualNormType $residualNormType --absoluteTolerance $absoluteTolerance --relativeTolerance $relativeTolerance
