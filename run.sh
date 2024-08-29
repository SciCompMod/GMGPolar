#!/bin/bash

verbose=1
maxOpenMPThreads=32
threadReductionFactor=1.0

R0=1e-5
Rmax=1.3
nr_exp=9
ntheta_exp=10
anisotropic_factor=0
divideBy2=0

write_grid_file=0
load_grid_file=0
file_grid_radii="_radii.txt"
file_grid_angles="_angles.txt"

DirBC_Interior=1 # Across-origin(0), u_D_Interior(1)

# Test Cases #
geometry=2 # Circular (0), Shafranov(1), Czarny(2), Culham (3)
problem=1 # CartesianR2(0), CartesianR6(1), PolarR6(2), RefinedRadius(3)
alpha_coeff=2 # Poisson(0), Sonnendrucker(1), Zoni(2), Zoni-Shifted(3)
beta_coeff=0 # Zero(0), Gyro - Alpha Inverse(1)
# Remark: For RefinedRadius choose alpha_coeff=3, beta_coeff=1
# Remark: For Culham Geometry choose geometry=3, problem=2,3, alpha_coeff=3, beta_coeff=1
# Remark: For Culham Geometry the provided exactSolution may be incorrect.

extrapolation=0
maxLevels=2
preSmoothingSteps=20
postSmoothingSteps=20
multigridCycle=0 # For multigridCycle > 0 we define reSmoothingSteps = preSmoothingSteps + postSmoothingSteps

maxIterations=150
residualNormType=0 # L2-Norm(0) = 0, Weighted L2-Norm(1), Infinity-Norm(2)
absoluteTolerance=1e-10
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

./build/gmgpolar --verbose $verbose --maxOpenMPThreads $maxOpenMPThreads --threadReductionFactor $threadReductionFactor --R0 $R0 --Rmax $Rmax --nr_exp $nr_exp --ntheta_exp $ntheta_exp --anisotropic_factor $anisotropic_factor --divideBy2 $divideBy2 --write_grid_file $write_grid_file --load_grid_file $load_grid_file --file_grid_radii "$file_grid_radii" --file_grid_angles "$file_grid_angles" --DirBC_Interior $DirBC_Interior --geometry $geometry --kappa_eps $kappa_eps --delta_e $delta_e --problem $problem --alpha_coeff $alpha_coeff --alpha_jump $alpha_jump --beta_coeff $beta_coeff --extrapolation $extrapolation --maxLevels $maxLevels --preSmoothingSteps $preSmoothingSteps --postSmoothingSteps $postSmoothingSteps --multigridCycle $multigridCycle --maxIterations $maxIterations --residualNormType $residualNormType --absoluteTolerance $absoluteTolerance --relativeTolerance $relativeTolerance
