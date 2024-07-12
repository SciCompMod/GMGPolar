#!/bin/bash

maxOpenMPThreads=1
finestLevelThreads=1
threadReductionFactor=1.0

# R0=1e-5
# Rmax=1.3
# nr_exp=4
# ntheta_exp=6
# anisotropic_factor=3
# divideBy2=5

R0=1e-5
Rmax=1.3
nr_exp=8
ntheta_exp=-1
anisotropic_factor=0
divideBy2=0

write_grid_file=1
load_grid_file=0
file_grid_radii="_radii.txt"
file_grid_angles="_angles.txt"

DirBC_Interior=0 # Across-origin(0), u_D_Interior(1)

# Test Cases #
geometry=2 # Circular 0), Shafranov(1), Czarny(2), Culham (3)
kappa_eps=0.3 # Shafranov: 0.3, Czarny = 0.3, else unused
delta_e=0.2 # Shafranov: 0.2, Czarny = 1.4, else unused
problem=0 # CartesianR2(0), CartesianR6(1), PolarR6(2), RefinedRadius(3)
alpha_coeff=1 # Poisson(0), Sonnendrucker(1), Zoni(2), Zoni-Shifted(3)
alpha_jump=0.5 # Poisson: 0.5, Sonnendrucker: 0.66, Zoni: 0.4837, Zoni-Shifted: 0.7081
beta_coeff=0 # Zero(0), Gyro - Alpha Inverse(1)
# Remark: For Culham Geometry choose geometry=3, problem=2,3, alpha_coeff=3, beta_coeff=1

extrapolation=0
maxLevels=-1
preSmoothingSteps=1
postSmoothingSteps=1
multigridCycle=0

maxIterations=150
residualNormType=1
absoluteTolerance=1e-12
relativeTolerance=1e-8

./build/gmgpolar --maxOpenMPThreads $maxOpenMPThreads --finestLevelThreads $finestLevelThreads --threadReductionFactor $threadReductionFactor --R0 $R0 --Rmax $Rmax --nr_exp $nr_exp --ntheta_exp $ntheta_exp --anisotropic_factor $anisotropic_factor --divideBy2 $divideBy2 --write_grid_file $write_grid_file --load_grid_file $load_grid_file --file_grid_radii "$file_grid_radii" --file_grid_angles "$file_grid_angles" --DirBC_Interior $DirBC_Interior --geometry $geometry --kappa_eps $kappa_eps --delta_e $delta_e --problem $problem --alpha_coeff $alpha_coeff --alpha_jump $alpha_jump --beta_coeff $beta_coeff --extrapolation $extrapolation --maxLevels $maxLevels --preSmoothingSteps $preSmoothingSteps --postSmoothingSteps $postSmoothingSteps --multigridCycle $multigridCycle --maxIterations $maxIterations --residualNormType $residualNormType --absoluteTolerance $absoluteTolerance --relativeTolerance $relativeTolerance







# maxOpenMPThreads=55
# finestLevelThreads=55
# threadReductionFactor=1.0

# R0=1e-5
# Rmax=1.3
# nr_exp=4
# ntheta_exp=6
# anisotropic_factor=3
# divideBy2=3

# write_grid_file=0
# load_grid_file=0
# file_grid_radii="_radii.txt"
# file_grid_angles="_angles.txt"

# DirBC_Interior=0 # Across-origin(0), u_D_Interior(1)

# # Test Cases #
# geometry=0 # Circular 0), Shafranov(1), Czarny(2), Culham (3)
# kappa_eps=0.3 # Shafranov: 0.3, Czarny = 0.3, else unused
# delta_e=0.2 # Shafranov: 0.2, Czarny = 1.4, else unused
# problem=0 # CartesianR2(0), CartesianR6(1), PolarR6(2), RefinedRadius(3)
# alpha_coeff=1 # Poisson(0), Sonnendrucker(1), Zoni(2), Zoni-Shifted(3)
# alpha_jump=0.5 # Poisson: 0.5, Sonnendrucker: 0.66, Zoni: 0.4837, Zoni-Shifted: 0.7081
# beta_coeff=0 # Zero(0), Gyro - Alpha Inverse(1)
# # Remark: For Culham Geometry choose geometry=3, problem=2,3, alpha_coeff=3, beta_coeff=1

# extrapolation=0
# maxLevels=-1
# preSmoothingSteps=1
# postSmoothingSteps=1
# multigridCycle=0

# maxIterations=150
# residualNormType=1
# absoluteTolerance=1e-8
# relativeTolerance=1e-8

# ./build/gmgpolar --maxOpenMPThreads $maxOpenMPThreads --finestLevelThreads $finestLevelThreads --threadReductionFactor $threadReductionFactor --R0 $R0 --Rmax $Rmax --nr_exp $nr_exp --ntheta_exp $ntheta_exp --anisotropic_factor $anisotropic_factor --divideBy2 $divideBy2 --write_grid_file $write_grid_file --load_grid_file $load_grid_file --file_grid_radii "$file_grid_radii" --file_grid_angles "$file_grid_angles" --DirBC_Interior $DirBC_Interior --geometry $geometry --kappa_eps $kappa_eps --delta_e $delta_e --problem $problem --alpha_coeff $alpha_coeff --alpha_jump $alpha_jump --beta_coeff $beta_coeff --extrapolation $extrapolation --maxLevels $maxLevels --preSmoothingSteps $preSmoothingSteps --postSmoothingSteps $postSmoothingSteps --multigridCycle $multigridCycle --maxIterations $maxIterations --residualNormType $residualNormType --absoluteTolerance $absoluteTolerance --relativeTolerance $relativeTolerance

