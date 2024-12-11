#!/bin/bash

# Verbosity level: 
# 0 - No output 
# 1 - Basic output 
verbose=1
# Set Paraview usage flag:
# 0 - Do not use Paraview 
# 1 - Enable Paraview to visualize the grid, solution and error
paraview=0

# OpenMP settings:
maxOpenMPThreads=16

# Finest grid parameters
R0=1e-8
Rmax=1.3
nr_exp=4
ntheta_exp=-1
anisotropic_factor=0
divideBy2=3

# Interior boundary condition: 
# 0: Across-origin
# 1: u_D_Interior
DirBC_Interior=0

# Full Multigrid Method:
# 0: Initial approximation is set to zero
# 1: Initial approximation obtained by nested iteration (recommended)
FMG=1
FMG_iterations=3
FMG_cycle=2 # V-Cycle(0), W-Cycle(1), F-Cycle(2)

# Extrapolation Method:
# 0: No extrapolation
# 1: Implicit extrapolation (recommended)
# 2: Implicit extrapolation with full grid smoothing (residuals cannot be used as convergence criteria)
# 3: Combination of both implicit extrapolation methods (May be usefull for FMG=0)
extrapolation=1
# Maximum number of multigrid levels:
maxLevels=5
gpuLevels=2
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
absoluteTolerance=1e-8
relativeTolerance=1e-8

./../build/gmgpolar --verbose $verbose --paraview $paraview --maxOpenMPThreads $maxOpenMPThreads --R0 $R0 --Rmax $Rmax --nr_exp $nr_exp --ntheta_exp $ntheta_exp --anisotropic_factor $anisotropic_factor --divideBy2 $divideBy2 --DirBC_Interior $DirBC_Interior --FMG $FMG --FMG_iterations $FMG_iterations --FMG_cycle $FMG_cycle --extrapolation $extrapolation --maxLevels $maxLevels --gpuLevels $gpuLevels --preSmoothingSteps $preSmoothingSteps --postSmoothingSteps $postSmoothingSteps --multigridCycle $multigridCycle --maxIterations $maxIterations --residualNormType $residualNormType --absoluteTolerance $absoluteTolerance --relativeTolerance $relativeTolerance