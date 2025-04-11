#!/bin/bash
#SBATCH --job-name=gmgpolar
#SBATCH --output=slurm-maxLevel13_GMGPolar_JEDI_profiles.out
#SBATCH --error=slurm-maxLevel13_GMGPolar_JEDI_profiles.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 72
#SBATCH --gres=gpu:1
#SBATCH --threads-per-core=1
#SBATCH --time=0:30:00
#SBATCH --exclusive
#SBATCH --account=training2508
#SBATCH --partition=all


nvidia-smi -L
echo $CUDA_VISIBLE_DEVICES 
# Verbosity level: 

# Set Paraview usage flag:
# 0 - Do not use Paraview 
# 1 - Enable Paraview to visualize the grid, solution and error
paraview=0

# OpenMP settings:
maxOpenMPThreads=72

# Finest grid parameters
R0=1e-8
Rmax=1.0
nr_exp=4
ntheta_exp=-1
anisotropic_factor=0

# Interior boundary condition: 
# 0: Across-origin
# 1: u_D_Interior
DirBC_Interior=1

# Full Multigrid Method:
# 0: Initial approximation is set to zero
# 1: Initial approximation obtained by nested iteration (recommended)
FMG=0
FMG_iterations=3
FMG_cycle=2 # V-Cycle(0), W-Cycle(1), F-Cycle(2)

# Extrapolation Method:
# 0: No extrapolation
# 1: Implicit extrapolation (recommended)
# 2: Implicit extrapolation with full grid smoothing (residuals cannot be used as convergence criteria)
# 3: Combination of both implicit extrapolation methods (May be usefull for FMG=0)
extrapolation=1
# Number of smoothing steps:
preSmoothingSteps=1
postSmoothingSteps=1
# Multigrid Cycle:
# 0: V-Cycle
# 1: W-Cycle
# 2: F-Cycle
multigridCycle=0

# Convergence criteria:
maxIterations=20
residualNormType=0 # L2-Norm(0) = 0, Weighted L2-Norm(1), Infinity-Norm(2)
absoluteTolerance=1e-100
relativeTolerance=1e-100

module load GCC CMake CUDA
spack load mumps metis
# MUMPS 5.5.1
# CUDA 12
gcc --version

# Maximum number of multigrid levels:
maxLevels=13
gpuLevels=-1

#### GPU ####
# Maximum number of multigrid levels:
maxLevels=13
gpuLevels=-1

# 0 - No output 
# 1 - Basic output
verbose=1

let divideBy2=4
while [ $divideBy2 -le 10 ]; do
srun --cpus-per-task=1 ./../build/gmgpolar --verbose $verbose --paraview $paraview --maxOpenMPThreads 1 --R0 $R0 --Rmax $Rmax --nr_exp $nr_exp --ntheta_exp $ntheta_exp --anisotropic_factor $anisotropic_factor --divideBy2 $divideBy2 --DirBC_Interior $DirBC_Interior --FMG $FMG --FMG_iterations $FMG_iterations --FMG_cycle $FMG_cycle --extrapolation $extrapolation --maxLevels $maxLevels --gpuLevels $gpuLevels --preSmoothingSteps $preSmoothingSteps --postSmoothingSteps $postSmoothingSteps --multigridCycle $multigridCycle --maxIterations $maxIterations --residualNormType $residualNormType --absoluteTolerance $absoluteTolerance --relativeTolerance $relativeTolerance
let divideBy2=divideBy2+1
done;


# Convergence criteria:
maxIterations=200
residualNormType=0 # L2-Norm(0) = 0, Weighted L2-Norm(1), Infinity-Norm(2)
absoluteTolerance=1e-14
relativeTolerance=1e-10

# 0 - No output 
# 1 - Basic output 
let divideBy2=4
while [ $divideBy2 -le 10 ]; do
srun --cpus-per-task=1 ./../build/gmgpolar --verbose $verbose --paraview $paraview --maxOpenMPThreads 1 --R0 $R0 --Rmax $Rmax --nr_exp $nr_exp --ntheta_exp $ntheta_exp --anisotropic_factor $anisotropic_factor --divideBy2 $divideBy2 --DirBC_Interior $DirBC_Interior --FMG $FMG --FMG_iterations $FMG_iterations --FMG_cycle $FMG_cycle --extrapolation $extrapolation --maxLevels $maxLevels --gpuLevels $gpuLevels --preSmoothingSteps $preSmoothingSteps --postSmoothingSteps $postSmoothingSteps --multigridCycle $multigridCycle --maxIterations $maxIterations --residualNormType $residualNormType --absoluteTolerance $absoluteTolerance --relativeTolerance $relativeTolerance
let divideBy2=divideBy2+1
done;

#### CPU #####
# Maximum number of multigrid levels:
maxLevels=13
gpuLevels=0

# Convergence criteria:
maxIterations=20
residualNormType=0 # L2-Norm(0) = 0, Weighted L2-Norm(1), Infinity-Norm(2)
absoluteTolerance=1e-100
relativeTolerance=1e-100

let divideBy2=4
while [ $divideBy2 -le 10 ]; do
srun --cpus-per-task=$maxOpenMPThreads ./../build/gmgpolar --verbose $verbose --paraview $paraview --maxOpenMPThreads $maxOpenMPThreads --R0 $R0 --Rmax $Rmax --nr_exp $nr_exp --ntheta_exp $ntheta_exp --anisotropic_factor $anisotropic_factor --divideBy2 $divideBy2 --DirBC_Interior $DirBC_Interior --FMG $FMG --FMG_iterations $FMG_iterations --FMG_cycle $FMG_cycle --extrapolation $extrapolation --maxLevels $maxLevels --gpuLevels $gpuLevels --preSmoothingSteps $preSmoothingSteps --postSmoothingSteps $postSmoothingSteps --multigridCycle $multigridCycle --maxIterations $maxIterations --residualNormType $residualNormType --absoluteTolerance $absoluteTolerance --relativeTolerance $relativeTolerance
let divideBy2=divideBy2+1
done;


# Convergence criteria:
maxIterations=200
residualNormType=0 # L2-Norm(0) = 0, Weighted L2-Norm(1), Infinity-Norm(2)
absoluteTolerance=1e-14
relativeTolerance=1e-10

# 0 - No output 
# 1 - Basic output 
let divideBy2=4
while [ $divideBy2 -le 10 ]; do
srun --cpus-per-task=$maxOpenMPThreads ./../build/gmgpolar --verbose $verbose --paraview $paraview --maxOpenMPThreads $maxOpenMPThreads --R0 $R0 --Rmax $Rmax --nr_exp $nr_exp --ntheta_exp $ntheta_exp --anisotropic_factor $anisotropic_factor --divideBy2 $divideBy2 --DirBC_Interior $DirBC_Interior --FMG $FMG --FMG_iterations $FMG_iterations --FMG_cycle $FMG_cycle --extrapolation $extrapolation --maxLevels $maxLevels --gpuLevels $gpuLevels --preSmoothingSteps $preSmoothingSteps --postSmoothingSteps $postSmoothingSteps --multigridCycle $multigridCycle --maxIterations $maxIterations --residualNormType $residualNormType --absoluteTolerance $absoluteTolerance --relativeTolerance $relativeTolerance
let divideBy2=divideBy2+1
done;

## profiles
divideBy2=10

maxLevels=13
gpuLevels=-1
srun --cpus-per-task=1 nsys profile --trace nvtx,cuda ./../build/gmgpolar --verbose $verbose --paraview $paraview --maxOpenMPThreads 1 --R0 $R0 --Rmax $Rmax --nr_exp $nr_exp --ntheta_exp $ntheta_exp --anisotropic_factor $anisotropic_factor --divideBy2 $divideBy2 --DirBC_Interior $DirBC_Interior --FMG $FMG --FMG_iterations $FMG_iterations --FMG_cycle $FMG_cycle --extrapolation $extrapolation --maxLevels $maxLevels --gpuLevels $gpuLevels --preSmoothingSteps $preSmoothingSteps --postSmoothingSteps $postSmoothingSteps --multigridCycle $multigridCycle --maxIterations $maxIterations --residualNormType $residualNormType --absoluteTolerance $absoluteTolerance --relativeTolerance $relativeTolerance

maxLevels=13
gpuLevels=0
srun --cpus-per-task=$maxOpenMPThreads nsys profile --trace nvtx,cuda ./../build/gmgpolar --verbose $verbose --paraview $paraview --maxOpenMPThreads $maxOpenMPThreads --R0 $R0 --Rmax $Rmax --nr_exp $nr_exp --ntheta_exp $ntheta_exp --anisotropic_factor $anisotropic_factor --divideBy2 $divideBy2 --DirBC_Interior $DirBC_Interior --FMG $FMG --FMG_iterations $FMG_iterations --FMG_cycle $FMG_cycle --extrapolation $extrapolation --maxLevels $maxLevels --gpuLevels $gpuLevels --preSmoothingSteps $preSmoothingSteps --postSmoothingSteps $postSmoothingSteps --multigridCycle $multigridCycle --maxIterations $maxIterations --residualNormType $residualNormType --absoluteTolerance $absoluteTolerance --relativeTolerance $relativeTolerance