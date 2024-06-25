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

m=55
srun ./build/gmgpolar --maxOpenMPThreads $m --finestLevelThreads $55 --threadReductionFactor 1.0 --R0 1e-5 --Rmax 1.3 --nr_exp 15 --ntheta_exp 15 --anisotropic_factor 0 --divideBy2 0 --write_grid_file 0 --load_grid_file 0 --file_grid_r "_radii.txt" --file_grid_theta "_angles.txt" --DirBC_Interior 0 --extrapolation 0 --maxLevels -1 --v1 1 --v2 1 --cycle 1