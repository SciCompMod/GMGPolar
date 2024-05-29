#!/bin/bash
#SBATCH --job-name=gmgpolar
#SBATCH --output=Output/slurm-%A-p6-r4-dbt7-mpk2-s3-e1--N1-R1-maxC128.out
#SBATCH --error=Output/slurm-%A-p6-r4-dbt7-mpk2-s3-e1--N1-R1-maxC128.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 60
#SBATCH --threads-per-core=1
#SBATCH -t 28800
#SBATCH --exclusive


# module purge
# module load likwid/likwid-5.2.1/likwid-5.2.1-gcc-9.4.0-szi3i3h

m=60
srun ./build/gmgpolar --nr_exp 10 --ntheta_exp 10 --maxOpenMPThreads $m --finestLevelThreads $m --problem 0 --geometry 0 --alpha_coeff 0 --beta_coeff 1 --DirBC_Interior 0

# mminus1=55
# srun --cpus-per-task=$m likwid-perfctr -f -m -C 0-$mminus1 -g FLOPS_DP ./build/gmgpolar --nr_exp 10 --ntheta_exp 10 --maxOpenMPThreads $m --finestLevelThreads $m --problem 0 --geometry 0 --alpha_coeff 0 --beta_coeff 1

# srun --cpus-per-task=$m likwid-pin -c N:0-$mminus1 ./build/gmgpolar --nr_exp 13 --ntheta_exp 14 --openmp $m --problem 0 --geometry 0 --alpha_coeff 0 --beta_coeff 1


# NODELIST: be-cpu05



# -m: number of regions
# -c: 0,...,m-1 pins cores used


# sacct -a --format="JOBID"

# module purge
# module load likwid/likwid-5.2.1/likwid-5.2.1-gcc-9.4.0-szi3i3h

# m=56
# mminus1=55
# # srun --cpus-per-task=$m likwid-perfctr -f -m -C 0-$mminus1 -g FLOPS_DP ./build/gmgpolar --nr_exp 13 --ntheta_exp 14 --openmp $m --problem 0 --geometry 0 --alpha_coeff 0 --beta_coeff 1

# srun --cpus-per-task=$m likwid-pin -c N:0-$mminus1 ./build/gmgpolar --nr_exp 13 --ntheta_exp 14 --openmp $m --problem 0 --geometry 0 --alpha_coeff 0 --beta_coeff 1









# potentially run benchmark on machine
# srun --cpus-per-task=128 likwid-bench -i 1000 -t stream_avx_fma -w S0:500MB:64


# let divideBy2=7
# let m=1
# while [ $m -le 128 ]; do
# let mminus1=m-1
# # for testing that pin works correctly, potentially use likwid-pin beforehand
# # srun --cpus-per-task=$m likwid-pin -c N:0-$mminus1 ./build/gmgpolar_simulation --openmp $m --matrix_free 1 -n 4 -a 3 --mod_pk 2 --DirBC_Interior 1 --divideBy2 0 -r 1e-8 --smoother 3 -E 1 --verbose 2 --debug 0 --optimized 1 --v1 1 --v2 1 -R 1.3 --prob 6 --maxiter 300 --alpha_coeff 2 --beta_coeff 1 --res_norm 3 --f_grid_r radii_files/Rmax1.3/aniso3/divide$divideBy2.txt --f_grid_theta angles_files/Rmax1.3/aniso3/divide$divideBy2.txt --rel_red_conv 1e-11
# srun --cpus-per-task=$m likwid-perfctr -f -m -C 0-$mminus1 -g FLOPS_DP ./build/gmgpolar_simulation --openmp $m --matrix_free 1 -n 4 -a 3 --mod_pk 2 --DirBC_Interior 1 --divideBy2 0 -r 1e-8 --smoother 3 -E 1 --verbose 2 --debug 0 --optimized 1 --v1 1 --v2 1 -R 1.3 --prob 6 --maxiter 300 --alpha_coeff 2 --beta_coeff 1 --res_norm 3 --f_grid_r radii_files/Rmax1.3/aniso3/divide$divideBy2.txt --f_grid_theta angles_files/Rmax1.3/aniso3/divide$divideBy2.txt --rel_red_conv 1e-11
# let m=m*2
# done;
# let m=1
# while [ $m -le 128 ]; do
# let mminus1=m-1
# srun --cpus-per-task=$m likwid-perfctr -f -m -C 0-$mminus1 -g MEM_DP ./build/gmgpolar_simulation --openmp $m --matrix_free 1 -n 4 -a 3 --mod_pk 2 --DirBC_Interior 1 --divideBy2 0 -r 1e-8 --smoother 3 -E 1 --verbose 2 --debug 0 --optimized 1 --v1 1 --v2 1 -R 1.3 --prob 6 --maxiter 300 --alpha_coeff 2 --beta_coeff 1 --res_norm 3 --f_grid_r radii_files/Rmax1.3/aniso3/divide$divideBy2.txt --f_grid_theta angles_files/Rmax1.3/aniso3/divide$divideBy2.txt --rel_red_conv 1e-11
# let m=m*2
# done;