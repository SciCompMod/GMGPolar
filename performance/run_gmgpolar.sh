#!/bin/bash
#SBATCH --job-name=gmgpolar-setup
#SBATCH --output=slurm-%A-setup.out
#SBATCH --error=slurm-%A-setup.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 5

#fixed variables
# If Origin is chosen, the node can be set as coarse or fine. Default: Coarse.
origin_NOT_coarse=0	# origin_NOT_coarse
# Choose anisotropy in angle-direction. Default: Off.
theta_aniso=0		# theta_aniso
# Smoother 3 is our default, 13 is used for some testing, should be not used
# for production. -> TODO: Rename smoother names
smoother=3		    # smoother (3,13)

# default variables
# If origin is not a particular node of the mesh, Dirichlet boundary conditions
# can be implemented on the most inner circle
DirBC_Interior=1	#  DirBC_Interior (0/1)
# Generalized radius of most inner circle. Defines if origin will be a particular node.
R0=1e-8		        # r (1e-8/1e-5/1e-2)
# Generalized radius of maximum outer circle.
R=1.3			    # R
# Anisotropy in radial direction.
fac_ani=3		    # a
# TODO: which nr_exp and divideby2 do we want to consider?
nr_exp=4		    # n

#changing variables
mod_pk=2		    # mod_pk=1: Shafranov geometry
prob=6			    # Prob=7: Simulate solution (23) of Bourne et al.
alpha_coeff=2
beta_coeff=1

# set to on
extrapolation=1		# E

debug=0
v1=1
v2=1
maxiter=300
res_norm=3
rel_red_conv=1e-11

nodes=1
ranks=1     # number of MPI Ranks
cores=128    # set OpenMP Num Threads to maximum number of cores requested

####################################
## create grids                   ##
####################################
create_grid=0
if [ $create_grid ]
then
cd ..
mkdir -p angles_files/Rmax"$R"/aniso"$fac_ani"/
mkdir -p radii_files/Rmax"$R"/aniso"$fac_ani"/
# Costly function as setup as expensive and sequential. Only run once.
for divideBy2 in 0 1 2 3 4 5 6 7 8     # create different grid sizes
do
	## ATTENTION / REMARK: 
	## Please note that these calls will abort/segfault as creation of grids and computation in one step
	## is not yet supported by GMGPolar. We will make this functionality available in a future commit.
	## Please ignore abort/segfault for the calls in this loop.
	# mod_pk has no effect on the creation of grids as the set of (r,theta) is
	# the same for all geometries, only the mapping F(r,theta) -> (x,y) changes.
	./build/gmgpolar_simulation -n $nr_exp -a $fac_ani --mod_pk 0 --DirBC_Interior $DirBC_Interior --divideBy2 $divideBy2 -r $R0  --smoother $smoother --verbose 2 --debug $debug --extrapolation $extrapolation --optimized 1 $ --v1 $v1 --v2 $v2 -R $R --prob $prob  --maxiter $maxiter --alpha_coeff $alpha_coeff --beta_coeff $beta_coeff --res_norm $res_norm --write_radii_angles 1 --f_grid_r "radii_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt" --f_grid_theta "angles_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt"
done
fi

echo "#!/bin/bash" > run_gmgpolar_sbatch.sh
# create a short name for your job
echo "#SBATCH --job-name=gmgpolar" >> run_gmgpolar_sbatch.sh
# Change fixed "dbt7" in file name if divideBy2 below is changed!
# stdout file %A=job id
echo "#SBATCH --output=slurm-%A-p$prob-r$nr_exp-dbt7-mpk$mod_pk-s$smoother-e$extrapolation--N$nodes-R$ranks-maxC$cores.out" >> run_gmgpolar_sbatch.sh
# stderr file
echo "#SBATCH --error=slurm-%A-p$prob-r$nr_exp-dbt7-mpk$mod_pk-s$smoother-e$extrapolation--N$nodes-R$ranks-maxC$cores.err" >> run_gmgpolar_sbatch.sh

echo "#SBATCH -N $nodes" >> run_gmgpolar_sbatch.sh
echo "#SBATCH -n $ranks" >> run_gmgpolar_sbatch.sh
echo "#SBATCH -c $cores" >> run_gmgpolar_sbatch.sh
# fix to one thread per core
echo "#SBATCH --threads-per-core=1" >> run_gmgpolar_sbatch.sh
# fix CPU frequency to 1.8 Mhz
echo "#SBATCH --cpu-freq=1800000" >> run_gmgpolar_sbatch.sh
## Estimation on 2x AMD EPYC 7702: nr=4, dbt=7, prob=7, mpk=1 takes 720 minutes for 1->2->...->128 cores ##
echo "#SBATCH -t 1600" >> run_gmgpolar_sbatch.sh
echo "#SBATCH --exclusive" >> run_gmgpolar_sbatch.sh
 
# remove potentially loaded and conflicting modules
echo "module purge" >> run_gmgpolar_sbatch.sh

# CARO
#echo "module load rev/23.05" >> run_gmgpolar_sbatch.sh
# spack install mumps@XXX+metis~mpi
echo "module load likwid/5.2.2" >> run_gmgpolar_sbatch.sh
# Local machine
# echo "module load PrgEnv/gcc10-openmpi" >> run_gmgpolar_sbatch.sh

# echo "cd .." >> run_gmgpolar_sbatch.sh
# echo "make -j16" >> run_gmgpolar_sbatch.sh

echo "# potentially run benchmark on machine" >> run_gmgpolar_sbatch.sh
echo "# srun --cpus-per-task=$((cores)) likwid-bench -i 1000 -t stream_avx_fma -w S0:500MB:64" >> run_gmgpolar_sbatch.sh

# to be defined for use case (3/4/5/6)
# Attention: divideBy is used as a dummy variable to access folders as grids are read in
echo "let divideBy2=7" >> run_gmgpolar_sbatch.sh

####################################
## solve system                   ##
####################################

# reduce cores as cores count from 0
max_threads=$((cores))
echo "let m=1" >> run_gmgpolar_sbatch.sh
# FLOPS-DP counter from 1 to cores many threads
echo "while [ \$m -le $max_threads ]; do" >> run_gmgpolar_sbatch.sh
echo "let mminus1=m-1" >> run_gmgpolar_sbatch.sh

echo "# for testing that pin works correctly, potentially use likwid-pin beforehand" >> run_gmgpolar_sbatch.sh
echo "# srun --cpus-per-task=\$m likwid-pin -c N:0-\$mminus1 ./build/gmgpolar_simulation --openmp \$m --matrix_free 1 -n $nr_exp -a $fac_ani --mod_pk $mod_pk --DirBC_Interior $DirBC_Interior --divideBy2 0 -r $R0 --smoother $smoother -E $extrapolation --verbose 2 --debug $debug --optimized 1 --v1 $v1 --v2 $v2 -R $R --prob $prob --maxiter $maxiter --alpha_coeff $alpha_coeff --beta_coeff $beta_coeff --res_norm $res_norm --f_grid_r "radii_files/Rmax"$R"/aniso"$fac_ani"/divide"\$divideBy2".txt" --f_grid_theta "angles_files/Rmax"$R"/aniso"$fac_ani"/divide"\$divideBy2".txt" --rel_red_conv $rel_red_conv" >> run_gmgpolar_sbatch.sh

echo "srun --cpus-per-task=\$m likwid-perfctr -f -m -C 0-\$mminus1 -g FLOPS_DP ./build/gmgpolar_simulation --openmp \$m --matrix_free 1 -n $nr_exp -a $fac_ani --mod_pk $mod_pk --DirBC_Interior $DirBC_Interior --divideBy2 0 -r $R0 --smoother $smoother -E $extrapolation --verbose 2 --debug $debug --optimized 1 --v1 $v1 --v2 $v2 -R $R --prob $prob --maxiter $maxiter --alpha_coeff $alpha_coeff --beta_coeff $beta_coeff --res_norm $res_norm --f_grid_r "radii_files/Rmax"$R"/aniso"$fac_ani"/divide"\$divideBy2".txt" --f_grid_theta "angles_files/Rmax"$R"/aniso"$fac_ani"/divide"\$divideBy2".txt" --rel_red_conv $rel_red_conv" >> run_gmgpolar_sbatch.sh
echo "let m=m*2" >> run_gmgpolar_sbatch.sh
echo "done;" >> run_gmgpolar_sbatch.sh

# # Memory (saturation) benchmarks
echo "let m=1" >> run_gmgpolar_sbatch.sh
echo "while [ \$m -le $max_threads ]; do" >> run_gmgpolar_sbatch.sh
echo "let mminus1=m-1" >> run_gmgpolar_sbatch.sh
echo "srun --cpus-per-task=\$m likwid-perfctr -f -m -C 0-\$mminus1 -g MEM_DP ./build/gmgpolar_simulation --openmp \$m --matrix_free 1 -n $nr_exp -a $fac_ani --mod_pk $mod_pk --DirBC_Interior $DirBC_Interior --divideBy2 0 -r $R0 --smoother $smoother -E $extrapolation --verbose 2 --debug $debug --optimized 1 --v1 $v1 --v2 $v2 -R $R --prob $prob --maxiter $maxiter --alpha_coeff $alpha_coeff --beta_coeff $beta_coeff --res_norm $res_norm --f_grid_r "radii_files/Rmax"$R"/aniso"$fac_ani"/divide"\$divideBy2".txt" --f_grid_theta "angles_files/Rmax"$R"/aniso"$fac_ani"/divide"\$divideBy2".txt" --rel_red_conv $rel_red_conv" >> run_gmgpolar_sbatch.sh
echo "let m=m*2" >> run_gmgpolar_sbatch.sh
echo "done;" >> run_gmgpolar_sbatch.sh

#submit the job
sbatch run_gmgpolar_sbatch.sh