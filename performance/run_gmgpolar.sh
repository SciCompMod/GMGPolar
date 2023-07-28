#!/bin/bash

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
mod_pk=1		    # mod_pk (0/1)
prob=5			    # prob
alpha_coeff=2		# Zoni shifted
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
cores=64    # set OpenMP Num Threads to maximum number of cores requested

echo "#!/bin/bash" > run_gmgpolar_sbatch.sh
# create a short name for your job
echo "#SBATCH --job-name=gmgpolar" >> run_gmgpolar_sbatch.sh
# stdout file %A=job id
echo "#SBATCH --output=slurm-%A-p$prob-r$nr_exp-mpk$mod_pk-s$smoother-e$extrapolation--N$nodes-R$ranks-maxC$cores.out" >> run_gmgpolar_sbatch.sh
# stderr file
echo "#SBATCH --error=slurm-%A-p$prob-r$nr_exp-mpk$mod_pk-s$smoother-e$extrapolation--N$nodes-R$ranks-maxC$cores.err" >> run_gmgpolar_sbatch.sh

echo "#SBATCH -N $nodes" >> run_gmgpolar_sbatch.sh
echo "#SBATCH -n $ranks" >> run_gmgpolar_sbatch.sh
echo "#SBATCH -c $cores" >> run_gmgpolar_sbatch.sh
echo "#SBATCH -t 600" >> run_gmgpolar_sbatch.sh
echo "#SBATCH --exclusive" >> run_gmgpolar_sbatch.sh
 
# remove potentially loaded and conflicting modules
echo "module purge" >> run_gmgpolar_sbatch.sh

# CARO
#echo "module load rev/23.05" >> run_gmgpolar_sbatch.sh
echo "module load mumps/5.5.1" >> run_gmgpolar_sbatch.sh
# Local machine
# echo "module load PrgEnv/gcc10-openmpi" >> run_gmgpolar_sbatch.sh

echo "cd .." >> run_gmgpolar_sbatch.sh
# echo "make -j16" >> run_gmgpolar_sbatch.sh

####################################
## create grids                   ##
####################################
mod_pk=0 # mod_pk has no effect on the creation of grids as the set of (r,theta) is
		 # the same for all geometries, only the mapping F(r,theta) -> (x,y) changes.

mkdir -p angles_files/Rmax"$R"/aniso"$fac_ani"/
mkdir -p radii_files/Rmax"$R"/aniso"$fac_ani"/

echo $prob $alpha_coeff $beta_coeff $fac_ani $extrapolation $mod_pk
for divideBy2 in 0 1 2 3 4 5 6          # create different grid sizes
do
	## ATTENTION / REMARK: 
	## Please note that these calls will abort/segfault as creation of grids and computation in one step
	## is not yet supported by GMGPolar. We will make this functionality available in a future commit.
	## Please ignore abort/segfault for the calls in this loop.
	./build/gmgpolar_simulation -n $nr_exp -a $fac_ani --mod_pk $mod_pk --DirBC_Interior $DirBC_Interior --divideBy2 $divideBy2 -r $R0  --smoother $smoother --verbose 2 --debug $debug --extrapolation $extrapolation --optimized 1 --openmp $openmp --v1 $v1 --v2 $v2 -R $R --prob $prob  --maxiter $maxiter --alpha_coeff $alpha_coeff --beta_coeff $beta_coeff --res_norm $res_norm --write_radii_angles 1 --f_grid_r "radii_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt" --f_grid_theta "angles_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt"
done

# to be defined for use case
divideBy2=4		    # divideBy2 (3/4/5/6)

####################################
## solve system                   ##
####################################

# reduce cores as cores count from 0
max_threads=$((cores))
echo "let m=1" >> run_gmgpolar_sbatch.sh
# FLOPS-DP counter from 1 to cores many threads
echo "while [ \$m -le $max_threads ]; do" >> run_gmgpolar_sbatch.sh
echo "let mminus1=m-1" >> run_gmgpolar_sbatch.sh
echo "likwid-perfctr -C 0-\$mminus1 -g FLOPS_DP ./build/gmgpolar_simulation --openmp \$m --matrix_free 1 -n $nr_exp -a $fac_ani --mod_pk $mod_pk --DirBC_Interior $DirBC_Interior --divideBy2 $divideBy2 -r $R0 --smoother $smoother -E $extrapolation --verbose 2 --debug $debug --optimized 1 --v1 $v1 --v2 $v2 -R $R --prob $prob --maxiter $maxiter --alpha_coeff $alpha_coeff --beta_coeff $beta_coeff --res_norm $res_norm --f_grid_r "radii_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt" --f_grid_theta "angles_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt" --rel_red_conv $rel_red_conv" >> run_gmgpolar_sbatch.sh
echo "let m=m*2" >> run_gmgpolar_sbatch.sh
echo "done;" >> run_gmgpolar_sbatch.sh

# # Memory (saturation) benchmarks
echo "let m=1" >> run_gmgpolar_sbatch.sh
echo "while [ \$m -le $max_threads ]; do" >> run_gmgpolar_sbatch.sh
echo "let mminus1=m-1" >> run_gmgpolar_sbatch.sh
echo "likwid-perfctr -C 0-\$mminus1 -g CACHES ./build/gmgpolar_simulation --openmp \$m --matrix_free 1 -n $nr_exp -a $fac_ani --mod_pk $mod_pk --DirBC_Interior $DirBC_Interior --divideBy2 $divideBy2 -r $R0 --smoother $smoother -E $extrapolation --verbose 2 --debug $debug --optimized 1 --v1 $v1 --v2 $v2 -R $R --prob $prob --maxiter $maxiter --alpha_coeff $alpha_coeff --beta_coeff $beta_coeff --res_norm $res_norm --f_grid_r "radii_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt" --f_grid_theta "angles_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt" --rel_red_conv $rel_red_conv" >> run_gmgpolar_sbatch.sh
echo "let m=m*2" >> run_gmgpolar_sbatch.sh
echo "done;" >> run_gmgpolar_sbatch.sh

#submit the job
sbatch run_gmgpolar_sbatch.sh