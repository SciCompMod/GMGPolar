#!/bin/bash
#SBATCH --job-name=gmgpolar     # create a short name for your job
#SBATCH --output=slurm-%A.out # stdout file %A=job id
#SBATCH --error=slurm-%A.err  # stderr file
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 6000
#SBATCH --exclusive

##
## This is the run configuration to reproduce GMGPolar results of Figure 8 a) - d) in
## Bourne et al. - Solver comparison for Poisson-like equations on tokamak geometries (2023)
## 
#fixed variables (some of these variables have no effect, variables have to be cleaned up for a better presentation)
debug=0
v1=1
v2=1
maxiter=300
nr_exp=4
res_norm=3
R0=1e-6
DirBC_Interior=0
smoother=3
rel_red_conv=1e-11

# Problem
# prob=7 # setting not regarded: modified later on in loop
alpha_coeff=2		# Zoni shifted
beta_coeff=1

# geometry/grid
R=1.0
fac_ani=0
# mod_pk=0 # setting not regarded: modified later on in loop
kappa_eps=0.3
delta_e=0.2
# divideBy2=0 # setting not regarded: modified later on in loop

# MG
openmp=4
# extrapolation=1 # setting not regarded: modified later on in loop

build_dir=build

#overwrite existing file with an empty file
break > output.txt

##########
## create grids
#########
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
	./${build_dir}/gmgpolar_simulation -n $nr_exp -a $fac_ani --mod_pk $mod_pk --DirBC_Interior $DirBC_Interior --divideBy2 $divideBy2 -r $R0  --smoother $smoother --verbose 2 --debug $debug --extrapolation $extrapolation --optimized 1 --openmp $openmp --v1 $v1 --v2 $v2 -R $R --prob $prob  --maxiter $maxiter --alpha_coeff $alpha_coeff --beta_coeff $beta_coeff --res_norm $res_norm --write_radii_angles 1 --f_grid_r "radii_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt" --f_grid_theta "angles_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt"
done


#################################################
##	TABLE 1-2: TEST CASES
#################################################
# module load mumps/mumps-5.4.1/mumps-5.4.1-gcc-11.3.0-openmpi-s7bj6hi
mkdir -p outputs

echo "prob alpha_coeff beta_coeff fac_ani extrapolation mod_pk"
# 1) Triangular/Czarny 2) Shafranov
for mod_pk in 2 1 # 2=Triangular/Czarny, 1=Shafranov
do
	# Cartesian + beta 0 + ani 0
	prob=7 # Solution (23) of Bourne et al.
		echo $prob $alpha_coeff $beta_coeff $fac_ani $extrapolation $mod_pk
		for extrapolation in 1
		do
			for divideBy2 in 0 1 2 3 4 5 6		#iterate over the different grid sizes
			do
			# note that the divideBy2 option here is only used as a dummy for looping. Grids need to be stored beforehand and are loaded here.
				echo "./${build_dir}/gmgpolar_simulation  -n "$nr_exp" -a "$fac_ani" --mod_pk "$mod_pk" --DirBC_Interior "$DirBC_Interior" --divideBy2 0 -r "$R0"  --smoother "$smoother" --verbose 2 --debug "$debug" --extrapolation "$extrapolation" --optimized 1 --openmp "$openmp" --v1 "$v1" --v2 "$v2" -R "$R" --prob "$prob"  --maxiter "$maxiter" --alpha_coeff "$alpha_coeff" --beta_coeff "$beta_coeff" --res_norm "$res_norm" --f_grid_r radii_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt --f_grid_theta angles_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt --rel_red_conv "$rel_red_conv" 1> outputs/job.out_"$fac_ani"_"$mod_pk"_"$prob"_"$beta_coeff"_"$extrapolation"_"$divideBy2"_"$rel_red_conv".txt 2> outputs/job.err_"$fac_ani"_"$mod_pk"_"$prob"_"$beta_coeff"_"$extrapolation"_"$divideBy2"_"$rel_red_conv".txt"
				./${build_dir}/gmgpolar_simulation  -n $nr_exp -a $fac_ani --mod_pk $mod_pk --DirBC_Interior $DirBC_Interior --divideBy2 0 -r $R0  --smoother $smoother --verbose 2 --debug $debug --extrapolation $extrapolation --optimized 1 --openmp $openmp --v1 $v1 --v2 $v2 -R $R --prob $prob  --maxiter $maxiter --alpha_coeff $alpha_coeff --beta_coeff $beta_coeff --res_norm $res_norm --f_grid_r "radii_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt" --f_grid_theta "angles_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt" --rel_red_conv $rel_red_conv 1> "outputs/job.out_"$fac_ani"_"$mod_pk"_"$prob"_"$beta_coeff"_"$extrapolation"_"$divideBy2"_"$rel_red_conv".txt" 2> "outputs/job.err_"$fac_ani"_"$mod_pk"_"$prob"_"$beta_coeff"_"$extrapolation"_"$divideBy2"_"$rel_red_conv".txt"
			done
		done

	# Polar + beta 0-1 + ani 0-1
	prob=6 # Solution (22) of Bourne et al.
		echo $prob $alpha_coeff $beta_coeff $fac_ani $extrapolation $mod_pk
		for extrapolation in 1
		do
			for divideBy2 in 0 1 2 3 4 5 6		#iterate over the different grid sizes
			do
			# note that the divideBy2 option here is only used as a dummy for looping. Grids need to be stored beforehand and are loaded here.
					echo "./${build_dir}/gmgpolar_simulation -n "$nr_exp" -a "$fac_ani" --mod_pk "$mod_pk" --DirBC_Interior "$DirBC_Interior" --divideBy2 0 -r "$R0"  --smoother "$smoother" --verbose 2 --debug "$debug" --extrapolation "$extrapolation" --optimized 1 --openmp "$openmp" --v1 "$v1" --v2 "$v2" -R "$R" --prob "$prob"  --maxiter "$maxiter" --alpha_coeff "$alpha_coeff" --beta_coeff "$beta_coeff" --res_norm "$res_norm" --f_grid_r radii_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt --f_grid_theta angles_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt --rel_red_conv "$rel_red_conv" 1> outputs/job.out_"$fac_ani"_"$mod_pk"_"$prob"_"$beta_coeff"_"$extrapolation"_"$divideBy2"_"$rel_red_conv".txt 2> outputs/job.err_"$fac_ani"_"$mod_pk"_"$prob"_"$beta_coeff"_"$extrapolation"_"$divideBy2"_"$rel_red_conv".txt"
				./${build_dir}/gmgpolar_simulation -n $nr_exp -a $fac_ani --mod_pk $mod_pk --DirBC_Interior $DirBC_Interior --divideBy2 0 -r $R0  --smoother $smoother --verbose 2 --debug $debug --extrapolation $extrapolation --optimized 1 --openmp $openmp --v1 $v1 --v2 $v2 -R $R --prob $prob  --maxiter $maxiter --alpha_coeff $alpha_coeff --beta_coeff $beta_coeff --res_norm $res_norm --f_grid_r "radii_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt" --f_grid_theta "angles_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt" --rel_red_conv $rel_red_conv 1> "outputs/job.out_"$fac_ani"_"$mod_pk"_"$prob"_"$beta_coeff"_"$extrapolation"_"$divideBy2"_"$rel_red_conv".txt" 2> "outputs/job.err_"$fac_ani"_"$mod_pk"_"$prob"_"$beta_coeff"_"$extrapolation"_"$divideBy2"_"$rel_red_conv".txt"
			done
		done
done
