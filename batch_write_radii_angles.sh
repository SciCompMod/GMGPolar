#!/bin/bash

#fixed variables
debug=0
v1=1
v2=1
cycle=1
compute_rho=0
level=-1
plotit=0
solveit=1
maxiter=300
periodic=1
origin_NOT_coarse=0
theta_aniso=0
paraview=0
discr=1
nr_exp=4
ntheta_exp=4
res_norm=3
R0=1e-5
DirBC_Interior=1
smoother=3
rel_red_conv=1e-11
write_radii_angles=1

# Problem
prob=7
alpha_coeff=2		# Zoni shifted
beta_coeff=1

# geometry/grid
R=1.0
fac_ani=1
mod_pk=0
kappa_eps=0.3
delta_e=0.2
divideBy2=0

# MG
openmp=4
extrapolation=1

#f_grid_r="radii_files/aniso3/divide0.txt"
#f_grid_theta="angles_files/aniso3/divide0.txt"
f_grid_r=""
f_grid_theta=""

build_dir=build_gnu

#overwrite existing file with an empty file
break > output.txt

#################################################
##	TABLE 1-2: TEST CASES
#################################################
mkdir -p outputs

for R in 1.0 1.3
do
	for fac_ani in 0 1 3
	do
		for divideBy2 in 0 1 2 3 4 5 6 7 8
		do
			./build/main -n $nr_exp -a $fac_ani --mod_pk $mod_pk --DirBC_Interior $DirBC_Interior --divideBy2 $divideBy2 -r $R0  --smoother $smoother --verbose 0 --debug $debug --extrapolation $extrapolation --maxiter 150 --optimized 1 --openmp $openmp --v1 $v1 --v2 $v2 -R $R --prob $prob  --maxiter $maxiter --alpha_coeff $alpha_coeff --beta_coeff $beta_coeff --res_norm $res_norm --write_radii_angles $write_radii_angles --f_grid_r "radii_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt" --f_grid_theta "angles_files/Rmax"$R"/aniso"$fac_ani"/divide"$divideBy2".txt"
		done
	done
done
