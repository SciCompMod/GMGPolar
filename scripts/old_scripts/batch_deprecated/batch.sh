#!/bin/bash

#fixed variables
debug=0			# D
v1=1			# v1
v2=1			# v2
cycle=1			# c
compute_rho=0		# compute_rho
level=-1		# l
plotit=0		# P
solveit=1		# S
maxiter=150		# maxiter
periodic=1		# periodic
origin_NOT_coarse=0	# origin_NOT_coarse
theta_aniso=0		# theta_aniso
paraview=0		# paraview
prob=5			# prob
R=1.3			# R
kappa_eps=0		# k
delta_e=0		# d
discr=3			# discr
fac_ani=3		# a
nr_exp=4		# n
ntheta_exp=4		# ntheta_exp


#changing variables
mod_pk=0		# mod_pk (0/1)
R0=0.1			# r (1e-8/1e-5/1e-2)
DirBC_Interior=0	#  DirBC_Interior (0/1)
divideBy2=0		# divideBy2 (3/4/5/6)
smoother=3		# smoother (3,13)
extrapolation=0		# E



#overwrite existing file with an empty file
break > output.txt



for extrapol in 0 1
do

	for smoother in 3 13			#iterate over old/new version of smoothing
	do

		for R0 in 1e-5 1e-2 1e-1		#iterate over inner radius
		do

			for mod_pk in 0 1			#iterate over geometry
			do

				for DirBC_Interior in 1 0		#iterate over treatment of origin (across-the-origin/ Diriclet BC)
				do

					for divideBy2 in 0 1 2 3		#iterate over the different grid sizes
					do
	
						#call the multigrid algorithm, and append the output to a file
						./build/main -n $nr_exp -a $fac_ani --mod_pk $mod_pk --DirBC_Interior $DirBC_Interior --divideBy2 $divideBy2 -r $R0  --smoother $smoother -E $extrapol >> output.txt	
				
					done
				done
			done
		done
	done
done

