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
R0=1e-5			# r (1e-8/1e-5/1e-2)
DirBC_Interior=0	#  DirBC_Interior (0/1)
divideBy2=0		# divideBy2 (3/4/5/6)
smoother=3		# smoother (3,13)
extrapol=0		# E

verb=0



#overwrite existing file with an empty file
break > output_omp.txt




for mod_pk in 1				#iterate over geometry
do

	for DirBC_Interior in 1 		#iterate over treatment of origin (across-the-origin/ Diriclet BC)
	do

		for divideBy2 in 3			#iterate over the different grid sizes
		do

			for omp in 1 2 4
			do

				for n in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20	#just do it several times to get an average value
				do
	
					#call the multigrid algorithm, and append the output to a file
					./build/main -n $nr_exp -a $fac_ani --mod_pk $mod_pk --DirBC_Interior $DirBC_Interior --divideBy2 $divideBy2 -r $R0  --smoother $smoother -E $extrapol --openmp $omp --verbose $verb >> output_omp.txt	

				done
			done
		done
	done
done

