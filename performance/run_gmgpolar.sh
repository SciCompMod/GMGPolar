#!/bin/bash

#fixed variables
origin_NOT_coarse=0	# origin_NOT_coarse  (not necessary!!! never used in C++ implementation)
theta_aniso=0		# theta_aniso (we never have an ansisotropy in theta-direction, so the value is always 0/false)

# changing variables ?!
prob=5			# prob (our test problem is always problem 5)
R=1.3			# R (outer radius, we always use R=1.3, which is also the default value)
kappa_eps=0		# k (Parameter kappa for the Shafranov geometry, parameters epsilon for the Czarny geometry, is set automatically depending on the geometry)
delta_e=0		# d (Parameterdelta for the Shafranov geometry, parameters e for the Czarny geometry, is set automatically depending on the geometry)
discr=3			# discr
fac_ani=3		# a (for the anisotropy in r-direction, number of refinements around r=2/3 R)
nr_exp=4		# n  (defines the number of nodes in r-direction: 
                #without any ansisotropy (fac_ani=0), we have nr=2^(nr_exp - 1) equally distributed nodes,
                #with anisotropy: we fisrt create nr=2^(nr_exp) - 1^(fac_ani) equally distributed nodes, and then, the grid is refined fac_ani-times,
                #Lastly, regardless of using any anisotropy or not, we refine again by splitting at the midpoint of all intervals 
                #(for a uniform grid refinement between the two finest levels for a successful extrapolation)
                #ntheta is then finally calculated as: ntheta= 2^(ceil(log_2(nr)))

#changing variables
mod_pk=0		# mod_pk (0/1) (Defines the geometry, 0: circular geometry, 1: Shafranov geometry=default, 2: Czarny geometry)
R0=0.1			# r (1e-8/1e-5/1e-2) (inner radius)
DirBC_Interior=0	#  DirBC_Interior (0/1) (0:across the origin discretization, 1: Dirichlet BC)
divideBy2=0		# divideBy2 (3/4/5/6) (defines how often to split the intervals of the grid at the midpoint)
smoother=3		# smoother (3,13) (3: coupled circle-radial version, 13: decoupled circle-radial version, we always use smoother=3 as the other version was not faster)
extrapolation=0		# E (0: no extrapolation, 1: implicit extrapolation, 2: implicit extrapolation with full grid smoothing)

#others
#cycle=1  #1:V-cycle, 2: W-cycle, we always use the V-cycle
#maxiter  #maximum number of iterations
#v1=1,v2=1    #number of pre-/post-smoothing steps
#optimised=1  #we always use the optimised version of the code
#matrix_free=1  #we always use the matrix-free version of the code



echo "#!/bin/bash" > run_gmgpolar_sbatch.sh
# create a short name for your job
echo "#SBATCH --job-name=gmgpolar" >> run_gmgpolar_sbatch.sh
# stdout file %A=job id
echo "#SBATCH --output=slurm-%A-p$prob-fa$fac_ani-r$nr_exp-mpk$mod_pk-s$smoother-e$extrapolation.out" >> run_gmgpolar_sbatch.sh
# stderr file
echo "#SBATCH --error=slurm-%A-p$prob-fa$fac_ani-r$nr_exp-mpk$mod_pk-s$smoother-e$extrapolation.err" >> run_gmgpolar_sbatch.sh

echo "#SBATCH -N 1" >> run_gmgpolar_sbatch.sh
echo "#SBATCH -n 1" >> run_gmgpolar_sbatch.sh
echo "#SBATCH -c 14" >> run_gmgpolar_sbatch.sh
echo "#SBATCH -t 6000" >> run_gmgpolar_sbatch.sh
echo "#SBATCH --exclusive" >> run_gmgpolar_sbatch.sh
 
# remove potentially loaded and conflicting modules
echo "module purge" >> run_gmgpolar_sbatch.sh

# gcc10
echo "module load PrgEnv/gcc10-openmpi" >> run_gmgpolar_sbatch.sh

echo "cd .." >> run_gmgpolar_sbatch.sh
echo "make -j16" >> run_gmgpolar_sbatch.sh

# FLOPS-DP counter
echo "for m in {0..1}; do" >> run_gmgpolar_sbatch.sh
echo "likwid-perfctr -C 0-$m -g FLOPS_DP ./build_gnu/main --matrix_free 1 -n $nr_exp -a $fac_ani --mod_pk $mod_pk --DirBC_Interior $DirBC_Interior --divideBy2 $divideBy2 -r $R0  --smoother $smoother -E $extrapolation" >> run_gmgpolar_sbatch.sh
echo "done;" >> run_gmgpolar_sbatch.sh

# memory (saturation) benchmarks
echo "for m in {0..1}; do" >> run_gmgpolar_sbatch.sh
echo "likwid-perfctr -C 0-$m -g CACHES ./build_gnu/main -n $nr_exp -a $fac_ani --mod_pk $mod_pk --DirBC_Interior $DirBC_Interior --divideBy2 $divideBy2 -r $R0  --smoother $smoother -E $extrapol" >> run_gmgpolar_sbatch.sh
echo "done;" >> run_gmgpolar_sbatch.sh

#submit the job
sbatch run_gmgpolar_sbatch.sh
