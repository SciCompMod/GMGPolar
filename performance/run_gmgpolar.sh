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
R0=1e-5		        # r (1e-8/1e-5/1e-2)
# Generalized radius of maximum outer circle.
R=1.3			    # R
# Anisotropy in radial direction.
fac_ani=3		    # a
# TODO: which nr_exp and divideby2 do we want to consider?
nr_exp=4		    # n

#changing variables
mod_pk=1		    # mod_pk (0/1)
prob=5			    # prob
# set to something other than 0 here? iteration over 0, 1, 2, 3, ... in another benchmark?
divideBy2=0		    # divideBy2 (3/4/5/6)
# set to on
extrapolation=1		# E

nodes=1
ranks=1     # number of MPI Ranks
cores=14    # set OpenMP Num Threads to maximum number of cores requested

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
echo "#SBATCH -t 6000" >> run_gmgpolar_sbatch.sh
# echo "#SBATCH --exclusive" >> run_gmgpolar_sbatch.sh
 
# remove potentially loaded and conflicting modules
echo "module purge" >> run_gmgpolar_sbatch.sh

# gcc10
echo "module load PrgEnv/gcc10-openmpi" >> run_gmgpolar_sbatch.sh

echo "cd .." >> run_gmgpolar_sbatch.sh
echo "make -j16" >> run_gmgpolar_sbatch.sh

# reduce cores as cores count from 0
max_threads=$((cores-1))
# FLOPS-DP counter from 1 to cores many threads
echo "for m in {0..$max_threads}; do" >> run_gmgpolar_sbatch.sh
echo "likwid-perfctr -C 0-\$m -g FLOPS_DP ./build_gnu/main --openmp \$((m+1)) --matrix_free 1 -n $nr_exp -a $fac_ani --mod_pk $mod_pk --DirBC_Interior $DirBC_Interior --divideBy2 $divideBy2 -r $R0 --smoother $smoother -E $extrapolation" >> run_gmgpolar_sbatch.sh
echo "done;" >> run_gmgpolar_sbatch.sh

# # Memory (saturation) benchmarks
echo "for m in {0..$max_threads}; do" >> run_gmgpolar_sbatch.sh
echo "likwid-perfctr -C 0-\$m -g CACHES ./build_gnu/main --openmp \$((m+1)) --matrix_free 1 -n $nr_exp -a $fac_ani --mod_pk $mod_pk --DirBC_Interior $DirBC_Interior --divideBy2 $divideBy2 -r $R0 --smoother $smoother -E $extrapolation" >> run_gmgpolar_sbatch.sh
echo "done;" >> run_gmgpolar_sbatch.sh

#submit the job
sbatch run_gmgpolar_sbatch.sh