#!/bin/bash

#fixed variables
origin_NOT_coarse=0	# origin_NOT_coarse
theta_aniso=0		# theta_aniso

# changing variables ?!
prob=5			# prob
R=1.3			# R
kappa_eps=0		# k
delta_e=0		# d
discr=3			# discr
fac_ani=3		# a
nr_exp=4		# n

#changing variables
mod_pk=0		# mod_pk (0/1)
R0=0.1			# r (1e-8/1e-5/1e-2)
DirBC_Interior=0	#  DirBC_Interior (0/1)
divideBy2=0		# divideBy2 (3/4/5/6)
smoother=3		# smoother (3,13)
extrapolation=0		# E

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
