#!/bin/bash
#SBATCH --job-name=Seq              # name of the job
#SBATCH --nodes=1                   # number of nodes
#SBATCH --ntasks-per-node=1         # number of MPI tasks per node
#SBATCH --time=00:03:00             # maximum execution time requested (HH:MM:SS)
#SBATCH --output=Seq%j.out          # name of output file
#SBATCH --error=Seq%j.out           # name of error file (here, in common with the output)
#SBATCH --account=emp@cpu
# #SBATCH --partition=gpu_p2

#cd ${SLURM_SUBMIT_DIR}

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
divideBy2=3		# divideBy2 (3/4/5/6)
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
echo "#SBATCH -t 25" >> run_gmgpolar_sbatch.sh
echo "#SBATCH --exclusive" >> run_gmgpolar_sbatch.sh
echo "#SBATCH --account=emp@cpu" >> run_gmgpolar_sbatch.sh
# remove potentially loaded and conflicting modules
echo "module purge" >> run_gmgpolar_sbatch.sh

# gcc10
echo "module load openmpi/4.1.1-cuda" >> run_gmgpolar_sbatch.sh
echo "module load cmake" >> run_gmgpolar_sbatch.sh

#echo "export OMP_PLACES=cores" >> run_gmgpolar_sbatch.sh
echo "export OMP_NUM_THREADS=14" >> run_gmgpolar_sbatch.sh

echo "cd $SCRATCH/GMGPolar/build" >> run_gmgpolar_sbatch.sh
echo "cmake --build . -- -j16" >> run_gmgpolar_sbatch.sh

# FLOPS-DP counter
echo "./gmgpolar_simulation --openmp 14  --matrix_free 1 -n $nr_exp -a $fac_ani --mod_pk $mod_pk --DirBC_Interior $DirBC_Interior --divideBy2 $divideBy2 -r $R0  --smoother $smoother -E $extrapolation" >> run_gmgpolar_sbatch.sh

#submit the job
sbatch run_gmgpolar_sbatch.sh

