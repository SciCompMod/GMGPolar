#SBATCH --job-name=gmgpolar
#SBATCH --output=slurm-%A-convergenceorder_GMGPolar_CARA.out
#SBATCH --error=slurm-%A-convergenceorder_GMGPolar_CARA.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH --threads-per-core=1
#SBATCH --time=3:30:00
#SBATCH --exclusive
#SBATCH --partition=naples128
#SBATCH --account=2476029

# Adjust parameters in src/weak_scaling.cpp

srun --cpus-per-task=56 ./../../../build/paper_weak_scaling
