#!/bin/sh

module purge



export MODULEPATH=/scratch/spack-23.2/share/spack/modules/linux-ubuntu20.04-x86_64_v3:/tools/modulesystem/modulefiles
module load mumps


# /scratch/spack-23.2/opt/spack/linux-ubuntu20.04-x86_64_v3/gcc-12.3.0/mumps-5.5.1-afnceqpp75o4zmfmqpbvr4whi2li2r4


srun ./build/gmgpolar --nr_exp 10 --ntheta_exp 10 --openmp 1 --problem 0 --geometry 0 --alpha_coeff 0 --beta_coeff 1