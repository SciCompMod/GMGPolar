#####################################
## DIRECTORIES
#####################################
#LIB_HOME = /home/leleux/Documents/CERFACS/lib
#LIB_HOME=/opt/homebrew/Cellar/scalapack/2.1.0_3/lib

#####################################
## COMPILER
#####################################
ERR = $(shell which icpc >/dev/null; echo $$?)

ifeq "$(ERR)" "0"
  COMPILER=intel
  CC = icc
  CXX = icpc
  LD = icpc
  #MT_FLAGS = -fopenmp -qoverride-limits -DUSE_MUMPS -I/scratch/ebourne/MUMPS_5.4.1/include
  MT_FLAGS = -fopenmp -qoverride-limits
  C_OTHERS = -std=c++14 -g
  L_OTHERS = -lz -Wwrite-strings
  #L_OTHERS = -lz -Wwrite-strings -L/scratch/ebourne/MUMPS_5.4.1/lib -L/scratch/ebourne/MUMPS_5.4.1/libseq/ -ldmumps -llapack -lblas -lifcore -lmumps_common -lpord -lmpiseq

else
  COMPILER=gnu
  CC = gcc
  CXX = g++
  LD = g++
  MT_FLAGS = -fopenmp
  C_OTHERS =
  L_OTHERS = -lgfortran -lz -Wwrite-strings
endif
