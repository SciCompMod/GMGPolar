#!/bin/bash

./build/gmgpolar \
\
--maxOpenMPThreads 40 \
--finestLevelThreads -1 \
--threadReductionFactor 1.0 \
\
--R0 1e-5 \
--Rmax 1.3 \
--nr_exp 5 \
--ntheta_exp 5 \
--anisotropic_factor 4 \
--divideBy2 5 \
--write_grid_file 0 \
--load_grid_file 0 \
--file_grid_r "_radii.txt" \
--file_grid_theta "_angles.txt" \
\
--DirBC_Interior 1 \
\
--extrapolation 0 \
--maxLevels -1 \
--v1 1 \
--v2 1 \
--cycle 1 \
# --help;