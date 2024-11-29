#pragma once

class LevelCache;
class Level;

#include <chrono>
#include <iostream>
#include <vector>

#include "dmumps_c.h"
#include "mpi.h"

#include "../InputFunctions/domainGeometry.h"
#include "../Level/level.h"
#include "../LinearAlgebra/matrix.h"
#include "../LinearAlgebra/vector.h"
#include "../LinearAlgebra/vector_operations.h"
#include "../PolarGrid/polargrid.h"
#include "../Stencil/stencil.h"
#include "../common/constants.h"

class DirectSolver
{
public:
    explicit DirectSolver(const PolarGrid& grid,
                          const LevelCache& level_cache,
                          const DomainGeometry& domain_geometry,
                          const DensityProfileCoefficients& density_profile_coefficients,
                          bool DirBC_Interior,
                          int num_omp_threads);

    virtual ~DirectSolver() = default;

    virtual void solveInPlace(Vector<double>& solution) = 0;

protected:
    const PolarGrid& grid_;
    const LevelCache& level_cache_;
    const DomainGeometry& domain_geometry_;
    const DensityProfileCoefficients& density_profile_coefficients_;
    const bool DirBC_Interior_;
    const int num_omp_threads_;
};
