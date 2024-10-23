#pragma once

class LevelCache;
class Level;

#include <chrono>
#include <iostream>
#include <vector>

#include "dmumps_c.h"
#include "mpi.h"

#include "../InputFunctions/domainGeometry.h"
#include "../LinearAlgebra/diagonalSolver.h"
#include "../LinearAlgebra/matrix.h"
#include "../LinearAlgebra/symmetricTridiagonalSolver.h"
#include "../LinearAlgebra/vector.h"
#include "../LinearAlgebra/vector_operations.h"
#include "../PolarGrid/polargrid.h"

#include "../Level/level.h"
#include "../Stencil/stencil.h"
#include "../common/constants.h"

class ExtrapolatedSmoother
{
public:
    explicit ExtrapolatedSmoother(const PolarGrid& grid,
                                  const LevelCache& level_cache,
                                  const DomainGeometry& domain_geometry,
                                  const DensityProfileCoefficients& density_profile_coefficients,
                                  bool DirBC_Interior,
                                  int num_omp_threads);
    virtual ~ExtrapolatedSmoother() = default;

    virtual void extrapolatedSmoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp) = 0;

protected:
    const PolarGrid& grid_;
    const LevelCache& level_cache_;
    const DomainGeometry& domain_geometry_;
    const DensityProfileCoefficients& density_profile_coefficients_;
    const bool DirBC_Interior_;
    const int num_omp_threads_;
};
