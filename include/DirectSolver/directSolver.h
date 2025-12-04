#pragma once

class LevelCache;
class Level;

#include <chrono>
#include <iostream>
#include <vector>

#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/densityProfileCoefficients.h"
#include "../Level/level.h"
#include "../PolarGrid/polargrid.h"
#include "../common/global_definitions.h"
#include "../LinearAlgebra/vector.h"
#include "../LinearAlgebra/vector_operations.h"
#include "../LinearAlgebra/coo_matrix.h"
#include "../LinearAlgebra/csr_matrix.h"
#include "../LinearAlgebra/sparseLUSolver.h"
#include "../Stencil/stencil.h"

class DirectSolver
{
public:
    explicit DirectSolver(const PolarGrid& grid, const LevelCache& level_cache, const DomainGeometry& domain_geometry,
                          const DensityProfileCoefficients& density_profile_coefficients, bool DirBC_Interior,
                          int num_omp_threads);

    virtual ~DirectSolver() = default;

    // Note: The rhs (right-hand side) vector gets overwritten during the solution process.
    virtual void solveInPlace(Vector<double> solution) const = 0;

protected:
    const PolarGrid& grid_;
    const LevelCache& level_cache_;
    const DomainGeometry& domain_geometry_;
    const DensityProfileCoefficients& density_profile_coefficients_;
    const bool DirBC_Interior_;
    const int num_omp_threads_;
};
