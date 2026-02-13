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
#include "../Definitions/global_definitions.h"
#include "../LinearAlgebra/Vector/vector.h"
#include "../LinearAlgebra/Vector/vector_operations.h"
#include "../LinearAlgebra/Solvers/tridiagonal_solver.h"
#include "../LinearAlgebra/Matrix/coo_matrix.h"
#include "../LinearAlgebra/Matrix/csr_matrix.h"
#include "../LinearAlgebra/Solvers/csr_lu_solver.h"
#include "../Stencil/stencil.h"

#ifdef GMGPOLAR_USE_MUMPS
    #include "dmumps_c.h"
    #include "mpi.h"s
#endif

#include "../LinearAlgebra/Solvers/symmetricTridiagonalSolver.h"

class Smoother
{
public:
    explicit Smoother(const PolarGrid& grid, const LevelCache& level_cache, const DomainGeometry& domain_geometry,
                      const DensityProfileCoefficients& density_profile_coefficients, bool DirBC_Interior,
                      int num_omp_threads);
    virtual ~Smoother() = default;

    virtual void smoothing(Vector<double> x, ConstVector<double> rhs, Vector<double> temp) = 0;

protected:
    const PolarGrid& grid_;
    const LevelCache& level_cache_;
    const DomainGeometry& domain_geometry_;
    const DensityProfileCoefficients& density_profile_coefficients_;
    const bool DirBC_Interior_;
    const int num_omp_threads_;
};
