#pragma once

#include <chrono>
#include <iostream>
#include <vector>

#include "../InputFunctions/domainGeometry.h"

template <concepts::DomainGeometry DomainGeometry>
class LevelCache;
template <concepts::DomainGeometry DomainGeometry>
class Level;

#include "../Level/level.h"
#include "../PolarGrid/polargrid.h"
#include "../Definitions/global_definitions.h"
#include "../LinearAlgebra/Vector/vector.h"
#include "../LinearAlgebra/Vector/vector_operations.h"
#include "../LinearAlgebra/Matrix/coo_matrix.h"
#include "../LinearAlgebra/Matrix/csr_matrix.h"
#include "../LinearAlgebra/Solvers/csr_lu_solver.h"
#include "../LinearAlgebra/Solvers/coo_mumps_solver.h"
#include "../Stencil/stencil.h"

#ifdef GMGPOLAR_USE_MUMPS
    #include "dmumps_c.h"
    #include "mpi.h"
#endif

template <concepts::DomainGeometry DomainGeometry>
class DirectSolver
{
public:
    explicit DirectSolver(const PolarGrid& grid, const LevelCache<DomainGeometry>& level_cache, bool DirBC_Interior,
                          int num_omp_threads)
        : grid_(grid)
        , level_cache_(level_cache)
        , DirBC_Interior_(DirBC_Interior)
        , num_omp_threads_(num_omp_threads)
    {
    }

    virtual ~DirectSolver() = default;

    // Note: The rhs (right-hand side) vector gets overwritten during the solution process.
    virtual void solveInPlace(Vector<double> solution) = 0;

protected:
    const PolarGrid& grid_;
    const LevelCache<DomainGeometry>& level_cache_;
    const bool DirBC_Interior_;
    const int num_omp_threads_;
};
