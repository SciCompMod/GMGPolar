#pragma once

#include <chrono>
#include <iostream>
#include <vector>

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

namespace gmgpolar
{

template <class LevelCacheType>
class DirectSolver
{
public:
    explicit DirectSolver(const PolarGrid& grid, const LevelCacheType& level_cache, bool DirBC_Interior)
        : grid_(grid)
        , level_cache_(level_cache)
        , DirBC_Interior_(DirBC_Interior)
    {
    }

    virtual ~DirectSolver() = default;

    // Note: The rhs (right-hand side) vector gets overwritten during the solution process.
    virtual void solveInPlace(Vector<double> solution) = 0;

    // Recomputes the solver matrix entries from the current LevelCache values
    // (arr, att, art, detDF, coeff_beta) and refreshes the factorization.
    // No reallocation: sparsity pattern and matrix storage are unchanged,
    // only the numeric values are overwritten.
    virtual void updateValues() = 0;

protected:
    const PolarGrid& grid_;
    const LevelCacheType& level_cache_;
    const bool DirBC_Interior_;
};
} // namespace gmgpolar
