#pragma once

#include <chrono>
#include <iostream>
#include <vector>

#include "../PolarGrid/polargrid.h"
#include "../Definitions/global_definitions.h"
#include "../LinearAlgebra/Vector/vector.h"
#include "../LinearAlgebra/Vector/vector_operations.h"
#include "../LinearAlgebra/Solvers/tridiagonal_solver.h"
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
class ExtrapolatedSmoother
{
public:
    explicit ExtrapolatedSmoother(const PolarGrid& grid, const LevelCacheType& level_cache, bool DirBC_Interior)
        : grid_(grid)
        , level_cache_(level_cache)
        , DirBC_Interior_(DirBC_Interior)
    {
    }
    virtual ~ExtrapolatedSmoother() = default;

    virtual void extrapolatedSmoothing(HostVector<double> x, HostConstVector<double> rhs, HostVector<double> temp) = 0;

protected:
    const PolarGrid& grid_;
    const LevelCacheType& level_cache_;
    const bool DirBC_Interior_;
};
} // namespace gmgpolar
