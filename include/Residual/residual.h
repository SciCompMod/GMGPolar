#pragma once

#include <chrono>
#include <iostream>
#include <vector>

#include "../PolarGrid/polargrid.h"
#include "../Definitions/global_definitions.h"
#include "../LinearAlgebra/Vector/vector.h"
#include "../LinearAlgebra/Vector/vector_operations.h"

namespace gmgpolar
{

template <class LevelCacheType>
class Residual
{
public:
    explicit Residual(const PolarGrid<DefaultMemorySpace>& grid, const LevelCacheType& level_cache,
                      const bool DirBC_Interior)
        : grid_(grid)
        , level_cache_(level_cache)
        , DirBC_Interior_(DirBC_Interior)
    {
    }
    virtual ~Residual() = default;

    virtual void applySystemOperator(Vector<double> result, ConstVector<double> x) const = 0;
    virtual void computeResidual(HostVector<double> result, HostConstVector<double> rhs,
                                 HostConstVector<double> x) const                                = 0;

protected:
    /* ------------------- */
    /* Constructor members */
    const PolarGrid<DefaultMemorySpace> grid_;
    const LevelCacheType level_cache_;
    const bool DirBC_Interior_;
};
} // namespace gmgpolar
