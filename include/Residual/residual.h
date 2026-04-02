#pragma once

#include <chrono>
#include <iostream>
#include <vector>

#include "../PolarGrid/polargrid.h"
#include "../Definitions/global_definitions.h"
#include "../LinearAlgebra/Vector/vector.h"
#include "../LinearAlgebra/Vector/vector_operations.h"

template <class LevelCacheType>
class Residual
{
public:
    explicit Residual(const PolarGrid& grid, const LevelCacheType& level_cache, const bool DirBC_Interior,
                      const int num_omp_threads)
        : grid_(grid)
        , level_cache_(level_cache)
        , DirBC_Interior_(DirBC_Interior)
        , num_omp_threads_(num_omp_threads)
    {
    }
    virtual ~Residual() = default;

    virtual void computeResidual(Vector<double> result, ConstVector<double> rhs, ConstVector<double> x) const = 0;

    virtual void applySystemOperator(Vector<double> result, ConstVector<double> x) const = 0;

protected:
    /* ------------------- */
    /* Constructor members */
    const PolarGrid& grid_;
    const LevelCacheType& level_cache_;
    const bool DirBC_Interior_;
    const int num_omp_threads_;
};
