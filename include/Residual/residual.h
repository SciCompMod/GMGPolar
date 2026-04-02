#pragma once

#include <chrono>
#include <iostream>
#include <vector>

#include "../InputFunctions/domainGeometry.h"

template <concepts::DomainGeometry DomainGeometry>
class LevelCache;

template <concepts::DomainGeometry DomainGeometry>
class Level;

#include "../InputFunctions/densityProfileCoefficients.h"
#include "../Level/level.h"
#include "../PolarGrid/polargrid.h"
#include "../Definitions/global_definitions.h"
#include "../LinearAlgebra/Vector/vector.h"
#include "../LinearAlgebra/Vector/vector_operations.h"

template <concepts::DomainGeometry DomainGeometry>
class Residual
{
public:
    explicit Residual(const PolarGrid& grid, const LevelCache<DomainGeometry>& level_cache,
                      const DensityProfileCoefficients& density_profile_coefficients, const bool DirBC_Interior,
                      const int num_omp_threads)
        : grid_(grid)
        , level_cache_(level_cache)
        , density_profile_coefficients_(density_profile_coefficients)
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
    const LevelCache<DomainGeometry>& level_cache_;
    const DensityProfileCoefficients& density_profile_coefficients_;
    const bool DirBC_Interior_;
    const int num_omp_threads_;
};
