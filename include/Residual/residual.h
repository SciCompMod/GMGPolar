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

class Residual
{
public:
    explicit Residual(const PolarGrid& grid, const LevelCache& level_cache, const DomainGeometry& domain_geometry,
                      const DensityProfileCoefficients& density_profile_coefficients, const bool DirBC_Interior,
                      const int num_omp_threads);
    virtual ~Residual() = default;

    virtual void computeResidual(Vector<double> result, ConstVector<double> rhs, ConstVector<double> x) const = 0;

protected:
    /* ------------------- */
    /* Constructor members */
    const PolarGrid& grid_;
    const LevelCache& level_cache_;
    const DomainGeometry& domain_geometry_;
    const DensityProfileCoefficients& density_profile_coefficients_;
    const bool DirBC_Interior_;
    const int num_omp_threads_;
};
