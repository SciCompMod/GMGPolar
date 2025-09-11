#pragma once

class LevelCache;
class Level;

#include <chrono>
#include <iostream>
#include <vector>

#include "../PolarGrid/polargrid.h"

#include "../InputFunctions/boundaryConditions.h"
#include "../InputFunctions/densityProfileCoefficients.h"
#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/sourceTerm.h"
#include "../Level/level.h"
#include "../LinearAlgebra/vector.h"
#include "../LinearAlgebra/vector_operations.h"
#include "../common/global_definitions.h"

class Residual
{
public:
    explicit Residual(const PolarGrid& grid, const LevelCache& level_cache, const DomainGeometry& domain_geometry,
                      const DensityProfileCoefficients& density_profile_coefficients, const bool DirBC_Interior,
                      const int num_omp_threads);
    virtual ~Residual() = default;

    virtual void computeResidual(Vector<double> result,
                                 const Vector<double> rhs,
                                 const Vector<double> x) const = 0;

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
