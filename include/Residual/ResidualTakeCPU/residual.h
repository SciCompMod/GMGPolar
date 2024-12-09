#pragma once

class LevelCache;
class Level;

#include <chrono>
#include <iostream>
#include <vector>

#include "../../PolarGrid/polargrid.h"

#include "../../InputFunctions/boundaryConditions.h"
#include "../../InputFunctions/densityProfileCoefficients.h"
#include "../../InputFunctions/domainGeometry.h"
#include "../../InputFunctions/sourceTerm.h"
#include "../../Level/level.h"
#include "../../LinearAlgebra/Vector/vector.h"
#include "../../LinearAlgebra/Vector/vector_operations.h"
#include "../../common/constants.h"

class ResidualTakeCPU
{
public:
    explicit ResidualTakeCPU(const Level& level, const DomainGeometry& domain_geometry,
                      const DensityProfileCoefficients& density_profile_coefficients, const bool DirBC_Interior);
    ~ResidualTakeCPU() = default;

    void computeResidual(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const;

private:
    /* ------------------- */
    /* Constructor members */
    const PolarGrid& grid_;
    const LevelCache& level_cache_;
    const DomainGeometry& domain_geometry_;
    const DensityProfileCoefficients& density_profile_coefficients_;
    const bool DirBC_Interior_;

    void applyCircleSection(const int i_r, Vector<double>& result, const Vector<double>& rhs,
                            const Vector<double>& x) const;
    void applyRadialSection(const int i_theta, Vector<double>& result, const Vector<double>& rhs,
                            const Vector<double>& x) const;
};
