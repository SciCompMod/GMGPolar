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
#include "../../LinearAlgebra/Vector/gpu_vector.h"
#include "../../common/constants.h"

class ResidualTakeGPU
{
public:
    explicit ResidualTakeGPU(const Level& level,
                      const DomainGeometry& domain_geometry,
                      const DensityProfileCoefficients& density_profile_coefficients,
                      const bool DirBC_Interior);
    ~ResidualTakeGPU() = default;

    void computeResidual(GPU_Vector<double>& result, const GPU_Vector<double>& rhs, const GPU_Vector<double>& x) const;

private:
    /* ------------------- */
    /* Constructor members */
    const Level& level_;
    const DomainGeometry& domain_geometry_;
    const DensityProfileCoefficients& density_profile_coefficients_;
    const bool DirBC_Interior_;
};