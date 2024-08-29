#pragma once

class LevelCache;
class Level;

#include <chrono>
#include <vector>
#include <iostream>

#include "../PolarGrid/polargrid.h"

#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/densityProfileCoefficients.h"
#include "../InputFunctions/boundaryConditions.h"
#include "../InputFunctions/sourceTerm.h"
#include "../LinearAlgebra/vector.h"
#include "../LinearAlgebra/matrix.h"
#include "../LinearAlgebra/vector_operations.h"
#include "../common/constants.h"
#include "../Level/level.h"

class Residual{
public:
    explicit Residual(const PolarGrid& grid, const LevelCache& level_cache, 
                      const DomainGeometry& domain_geometry, const bool DirBC_Interior, const int num_omp_threads);
    
    void computeResidualTake0(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const;
    void computeResidual(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const;

private:
    /* ------------------- */
    /* Constructor members */
    const PolarGrid& grid_;
    const std::vector<double>& sin_theta_cache_;
    const std::vector<double>& cos_theta_cache_;
    const std::vector<double>& coeff_alpha_cache_;
    const std::vector<double>& coeff_beta_cache_;
    const DomainGeometry& domain_geometry_;
    const bool DirBC_Interior_;
    const int num_omp_threads_;

    void applyAGiveCircleSection(const int i_r, Vector<double>& result, const Vector<double>& x, const double& factor) const;
    void applyAGiveRadialSection(const int i_theta, Vector<double>& result, const Vector<double>& x, const double& factor) const;
};