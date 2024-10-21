#pragma once

class LevelCache;
class Level;

#include <chrono>
#include <vector>
#include <iostream>

#include "mpi.h" 
#include "dmumps_c.h"   

#include "../PolarGrid/polargrid.h"
#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/densityProfileCoefficients.h"
#include "../InputFunctions/boundaryConditions.h"
#include "../InputFunctions/sourceTerm.h"
#include "../LinearAlgebra/vector.h"
#include "../LinearAlgebra/matrix.h"
#include "../LinearAlgebra/vector_operations.h"
#include "../LinearAlgebra/symmetricTridiagonalSolver.h"
#include "../common/constants.h"
#include "../Level/level.h"
#include "../Stencil/stencil.h"

class Smoother {
public:
    explicit Smoother(
        const PolarGrid& grid, const LevelCache& level_cache, 
        const DomainGeometry& domain_geometry, const DensityProfileCoefficients& density_profile_coefficients,
        bool DirBC_Interior, int num_omp_threads
    );
    virtual ~Smoother() = default;

    virtual void smoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp) = 0;

protected:
    const PolarGrid& grid_;
    const LevelCache& level_cache_;
    const DomainGeometry& domain_geometry_;
    const DensityProfileCoefficients& density_profile_coefficients_;
    const bool DirBC_Interior_;
    const int num_omp_threads_;
};

