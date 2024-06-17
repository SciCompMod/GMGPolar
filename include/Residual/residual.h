#pragma once

class LevelCache;
class Level;

#include "../PolarGrid/polargrid.h"

#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/systemParameters.h"

#include "../LinearAlgebra/vector.h"
#include "../LinearAlgebra/matrix.h"
#include "../LinearAlgebra/operations.h"

#include "../common/constants.h"

#include "../Level/level.h"

#include "../TaskDistribution/taskDistribution.h"

#include <chrono>
#include <vector>
#include <iostream>

class Residual{
public:
    explicit Residual(const PolarGrid& grid, const LevelCache& level_data, 
        const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, 
        const int maxOpenMPThreads, const int openMPTaskThreads
    );
    
    void computeResidual_V1(Vector<double>& result, const Vector<double>& x);
    void computeResidual_V2(Vector<double>& result, const Vector<double>& x);
    void computeResidual_V3(Vector<double>& result, const Vector<double>& x);

private:
    const PolarGrid& grid_;
    const std::vector<double>& sin_theta_;
    const std::vector<double>& cos_theta_;
    const DomainGeometry& domain_geometry_;
    const SystemParameters& system_parameters_;
    const bool DirBC_Interior_;
    const int maxOpenMPThreads_;
    const int openMPTaskThreads_;
};