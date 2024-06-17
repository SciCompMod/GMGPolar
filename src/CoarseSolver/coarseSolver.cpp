#include "../../include/CoarseSolver/coarseSolver.h"

CoarseSolver::CoarseSolver(const PolarGrid& grid, const LevelCache& level_data, 
    const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, 
    const int maxOpenMPThreads, const int openMPTaskThreads
) :
    grid_(grid), 
    sin_theta_(level_data.sin_theta()),
    cos_theta_(level_data.cos_theta()),
    domain_geometry_(domain_geometry),
    system_parameters_(system_parameters),
    DirBC_Interior_(DirBC_Interior),
    maxOpenMPThreads_(maxOpenMPThreads),
    openMPTaskThreads_(openMPTaskThreads)
{
    buildMatrixA(matrixA_);
    initializeMumps(mumps_, matrixA_);
}

CoarseSolver::~CoarseSolver() {
    deleteMumps(mumps_);
}

void CoarseSolver::solveInPlace(Vector<double>& x) {
    subtractSymmetryShift(x);
    solveMumps(x);
}