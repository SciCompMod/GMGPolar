#include "../../include/DirectSolver/directSolver.h"

DirectSolver::DirectSolver(const PolarGrid& grid, const LevelCache& level_cache, 
    const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, 
    const int maxOpenMPThreads, const int openMPTaskThreads
) :
    grid_(grid), 
    sin_theta_(level_cache.sin_theta()),
    cos_theta_(level_cache.cos_theta()),
    domain_geometry_(domain_geometry),
    system_parameters_(system_parameters),
    DirBC_Interior_(DirBC_Interior),
    maxOpenMPThreads_(maxOpenMPThreads),
    openMPTaskThreads_(openMPTaskThreads)
{
    buildMatrixA(matrixA_);
    initializeMumps(mumps_, matrixA_);
}

DirectSolver::~DirectSolver() {
    deleteMumps(mumps_);
}

void DirectSolver::solveInPlace(Vector<double>& x) {
    subtractSymmetryShift(x);
    solveMumps(x);
}
