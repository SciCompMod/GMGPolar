#include "../../../include/Smoother/SmootherGive/smootherGive.h"

// clang-format off
SmootherGive::SmootherGive(
    const PolarGrid& grid, const LevelCache& level_cache, 
    const DomainGeometry& domain_geometry, const DensityProfileCoefficients& density_profile_coefficients,
    bool DirBC_Interior, int num_omp_threads
) :
    Smoother(grid, level_cache, domain_geometry, density_profile_coefficients, DirBC_Interior, num_omp_threads) 
{
    buildAscMatrices();
    initializeMumpsSolver(inner_boundary_mumps_solver_, inner_boundary_circle_matrix_);    
}

SmootherGive::~SmootherGive(){
    finalizeMumpsSolver(inner_boundary_mumps_solver_);
}

void SmootherGive::smoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp){
    smoothingInPlaceForLoop(x, rhs, temp); /* This is the fastest option */
    // smoothingInPlaceTaskLoop(x, rhs, temp);
    // smoothingInPlaceTaskDependencies(x, rhs, temp);
}
// clang-format on