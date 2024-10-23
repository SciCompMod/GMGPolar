#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherGive/extrapolatedSmootherGive.h"

// clang-format off
ExtrapolatedSmootherGive::ExtrapolatedSmootherGive(
    const PolarGrid& grid, const LevelCache& level_cache, 
    const DomainGeometry& domain_geometry, const DensityProfileCoefficients& density_profile_coefficients,
    const bool DirBC_Interior, const int num_omp_threads
) :
    ExtrapolatedSmoother(grid, level_cache, domain_geometry, density_profile_coefficients, DirBC_Interior, num_omp_threads) 
{
    buildAscMatrices();
    initializeMumpsSolver(inner_boundary_mumps_solver_, inner_boundary_circle_matrix_);
}

ExtrapolatedSmootherGive::~ExtrapolatedSmootherGive()
{
    finalizeMumpsSolver(inner_boundary_mumps_solver_);
}

void ExtrapolatedSmootherGive::extrapolatedSmoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp)
{
    extrapolatedSmoothingInPlaceForLoop(x, rhs, temp); /* This is the fastest option */
    // extrapolatedSmoothingInPlaceTaskLoop(x, rhs, temp);
    // extrapolatedSmoothingInPlaceTaskDependencies(x, rhs, temp);
}