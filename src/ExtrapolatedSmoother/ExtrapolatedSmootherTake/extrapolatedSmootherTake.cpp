#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTake/extrapolatedSmootherTake.h"

ExtrapolatedSmootherTake::ExtrapolatedSmootherTake(const PolarGrid& grid, const LevelCache& level_cache,
                                                   const DomainGeometry& domain_geometry,
                                                   const DensityProfileCoefficients& density_profile_coefficients,
                                                   const bool DirBC_Interior, const int num_omp_threads)
    : ExtrapolatedSmoother(grid, level_cache, domain_geometry, density_profile_coefficients, DirBC_Interior,
                           num_omp_threads)
{
    buildAscMatrices();
#ifdef GMGPOLAR_USE_MUMPS
    initializeMumpsSolver(inner_boundary_mumps_solver_, inner_boundary_circle_matrix_);
#else
    inner_boundary_lu_solver_ = SparseLUSolver<double>(inner_boundary_circle_matrix_);
#endif
}

ExtrapolatedSmootherTake::~ExtrapolatedSmootherTake()
{
#ifdef GMGPOLAR_USE_MUMPS
    finalizeMumpsSolver(inner_boundary_mumps_solver_);
#endif
}
