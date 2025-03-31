#include "../../../include/Smoother/SmootherTake/smootherTake.h"

SmootherTake::SmootherTake(const PolarGrid& grid, const LevelCache& level_cache, const DomainGeometry& domain_geometry,
                           const DensityProfileCoefficients& density_profile_coefficients, bool DirBC_Interior,
                           int num_omp_threads)
    : Smoother(grid, level_cache, domain_geometry, density_profile_coefficients, DirBC_Interior, num_omp_threads)
{
    buildAscMatrices();
#ifdef GMGPOLAR_USE_MUMPS
    initializeMumpsSolver(inner_boundary_mumps_solver_, inner_boundary_circle_matrix_);
#else
    inner_boundary_lu_solver_ = SparseLUSolver<double>(inner_boundary_circle_matrix_);
#endif
}

SmootherTake::~SmootherTake()
{
#ifdef GMGPOLAR_USE_MUMPS
    finalizeMumpsSolver(inner_boundary_mumps_solver_);
#endif
}
