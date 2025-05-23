#include "../../include/ExtrapolatedSmoother/extrapolatedSmoother.h"

ExtrapolatedSmoother::ExtrapolatedSmoother(const PolarGrid& grid, const LevelCache& level_cache,
                                           const DomainGeometry& domain_geometry,
                                           const DensityProfileCoefficients& density_profile_coefficients,
                                           const bool DirBC_Interior, const int num_omp_threads)
    : grid_(grid)
    , level_cache_(level_cache)
    , domain_geometry_(domain_geometry)
    , density_profile_coefficients_(density_profile_coefficients)
    , DirBC_Interior_(DirBC_Interior)
    , num_omp_threads_(num_omp_threads)
{
}
