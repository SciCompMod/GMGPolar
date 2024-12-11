#include "../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTakeCPU/extrapolatedSmoother.h"

ExtrapolatedSmootherTakeCPU::ExtrapolatedSmootherTakeCPU(const Level& level, const DomainGeometry& domain_geometry,
        const DensityProfileCoefficients& density_profile_coefficients, const bool DirBC_Interior)
    : grid_(level.grid())
    , level_cache_(level.levelCache())
    , domain_geometry_(domain_geometry)
    , density_profile_coefficients_(density_profile_coefficients)
    , DirBC_Interior_(DirBC_Interior)
{
    buildAscMatrices();
    initializeMumpsSolver(inner_boundary_mumps_solver_, inner_boundary_circle_matrix_);
}

ExtrapolatedSmootherTakeCPU::~ExtrapolatedSmootherTakeCPU()
{
    finalizeMumpsSolver(inner_boundary_mumps_solver_);
}
