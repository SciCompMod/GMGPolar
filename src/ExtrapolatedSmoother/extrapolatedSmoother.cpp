#include "../../include/ExtrapolatedSmoother/extrapolated_smoother.h"

ExtrapolatedSmoother::ExtrapolatedSmoother(const PolarGrid& grid, const LevelCache& level_cache, 
                                           const DomainGeometry& domain_geometry, 
                                           const bool DirBC_Interior, const int num_omp_threads
) :
    grid_(grid), 
    sin_theta_cache_(level_cache.sin_theta()),
    cos_theta_cache_(level_cache.cos_theta()),
    coeff_alpha_cache_(level_cache.coeff_alpha()),
    coeff_beta_cache_(level_cache.coeff_beta()),
    domain_geometry_(domain_geometry),
    DirBC_Interior_(DirBC_Interior),
    num_omp_threads_(num_omp_threads)
{
    buildAscMatrices();
    initializeMumpsSolver(inner_boundary_mumps_solver_, inner_boundary_circle_matrix_);    
}

ExtrapolatedSmoother::~ExtrapolatedSmoother(){
    finalizeMumpsSolver(inner_boundary_mumps_solver_);
}

