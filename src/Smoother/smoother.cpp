#include "../../include/Smoother/smoother.h"

Smoother::Smoother(const PolarGrid& grid, const LevelCache& level_cache, 
                   const DomainGeometry& domain_geometry,
                   bool DirBC_Interior, int num_omp_threads
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

Smoother::~Smoother(){
    finalizeMumpsSolver(inner_boundary_mumps_solver_);
}

void Smoother::smoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp){
    smoothingInPlaceForLoop(x, rhs, temp); /* This is the fastest option */
    // smoothingInPlaceTaskLoop(x, rhs, temp);
    // smoothingInPlaceTaskDependencies(x, rhs, temp);
}