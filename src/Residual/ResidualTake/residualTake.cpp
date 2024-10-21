#include "../../../include/Residual/ResidualTake/residualTake.h"


ResidualTake::ResidualTake(
    const PolarGrid& grid, const LevelCache& level_cache, 
    const DomainGeometry& domain_geometry, const DensityProfileCoefficients& density_profile_coefficients,
    bool DirBC_Interior, int num_omp_threads
) :
    Residual(grid, level_cache, domain_geometry, density_profile_coefficients, DirBC_Interior, num_omp_threads)
{}


/* ------------------ */
/* result = rhs - A*x */
void ResidualTake::computeResidual(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const {
    assert(result.size() == x.size());
    omp_set_num_threads(num_omp_threads_);

    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    if(omp_get_max_threads() == 1){
        /* Single-threaded execution */
        for(int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
            applyCircleSection(i_r, result, rhs, x);
        }
        for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
            applyRadialSection(i_theta, result, rhs, x);
        }
    }
    else{
        /* Multi-threaded execution */
        #pragma omp parallel
        {
            /* Circle Section */
            #pragma omp for nowait
            for(int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
                applyCircleSection(i_r, result, rhs, x);
            }
            /* Radial Section */
            #pragma omp for nowait
            for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
                applyRadialSection(i_theta, result, rhs, x);
            }
        }
    }
}