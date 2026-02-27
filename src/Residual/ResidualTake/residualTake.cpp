#include "../../../include/Residual/ResidualTake/residualTake.h"

ResidualTake::ResidualTake(const PolarGrid& grid, const LevelCache& level_cache, const DomainGeometry& domain_geometry,
                           const DensityProfileCoefficients& density_profile_coefficients, bool DirBC_Interior,
                           int num_omp_threads)
    : Residual(grid, level_cache, domain_geometry, density_profile_coefficients, DirBC_Interior, num_omp_threads)
{
}

/* ------------ */
/* result = A*x */
void ResidualTake::applySystemOperator(Vector<double> result, ConstVector<double> x) const
{
    assert(result.size() == x.size());

    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

#pragma omp parallel num_threads(num_omp_threads_)
    {
/* Circle Section */
#pragma omp for nowait
        for (int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
            applyCircleSection(i_r, result, x);
        }

/* Radial Section */
#pragma omp for nowait
        for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
            applyRadialSection(i_theta, result, x);
        }
    }
}

/* ------------------ */
/* result = rhs - A*x */
void ResidualTake::computeResidual(Vector<double> result, ConstVector<double> rhs, ConstVector<double> x) const
{
    assert(result.size() == x.size());

    applySystemOperator(result, x);

    // Subtract A*x from rhs to get the residual.
    const int n = result.size();
#pragma omp parallel for num_threads(num_omp_threads_)
    for (int i = 0; i < n; i++) {
        result[i] = rhs[i] - result[i];
    }
}