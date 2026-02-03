#include "../../include/Interpolation/interpolation.h"

/* Remark: This injection is not scaled. */

inline void coarseNodeInjection(int i_r_coarse, int i_theta_coarse, const PolarGrid& fine_grid,
                                const PolarGrid& coarse_grid, Vector<double>& coarse_result,
                                ConstVector<double>& fine_values)
{
    int i_r     = i_r_coarse * 2;
    int i_theta = i_theta_coarse * 2;

    coarse_result[coarse_grid.index(i_r_coarse, i_theta_coarse)] = fine_values[fine_grid.index(i_r, i_theta)];
}

void Interpolation::applyInjection(const PolarGrid& fine_grid, const PolarGrid& coarse_grid,
                                   Vector<double> coarse_result, ConstVector<double> fine_values) const
{
    assert(std::ssize(fine_values) == fine_grid.numberOfNodes());
    assert(std::ssize(coarse_result) == coarse_grid.numberOfNodes());

#pragma omp parallel num_threads(max_omp_threads_)
    {
/* For loop matches circular access pattern */
#pragma omp for nowait
        for (int i_r_coarse = 0; i_r_coarse < coarse_grid.numberSmootherCircles(); i_r_coarse++) {
            for (int i_theta_coarse = 0; i_theta_coarse < coarse_grid.ntheta(); i_theta_coarse++) {
                coarseNodeInjection(i_r_coarse, i_theta_coarse, fine_grid, coarse_grid, coarse_result, fine_values);
            }
        }

/* For loop matches radial access pattern */
#pragma omp for nowait
        for (int i_theta_coarse = 0; i_theta_coarse < coarse_grid.ntheta(); i_theta_coarse++) {
            for (int i_r_coarse = coarse_grid.numberSmootherCircles(); i_r_coarse < coarse_grid.nr(); i_r_coarse++) {
                coarseNodeInjection(i_r_coarse, i_theta_coarse, fine_grid, coarse_grid, coarse_result, fine_values);
            }
        }
    }
}
