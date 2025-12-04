#include "../../include/Interpolation/interpolation.h"

/* Remark: This injection is not scaled. */

void Interpolation::applyInjection_(const PolarGrid& coarseGrid, const PolarGrid& fineGrid, Vector<double> result,
                                   ConstVector<double> x, int nthreads) const
{
    assert(x.size() == static_cast<uint>(fineGrid.numberOfNodes()));
    assert(result.size() == static_cast<uint>(coarseGrid.numberOfNodes()));

#pragma omp parallel num_threads(nthreads) if (fineGrid.numberOfNodes() > 10'000)
    {
/* For loop matches circular access pattern */
#pragma omp for nowait
        for (int i_r_coarse = 0; i_r_coarse < coarseGrid.numberSmootherCircles(); i_r_coarse++) {
            int i_r = i_r_coarse * 2;
            for (int i_theta_coarse = 0; i_theta_coarse < coarseGrid.ntheta(); i_theta_coarse++) {
                int i_theta                                          = i_theta_coarse * 2;
                result[coarseGrid.index(i_r_coarse, i_theta_coarse)] = x[fineGrid.index(i_r, i_theta)];
            }
        }

/* For loop matches circular access pattern */
#pragma omp for nowait
        for (int i_theta_coarse = 0; i_theta_coarse < coarseGrid.ntheta(); i_theta_coarse++) {
            int i_theta = i_theta_coarse * 2;
            for (int i_r_coarse = coarseGrid.numberSmootherCircles(); i_r_coarse < coarseGrid.nr(); i_r_coarse++) {
                int i_r                                              = i_r_coarse * 2;
                result[coarseGrid.index(i_r_coarse, i_theta_coarse)] = x[fineGrid.index(i_r, i_theta)];
            }
        }
    }
}
