#include "../../include/Interpolation/interpolation.h"

void Interpolation::applyInjection(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const{
    assert(toLevel.level() == fromLevel.level() + 1);

    omp_set_num_threads(maxOpenMPThreads_);

    const PolarGrid& fineGrid = fromLevel.grid();
    const PolarGrid& coarseGrid = toLevel.grid();

    assert(x.size() == fineGrid.number_of_nodes());
    assert(result.size() == coarseGrid.number_of_nodes());

    const int coarseNumberSmootherCircles = coarseGrid.numberSmootherCircles();

    #pragma omp parallel num_threads(maxOpenMPThreads_) if(fineGrid.number_of_nodes() > 10'000)
    {
        /* For loop matches circular access pattern */
        #pragma omp for nowait
        for (int i_r_coarse = 0; i_r_coarse < coarseNumberSmootherCircles; i_r_coarse++){
            int i_r = i_r_coarse << 1;
            for (int i_theta_coarse = 0; i_theta_coarse < coarseGrid.ntheta(); i_theta_coarse++){
                int i_theta = i_theta_coarse << 1;
                result[coarseGrid.index(i_r_coarse, i_theta_coarse)] = x[fineGrid.index(i_r, i_theta)];
            }
        }

        /* For loop matches circular access pattern */
        #pragma omp for nowait
        for (int i_theta_coarse = 0; i_theta_coarse < coarseGrid.ntheta(); i_theta_coarse++){
            int i_theta = i_theta_coarse << 1;
            for (int i_r_coarse = coarseNumberSmootherCircles; i_r_coarse < coarseGrid.nr(); i_r_coarse++){
                int i_r = i_r_coarse << 1;
                result[coarseGrid.index(i_r_coarse, i_theta_coarse)] = x[fineGrid.index(i_r, i_theta)];
            }
        }
    }
}
