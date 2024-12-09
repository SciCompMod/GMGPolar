#include "../../include/Interpolation/interpolation.h"

#include "../../include/LinearAlgebra/Vector/vector_operations.h"

/* Remark: This injection is not scaled. */

void Interpolation::applyInjection(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const
{
    assert(toLevel.level() == fromLevel.level() + 1);

    const PolarGrid& fineGrid = fromLevel.grid();
    const PolarGrid& coarseGrid = toLevel.grid();

    assert(x.size() == fineGrid.numberOfNodes());
    assert(result.size() == coarseGrid.numberOfNodes());

    #pragma omp parallel if (fineGrid.numberOfNodes() > 10'000)
    {
        /* For loop matches circular access pattern */
        #pragma omp for nowait
        for (int i_r_coarse = 0; i_r_coarse < coarseGrid.numberSmootherCircles(); i_r_coarse++)
        {
            int i_r = i_r_coarse << 1;
            for (int i_theta_coarse = 0; i_theta_coarse < coarseGrid.ntheta(); i_theta_coarse++)
            {
                int i_theta = i_theta_coarse << 1;
                result[coarseGrid.index(i_r_coarse, i_theta_coarse)] = x[fineGrid.index(i_r, i_theta)];
            }
        }

        /* For loop matches circular access pattern */
        #pragma omp for nowait
        for (int i_theta_coarse = 0; i_theta_coarse < coarseGrid.ntheta(); i_theta_coarse++)
        {
            int i_theta = i_theta_coarse << 1;
            for (int i_r_coarse = coarseGrid.numberSmootherCircles(); i_r_coarse < coarseGrid.nr(); i_r_coarse++)
            {
                int i_r = i_r_coarse << 1;
                result[coarseGrid.index(i_r_coarse, i_theta_coarse)] = x[fineGrid.index(i_r, i_theta)];
            }
        }
    }
}
