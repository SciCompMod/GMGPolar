#include "../../include/Interpolation/interpolation.h"

#include "../../include/LinearAlgebra/Vector/vector_operations.h"

#define FINE_NODE_EXTRAPOLATED_PROLONGATION() \
do { \
    if(i_r & 1) { \
        if(i_theta & 1) { \
            /* i_r % 2 == 1, i_theta % 2 == 1 */ \
            /* Fine node in the center of four coarse nodes */ \
            result[fineGrid.index(i_r, i_theta)] = 0.5 * ( \
                x[coarseGrid.index(i_r_coarse+1, i_theta_coarse)] + /* Bottom right */ \
                x[coarseGrid.index(i_r_coarse, i_theta_coarse+1)] /* Top left */ \
            ); \
        } \
        else { \
            /* i_r % 2 == 1, i_theta % 2 == 0 */ \
            /* Fine node between coarse nodes in radial direction */ \
            result[fineGrid.index(i_r, i_theta)] = 0.5 * ( \
                x[coarseGrid.index(i_r_coarse, i_theta_coarse)] + /* left */ \
                x[coarseGrid.index(i_r_coarse+1, i_theta_coarse)] /* right */ \
            ); \
        } \
    } \
    else { \
        if(i_theta & 1) { \
            /* i_r % 2 == 0, i_theta % 2 == 1 */ \
            /* Fine node between coarse nodes in theta direction */ \
            result[fineGrid.index(i_r, i_theta)] = 0.5 * ( \
                x[coarseGrid.index(i_r_coarse, i_theta_coarse)] + /* bottom */ \
                x[coarseGrid.index(i_r_coarse, i_theta_coarse+1)] /* top */ \
            ); \
        } \
        else { \
            /* i_r % 2 == 0, i_theta % 2 == 0 */ \
            /* Fine node appears in coarse grid */ \
            result[fineGrid.index(i_r, i_theta)] = x[coarseGrid.index(i_r_coarse, i_theta_coarse)]; /* center */ \
        } \
    } \
} while(0)

void Interpolation::applyExtrapolatedProlongation(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const
{
    assert(toLevel.level() == fromLevel.level() - 1);

    const PolarGrid& coarseGrid = fromLevel.grid();
    const PolarGrid& fineGrid = toLevel.grid();

    assert(x.size() == coarseGrid.numberOfNodes());
    assert(result.size() == fineGrid.numberOfNodes());

    #pragma omp parallel if (fineGrid.numberOfNodes() > 10'000)
    {
        /* Circluar Indexing Section */
        /* For loop matches circular access pattern */
        #pragma omp for nowait
        for (int i_r = 0; i_r < fineGrid.numberSmootherCircles(); i_r++)
        {
            int i_r_coarse = i_r >> 1;
            for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++)
            {
                int i_theta_coarse = i_theta >> 1;
                FINE_NODE_EXTRAPOLATED_PROLONGATION();
            }
        }

        /* Radial Indexing Section */
        /* For loop matches radial access pattern */
        #pragma omp for nowait
        for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++)
        {
            int i_theta_coarse = i_theta >> 1;
            for (int i_r = fineGrid.numberSmootherCircles(); i_r < fineGrid.nr(); i_r++)
            {
                int i_r_coarse = i_r >> 1;
                FINE_NODE_EXTRAPOLATED_PROLONGATION();
            }
        }
    }
}
