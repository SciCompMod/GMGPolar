#include "../../include/Interpolation/interpolation.h"

#include "../../include/LinearAlgebra/Vector/vector_operations.h"

// We use the anisotropic bilinear interpolation stencil.
// For an isotropic mesh, this stencil reduces to
//           |1  2  1|
// P = 1/4 * |2  4  2|
//           |1  2  1|

#define FINE_NODE_PROLONGATION() \
do { \
    if(i_r & 1) { \
        if(i_theta & 1) { \
            /* i_r % 2 == 1, i_theta % 2 == 1 */ \
            /* Fine node in the center of four coarse nodes */ \
            double h1 = fineGrid.radialSpacing(i_r-1); \
            double h2 = fineGrid.radialSpacing(i_r); \
            double k1 = fineGrid.angularSpacing(i_theta-1); \
            double k2 = fineGrid.angularSpacing(i_theta); \
            int i_theta_coarse_P1 = coarseGrid.wrapThetaIndex(i_theta_coarse+1); \
            double divisor = (h1+h2) * (k1+k2); \
            double value = ( \
                h1*k1*x[coarseGrid.index(i_r_coarse, i_theta_coarse)] +  /* Bottom left */ \
                h2*k1*x[coarseGrid.index(i_r_coarse+1, i_theta_coarse)] +  /* Bottom right */ \
                h1*k2*x[coarseGrid.index(i_r_coarse, i_theta_coarse_P1)] + /* Top left */ \
                h2*k2*x[coarseGrid.index(i_r_coarse+1, i_theta_coarse_P1)]   /* Top right */ \
            ); \
            result[fineGrid.index(i_r, i_theta)] = value / divisor; \
        } \
        else { \
            /* i_r % 2 == 1, i_theta % 2 == 0 */ \
            /* Fine node between coarse nodes in radial direction */ \
            double h1 = fineGrid.radialSpacing(i_r-1); \
            double h2 = fineGrid.radialSpacing(i_r); \
            double divisor = (h1+h2); \
            double value = ( \
                h1*x[coarseGrid.index(i_r_coarse, i_theta_coarse)] +  /* left */ \
                h2*x[coarseGrid.index(i_r_coarse+1, i_theta_coarse)]    /* right */ \
            ); \
            result[fineGrid.index(i_r, i_theta)] = value / divisor; \
        } \
    } \
    else { \
        if(i_theta & 1) { \
            /* i_r % 2 == 0, i_theta % 2 == 1 */ \
            /* Fine node between coarse nodes in theta direction */ \
            double k1 = fineGrid.angularSpacing(i_theta-1); \
            double k2 = fineGrid.angularSpacing(i_theta); \
            int i_theta_coarse_P1 = coarseGrid.wrapThetaIndex(i_theta_coarse+1); \
            double divisor = (k1+k2); \
            double value = ( \
                k1*x[coarseGrid.index(i_r_coarse, i_theta_coarse)] + /* bottom */ \
                k2*x[coarseGrid.index(i_r_coarse, i_theta_coarse_P1)] /* top */ \
            ); \
            result[fineGrid.index(i_r, i_theta)] = value / divisor; \
        } \
        else { \
            /* i_r % 2 == 0, i_theta % 2 == 0 */ \
            /* Fine node appears in coarse grid */ \
            result[fineGrid.index(i_r, i_theta)] = x[coarseGrid.index(i_r_coarse, i_theta_coarse)]; /* center */ \
        } \
    } \
} while(0)

void Interpolation::applyProlongation(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const
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
                FINE_NODE_PROLONGATION();
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
                FINE_NODE_PROLONGATION();
            }
        }
    }
}
