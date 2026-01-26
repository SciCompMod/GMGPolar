#include "../../include/Interpolation/interpolation.h"

/*
 * Restriction Operator
 * --------------------
 *
 * We use the transpose of the anisotropic bilinear interpolation stencil: R = P^T
 * For an isotropic (uniform) mesh, this stencil reduces to:
 *
 *           |1  2  1|
 *   R = 1/4 |2  4  2|
 *           |1  2  1|
 *
 * Each coarse node accumulates contributions from its corresponding 9 surrounding fine nodes.
 *
 * Weights are determined by the anisotropic mesh spacings:
 *  - h1, h2, h3, h4 in radial direction
 *  - k1, k2, k3, k4 in angular direction
 *
 * Boundary handling:
 *  - Angular direction: periodic (wrapping)
 *  - Radial direction: check domain boundaries (inner/outer radius)
 */

#define COARSE_NODE_RESTRICTION()                                                                                      \
    do {                                                                                                               \
        int i_r     = i_r_coarse * 2;                                                                                  \
        int i_theta = i_theta_coarse * 2;                                                                              \
                                                                                                                       \
        /* Angular indices with periodic wrapping */                                                                   \
        int i_theta_M2 = fine_grid.wrapThetaIndex(i_theta - 2);                                                        \
        int i_theta_M1 = fine_grid.wrapThetaIndex(i_theta - 1);                                                        \
        int i_theta_P1 = fine_grid.wrapThetaIndex(i_theta + 1);                                                        \
                                                                                                                       \
        /* Angular spacings */                                                                                         \
        double k1 = fine_grid.angularSpacing(i_theta_M2);                                                              \
        double k2 = fine_grid.angularSpacing(i_theta_M1);                                                              \
        double k3 = fine_grid.angularSpacing(i_theta);                                                                 \
        double k4 = fine_grid.angularSpacing(i_theta_P1);                                                              \
                                                                                                                       \
        /* Center + Angular contributions (always present) */                                                          \
        double value = fine_values[fine_grid.index(i_r, i_theta)] +                                                    \
                       k2 / (k1 + k2) * fine_values[fine_grid.index(i_r, i_theta_M1)] +                                \
                       k3 / (k3 + k4) * fine_values[fine_grid.index(i_r, i_theta_P1)];                                 \
                                                                                                                       \
        /* Left contributions (if not at inner boundary) */                                                            \
        if (i_r_coarse > 0) {                                                                                          \
            double h1 = fine_grid.radialSpacing(i_r - 2);                                                              \
            double h2 = fine_grid.radialSpacing(i_r - 1);                                                              \
            value += h2 / (h1 + h2) * fine_values[fine_grid.index(i_r - 1, i_theta)] +                                 \
                     h2 * k2 / ((h1 + h2) * (k1 + k2)) * fine_values[fine_grid.index(i_r - 1, i_theta_M1)] +           \
                     h2 * k3 / ((h1 + h2) * (k3 + k4)) * fine_values[fine_grid.index(i_r - 1, i_theta_P1)];            \
        }                                                                                                              \
                                                                                                                       \
        /* Right contributions (if not at outer boundary) */                                                           \
        if (i_r_coarse < coarse_grid.nr() - 1) {                                                                       \
            double h3 = fine_grid.radialSpacing(i_r);                                                                  \
            double h4 = fine_grid.radialSpacing(i_r + 1);                                                              \
            value += h3 / (h3 + h4) * fine_values[fine_grid.index(i_r + 1, i_theta)] +                                 \
                     h3 * k2 / ((h3 + h4) * (k1 + k2)) * fine_values[fine_grid.index(i_r + 1, i_theta_M1)] +           \
                     h3 * k3 / ((h3 + h4) * (k3 + k4)) * fine_values[fine_grid.index(i_r + 1, i_theta_P1)];            \
        }                                                                                                              \
                                                                                                                       \
        coarse_result[coarse_grid.index(i_r_coarse, i_theta_coarse)] = value;                                          \
    } while (0)

void Interpolation::applyRestriction(const PolarGrid& fine_grid, const PolarGrid& coarse_grid,
                                     Vector<double> coarse_result, ConstVector<double> fine_values) const
{
    assert(std::ssize(fine_values) == fine_grid.numberOfNodes());
    assert(std::ssize(coarse_result) == coarse_grid.numberOfNodes());

    /* We split the loops into two regions to better respect the */
    /* access patterns of the smoother and improve cache locality. */

#pragma omp parallel num_threads(max_omp_threads_)
    {
        /* Circular Indexing Section */
#pragma omp for nowait
        for (int i_r_coarse = 0; i_r_coarse < coarse_grid.numberSmootherCircles(); i_r_coarse++) {
            for (int i_theta_coarse = 0; i_theta_coarse < coarse_grid.ntheta(); i_theta_coarse++) {
                COARSE_NODE_RESTRICTION();
            }
        }

        /* Radial Indexing Section */
#pragma omp for nowait
        for (int i_theta_coarse = 0; i_theta_coarse < coarse_grid.ntheta(); i_theta_coarse++) {
            for (int i_r_coarse = coarse_grid.numberSmootherCircles(); i_r_coarse < coarse_grid.nr(); i_r_coarse++) {
                COARSE_NODE_RESTRICTION();
            }
        }
    }
}
