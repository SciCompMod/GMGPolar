#include "../../include/Interpolation/interpolation.h"

/*
 * Extrapolated Restriction Operator
 * ----------------------------------
 *
 * This is the transpose of the extrapolated prolongation operator: R = P^T
 * Used between the finest grids in the multigrid hierarchy where uniform refinement is assumed.
 * All spacings are equal, so weights are simple factors of 0.5.
 * Each coarse node accumulates contributions from its corresponding 9 surrounding fine nodes.
 *
 * The stencil accumulates:
 *  - Center (weight 1.0)
 *  - Left, Right, Bottom, Top (weight 0.5 each)
 *  - Bottom-Right and Top-Left diagonal (weight 0.5 each)
 *
 * Note: Bottom-Left and Top-Right are NOT included (consistent with diagonal averaging in prolongation)
 *
 * Boundary handling:
 *  - Angular direction: periodic (wrapping)
 *  - Radial direction: check domain boundaries
 */

#define COARSE_NODE_EXTRAPOLATED_RESTRICTION()                                                                         \
    do {                                                                                                               \
        int i_r     = i_r_coarse * 2;                                                                                  \
        int i_theta = i_theta_coarse * 2;                                                                              \
                                                                                                                       \
        /* Center + Angular contributions (always present) */                                                          \
        double value = fine_values[fine_grid.index(i_r, i_theta)] +                                                    \
                       0.5 * fine_values[fine_grid.index(i_r, i_theta - 1)] +                                          \
                       0.5 * fine_values[fine_grid.index(i_r, i_theta + 1)];                                           \
                                                                                                                       \
        /* Left contributions (if not at inner boundary) */                                                            \
        if (i_r_coarse > 0) {                                                                                          \
            value += 0.5 * fine_values[fine_grid.index(i_r - 1, i_theta)] +                                            \
                     0.5 * fine_values[fine_grid.index(i_r - 1, i_theta + 1)]; /* Top-Left diagonal */                 \
        }                                                                                                              \
                                                                                                                       \
        /* Right contributions (if not at outer boundary) */                                                           \
        if (i_r_coarse < coarse_grid.nr() - 1) {                                                                       \
            value += 0.5 * fine_values[fine_grid.index(i_r + 1, i_theta)] +                                            \
                     0.5 * fine_values[fine_grid.index(i_r + 1, i_theta - 1)]; /* Bottom-Right diagonal */             \
        }                                                                                                              \
                                                                                                                       \
        coarse_result[coarse_grid.index(i_r_coarse, i_theta_coarse)] = value;                                          \
    } while (0)

void Interpolation::applyExtrapolatedRestriction(const PolarGrid& fine_grid, const PolarGrid& coarse_grid,
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
                COARSE_NODE_EXTRAPOLATED_RESTRICTION();
            }
        }

        /* Radial Indexing Section */
#pragma omp for nowait
        for (int i_theta_coarse = 0; i_theta_coarse < coarse_grid.ntheta(); i_theta_coarse++) {
            for (int i_r_coarse = coarse_grid.numberSmootherCircles(); i_r_coarse < coarse_grid.nr(); i_r_coarse++) {
                COARSE_NODE_EXTRAPOLATED_RESTRICTION();
            }
        }
    }
}
