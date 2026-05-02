#include "../../include/Interpolation/interpolation.h"
using namespace gmgpolar;

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

static inline void coarseNodeRestriction(const int i_r_coarse, const int i_theta_coarse, const PolarGrid& fine_grid,
                                         const PolarGrid& coarse_grid, Vector<double>& coarse_result,
                                         ConstVector<double>& fine_values)
{
    const int i_r     = i_r_coarse * 2;
    const int i_theta = i_theta_coarse * 2;

    /* Angular indices with periodic wrapping */
    const int i_theta_M2 = fine_grid.wrapThetaIndex(i_theta - 2);
    const int i_theta_M1 = fine_grid.wrapThetaIndex(i_theta - 1);
    const int i_theta_P1 = fine_grid.wrapThetaIndex(i_theta + 1);

    /* Angular spacings */
    const double k1 = fine_grid.angularSpacing(i_theta_M2);
    const double k2 = fine_grid.angularSpacing(i_theta_M1);
    const double k3 = fine_grid.angularSpacing(i_theta);
    const double k4 = fine_grid.angularSpacing(i_theta_P1);

    /* Center + Angular contributions (always present) */
    double value = fine_values[fine_grid.index(i_r, i_theta)] +
                   k2 / (k1 + k2) * fine_values[fine_grid.index(i_r, i_theta_M1)] +
                   k3 / (k3 + k4) * fine_values[fine_grid.index(i_r, i_theta_P1)];

    /* Left contributions (if not at inner boundary) */
    if (i_r_coarse > 0) {
        const double h1 = fine_grid.radialSpacing(i_r - 2);
        const double h2 = fine_grid.radialSpacing(i_r - 1);
        value += h2 / (h1 + h2) * fine_values[fine_grid.index(i_r - 1, i_theta)] +
                 h2 * k2 / ((h1 + h2) * (k1 + k2)) * fine_values[fine_grid.index(i_r - 1, i_theta_M1)] +
                 h2 * k3 / ((h1 + h2) * (k3 + k4)) * fine_values[fine_grid.index(i_r - 1, i_theta_P1)];
    }

    /* Right contributions (if not at outer boundary) */
    if (i_r_coarse < coarse_grid.nr() - 1) {
        const double h3 = fine_grid.radialSpacing(i_r);
        const double h4 = fine_grid.radialSpacing(i_r + 1);
        value += h3 / (h3 + h4) * fine_values[fine_grid.index(i_r + 1, i_theta)] +
                 h3 * k2 / ((h3 + h4) * (k1 + k2)) * fine_values[fine_grid.index(i_r + 1, i_theta_M1)] +
                 h3 * k3 / ((h3 + h4) * (k3 + k4)) * fine_values[fine_grid.index(i_r + 1, i_theta_P1)];
    }

    coarse_result[coarse_grid.index(i_r_coarse, i_theta_coarse)] = value;
}

void Interpolation::applyRestriction(const PolarGrid& fine_grid, const PolarGrid& coarse_grid,
                                     Vector<double> coarse_result, ConstVector<double> fine_values) const
{
    assert(std::ssize(fine_values) == fine_grid.numberOfNodes());
    assert(std::ssize(coarse_result) == coarse_grid.numberOfNodes());

    const PolarGrid* fine_grid_ptr   = &fine_grid;
    const PolarGrid* coarse_grid_ptr = &coarse_grid;

    /* We split the loops into two regions to better respect the */
    /* access patterns of the smoother and improve cache locality. */

    // The For loop matches circular access pattern */
    Kokkos::parallel_for(
        "Interpolation: Restriction (Circular)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>( // Rank of the index space
            {0, 0}, // Starting point of the index space
            {coarse_grid.numberSmootherCircles(), coarse_grid.ntheta()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_r_coarse, const int i_theta_coarse) {
            coarseNodeRestriction(i_r_coarse, i_theta_coarse, *fine_grid_ptr, *coarse_grid_ptr, coarse_result,
                                  fine_values);
        });

    /* For loop matches radial access pattern */
    Kokkos::parallel_for(
        "Interpolation: Restriction (Radial)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>( // Rank of the index space
            {0, coarse_grid.numberSmootherCircles()}, // Starting point of the index space
            {coarse_grid.ntheta(), coarse_grid.nr()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_theta_coarse, const int i_r_coarse) {
            coarseNodeRestriction(i_r_coarse, i_theta_coarse, *fine_grid_ptr, *coarse_grid_ptr, coarse_result,
                                  fine_values);
        });

    Kokkos::fence();
}
