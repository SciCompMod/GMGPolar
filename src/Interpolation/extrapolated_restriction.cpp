#include "../../include/Interpolation/interpolation.h"
using namespace gmgpolar;

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

static KOKKOS_INLINE_FUNCTION void coarseNodeExtrapolatedRestriction(const int i_r_coarse, const int i_theta_coarse,
                                                                     const PolarGrid& fine_grid,
                                                                     const PolarGrid& coarse_grid,
                                                                     Vector<double>& coarse_result,
                                                                     ConstVector<double>& fine_values)
{
    const int i_r     = i_r_coarse * 2;
    const int i_theta = i_theta_coarse * 2;

    /* Center + Angular contributions (always present) */
    double value = fine_values[fine_grid.index(i_r, i_theta)] + 0.5 * fine_values[fine_grid.index(i_r, i_theta - 1)] +
                   0.5 * fine_values[fine_grid.index(i_r, i_theta + 1)];

    /* Left contributions (if not at inner boundary) */
    if (i_r_coarse > 0) {
        value += 0.5 * fine_values[fine_grid.index(i_r - 1, i_theta)] +
                 0.5 * fine_values[fine_grid.index(i_r - 1, i_theta + 1)]; /* Top-Left diagonal */
    }

    /* Right contributions (if not at outer boundary) */
    if (i_r_coarse < coarse_grid.nr() - 1) {
        value += 0.5 * fine_values[fine_grid.index(i_r + 1, i_theta)] +
                 0.5 * fine_values[fine_grid.index(i_r + 1, i_theta - 1)]; /* Bottom-Right diagonal */
    }

    coarse_result[coarse_grid.index(i_r_coarse, i_theta_coarse)] = value;
}

void Interpolation::applyExtrapolatedRestriction(const PolarGrid& fine_grid, const PolarGrid& coarse_grid,
                                                 Vector<double> coarse_result, ConstVector<double> fine_values) const
{
    assert(std::ssize(fine_values) == fine_grid.numberOfNodes());
    assert(std::ssize(coarse_result) == coarse_grid.numberOfNodes());

    /* We split the loops into two regions to better respect the */
    /* access patterns of the smoother and improve cache locality. */

    // The For loop matches circular access pattern */
    Kokkos::parallel_for(
        "Interpolation: Extrapolated Restriction (Circular)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>( // Rank of the index space
            {0, 0}, // Starting point of the index space
            {coarse_grid.numberSmootherCircles(), coarse_grid.ntheta()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_r_coarse, const int i_theta_coarse) {
            coarseNodeExtrapolatedRestriction(i_r_coarse, i_theta_coarse, fine_grid, coarse_grid, coarse_result,
                                              fine_values);
        });

    /* For loop matches radial access pattern */
    Kokkos::parallel_for(
        "Interpolation: Extrapolated Restriction (Radial)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>( // Rank of the index space
            {0, coarse_grid.numberSmootherCircles()}, // Starting point of the index space
            {coarse_grid.ntheta(), coarse_grid.nr()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_theta_coarse, const int i_r_coarse) {
            coarseNodeExtrapolatedRestriction(i_r_coarse, i_theta_coarse, fine_grid, coarse_grid, coarse_result,
                                              fine_values);
        });

    Kokkos::fence();
}
