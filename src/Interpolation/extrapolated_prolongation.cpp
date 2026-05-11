#include "../../include/Interpolation/interpolation.h"
using namespace gmgpolar;

/*
 * Extrapolated Prolongation Operator
 * ----------------------------------
 *
 * Extrapolated prolongation is used between the finest most grids in the multigrid hierarchy.
 * It assumes the fine grid comes from a uniform refinement of the coarse grid, so all spacings are equal.
 * Thus fine values can be computed simply by averaging the neighboring coarse nodes.
 *
 * A fine node is classified by the parity of its (i_r, i_theta) indices:
 *
 * 1) (even, even)
 *    Fine node coincides with a coarse node
 *    -> copy value
 *
 * 2) (odd, even)
 *    Node lies between two coarse nodes in radial direction
 *
 *        X ---- O ---- X
 *
 *    -> arithmetic mean of left + right coarse node
 *
 * 3) (even, odd)
 *    Node lies between two coarse nodes in angular direction
 *
 *        X
 *        |
 *        O
 *        |
 *        X
 *
 *    -> arithmetic mean of bottom + top coarse node
 *
 * 4) (odd, odd)
 *    Node lies inside a coarse cell
 *    We extrapolate/average across the diagonal:
 *
 *          X      
 *            \   
 *              O
 *                \
 *                 X
 * 
 *   -> arithmetic mean of bottom right + top left coarse node
 *
 */

static KOKKOS_INLINE_FUNCTION void fineNodeExtrapolatedProlongation(const int i_r, const int i_theta,
                                                                    const PolarGrid& coarse_grid,
                                                                    const PolarGrid& fine_grid,
                                                                    Vector<double>& fine_result,
                                                                    ConstVector<double>& coarse_values)
{
    const int i_r_coarse     = i_r / 2;
    const int i_theta_coarse = i_theta / 2;

    if (i_r & 1) {
        if (i_theta & 1) { /* (odd, odd) -> node in center of coarse cell */
            const double value =
                0.5 * (coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse)] + /* Bottom right */
                       coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse + 1)] /* Top left */
                      );
            fine_result[fine_grid.index(i_r, i_theta)] = value;
        }
        else { /* (odd, even) -> between coarse nodes in radial direction */
            const double value = 0.5 * (coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)] + /* Left */
                                        coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse)] /* Right */
                                       );
            fine_result[fine_grid.index(i_r, i_theta)] = value;
        }
    }
    else {
        if (i_theta & 1) { /* (even, odd) -> between coarse nodes in angular direction */
            const double value = 0.5 * (coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)] + /* Bottom */
                                        coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse + 1)] /* Top */
                                       );
            fine_result[fine_grid.index(i_r, i_theta)] = value;
        }
        else { /* (even, even) -> node lies exactly on coarse grid */
            fine_result[fine_grid.index(i_r, i_theta)] =
                coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)]; /* Center */
        }
    }
}

void Interpolation::applyExtrapolatedProlongation(const PolarGrid& coarse_grid, const PolarGrid& fine_grid,
                                                  Vector<double> fine_result, ConstVector<double> coarse_values) const
{
    assert(std::ssize(coarse_values) == coarse_grid.numberOfNodes());
    assert(std::ssize(fine_result) == fine_grid.numberOfNodes());

    /* We split the loops into two regions to better respect the */
    /* access patterns of the smoother and improve cache locality. */

    // The For loop matches circular access pattern */
    Kokkos::parallel_for(
        "Interpolation: Extrapolated Prolongation (Circular)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>( // Rank of the index space
            {0, 0}, // Starting point of the index space
            {fine_grid.numberSmootherCircles(), fine_grid.ntheta()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_r, const int i_theta) {
            fineNodeExtrapolatedProlongation(i_r, i_theta, coarse_grid, fine_grid, fine_result, coarse_values);
        });

    /* For loop matches radial access pattern */
    Kokkos::parallel_for(
        "Interpolation: Extrapolated Prolongation (Radial)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>( // Rank of the index space
            {0, fine_grid.numberSmootherCircles()}, // Starting point of the index space
            {fine_grid.ntheta(), fine_grid.nr()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_theta, const int i_r) {
            fineNodeExtrapolatedProlongation(i_r, i_theta, coarse_grid, fine_grid, fine_result, coarse_values);
        });

    Kokkos::fence();
}
