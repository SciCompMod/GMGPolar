#include "../../include/Interpolation/interpolation.h"
using namespace gmgpolar;

/*
 * Prolongation Operator
 * ---------------------
 *
 * We use an anisotropic bilinear interpolation stencil.
 * For an isotropic (uniform) mesh, this stencil reduces to:
 *
 *           |1  2  1|
 *   P = 1/4 |2  4  2|
 *           |1  2  1|
 *
 * A fine node is classified by the parity of its (i_r, i_theta) indices:
 *
 * 1) (even, even)
 *    The fine node coincides with a coarse node: 
 *    -> Copy value
 *
 * 2) (odd, even)
 *    Node lies between two coarse nodes in radial direction:
 * 
 *        X ---- O ---- X
 *
 *    -> 1D radial interpolation with anisotropic weights
 *
 * 3) (even, odd)
 *    Node lies between two coarse nodes in angular direction:
 *
 *        X
 *        |
 *        O
 *        |
 *        X
 *
 *    -> 1D angular interpolation with anisotropic weights
 *
 * 4) (odd, odd)
 *    Node lies inside a coarse cell, surrounded by 4 coarse nodes:
 *
 *          X      X
 *            \   /
 *              O
 *            /   \
 *          X      X
 *
 *    -> full anisotropic bilinear interpolation
 *
 * All weights are defined via local mesh spacings:
 *  - h1, h2 in radial direction
 *  - k1, k2 in angular direction
 */

static inline void fineNodeProlongation(const int i_r, const int i_theta, const PolarGrid& coarse_grid,
                                        const PolarGrid& fine_grid, Vector<double>& fine_result,
                                        ConstVector<double>& coarse_values)
{
    const int i_r_coarse     = i_r / 2;
    const int i_theta_coarse = i_theta / 2;

    if (i_r & 1) {
        if (i_theta & 1) { /* (odd, odd) -> fine node in center of coarse cell */
            const double h1 = fine_grid.radialSpacing(i_r - 1);
            const double h2 = fine_grid.radialSpacing(i_r);
            const double k1 = fine_grid.angularSpacing(i_theta - 1);
            const double k2 = fine_grid.angularSpacing(i_theta);

            const double value =
                (h1 * k1 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)] + /* Bottom left */
                 h2 * k1 * coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse)] + /* Bottom right */
                 h1 * k2 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse + 1)] + /* Top left */
                 h2 * k2 * coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse + 1)] /* Top right */
                 ) /
                ((h1 + h2) * (k1 + k2));

            fine_result[fine_grid.index(i_r, i_theta)] = value;
        }
        else { /* (odd, even) -> between coarse nodes in radial direction */
            const double h1 = fine_grid.radialSpacing(i_r - 1);
            const double h2 = fine_grid.radialSpacing(i_r);

            const double value = (h1 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)] + /* Left */
                                  h2 * coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse)] /* Right */
                                  ) /
                                 (h1 + h2);

            fine_result[fine_grid.index(i_r, i_theta)] = value;
        }
    }
    else {
        if (i_theta & 1) { /* (even, odd) -> between coarse nodes in angular direction */
            const double k1 = fine_grid.angularSpacing(i_theta - 1);
            const double k2 = fine_grid.angularSpacing(i_theta);

            const double value = (k1 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)] + /* Bottom */
                                  k2 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse + 1)] /* Top */
                                  ) /
                                 (k1 + k2);

            fine_result[fine_grid.index(i_r, i_theta)] = value;
        }
        else { /* (even, even) -> node lies on coarse grid */
            fine_result[fine_grid.index(i_r, i_theta)] =
                coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)]; /* Center */
        }
    }
}

void Interpolation::applyProlongation(const PolarGrid& coarse_grid, const PolarGrid& fine_grid,
                                      Vector<double> fine_result, ConstVector<double> coarse_values) const
{
    assert(std::ssize(coarse_values) == coarse_grid.numberOfNodes());
    assert(std::ssize(fine_result) == fine_grid.numberOfNodes());

    /* We split the loops into two regions to better respect the */
    /* access patterns of the smoother and improve cache locality. */

    // The For loop matches circular access pattern */
    Kokkos::parallel_for(
        "Interpolation: Prolongation (Circular)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>( // Rank of the index space
            {0, 0}, // Starting point of the index space
            {fine_grid.numberSmootherCircles(), fine_grid.ntheta()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_r, const int i_theta) {
            fineNodeProlongation(i_r, i_theta, coarse_grid, fine_grid, fine_result, coarse_values);
        });

    /* For loop matches radial access pattern */
    Kokkos::parallel_for(
        "Interpolation: Prolongation (Radial)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>( // Rank of the index space
            {0, fine_grid.numberSmootherCircles()}, // Starting point of the index space
            {fine_grid.ntheta(), fine_grid.nr()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_theta, const int i_r) {
            fineNodeProlongation(i_r, i_theta, coarse_grid, fine_grid, fine_result, coarse_values);
        });

    Kokkos::fence();
}
