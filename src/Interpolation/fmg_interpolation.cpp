#include "../../include/Interpolation/interpolation.h"
using namespace gmgpolar;

/*
 * Bicubic FMG Interpolator using the Lagrange Interpolating Polynomial
 * --------------------------------------------------------------------
 * In contrast to the interpolation of the multigrid cycle, which is applied to corrections,
 * the FMG interpolation transfers approximations of the solution to the fine grid.
 * The FMG scheme requires interpolation of higher order than the discretization order 
 * of the PDE. For second-order discretizations, bicubic interpolation is commonly used.
 * For a uniform grid, this reduces to the 1D stencil: 1/16 * [-1, 9, 9, -1]
 *
 * The method uses Lagrange interpolation with 4 coarse grid nodes to calculate a single interpolated value.
 * Nodes with an X are the 4 coarse interpolation points used for the interpolation.
 * The node marked with an O is the interpolated value which we aim to obtain.
 *
 * |                     |                     |                     |
 * |                     |                     |                     |
 * v                     v                     v                     v
 * X ---------h0-------- X ---h1--- O ---h2--- X ---------h3-------- X
 *                                  ^
 *                                  |
 *                                  |
 * Lagrange Interpolation:
 * -----------------------
 * The interpolated value y is determined by forming a weighted sum of the function values:
 *
 *     y = w0 * f0 + w1 * f1 + w2 * f2 + w3 * f3
 *
 * Each weight w_i is calculated using the Lagrange interpolation formula. For example:
 *
 *     w0 = ((z - x1) / (x0 - x1)) * ((z - x2) / (x0 - x2)) * ((z - x3) / (x0 - x3))
 *
 * Similar expressions are used to compute weights w1, w2, and w3, where each weight excludes its corresponding node
 * from the product.
 * 
 * First, the interpolation is performed along the angular direction. Thanks to the periodic boundary conditions,
 * no additional linear interpolation is required near the edges. In the subsequent step, the angular interpolation
 * results are further refined through radial interpolation to produce the final interpolated value.
 * 
 * Step 1: Interpolate the 4 fine nodes marked with an O.
 * Step 2: Interpolate the node Z using the results from Step 1.
 * 
 * X         X          X         X
 * 
 * 
 * 
 * X         X          X         X
 * 
 * O         O     Z    O         O
 * 
 * X         X          X         X
 * 
 * 
 * 
 * X         X          X         X
 * 
 */

static KOKKOS_INLINE_FUNCTION void fineNodeFMGInterpolation(const int i_r, const int i_theta,
                                                            const PolarGrid& coarse_grid, const PolarGrid& fine_grid,
                                                            Vector<double>& fine_result,
                                                            ConstVector<double>& coarse_values)
{
    const int i_r_coarse     = i_r / 2;
    const int i_theta_coarse = i_theta / 2;

    /* Case 1: On the boundary */
    if (i_r == 0 || i_r == fine_grid.nr() - 1) {
        if (i_theta & 1) {
            const double k0 = coarse_grid.angularSpacing(i_theta_coarse - 1);
            const double k1 = fine_grid.angularSpacing(i_theta - 1);
            const double k2 = fine_grid.angularSpacing(i_theta);
            const double k3 = coarse_grid.angularSpacing(i_theta_coarse + 1);

            const double w_theta0 = -k1 / k0 * k2 / (k0 + k1 + k2) * (k2 + k3) / (k0 + k1 + k2 + k3);
            const double w_theta1 = (k0 + k1) / k0 * k2 / (k1 + k2) * (k2 + k3) / (k1 + k2 + k3);
            const double w_theta2 = (k0 + k1) / (k0 + k1 + k2) * k1 / (k1 + k2) * (k2 + k3) / k3;
            const double w_theta3 = -(k0 + k1) / (k0 + k1 + k2 + k3) * k1 / (k1 + k2 + k3) * k2 / k3;

            fine_result[fine_grid.index(i_r, i_theta)] =
                (w_theta0 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse - 1)] + /* (0, -3) */
                 w_theta1 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)] + /* (0, -1) */
                 w_theta2 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse + 1)] + /* (0, +1) */
                 w_theta3 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse + 2)] /* (0, +3) */
                );
        }
        else {
            fine_result[fine_grid.index(i_r, i_theta)] =
                coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)]; /* center */
        }
    } /* Case 2: Next to the boundary */
    else if (i_r == 1 || i_r == fine_grid.nr() - 2) {
        if (i_theta & 1) {
            const double k0 = coarse_grid.angularSpacing(i_theta_coarse - 1);
            const double k1 = fine_grid.angularSpacing(i_theta - 1);
            const double k2 = fine_grid.angularSpacing(i_theta);
            const double k3 = coarse_grid.angularSpacing(i_theta_coarse + 1);

            const double w_theta0 = -k1 / k0 * k2 / (k0 + k1 + k2) * (k2 + k3) / (k0 + k1 + k2 + k3);
            const double w_theta1 = (k0 + k1) / k0 * k2 / (k1 + k2) * (k2 + k3) / (k1 + k2 + k3);
            const double w_theta2 = (k0 + k1) / (k0 + k1 + k2) * k1 / (k1 + k2) * (k2 + k3) / k3;
            const double w_theta3 = -(k0 + k1) / (k0 + k1 + k2 + k3) * k1 / (k1 + k2 + k3) * k2 / k3;

            const double left_value =
                (w_theta0 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse - 1)] + /* (-1, -3) */
                 w_theta1 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)] + /* (-1, -1) */
                 w_theta2 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse + 1)] + /* (-1, +1) */
                 w_theta3 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse + 2)] /* (-1, +3) */
                );
            const double right_value =
                (w_theta0 * coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse - 1)] + /* (+1, -3) */
                 w_theta1 * coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse)] + /* (+1, -1) */
                 w_theta2 * coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse + 1)] + /* (+1, +1) */
                 w_theta3 * coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse + 2)] /* (+1, +3) */
                );

            const double h1                            = fine_grid.radialSpacing(i_r - 1);
            const double h2                            = fine_grid.radialSpacing(i_r);
            fine_result[fine_grid.index(i_r, i_theta)] = (h1 * left_value + h2 * right_value) / (h1 + h2);
        }
        else {
            const double h1 = fine_grid.radialSpacing(i_r - 1);
            const double h2 = fine_grid.radialSpacing(i_r);
            fine_result[fine_grid.index(i_r, i_theta)] =
                (h1 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)] + /* left */
                 h2 * coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse)] /* right */
                 ) /
                (h1 + h2);
        }
    }
    else { /* Case 3: In the interior */
        if (i_r & 1) {
            if (i_theta & 1) {
                const double k0 = coarse_grid.angularSpacing(i_theta_coarse - 1);
                const double k1 = fine_grid.angularSpacing(i_theta - 1);
                const double k2 = fine_grid.angularSpacing(i_theta);
                const double k3 = coarse_grid.angularSpacing(i_theta_coarse + 1);

                const double w_theta0 = -k1 / k0 * k2 / (k0 + k1 + k2) * (k2 + k3) / (k0 + k1 + k2 + k3);
                const double w_theta1 = (k0 + k1) / k0 * k2 / (k1 + k2) * (k2 + k3) / (k1 + k2 + k3);
                const double w_theta2 = (k0 + k1) / (k0 + k1 + k2) * k1 / (k1 + k2) * (k2 + k3) / k3;
                const double w_theta3 = -(k0 + k1) / (k0 + k1 + k2 + k3) * k1 / (k1 + k2 + k3) * k2 / k3;

                const double outer_left_value =
                    (w_theta0 * coarse_values[coarse_grid.index(i_r_coarse - 1, i_theta_coarse - 1)] + /* (-3, -3) */
                     w_theta1 * coarse_values[coarse_grid.index(i_r_coarse - 1, i_theta_coarse)] + /* (-3, -1) */
                     w_theta2 * coarse_values[coarse_grid.index(i_r_coarse - 1, i_theta_coarse + 1)] + /* (-3, +1) */
                     w_theta3 * coarse_values[coarse_grid.index(i_r_coarse - 1, i_theta_coarse + 2)] /* (-3, +3) */
                    );
                const double inner_left_value =
                    (w_theta0 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse - 1)] + /* (-1, -3) */
                     w_theta1 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)] + /* (-1, -1) */
                     w_theta2 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse + 1)] + /* (-1, +1) */
                     w_theta3 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse + 2)] /* (-1, +3) */
                    );
                const double inner_right_value =
                    (w_theta0 * coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse - 1)] + /* (+1, -3) */
                     w_theta1 * coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse)] + /* (+1, -1) */
                     w_theta2 * coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse + 1)] + /* (+1, +1) */
                     w_theta3 * coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse + 2)] /* (+1, +3) */
                    );
                const double outer_right_value =
                    (w_theta0 * coarse_values[coarse_grid.index(i_r_coarse + 2, i_theta_coarse - 1)] + /* (+3, -3) */
                     w_theta1 * coarse_values[coarse_grid.index(i_r_coarse + 2, i_theta_coarse)] + /* (+3, -1) */
                     w_theta2 * coarse_values[coarse_grid.index(i_r_coarse + 2, i_theta_coarse + 1)] + /* (+3, +1) */
                     w_theta3 * coarse_values[coarse_grid.index(i_r_coarse + 2, i_theta_coarse + 2)] /* (+3, +3) */
                    );

                const double h0 = coarse_grid.radialSpacing(i_r_coarse - 1);
                const double h1 = fine_grid.radialSpacing(i_r - 1);
                const double h2 = fine_grid.radialSpacing(i_r);
                const double h3 = coarse_grid.radialSpacing(i_r_coarse + 1);

                const double w_r0 = -h1 / h0 * h2 / (h0 + h1 + h2) * (h2 + h3) / (h0 + h1 + h2 + h3);
                const double w_r1 = (h0 + h1) / h0 * h2 / (h1 + h2) * (h2 + h3) / (h1 + h2 + h3);
                const double w_r2 = (h0 + h1) / (h0 + h1 + h2) * h1 / (h1 + h2) * (h2 + h3) / h3;
                const double w_r3 = -(h0 + h1) / (h0 + h1 + h2 + h3) * h1 / (h1 + h2 + h3) * h2 / h3;

                fine_result[fine_grid.index(i_r, i_theta)] = (w_r0 * outer_left_value + w_r1 * inner_left_value +
                                                              w_r2 * inner_right_value + w_r3 * outer_right_value);
            }
            else {
                const double h0 = coarse_grid.radialSpacing(i_r_coarse - 1);
                const double h1 = fine_grid.radialSpacing(i_r - 1);
                const double h2 = fine_grid.radialSpacing(i_r);
                const double h3 = coarse_grid.radialSpacing(i_r_coarse + 1);

                const double w_r0 = -h1 / h0 * h2 / (h0 + h1 + h2) * (h2 + h3) / (h0 + h1 + h2 + h3);
                const double w_r1 = (h0 + h1) / h0 * h2 / (h1 + h2) * (h2 + h3) / (h1 + h2 + h3);
                const double w_r2 = (h0 + h1) / (h0 + h1 + h2) * h1 / (h1 + h2) * (h2 + h3) / h3;
                const double w_r3 = -(h0 + h1) / (h0 + h1 + h2 + h3) * h1 / (h1 + h2 + h3) * h2 / h3;

                fine_result[fine_grid.index(i_r, i_theta)] =
                    (w_r0 * coarse_values[coarse_grid.index(i_r_coarse - 1, i_theta_coarse)] + /* (-3, 0) */
                     w_r1 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)] + /* (-1, 0) */
                     w_r2 * coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse)] + /* (+1, 0) */
                     w_r3 * coarse_values[coarse_grid.index(i_r_coarse + 2, i_theta_coarse)] /* (+3, 0) */
                    );
            }
        }
        else {
            if (i_theta & 1) {
                const double k0 = coarse_grid.angularSpacing(i_theta_coarse - 1);
                const double k1 = fine_grid.angularSpacing(i_theta - 1);
                const double k2 = fine_grid.angularSpacing(i_theta);
                const double k3 = coarse_grid.angularSpacing(i_theta_coarse + 1);

                const double w_theta0 = -k1 / k0 * k2 / (k0 + k1 + k2) * (k2 + k3) / (k0 + k1 + k2 + k3);
                const double w_theta1 = (k0 + k1) / k0 * k2 / (k1 + k2) * (k2 + k3) / (k1 + k2 + k3);
                const double w_theta2 = (k0 + k1) / (k0 + k1 + k2) * k1 / (k1 + k2) * (k2 + k3) / k3;
                const double w_theta3 = -(k0 + k1) / (k0 + k1 + k2 + k3) * k1 / (k1 + k2 + k3) * k2 / k3;

                fine_result[fine_grid.index(i_r, i_theta)] =
                    (w_theta0 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse - 1)] + /* (0, -3) */
                     w_theta1 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)] + /* (0, -1) */
                     w_theta2 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse + 1)] + /* (0, +1) */
                     w_theta3 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse + 2)] /* (0, +3) */
                    );
            }
            else {
                fine_result[fine_grid.index(i_r, i_theta)] =
                    coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)]; /* center */
            }
        }
    }
}

void Interpolation::applyFMGInterpolation(const PolarGrid& coarse_grid, const PolarGrid& fine_grid,
                                          Vector<double> fine_result, ConstVector<double> coarse_values) const
{
    assert(std::ssize(coarse_values) == coarse_grid.numberOfNodes());
    assert(std::ssize(fine_result) == fine_grid.numberOfNodes());

    /* We split the loops into two regions to better respect the */
    /* access patterns of the smoother and improve cache locality. */

    // The For loop matches circular access pattern */
    Kokkos::parallel_for(
        "Interpolation: FMG-Interpolation (Circular)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>( // Rank of the index space
            {0, 0}, // Starting point of the index space
            {fine_grid.numberSmootherCircles(), fine_grid.ntheta()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_r, const int i_theta) {
            fineNodeFMGInterpolation(i_r, i_theta, coarse_grid, fine_grid, fine_result, coarse_values);
        });

    /* For loop matches radial access pattern */
    Kokkos::parallel_for(
        "Interpolation: FMG-Interpolation (Radial)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>( // Rank of the index space
            {0, fine_grid.numberSmootherCircles()}, // Starting point of the index space
            {fine_grid.ntheta(), fine_grid.nr()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_theta, const int i_r) {
            fineNodeFMGInterpolation(i_r, i_theta, coarse_grid, fine_grid, fine_result, coarse_values);
        });

    Kokkos::fence();
}
