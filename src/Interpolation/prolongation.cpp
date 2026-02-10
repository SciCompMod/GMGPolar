#include "../../include/Interpolation/interpolation.h"

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

inline void fineNodeProlongation(int i_r, int i_theta, int i_r_coarse, int i_theta_coarse, const PolarGrid& coarse_grid,
                                 const PolarGrid& fine_grid, Vector<double>& fine_result,
                                 ConstVector<double>& coarse_values)
{
    if (i_r & 1) {
        if (i_theta & 1) { /* (odd, odd) -> fine node in center of coarse cell */
            double h1 = fine_grid.radialSpacing(i_r - 1);
            double h2 = fine_grid.radialSpacing(i_r);
            double k1 = fine_grid.angularSpacing(i_theta - 1);
            double k2 = fine_grid.angularSpacing(i_theta);

            double value =
                (h1 * k1 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)] + /* Bottom left */
                 h2 * k1 * coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse)] + /* Bottom right */
                 h1 * k2 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse + 1)] + /* Top left */
                 h2 * k2 * coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse + 1)] /* Top right */
                 ) /
                ((h1 + h2) * (k1 + k2));

            fine_result[fine_grid.index(i_r, i_theta)] = value;
        }
        else { /* (odd, even) -> between coarse nodes in radial direction */
            double h1 = fine_grid.radialSpacing(i_r - 1);
            double h2 = fine_grid.radialSpacing(i_r);

            double value = (h1 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)] + /* Left */
                            h2 * coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse)] /* Right */
                            ) /
                           (h1 + h2);

            fine_result[fine_grid.index(i_r, i_theta)] = value;
        }
    }
    else {
        if (i_theta & 1) { /* (even, odd) -> between coarse nodes in angular direction */
            double k1 = fine_grid.angularSpacing(i_theta - 1);
            double k2 = fine_grid.angularSpacing(i_theta);

            double value = (k1 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)] + /* Bottom */
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

#pragma omp parallel num_threads(max_omp_threads_)
    {
        /* Circular Indexing Section */
#pragma omp for nowait
        for (int i_r = 0; i_r < fine_grid.numberSmootherCircles(); i_r++) {
            int i_r_coarse = i_r / 2;
            for (int i_theta = 0; i_theta < fine_grid.ntheta(); i_theta++) {
                int i_theta_coarse = i_theta / 2;
                fineNodeProlongation(i_r, i_theta, i_r_coarse, i_theta_coarse, coarse_grid, fine_grid, fine_result,
                                     coarse_values);
            }
        }

        /* Radial Indexing Section */
#pragma omp for nowait
        for (int i_theta = 0; i_theta < fine_grid.ntheta(); i_theta++) {
            int i_theta_coarse = i_theta / 2;
            for (int i_r = fine_grid.numberSmootherCircles(); i_r < fine_grid.nr(); i_r++) {
                int i_r_coarse = i_r / 2;
                fineNodeProlongation(i_r, i_theta, i_r_coarse, i_theta_coarse, coarse_grid, fine_grid, fine_result,
                                     coarse_values);
            }
        }
    }
}
