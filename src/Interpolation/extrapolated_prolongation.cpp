#include "../../include/Interpolation/interpolation.h"

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

inline void fineNodeExtrapolatedProlongation(int i_r, int i_theta, int i_r_coarse, int i_theta_coarse,
                                             const PolarGrid& coarse_grid, const PolarGrid& fine_grid,
                                             Vector<double>& fine_result, ConstVector<double>& coarse_values)
{
    if (i_r & 1) {
        if (i_theta & 1) { /* (odd, odd) -> node in center of coarse cell */
            double value = 0.5 * (coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse)] + /* Bottom right */
                                  coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse + 1)] /* Top left */
                                 );
            fine_result[fine_grid.index(i_r, i_theta)] = value;
        }
        else { /* (odd, even) -> between coarse nodes in radial direction */
            double value = 0.5 * (coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)] + /* Left */
                                  coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse)] /* Right */
                                 );
            fine_result[fine_grid.index(i_r, i_theta)] = value;
        }
    }
    else {
        if (i_theta & 1) { /* (even, odd) -> between coarse nodes in angular direction */
            double value = 0.5 * (coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)] + /* Bottom */
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

#pragma omp parallel num_threads(max_omp_threads_)
    {
        /* Circular Indexing Section */
#pragma omp for nowait
        for (int i_r = 0; i_r < fine_grid.numberSmootherCircles(); i_r++) {
            int i_r_coarse = i_r / 2;
            for (int i_theta = 0; i_theta < fine_grid.ntheta(); i_theta++) {
                int i_theta_coarse = i_theta / 2;
                fineNodeExtrapolatedProlongation(i_r, i_theta, i_r_coarse, i_theta_coarse, coarse_grid, fine_grid,
                                                 fine_result, coarse_values);
            }
        }

        /* Radial Indexing Section */
#pragma omp for nowait
        for (int i_theta = 0; i_theta < fine_grid.ntheta(); i_theta++) {
            int i_theta_coarse = i_theta / 2;
            for (int i_r = fine_grid.numberSmootherCircles(); i_r < fine_grid.nr(); i_r++) {
                int i_r_coarse = i_r / 2;
                fineNodeExtrapolatedProlongation(i_r, i_theta, i_r_coarse, i_theta_coarse, coarse_grid, fine_grid,
                                                 fine_result, coarse_values);
            }
        }
    }
}
