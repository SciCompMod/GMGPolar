#include "../../include/Interpolation/interpolation.h"

/* For the restriction we use R_ex = P_ex^T */

void Interpolation::applyExtrapolatedRestriction0(const PolarGrid& fine_grid, const PolarGrid& coarse_grid,
                                                  Vector<double> coarse_result, ConstVector<double> fine_values) const
{
    assert(fine_values.size() == static_cast<uint>(fine_grid.numberOfNodes()));
    assert(coarse_result.size() == static_cast<uint>(coarse_grid.numberOfNodes()));

#pragma omp parallel for num_threads(max_omp_threads_)
    for (int index = 0; index < coarse_grid.numberOfNodes(); index++) {
        MultiIndex coarse_node = coarse_grid.multiIndex(index);
        MultiIndex fine_node(2 * coarse_node[0], 2 * coarse_node[1]);

        std::array<std::pair<int, int>, space_dimension> neighbors;

        // Center
        double value = fine_values[fine_grid.index(fine_node)];

        fine_grid.adjacentNeighborsOf(fine_node, neighbors);

        // Left
        if (neighbors[0].first != -1) {
            value += 0.5 * fine_values[neighbors[0].first];
        }

        // Right
        if (neighbors[0].second != -1) {
            value += 0.5 * fine_values[neighbors[0].second];
        }

        // Bottom
        if (neighbors[1].first != -1) {
            value += 0.5 * fine_values[neighbors[1].first];
        }

        // Top
        if (neighbors[1].second != -1) {
            value += 0.5 * fine_values[neighbors[1].second];
        }

        fine_grid.diagonalNeighborsOf(fine_node, neighbors);

        // Bottom Right
        if (neighbors[0].second != -1) {
            value += 0.5 * fine_values[neighbors[0].second];
        }

        // Top Left
        if (neighbors[1].first != -1) {
            value += 0.5 * fine_values[neighbors[1].first];
        }

        coarse_result[index] = value;
    }
}

// -------------------------------------- //
// Optimized version of applyRestriction0 //
// -------------------------------------- //

void Interpolation::applyExtrapolatedRestriction(const PolarGrid& fine_grid, const PolarGrid& coarse_grid,
                                                 Vector<double> coarse_result, ConstVector<double> fine_values) const
{
    assert(fine_values.size() == static_cast<uint>(fine_grid.numberOfNodes()));
    assert(coarse_result.size() == static_cast<uint>(coarse_grid.numberOfNodes()));

    const int coarseNumberSmootherCircles = coarse_grid.numberSmootherCircles();

#pragma omp parallel num_threads(max_omp_threads_)
    {
/* For loop matches circular access pattern */
#pragma omp for nowait
        for (int i_r_coarse = 0; i_r_coarse < coarseNumberSmootherCircles; i_r_coarse++) {
            int i_r = i_r_coarse * 2;
            for (int i_theta_coarse = 0; i_theta_coarse < coarse_grid.ntheta(); i_theta_coarse++) {
                int i_theta = i_theta_coarse * 2;

                if (0 < i_r_coarse && i_r_coarse < coarseNumberSmootherCircles - 1) {
                    int i_theta_M1 = fine_grid.wrapThetaIndex(i_theta - 1);
                    int i_theta_P1 = fine_grid.wrapThetaIndex(i_theta + 1);

                    coarse_result[coarse_grid.index(i_r_coarse, i_theta_coarse)] =
                        // Center
                        fine_values[fine_grid.index(i_r, i_theta)] +
                        // Left, Right, Bottom, Top
                        0.5 * fine_values[fine_grid.index(i_r - 1, i_theta)] +
                        0.5 * fine_values[fine_grid.index(i_r + 1, i_theta)] +
                        0.5 * fine_values[fine_grid.index(i_r, i_theta_M1)] +
                        0.5 * fine_values[fine_grid.index(i_r, i_theta_P1)] +
                        // Bottom Right, Top Left
                        0.5 * fine_values[fine_grid.index(i_r + 1, i_theta_M1)] +
                        0.5 * fine_values[fine_grid.index(i_r - 1, i_theta_P1)];
                }
                else {
                    /* First and Last Circle have to be checked for domain boundary */
                    int i_theta_M1 = fine_grid.wrapThetaIndex(i_theta - 1);
                    int i_theta_P1 = fine_grid.wrapThetaIndex(i_theta + 1);
                    // Center, Bottom, Top
                    double value = fine_values[fine_grid.index(i_r, i_theta)] +
                                   0.5 * fine_values[fine_grid.index(i_r, i_theta_M1)] +
                                   0.5 * fine_values[fine_grid.index(i_r, i_theta_P1)];

                    if (i_r_coarse > 0) {
                        // Left, Top Left
                        value += 0.5 * fine_values[fine_grid.index(i_r - 1, i_theta)] +
                                 0.5 * fine_values[fine_grid.index(i_r - 1, i_theta_P1)];
                    }
                    if (i_r_coarse < coarse_grid.nr() - 1) {
                        // Right, Bottom Right
                        value += 0.5 * fine_values[fine_grid.index(i_r + 1, i_theta)] +
                                 0.5 * fine_values[fine_grid.index(i_r + 1, i_theta_M1)];
                    }
                    coarse_result[coarse_grid.index(i_r_coarse, i_theta_coarse)] = value;
                }
            }
        }

/* For loop matches radial access pattern */
#pragma omp for nowait
        for (int i_theta_coarse = 0; i_theta_coarse < coarse_grid.ntheta(); i_theta_coarse++) {
            int i_theta = i_theta_coarse * 2;
            for (int i_r_coarse = coarseNumberSmootherCircles; i_r_coarse < coarse_grid.nr(); i_r_coarse++) {
                int i_r = i_r_coarse * 2;

                if (coarse_grid.numberSmootherCircles() < i_r_coarse && i_r_coarse < coarse_grid.nr() - 1) {
                    int i_theta_M1 = fine_grid.wrapThetaIndex(i_theta - 1);
                    int i_theta_P1 = fine_grid.wrapThetaIndex(i_theta + 1);

                    coarse_result[coarse_grid.index(i_r_coarse, i_theta_coarse)] =
                        // Center
                        fine_values[fine_grid.index(i_r, i_theta)] +
                        // Left, Right, Bottom, Top
                        0.5 * fine_values[fine_grid.index(i_r - 1, i_theta)] +
                        0.5 * fine_values[fine_grid.index(i_r + 1, i_theta)] +
                        0.5 * fine_values[fine_grid.index(i_r, i_theta_M1)] +
                        0.5 * fine_values[fine_grid.index(i_r, i_theta_P1)] +
                        // Bottom Right, Top Left
                        0.5 * fine_values[fine_grid.index(i_r + 1, i_theta_M1)] +
                        0.5 * fine_values[fine_grid.index(i_r - 1, i_theta_P1)];
                }
                else {
                    /* First and Last radial nodes have to be checked for domain boundary */
                    int i_theta_M1 = fine_grid.wrapThetaIndex(i_theta - 1);
                    int i_theta_P1 = fine_grid.wrapThetaIndex(i_theta + 1);
                    // Center, Bottom, Top
                    double value = fine_values[fine_grid.index(i_r, i_theta)] +
                                   0.5 * fine_values[fine_grid.index(i_r, i_theta_M1)] +
                                   0.5 * fine_values[fine_grid.index(i_r, i_theta_P1)];
                    if (i_r_coarse > 0) {
                        // Left, Top Left
                        value += 0.5 * fine_values[fine_grid.index(i_r - 1, i_theta)] +
                                 0.5 * fine_values[fine_grid.index(i_r - 1, i_theta_P1)];
                    }
                    if (i_r_coarse < coarse_grid.nr() - 1) {
                        // Right, Bottom Right
                        value += 0.5 * fine_values[fine_grid.index(i_r + 1, i_theta)] +
                                 0.5 * fine_values[fine_grid.index(i_r + 1, i_theta_M1)];
                    }
                    coarse_result[coarse_grid.index(i_r_coarse, i_theta_coarse)] = value;
                }
            }
        }
    }
}
