#include "../../include/Interpolation/interpolation.h"

// For the restriction we use R = P^T.
// The restriction for an siotropic mesh reduces to
//           |1  2  1|
// R = 1/4 * |2  4  2| = P^T
//           |1  2  1|

void Interpolation::applyRestriction0(const PolarGrid& fine_grid, const PolarGrid& coarse_grid,
                                      Vector<double> coarse_result, ConstVector<double> fine_values) const
{
    assert(fine_values.size() == static_cast<uint>(fine_grid.numberOfNodes()));
    assert(coarse_result.size() == static_cast<uint>(coarse_grid.numberOfNodes()));

#pragma omp parallel for num_threads(max_omp_threads_)
    for (int index = 0; index < coarse_grid.numberOfNodes(); index++) {
        MultiIndex coarse_node = coarse_grid.multiIndex(index);
        MultiIndex fine_node(2 * coarse_node[0], 2 * coarse_node[1]);

        std::array<std::pair<double, double>, space_dimension> neighbor_distance;
        std::array<std::pair<int, int>, space_dimension> neighbors;

        // Center
        double value = fine_values[fine_grid.index(fine_node)];

        fine_grid.adjacentNeighborsOf(fine_node, neighbors);

        // Left
        if (neighbors[0].first != -1) {
            fine_grid.adjacentNeighborDistances(fine_grid.multiIndex(neighbors[0].first), neighbor_distance);
            value += neighbor_distance[0].second * fine_values[neighbors[0].first] /
                     (neighbor_distance[0].first + neighbor_distance[0].second);
        }

        // Right
        if (neighbors[0].second != -1) {
            fine_grid.adjacentNeighborDistances(fine_grid.multiIndex(neighbors[0].second), neighbor_distance);
            value += neighbor_distance[0].first * fine_values[neighbors[0].second] /
                     (neighbor_distance[0].first + neighbor_distance[0].second);
        }

        // Bottom
        if (neighbors[1].first != -1) {
            fine_grid.adjacentNeighborDistances(fine_grid.multiIndex(neighbors[1].first), neighbor_distance);
            value += neighbor_distance[1].second * fine_values[neighbors[1].first] /
                     (neighbor_distance[1].first + neighbor_distance[1].second);
        }

        // Top
        if (neighbors[1].second != -1) {
            fine_grid.adjacentNeighborDistances(fine_grid.multiIndex(neighbors[1].second), neighbor_distance);
            value += neighbor_distance[1].first * fine_values[neighbors[1].second] /
                     (neighbor_distance[1].first + neighbor_distance[1].second);
        }

        fine_grid.diagonalNeighborsOf(fine_node, neighbors);

        // Bottom Left
        if (neighbors[0].first != -1) {
            fine_grid.adjacentNeighborDistances(fine_grid.multiIndex(neighbors[0].first), neighbor_distance);
            value += neighbor_distance[0].second * neighbor_distance[1].second * fine_values[neighbors[0].first] /
                     ((neighbor_distance[0].first + neighbor_distance[0].second) *
                      (neighbor_distance[1].first + neighbor_distance[1].second));
        }

        // Bottom Right
        if (neighbors[0].second != -1) {
            fine_grid.adjacentNeighborDistances(fine_grid.multiIndex(neighbors[0].second), neighbor_distance);
            value += neighbor_distance[0].first * neighbor_distance[1].second * fine_values[neighbors[0].second] /
                     ((neighbor_distance[0].first + neighbor_distance[0].second) *
                      (neighbor_distance[1].first + neighbor_distance[1].second));
        }

        // Top Left
        if (neighbors[1].first != -1) {
            fine_grid.adjacentNeighborDistances(fine_grid.multiIndex(neighbors[1].first), neighbor_distance);
            value += neighbor_distance[0].second * neighbor_distance[1].first * fine_values[neighbors[1].first] /
                     ((neighbor_distance[0].first + neighbor_distance[0].second) *
                      (neighbor_distance[1].first + neighbor_distance[1].second));
        }

        // Top Right
        if (neighbors[1].second != -1) {
            fine_grid.adjacentNeighborDistances(fine_grid.multiIndex(neighbors[1].second), neighbor_distance);
            value += neighbor_distance[0].first * neighbor_distance[1].first * fine_values[neighbors[1].second] /
                     ((neighbor_distance[0].first + neighbor_distance[0].second) *
                      (neighbor_distance[1].first + neighbor_distance[1].second));
        }

        coarse_result[index] = value;
    }
}

// -------------------------------------- //
// Optimized version of applyRestriction0 //
// -------------------------------------- //

void Interpolation::applyRestriction(const PolarGrid& fine_grid, const PolarGrid& coarse_grid,
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
                    double h1 = fine_grid.radialSpacing(i_r - 2);
                    double h2 = fine_grid.radialSpacing(i_r - 1);
                    double h3 = fine_grid.radialSpacing(i_r);
                    double h4 = fine_grid.radialSpacing(i_r + 1);

                    int i_theta_M2 = fine_grid.wrapThetaIndex(i_theta - 2);
                    int i_theta_M1 = fine_grid.wrapThetaIndex(i_theta - 1);
                    int i_theta_P1 = fine_grid.wrapThetaIndex(i_theta + 1);

                    double k1 = fine_grid.angularSpacing(i_theta_M2);
                    double k2 = fine_grid.angularSpacing(i_theta_M1);
                    double k3 = fine_grid.angularSpacing(i_theta);
                    double k4 = fine_grid.angularSpacing(i_theta_P1);

                    coarse_result[coarse_grid.index(i_r_coarse, i_theta_coarse)] =
                        // Center
                        fine_values[fine_grid.index(i_r, i_theta)] +
                        // Left, Right, Bottom, Top
                        h2 * fine_values[fine_grid.index(i_r - 1, i_theta)] / (h1 + h2) +
                        h3 * fine_values[fine_grid.index(i_r + 1, i_theta)] / (h3 + h4) +
                        k2 * fine_values[fine_grid.index(i_r, i_theta_M1)] / (k1 + k2) +
                        k3 * fine_values[fine_grid.index(i_r, i_theta_P1)] / (k3 + k4) +
                        // Bottom Left, Bottom Right, Top Left, Top Right
                        h2 * k2 * fine_values[fine_grid.index(i_r - 1, i_theta_M1)] / ((h1 + h2) * (k1 + k2)) +
                        h3 * k2 * fine_values[fine_grid.index(i_r + 1, i_theta_M1)] / ((h3 + h4) * (k1 + k2)) +
                        h2 * k3 * fine_values[fine_grid.index(i_r - 1, i_theta_P1)] / ((h1 + h2) * (k3 + k4)) +
                        h3 * k3 * fine_values[fine_grid.index(i_r + 1, i_theta_P1)] / ((h3 + h4) * (k3 + k4));
                }
                else {
                    /* First and Last Circle have to be checked for domain boundary */
                    // Middle Part
                    int i_theta_M2 = fine_grid.wrapThetaIndex(i_theta - 2);
                    int i_theta_M1 = fine_grid.wrapThetaIndex(i_theta - 1);
                    int i_theta_P1 = fine_grid.wrapThetaIndex(i_theta + 1);
                    double k1      = fine_grid.angularSpacing(i_theta_M2);
                    double k2      = fine_grid.angularSpacing(i_theta_M1);
                    double k3      = fine_grid.angularSpacing(i_theta);
                    double k4      = fine_grid.angularSpacing(i_theta_P1);
                    // Center, Bottom, Top
                    double value = fine_values[fine_grid.index(i_r, i_theta)] +
                                   k2 * fine_values[fine_grid.index(i_r, i_theta_M1)] / (k1 + k2) +
                                   k3 * fine_values[fine_grid.index(i_r, i_theta_P1)] / (k3 + k4);

                    if (i_r_coarse > 0) {
                        // Left Part
                        double h1 = fine_grid.radialSpacing(i_r - 2);
                        double h2 = fine_grid.radialSpacing(i_r - 1);
                        // Left, Bottom Left, Top Left
                        value += h2 * fine_values[fine_grid.index(i_r - 1, i_theta)] / (h1 + h2) +
                                 h2 * k2 * fine_values[fine_grid.index(i_r - 1, i_theta_M1)] / ((h1 + h2) * (k1 + k2)) +
                                 h2 * k3 * fine_values[fine_grid.index(i_r - 1, i_theta_P1)] / ((h1 + h2) * (k3 + k4));
                    }
                    if (i_r_coarse < coarse_grid.nr() - 1) {
                        // Right Part
                        double h3 = fine_grid.radialSpacing(i_r);
                        double h4 = fine_grid.radialSpacing(i_r + 1);
                        // Right, Bottom Right, Top Right
                        value += h3 * fine_values[fine_grid.index(i_r + 1, i_theta)] / (h3 + h4) +
                                 h3 * k2 * fine_values[fine_grid.index(i_r + 1, i_theta_M1)] / ((h3 + h4) * (k1 + k2)) +
                                 h3 * k3 * fine_values[fine_grid.index(i_r + 1, i_theta_P1)] / ((h3 + h4) * (k3 + k4));
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
                    double h1 = fine_grid.radialSpacing(i_r - 2);
                    double h2 = fine_grid.radialSpacing(i_r - 1);
                    double h3 = fine_grid.radialSpacing(i_r);
                    double h4 = fine_grid.radialSpacing(i_r + 1);

                    int i_theta_M2 = fine_grid.wrapThetaIndex(i_theta - 2);
                    int i_theta_M1 = fine_grid.wrapThetaIndex(i_theta - 1);
                    int i_theta_P1 = fine_grid.wrapThetaIndex(i_theta + 1);
                    double k1      = fine_grid.angularSpacing(i_theta_M2);
                    double k2      = fine_grid.angularSpacing(i_theta_M1);
                    double k3      = fine_grid.angularSpacing(i_theta);
                    double k4      = fine_grid.angularSpacing(i_theta_P1);

                    coarse_result[coarse_grid.index(i_r_coarse, i_theta_coarse)] =
                        // Center
                        fine_values[fine_grid.index(i_r, i_theta)] +
                        // Left, Right, Bottom, Top
                        h2 * fine_values[fine_grid.index(i_r - 1, i_theta)] / (h1 + h2) +
                        h3 * fine_values[fine_grid.index(i_r + 1, i_theta)] / (h3 + h4) +
                        k2 * fine_values[fine_grid.index(i_r, i_theta_M1)] / (k1 + k2) +
                        k3 * fine_values[fine_grid.index(i_r, i_theta_P1)] / (k3 + k4) +
                        // Bottom Left, Bottom Right, Top Left, Top Right
                        h2 * k2 * fine_values[fine_grid.index(i_r - 1, i_theta_M1)] / ((h1 + h2) * (k1 + k2)) +
                        h3 * k2 * fine_values[fine_grid.index(i_r + 1, i_theta_M1)] / ((h3 + h4) * (k1 + k2)) +
                        h2 * k3 * fine_values[fine_grid.index(i_r - 1, i_theta_P1)] / ((h1 + h2) * (k3 + k4)) +
                        h3 * k3 * fine_values[fine_grid.index(i_r + 1, i_theta_P1)] / ((h3 + h4) * (k3 + k4));
                }
                else {
                    /* First and Last radial nodes have to be checked for domain boundary */
                    // Middle Part
                    int i_theta_M2 = fine_grid.wrapThetaIndex(i_theta - 2);
                    int i_theta_M1 = fine_grid.wrapThetaIndex(i_theta - 1);
                    int i_theta_P1 = fine_grid.wrapThetaIndex(i_theta + 1);
                    double k1      = fine_grid.angularSpacing(i_theta_M2);
                    double k2      = fine_grid.angularSpacing(i_theta_M1);
                    double k3      = fine_grid.angularSpacing(i_theta);
                    double k4      = fine_grid.angularSpacing(i_theta_P1);
                    // Center, Bottom, Top
                    double value = fine_values[fine_grid.index(i_r, i_theta)] +
                                   k2 * fine_values[fine_grid.index(i_r, i_theta_M1)] / (k1 + k2) +
                                   k3 * fine_values[fine_grid.index(i_r, i_theta_P1)] / (k3 + k4);
                    if (i_r_coarse > 0) {
                        // Left Part
                        double h1 = fine_grid.radialSpacing(i_r - 2);
                        double h2 = fine_grid.radialSpacing(i_r - 1);
                        // Left, Bottom Left, Top Left
                        value += h2 * fine_values[fine_grid.index(i_r - 1, i_theta)] / (h1 + h2) +
                                 h2 * k2 * fine_values[fine_grid.index(i_r - 1, i_theta_M1)] / ((h1 + h2) * (k1 + k2)) +
                                 h2 * k3 * fine_values[fine_grid.index(i_r - 1, i_theta_P1)] / ((h1 + h2) * (k3 + k4));
                    }
                    if (i_r_coarse < coarse_grid.nr() - 1) {
                        // Right Part
                        double h3 = fine_grid.radialSpacing(i_r);
                        double h4 = fine_grid.radialSpacing(i_r + 1);
                        // Right, Bottom Right, Top Right
                        value += h3 * fine_values[fine_grid.index(i_r + 1, i_theta)] / (h3 + h4) +
                                 h3 * k2 * fine_values[fine_grid.index(i_r + 1, i_theta_M1)] / ((h3 + h4) * (k1 + k2)) +
                                 h3 * k3 * fine_values[fine_grid.index(i_r + 1, i_theta_P1)] / ((h3 + h4) * (k3 + k4));
                    }
                    coarse_result[coarse_grid.index(i_r_coarse, i_theta_coarse)] = value;
                }
            }
        }
    }
}
//clang-format on
