#include "../../include/Interpolation/interpolation.h"

// For the restriction we use R = P^T.
// The restriction for an siotropic mesh reduces to
//           |1  2  1|
// R = 1/4 * |2  4  2| = P^T
//           |1  2  1|

void Interpolation::applyRestriction0(const Level& fromLevel, const Level& toLevel, Vector<double>& result,
                                      const Vector<double>& x) const
{
    assert(toLevel.level() == fromLevel.level() + 1);

    omp_set_num_threads(threads_per_level_[toLevel.level()]);

    const PolarGrid& fineGrid   = fromLevel.grid();
    const PolarGrid& coarseGrid = toLevel.grid();

    assert(x.size() == fineGrid.numberOfNodes());
    assert(result.size() == coarseGrid.numberOfNodes());

#pragma omp parallel for
    for (int index = 0; index < coarseGrid.numberOfNodes(); index++) {
        MultiIndex coarse_node = coarseGrid.multiIndex(index);
        MultiIndex fine_node(2 * coarse_node[0], 2 * coarse_node[1]);

        std::array<std::pair<double, double>, space_dimension> neighbor_distance;
        std::array<std::pair<int, int>, space_dimension> neighbors;

        // Center
        double value = x[fineGrid.index(fine_node)];

        fineGrid.adjacentNeighborsOf(fine_node, neighbors);

        // Left
        if (neighbors[0].first != -1) {
            fineGrid.adjacentNeighborDistances(fineGrid.multiIndex(neighbors[0].first), neighbor_distance);
            value += neighbor_distance[0].second * x[neighbors[0].first] /
                     (neighbor_distance[0].first + neighbor_distance[0].second);
        }

        // Right
        if (neighbors[0].second != -1) {
            fineGrid.adjacentNeighborDistances(fineGrid.multiIndex(neighbors[0].second), neighbor_distance);
            value += neighbor_distance[0].first * x[neighbors[0].second] /
                     (neighbor_distance[0].first + neighbor_distance[0].second);
        }

        // Bottom
        if (neighbors[1].first != -1) {
            fineGrid.adjacentNeighborDistances(fineGrid.multiIndex(neighbors[1].first), neighbor_distance);
            value += neighbor_distance[1].second * x[neighbors[1].first] /
                     (neighbor_distance[1].first + neighbor_distance[1].second);
        }

        // Top
        if (neighbors[1].second != -1) {
            fineGrid.adjacentNeighborDistances(fineGrid.multiIndex(neighbors[1].second), neighbor_distance);
            value += neighbor_distance[1].first * x[neighbors[1].second] /
                     (neighbor_distance[1].first + neighbor_distance[1].second);
        }

        fineGrid.diagonalNeighborsOf(fine_node, neighbors);

        // Bottom Left
        if (neighbors[0].first != -1) {
            fineGrid.adjacentNeighborDistances(fineGrid.multiIndex(neighbors[0].first), neighbor_distance);
            value += neighbor_distance[0].second * neighbor_distance[1].second * x[neighbors[0].first] /
                     ((neighbor_distance[0].first + neighbor_distance[0].second) *
                      (neighbor_distance[1].first + neighbor_distance[1].second));
        }

        // Bottom Right
        if (neighbors[0].second != -1) {
            fineGrid.adjacentNeighborDistances(fineGrid.multiIndex(neighbors[0].second), neighbor_distance);
            value += neighbor_distance[0].first * neighbor_distance[1].second * x[neighbors[0].second] /
                     ((neighbor_distance[0].first + neighbor_distance[0].second) *
                      (neighbor_distance[1].first + neighbor_distance[1].second));
        }

        // Top Left
        if (neighbors[1].first != -1) {
            fineGrid.adjacentNeighborDistances(fineGrid.multiIndex(neighbors[1].first), neighbor_distance);
            value += neighbor_distance[0].second * neighbor_distance[1].first * x[neighbors[1].first] /
                     ((neighbor_distance[0].first + neighbor_distance[0].second) *
                      (neighbor_distance[1].first + neighbor_distance[1].second));
        }

        // Top Right
        if (neighbors[1].second != -1) {
            fineGrid.adjacentNeighborDistances(fineGrid.multiIndex(neighbors[1].second), neighbor_distance);
            value += neighbor_distance[0].first * neighbor_distance[1].first * x[neighbors[1].second] /
                     ((neighbor_distance[0].first + neighbor_distance[0].second) *
                      (neighbor_distance[1].first + neighbor_distance[1].second));
        }

        result[index] = value;
    }
}

// -------------------------------------- //
// Optimized version of applyRestriction0 //
// -------------------------------------- //

void Interpolation::applyRestriction(const Level& fromLevel, const Level& toLevel, Vector<double>& result,
                                     const Vector<double>& x) const
{
    assert(toLevel.level() == fromLevel.level() + 1);

    omp_set_num_threads(threads_per_level_[toLevel.level()]);

    const PolarGrid& fineGrid   = fromLevel.grid();
    const PolarGrid& coarseGrid = toLevel.grid();

    assert(x.size() == fineGrid.numberOfNodes());
    assert(result.size() == coarseGrid.numberOfNodes());

    const int coarseNumberSmootherCircles = coarseGrid.numberSmootherCircles();

#pragma omp parallel if (fineGrid.numberOfNodes() > 10'000)
    {
/* For loop matches circular access pattern */
#pragma omp for nowait
        for (int i_r_coarse = 0; i_r_coarse < coarseNumberSmootherCircles; i_r_coarse++) {
            int i_r = i_r_coarse << 1;
            for (int i_theta_coarse = 0; i_theta_coarse < coarseGrid.ntheta(); i_theta_coarse++) {
                int i_theta = i_theta_coarse << 1;

                if (0 < i_r_coarse && i_r_coarse < coarseNumberSmootherCircles - 1) {
                    double h1 = fineGrid.radialSpacing(i_r - 2);
                    double h2 = fineGrid.radialSpacing(i_r - 1);
                    double h3 = fineGrid.radialSpacing(i_r);
                    double h4 = fineGrid.radialSpacing(i_r + 1);

                    int i_theta_M2 = fineGrid.wrapThetaIndex(i_theta - 2);
                    int i_theta_M1 = fineGrid.wrapThetaIndex(i_theta - 1);
                    int i_theta_P1 = fineGrid.wrapThetaIndex(i_theta + 1);

                    double k1 = fineGrid.angularSpacing(i_theta_M2);
                    double k2 = fineGrid.angularSpacing(i_theta_M1);
                    double k3 = fineGrid.angularSpacing(i_theta);
                    double k4 = fineGrid.angularSpacing(i_theta_P1);

                    result[coarseGrid.index(i_r_coarse, i_theta_coarse)] =
                        // Center
                        x[fineGrid.index(i_r, i_theta)] +
                        // Left, Right, Bottom, Top
                        h2 * x[fineGrid.index(i_r - 1, i_theta)] / (h1 + h2) +
                        h3 * x[fineGrid.index(i_r + 1, i_theta)] / (h3 + h4) +
                        k2 * x[fineGrid.index(i_r, i_theta_M1)] / (k1 + k2) +
                        k3 * x[fineGrid.index(i_r, i_theta_P1)] / (k3 + k4) +
                        // Bottom Left, Bottom Right, Top Left, Top Right
                        h2 * k2 * x[fineGrid.index(i_r - 1, i_theta_M1)] / ((h1 + h2) * (k1 + k2)) +
                        h3 * k2 * x[fineGrid.index(i_r + 1, i_theta_M1)] / ((h3 + h4) * (k1 + k2)) +
                        h2 * k3 * x[fineGrid.index(i_r - 1, i_theta_P1)] / ((h1 + h2) * (k3 + k4)) +
                        h3 * k3 * x[fineGrid.index(i_r + 1, i_theta_P1)] / ((h3 + h4) * (k3 + k4));
                }
                else {
                    /* First and Last Circle have to be checked for domain boundary */
                    // Middle Part
                    int i_theta_M2 = fineGrid.wrapThetaIndex(i_theta - 2);
                    int i_theta_M1 = fineGrid.wrapThetaIndex(i_theta - 1);
                    int i_theta_P1 = fineGrid.wrapThetaIndex(i_theta + 1);
                    double k1      = fineGrid.angularSpacing(i_theta_M2);
                    double k2      = fineGrid.angularSpacing(i_theta_M1);
                    double k3      = fineGrid.angularSpacing(i_theta);
                    double k4      = fineGrid.angularSpacing(i_theta_P1);
                    // Center, Bottom, Top
                    double value = x[fineGrid.index(i_r, i_theta)] +
                                   k2 * x[fineGrid.index(i_r, i_theta_M1)] / (k1 + k2) +
                                   k3 * x[fineGrid.index(i_r, i_theta_P1)] / (k3 + k4);

                    if (i_r_coarse > 0) {
                        // Left Part
                        double h1 = fineGrid.radialSpacing(i_r - 2);
                        double h2 = fineGrid.radialSpacing(i_r - 1);
                        // Left, Bottom Left, Top Left
                        value += h2 * x[fineGrid.index(i_r - 1, i_theta)] / (h1 + h2) +
                                 h2 * k2 * x[fineGrid.index(i_r - 1, i_theta_M1)] / ((h1 + h2) * (k1 + k2)) +
                                 h2 * k3 * x[fineGrid.index(i_r - 1, i_theta_P1)] / ((h1 + h2) * (k3 + k4));
                    }
                    if (i_r_coarse < coarseGrid.nr() - 1) {
                        // Right Part
                        double h3 = fineGrid.radialSpacing(i_r);
                        double h4 = fineGrid.radialSpacing(i_r + 1);
                        // Right, Bottom Right, Top Right
                        value += h3 * x[fineGrid.index(i_r + 1, i_theta)] / (h3 + h4) +
                                 h3 * k2 * x[fineGrid.index(i_r + 1, i_theta_M1)] / ((h3 + h4) * (k1 + k2)) +
                                 h3 * k3 * x[fineGrid.index(i_r + 1, i_theta_P1)] / ((h3 + h4) * (k3 + k4));
                    }
                    result[coarseGrid.index(i_r_coarse, i_theta_coarse)] = value;
                }
            }
        }

/* For loop matches circular access pattern */
#pragma omp for nowait
        for (int i_theta_coarse = 0; i_theta_coarse < coarseGrid.ntheta(); i_theta_coarse++) {
            int i_theta = i_theta_coarse << 1;
            for (int i_r_coarse = coarseNumberSmootherCircles; i_r_coarse < coarseGrid.nr(); i_r_coarse++) {
                int i_r = i_r_coarse << 1;

                if (coarseGrid.numberSmootherCircles() < i_r_coarse && i_r_coarse < coarseGrid.nr() - 1) {
                    double h1 = fineGrid.radialSpacing(i_r - 2);
                    double h2 = fineGrid.radialSpacing(i_r - 1);
                    double h3 = fineGrid.radialSpacing(i_r);
                    double h4 = fineGrid.radialSpacing(i_r + 1);

                    int i_theta_M2 = fineGrid.wrapThetaIndex(i_theta - 2);
                    int i_theta_M1 = fineGrid.wrapThetaIndex(i_theta - 1);
                    int i_theta_P1 = fineGrid.wrapThetaIndex(i_theta + 1);
                    double k1      = fineGrid.angularSpacing(i_theta_M2);
                    double k2      = fineGrid.angularSpacing(i_theta_M1);
                    double k3      = fineGrid.angularSpacing(i_theta);
                    double k4      = fineGrid.angularSpacing(i_theta_P1);

                    result[coarseGrid.index(i_r_coarse, i_theta_coarse)] =
                        // Center
                        x[fineGrid.index(i_r, i_theta)] +
                        // Left, Right, Bottom, Top
                        h2 * x[fineGrid.index(i_r - 1, i_theta)] / (h1 + h2) +
                        h3 * x[fineGrid.index(i_r + 1, i_theta)] / (h3 + h4) +
                        k2 * x[fineGrid.index(i_r, i_theta_M1)] / (k1 + k2) +
                        k3 * x[fineGrid.index(i_r, i_theta_P1)] / (k3 + k4) +
                        // Bottom Left, Bottom Right, Top Left, Top Right
                        h2 * k2 * x[fineGrid.index(i_r - 1, i_theta_M1)] / ((h1 + h2) * (k1 + k2)) +
                        h3 * k2 * x[fineGrid.index(i_r + 1, i_theta_M1)] / ((h3 + h4) * (k1 + k2)) +
                        h2 * k3 * x[fineGrid.index(i_r - 1, i_theta_P1)] / ((h1 + h2) * (k3 + k4)) +
                        h3 * k3 * x[fineGrid.index(i_r + 1, i_theta_P1)] / ((h3 + h4) * (k3 + k4));
                }
                else {
                    /* First and Last radial nodes have to be checked for domain boundary */
                    // Middle Part
                    int i_theta_M2 = fineGrid.wrapThetaIndex(i_theta - 2);
                    int i_theta_M1 = fineGrid.wrapThetaIndex(i_theta - 1);
                    int i_theta_P1 = fineGrid.wrapThetaIndex(i_theta + 1);
                    double k1      = fineGrid.angularSpacing(i_theta_M2);
                    double k2      = fineGrid.angularSpacing(i_theta_M1);
                    double k3      = fineGrid.angularSpacing(i_theta);
                    double k4      = fineGrid.angularSpacing(i_theta_P1);
                    // Center, Bottom, Top
                    double value = x[fineGrid.index(i_r, i_theta)] +
                                   k2 * x[fineGrid.index(i_r, i_theta_M1)] / (k1 + k2) +
                                   k3 * x[fineGrid.index(i_r, i_theta_P1)] / (k3 + k4);
                    if (i_r_coarse > 0) {
                        // Left Part
                        double h1 = fineGrid.radialSpacing(i_r - 2);
                        double h2 = fineGrid.radialSpacing(i_r - 1);
                        // Left, Bottom Left, Top Left
                        value += h2 * x[fineGrid.index(i_r - 1, i_theta)] / (h1 + h2) +
                                 h2 * k2 * x[fineGrid.index(i_r - 1, i_theta_M1)] / ((h1 + h2) * (k1 + k2)) +
                                 h2 * k3 * x[fineGrid.index(i_r - 1, i_theta_P1)] / ((h1 + h2) * (k3 + k4));
                    }
                    if (i_r_coarse < coarseGrid.nr() - 1) {
                        // Right Part
                        double h3 = fineGrid.radialSpacing(i_r);
                        double h4 = fineGrid.radialSpacing(i_r + 1);
                        // Right, Bottom Right, Top Right
                        value += h3 * x[fineGrid.index(i_r + 1, i_theta)] / (h3 + h4) +
                                 h3 * k2 * x[fineGrid.index(i_r + 1, i_theta_M1)] / ((h3 + h4) * (k1 + k2)) +
                                 h3 * k3 * x[fineGrid.index(i_r + 1, i_theta_P1)] / ((h3 + h4) * (k3 + k4));
                    }
                    result[coarseGrid.index(i_r_coarse, i_theta_coarse)] = value;
                }
            }
        }
    }
}
//clang-format on