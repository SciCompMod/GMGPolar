#include "../../include/Interpolation/interpolation.h"

/* For the restriction we use R_ex = P_ex^T */

void Interpolation::applyExtrapolatedRestriction0(
    const Level& fromLevel, const Level& toLevel, Vector<double> result,
    const Vector<double> x) const
{
    assert(toLevel.level_depth() == fromLevel.level_depth() + 1);

    const PolarGrid& fineGrid   = fromLevel.grid();
    const PolarGrid& coarseGrid = toLevel.grid();

    assert(x.size() == fineGrid.numberOfNodes());
    assert(result.size() == coarseGrid.numberOfNodes());

#pragma omp parallel for num_threads(threads_per_level_[toLevel.level_depth()])
    for (int index = 0; index < coarseGrid.numberOfNodes(); index++) {
        MultiIndex coarse_node = coarseGrid.multiIndex(index);
        MultiIndex fine_node(2 * coarse_node[0], 2 * coarse_node[1]);

        std::array<std::pair<int, int>, space_dimension> neighbors;

        // Center
        double value = x[fineGrid.index(fine_node)];

        fineGrid.adjacentNeighborsOf(fine_node, neighbors);

        // Left
        if (neighbors[0].first != -1) {
            value += 0.5 * x[neighbors[0].first];
        }

        // Right
        if (neighbors[0].second != -1) {
            value += 0.5 * x[neighbors[0].second];
        }

        // Bottom
        if (neighbors[1].first != -1) {
            value += 0.5 * x[neighbors[1].first];
        }

        // Top
        if (neighbors[1].second != -1) {
            value += 0.5 * x[neighbors[1].second];
        }

        fineGrid.diagonalNeighborsOf(fine_node, neighbors);

        // Bottom Right
        if (neighbors[0].second != -1) {
            value += 0.5 * x[neighbors[0].second];
        }

        // Top Left
        if (neighbors[1].first != -1) {
            value += 0.5 * x[neighbors[1].first];
        }

        result[index] = value;
    }
}

// -------------------------------------- //
// Optimized version of applyRestriction0 //
// -------------------------------------- //

void Interpolation::applyExtrapolatedRestriction(
    const Level& fromLevel, const Level& toLevel, Vector<double> result,
    const Vector<double> x) const
{
    assert(toLevel.level_depth() == fromLevel.level_depth() + 1);

    const PolarGrid& fineGrid   = fromLevel.grid();
    const PolarGrid& coarseGrid = toLevel.grid();

    assert(x.size() == fineGrid.numberOfNodes());
    assert(result.size() == coarseGrid.numberOfNodes());

    const int coarseNumberSmootherCircles = coarseGrid.numberSmootherCircles();

#pragma omp parallel num_threads(threads_per_level_[toLevel.level_depth()]) if (fineGrid.numberOfNodes() > 10'000)
    {
/* For loop matches circular access pattern */
#pragma omp for nowait
        for (int i_r_coarse = 0; i_r_coarse < coarseNumberSmootherCircles; i_r_coarse++) {
            int i_r = i_r_coarse * 2;
            for (int i_theta_coarse = 0; i_theta_coarse < coarseGrid.ntheta(); i_theta_coarse++) {
                int i_theta = i_theta_coarse * 2;

                if (0 < i_r_coarse && i_r_coarse < coarseNumberSmootherCircles - 1) {
                    int i_theta_M1 = fineGrid.wrapThetaIndex(i_theta - 1);
                    int i_theta_P1 = fineGrid.wrapThetaIndex(i_theta + 1);

                    result[coarseGrid.index(i_r_coarse, i_theta_coarse)] =
                        // Center
                        x[fineGrid.index(i_r, i_theta)] +
                        // Left, Right, Bottom, Top
                        0.5 * x[fineGrid.index(i_r - 1, i_theta)] + 0.5 * x[fineGrid.index(i_r + 1, i_theta)] +
                        0.5 * x[fineGrid.index(i_r, i_theta_M1)] + 0.5 * x[fineGrid.index(i_r, i_theta_P1)] +
                        // Bottom Right, Top Left
                        0.5 * x[fineGrid.index(i_r + 1, i_theta_M1)] + 0.5 * x[fineGrid.index(i_r - 1, i_theta_P1)];
                }
                else {
                    /* First and Last Circle have to be checked for domain boundary */
                    int i_theta_M1 = fineGrid.wrapThetaIndex(i_theta - 1);
                    int i_theta_P1 = fineGrid.wrapThetaIndex(i_theta + 1);
                    // Center, Bottom, Top
                    double value = x[fineGrid.index(i_r, i_theta)] + 0.5 * x[fineGrid.index(i_r, i_theta_M1)] +
                                   0.5 * x[fineGrid.index(i_r, i_theta_P1)];

                    if (i_r_coarse > 0) {
                        // Left, Top Left
                        value +=
                            0.5 * x[fineGrid.index(i_r - 1, i_theta)] + 0.5 * x[fineGrid.index(i_r - 1, i_theta_P1)];
                    }
                    if (i_r_coarse < coarseGrid.nr() - 1) {
                        // Right, Bottom Right
                        value +=
                            0.5 * x[fineGrid.index(i_r + 1, i_theta)] + 0.5 * x[fineGrid.index(i_r + 1, i_theta_M1)];
                    }
                    result[coarseGrid.index(i_r_coarse, i_theta_coarse)] = value;
                }
            }
        }

/* For loop matches circular access pattern */
#pragma omp for nowait
        for (int i_theta_coarse = 0; i_theta_coarse < coarseGrid.ntheta(); i_theta_coarse++) {
            int i_theta = i_theta_coarse * 2;
            for (int i_r_coarse = coarseNumberSmootherCircles; i_r_coarse < coarseGrid.nr(); i_r_coarse++) {
                int i_r = i_r_coarse * 2;

                if (coarseGrid.numberSmootherCircles() < i_r_coarse && i_r_coarse < coarseGrid.nr() - 1) {
                    int i_theta_M1 = fineGrid.wrapThetaIndex(i_theta - 1);
                    int i_theta_P1 = fineGrid.wrapThetaIndex(i_theta + 1);

                    result[coarseGrid.index(i_r_coarse, i_theta_coarse)] =
                        // Center
                        x[fineGrid.index(i_r, i_theta)] +
                        // Left, Right, Bottom, Top
                        0.5 * x[fineGrid.index(i_r - 1, i_theta)] + 0.5 * x[fineGrid.index(i_r + 1, i_theta)] +
                        0.5 * x[fineGrid.index(i_r, i_theta_M1)] + 0.5 * x[fineGrid.index(i_r, i_theta_P1)] +
                        // Bottom Right, Top Left
                        0.5 * x[fineGrid.index(i_r + 1, i_theta_M1)] + 0.5 * x[fineGrid.index(i_r - 1, i_theta_P1)];
                }
                else {
                    /* First and Last radial nodes have to be checked for domain boundary */
                    int i_theta_M1 = fineGrid.wrapThetaIndex(i_theta - 1);
                    int i_theta_P1 = fineGrid.wrapThetaIndex(i_theta + 1);
                    // Center, Bottom, Top
                    double value = x[fineGrid.index(i_r, i_theta)] + 0.5 * x[fineGrid.index(i_r, i_theta_M1)] +
                                   0.5 * x[fineGrid.index(i_r, i_theta_P1)];
                    if (i_r_coarse > 0) {
                        // Left, Top Left
                        value +=
                            0.5 * x[fineGrid.index(i_r - 1, i_theta)] + 0.5 * x[fineGrid.index(i_r - 1, i_theta_P1)];
                    }
                    if (i_r_coarse < coarseGrid.nr() - 1) {
                        // Right, Bottom Right
                        value +=
                            0.5 * x[fineGrid.index(i_r + 1, i_theta)] + 0.5 * x[fineGrid.index(i_r + 1, i_theta_M1)];
                    }
                    result[coarseGrid.index(i_r_coarse, i_theta_coarse)] = value;
                }
            }
        }
    }
}
