#include <gtest/gtest.h>

#include <random>

#include "../test_tools.h"

#include "../../include/GMGPolar/gmgpolar.h"
#include "../../include/Interpolation/interpolation.h"
#include "../../include/InputFunctions/domainGeometry.h"
#include "../../include/InputFunctions/densityProfileCoefficients.h"

#include "../../include/InputFunctions/DomainGeometry/circularGeometry.h"
#include "../../include/InputFunctions/DensityProfileCoefficients/poissonCoefficients.h"

namespace ExtrapolatedRestrictionTest
{
/* In src/Interpolation/restriction.cpp the Restriction Operator is implemented with "Take". */
/* Here we test against the "Give" version. */

void applyExtrapolatedRestrictionGive0(const PolarGrid& fine_grid, const PolarGrid& coarse_grid,
                                       Vector<double> coarse_result, ConstVector<double> fine_values)
{
    assert(fine_values.size() == static_cast<uint>(fine_grid.numberOfNodes()));
    assert(coarse_result.size() == static_cast<uint>(coarse_grid.numberOfNodes()));

    assign(coarse_result, 0.0);

    for (int index = 0; index < fine_grid.numberOfNodes(); index++) {
        std::array<std::pair<double, double>, space_dimension> neighbor_distance;
        std::array<std::pair<int, int>, space_dimension> neighbors;

        MultiIndex fine_node = fine_grid.multiIndex(index);

        // Fine node appears in coarse grid
        if (fine_node[0] % 2 == 0 && fine_node[1] % 2 == 0) {
            // Input x needs a fine grid index: fine_values[FINE_INDEX]
            // Result needs a coarse grid index: coarse_result[COARSE_INDEX]
            MultiIndex coarse_node(fine_node[0] / 2, fine_node[1] / 2);
            coarse_result[coarse_grid.index(coarse_node)] += fine_values[index];
        }

        // Fine node between coarse nodes in theta direction
        if (fine_node[0] % 2 == 0 && fine_node[1] % 2 == 1) {
            fine_grid.adjacentNeighborDistances(fine_node, neighbor_distance);
            double k1 = neighbor_distance[1].first;
            double k2 = neighbor_distance[1].second;

            fine_grid.adjacentNeighborsOf(fine_node, neighbors);

            MultiIndex bottom_coarse_node(fine_node[0] / 2, fine_node[1] / 2);
            MultiIndex top_coarse_node(fine_node[0] / 2, (fine_node[1] / 2 + 1) % coarse_grid.ntheta());

            coarse_result[coarse_grid.index(bottom_coarse_node)] += fine_values[index] / 2.0;
            coarse_result[coarse_grid.index(top_coarse_node)] += fine_values[index] / 2.0;
        }

        // Fine node between coarse nodes in radial direction
        if (fine_node[0] % 2 == 1 && fine_node[1] % 2 == 0) {
            fine_grid.adjacentNeighborDistances(fine_node, neighbor_distance);
            double h1 = neighbor_distance[0].first;
            double h2 = neighbor_distance[0].second;

            fine_grid.adjacentNeighborsOf(fine_node, neighbors);

            MultiIndex left_coarse_node(fine_node[0] / 2, fine_node[1] / 2);
            MultiIndex right_coarse_node(fine_node[0] / 2 + 1, fine_node[1] / 2);

            coarse_result[coarse_grid.index(left_coarse_node)] += fine_values[index] / 2.0;
            coarse_result[coarse_grid.index(right_coarse_node)] += fine_values[index] / 2.0;
        }

        //Fine node in the center of four coarse nodes
        if (fine_node[0] % 2 == 1 && fine_node[1] % 2 == 1) {

            fine_grid.adjacentNeighborDistances(fine_node, neighbor_distance);
            double h1 = neighbor_distance[0].first;
            double h2 = neighbor_distance[0].second;
            double k1 = neighbor_distance[1].first;
            double k2 = neighbor_distance[1].second;

            fine_grid.adjacentNeighborsOf(fine_node, neighbors);

            MultiIndex bottom_right_coarse_node(fine_node[0] / 2 + 1, fine_node[1] / 2);
            MultiIndex top_left_coarse_node(fine_node[0] / 2, (fine_node[1] / 2 + 1) % coarse_grid.ntheta());

            coarse_result[coarse_grid.index(bottom_right_coarse_node)] += fine_values[index] / 2.0;
            coarse_result[coarse_grid.index(top_left_coarse_node)] += fine_values[index] / 2.0;
        }
    }
}
} // namespace ExtrapolatedRestrictionTest

using namespace ExtrapolatedRestrictionTest;

TEST(ExtrapolatedRestrictionTest, applyExtrapolatedRestriction)
{
    std::vector<double> fine_radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> fine_angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    int maxOpenMPThreads = 16;
    bool DirBC_Interior  = true;

    PolarGrid finest_grid(fine_radii, fine_angles);
    PolarGrid coarse_grid = coarseningGrid(finest_grid);

    Interpolation interpolation_operator(maxOpenMPThreads, DirBC_Interior);

    unsigned int seed = 42;
    Vector<double> x  = generate_random_sample_data(finest_grid, seed, 0.0, 1.0);

    // Apply prolongation to both functions
    Vector<double> result1("result1", coarse_grid.numberOfNodes());
    Vector<double> result2("result2", coarse_grid.numberOfNodes());
    Vector<double> result3("result3", coarse_grid.numberOfNodes());

    interpolation_operator.applyExtrapolatedRestriction0(finest_grid, coarse_grid, result1, x);
    interpolation_operator.applyExtrapolatedRestriction(finest_grid, coarse_grid, result2, x);

    applyExtrapolatedRestrictionGive0(finest_grid, coarse_grid, result3, x);

    ASSERT_EQ(result1.size(), result2.size());
    for (uint i = 0; i < result1.size(); ++i) {
        ASSERT_DOUBLE_EQ(result1[i], result2[i]);
    }
    ASSERT_EQ(result2.size(), result3.size());
    for (uint i = 0; i < result2.size(); ++i) {
        ASSERT_DOUBLE_EQ(result2[i], result3[i]);
    }
}
