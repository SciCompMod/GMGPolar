#include <gtest/gtest.h>

#include <random>

#include "../../include/GMGPolar/gmgpolar.h"
#include "../../include/Interpolation/interpolation.h"

// Function to generate sample data for vector x using random values with seed
Vector<double> generate_random_sample_data(const PolarGrid& grid, unsigned int seed) {
    Vector<double> x(grid.number_of_nodes());
    std::mt19937 gen(seed);  // Standard mersenne_twister_engine seeded with seed
    std::uniform_real_distribution<double> dist(0.0, 1.0);  // Generate random double between 0 and 1
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = dist(gen);
    }
    return x;
}

/* In src/Interpolation/restriction.cpp the Restriction Operator is implemented with "Take". */
/* Here we test against the "Give" version. */

void applyExtrapolatedRestrictionGive0(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x){
    assert(toLevel.level() == fromLevel.level() + 1);

    const PolarGrid& fineGrid = fromLevel.grid();
    const PolarGrid& coarseGrid = toLevel.grid();

    assert(x.size() == fineGrid.number_of_nodes());
    assert(result.size() == coarseGrid.number_of_nodes());
    
    assign(result, 0.0);

    for(int index = 0; index < fineGrid.number_of_nodes(); index ++){
        std::array<std::pair<double,double>, space_dimension> neighbor_distance;
        std::array<std::pair<int,int>, space_dimension> neighbors;

        MultiIndex fine_node = fineGrid.multiindex(index);

        // Fine node appears in coarse grid
        if(fine_node[0] % 2 == 0 && fine_node[1] % 2 == 0){
            // Input x needs a fine grid index: x[FINE_INDEX]
            // Result needs a coarse grid index: result[COARSE_INDEX]
            MultiIndex coarse_node(fine_node[0] / 2, fine_node[1] / 2);
            result[coarseGrid.index(coarse_node)] += x[index];
        }

        // Fine node between coarse nodes in theta direction
        if(fine_node[0] % 2 == 0 && fine_node[1] % 2 == 1){
            fineGrid.adjacent_neighbor_distances(fine_node, neighbor_distance);
            double k1 = neighbor_distance[1].first;
            double k2 = neighbor_distance[1].second;

            fineGrid.adjacent_neighbors_of(fine_node, neighbors);

            MultiIndex bottom_coarse_node(fine_node[0] / 2, fine_node[1] / 2);
            MultiIndex top_coarse_node(fine_node[0] / 2, (fine_node[1] / 2 + 1) % coarseGrid.ntheta());

            result[coarseGrid.index(bottom_coarse_node)] += x[index] / 2.0;
            result[coarseGrid.index(top_coarse_node)] += x[index] / 2.0;
        }


        // Fine node between coarse nodes in radial direction
        if(fine_node[0] % 2 == 1 && fine_node[1] % 2 == 0){
            fineGrid.adjacent_neighbor_distances(fine_node, neighbor_distance);
            double h1 = neighbor_distance[0].first;
            double h2 = neighbor_distance[0].second;

            fineGrid.adjacent_neighbors_of(fine_node, neighbors);

            MultiIndex left_coarse_node(fine_node[0] / 2, fine_node[1] / 2);
            MultiIndex right_coarse_node(fine_node[0] / 2 + 1, fine_node[1] / 2);

            result[coarseGrid.index(left_coarse_node)] += x[index] / 2.0;
            result[coarseGrid.index(right_coarse_node)] += x[index] / 2.0;
        }


        //Fine node in the center of four coarse nodes
        if(fine_node[0] % 2 == 1 && fine_node[1] % 2 == 1){

            fineGrid.adjacent_neighbor_distances(fine_node, neighbor_distance);
            double h1 = neighbor_distance[0].first;
            double h2 = neighbor_distance[0].second;
            double k1 = neighbor_distance[1].first;
            double k2 = neighbor_distance[1].second;

            fineGrid.adjacent_neighbors_of(fine_node, neighbors);

            MultiIndex bottom_right_coarse_node(fine_node[0] / 2 + 1, fine_node[1] / 2);
            MultiIndex top_left_coarse_node(fine_node[0] / 2, (fine_node[1] / 2 + 1) % coarseGrid.ntheta());

            result[coarseGrid.index(bottom_right_coarse_node)] += x[index] / 2.0;
            result[coarseGrid.index(top_left_coarse_node)] += x[index] / 2.0;
        }
    }
}



TEST(ExtrapolatedRestrictionTest, applyExtrapolatedRestriction) {
    std::vector<double> fine_radii = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> fine_angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    auto finest_grid = std::make_unique<PolarGrid>(fine_radii, fine_angles);
    auto coarse_grid = std::make_unique<PolarGrid>(coarseningGrid(*finest_grid));

    auto finest_levelCache = std::make_unique<LevelCache>(*finest_grid);
    auto coarse_levelCache = std::make_unique<LevelCache>(*coarse_grid);

    Level finest_level(0, std::move(finest_grid), std::move(finest_levelCache));
    Level coarse_level(1, std::move(coarse_grid), std::move(coarse_levelCache));

    const int maxOpenMPThreads = 1;
    const std::vector<int> taskingThreads = {1, 1};
    
    Interpolation interpolation_operator(maxOpenMPThreads, taskingThreads);

    unsigned int seed = 42;
    Vector<double> x = generate_random_sample_data(finest_level.grid(), seed);

    // Apply prolongation to both functions
    Vector<double> result1(coarse_level.grid().number_of_nodes());
    Vector<double> result2(coarse_level.grid().number_of_nodes());
    Vector<double> result3(coarse_level.grid().number_of_nodes());

    interpolation_operator.applyExtrapolatedRestriction0(finest_level, coarse_level, result1, x);
    interpolation_operator.applyExtrapolatedRestriction(finest_level, coarse_level, result2, x);

    applyExtrapolatedRestrictionGive0(finest_level, coarse_level, result3, x);

    ASSERT_EQ(result1.size(), result2.size());
    for (size_t i = 0; i < result1.size(); ++i) {
        ASSERT_DOUBLE_EQ(result1[i], result2[i]);
    }
    ASSERT_EQ(result2.size(), result3.size());
    for (size_t i = 0; i < result2.size(); ++i) {
        ASSERT_DOUBLE_EQ(result2[i], result3[i]);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}