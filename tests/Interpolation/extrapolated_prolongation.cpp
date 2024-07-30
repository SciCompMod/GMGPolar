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

TEST(ExtrapolatedProlongationTest, ExtrapolatedProlongationSmoothingRadius) {
    std::vector<double> fine_radii = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> fine_angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    auto finest_grid = std::make_unique<PolarGrid>(fine_radii, fine_angles);
    auto coarse_grid = std::make_unique<PolarGrid>(coarseningGrid(*finest_grid));

    auto finest_levelCache = std::make_unique<LevelCache>(*finest_grid);
    auto coarse_levelCache = std::make_unique<LevelCache>(*coarse_grid);

    Level finest_level(0, std::move(finest_grid), std::move(finest_levelCache));
    Level coarse_level(1, std::move(coarse_grid), std::move(coarse_levelCache));

    const int maxOpenMPThreads = 5;
    const std::vector<int> taskingThreads = {5,5};
    
    Interpolation interpolation_operator(maxOpenMPThreads, taskingThreads);

    unsigned int seed = 42;
    Vector<double> x = generate_random_sample_data(coarse_level.grid(), seed);

    // Apply prolongation to both functions
    Vector<double> result1(finest_level.grid().number_of_nodes());
    Vector<double> result2(finest_level.grid().number_of_nodes());

    interpolation_operator.applyExtrapolatedProlongation0(coarse_level, finest_level, result1, x);
    interpolation_operator.applyExtrapolatedProlongation(coarse_level, finest_level, result2, x);

    ASSERT_EQ(result1.size(), result2.size());
    for (size_t i = 0; i < result1.size(); ++i) {
        ASSERT_DOUBLE_EQ(result1[i], result2[i]);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}