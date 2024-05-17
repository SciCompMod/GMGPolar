#include <gtest/gtest.h>
#include "../include/GMGPolar/gmgpolar.h"
#include "../include/Operator/operator.h"

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


TEST(RestrictionTest, RestrictionSmoothingRadius) {
    std::vector<double> fine_radii = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> fine_angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};
    double smoothing_radius = 0.9;
    PolarGrid fine_grid(fine_radii, fine_angles, smoothing_radius);
    auto coarse_grid = coarseningGrid(fine_grid);

    Level fromLevel(0);
    fromLevel.setGrid(std::make_unique<PolarGrid>(std::move(fine_grid)));

    Level toLevel(1);
    toLevel.setGrid(std::move(coarse_grid));

    unsigned int seed = 42;
    Vector<double> x = generate_random_sample_data(fromLevel.grid(), seed);

    // Apply prolongation to both functions
    Vector<double> result1(toLevel.grid().number_of_nodes());
    Vector<double> result2(toLevel.grid().number_of_nodes());
    Vector<double> result3(toLevel.grid().number_of_nodes());
    Vector<double> result4(toLevel.grid().number_of_nodes());
    Vector<double> result5(toLevel.grid().number_of_nodes());
    Vector<double> result6(toLevel.grid().number_of_nodes());

    std::unique_ptr<Operator> op = std::make_unique<Operator>();

    op->applyRestrictionTake0(fromLevel, toLevel, result1, x);
    op->applyRestrictionTake(fromLevel, toLevel, result2, x);
    op->applyRestrictionTakeTasks(fromLevel, toLevel, result3, x);
    op->applyRestrictionGive0(fromLevel, toLevel, result4, x);
    op->applyRestrictionGive(fromLevel, toLevel, result5, x);
    op->applyRestrictionGiveTasks(fromLevel, toLevel, result6, x);

    ASSERT_EQ(result1.size(), result2.size());
    ASSERT_EQ(result1.size(), result3.size());
    ASSERT_EQ(result1.size(), result4.size());
    ASSERT_EQ(result1.size(), result5.size());
    ASSERT_EQ(result1.size(), result6.size());

    for (size_t i = 0; i < result1.size(); ++i) {
        ASSERT_DOUBLE_EQ(result1[i], result2[i]);
        ASSERT_DOUBLE_EQ(result1[i], result3[i]);
        ASSERT_DOUBLE_EQ(result1[i], result4[i]);
        ASSERT_DOUBLE_EQ(result1[i], result5[i]);
        ASSERT_DOUBLE_EQ(result1[i], result6[i]);
    }
}

TEST(RestrictionTest, RestrictionRadialIndexing) {
    std::vector<double> fine_radii = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> fine_angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};
    double smoothing_radius = -1.0;
    PolarGrid fine_grid(fine_radii, fine_angles, smoothing_radius);
    auto coarse_grid = coarseningGrid(fine_grid);

    Level fromLevel(0);
    fromLevel.setGrid(std::make_unique<PolarGrid>(std::move(fine_grid)));

    Level toLevel(1);
    toLevel.setGrid(std::move(coarse_grid));

    unsigned int seed = 42;
    Vector<double> x = generate_random_sample_data(fromLevel.grid(), seed);

    // Apply prolongation to both functions
    Vector<double> result1(toLevel.grid().number_of_nodes());
    Vector<double> result2(toLevel.grid().number_of_nodes());
    Vector<double> result3(toLevel.grid().number_of_nodes());
    Vector<double> result4(toLevel.grid().number_of_nodes());
    Vector<double> result5(toLevel.grid().number_of_nodes());
    Vector<double> result6(toLevel.grid().number_of_nodes());

    std::unique_ptr<Operator> op = std::make_unique<Operator>();

    op->applyRestrictionTake0(fromLevel, toLevel, result1, x);
    op->applyRestrictionTake(fromLevel, toLevel, result2, x);
    op->applyRestrictionTakeTasks(fromLevel, toLevel, result3, x);
    op->applyRestrictionGive0(fromLevel, toLevel, result4, x);
    op->applyRestrictionGive(fromLevel, toLevel, result5, x);
    op->applyRestrictionGiveTasks(fromLevel, toLevel, result6, x);

    ASSERT_EQ(result1.size(), result2.size());
    ASSERT_EQ(result1.size(), result3.size());
    ASSERT_EQ(result1.size(), result4.size());
    ASSERT_EQ(result1.size(), result5.size());
    ASSERT_EQ(result1.size(), result6.size());

    for (size_t i = 0; i < result1.size(); ++i) {
        ASSERT_DOUBLE_EQ(result1[i], result2[i]);
        ASSERT_DOUBLE_EQ(result1[i], result3[i]);
        ASSERT_DOUBLE_EQ(result1[i], result4[i]);
        ASSERT_DOUBLE_EQ(result1[i], result5[i]);
        ASSERT_DOUBLE_EQ(result1[i], result6[i]);
    }
}

TEST(RestrictionTest, RestrictionCircularIndexing) {
    std::vector<double> fine_radii = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> fine_angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};
    double smoothing_radius = fine_radii.back() + 1.0;
    PolarGrid fine_grid(fine_radii, fine_angles, smoothing_radius);
    auto coarse_grid = coarseningGrid(fine_grid);

    Level fromLevel(0);
    fromLevel.setGrid(std::make_unique<PolarGrid>(std::move(fine_grid)));

    Level toLevel(1);
    toLevel.setGrid(std::move(coarse_grid));

    unsigned int seed = 42;
    Vector<double> x = generate_random_sample_data(fromLevel.grid(), seed);

    // Apply prolongation to both functions
    Vector<double> result1(toLevel.grid().number_of_nodes());
    Vector<double> result2(toLevel.grid().number_of_nodes());
    Vector<double> result3(toLevel.grid().number_of_nodes());
    Vector<double> result4(toLevel.grid().number_of_nodes());
    Vector<double> result5(toLevel.grid().number_of_nodes());
    Vector<double> result6(toLevel.grid().number_of_nodes());

    std::unique_ptr<Operator> op = std::make_unique<Operator>();

    op->applyRestrictionTake0(fromLevel, toLevel, result1, x);
    op->applyRestrictionTake(fromLevel, toLevel, result2, x);
    op->applyRestrictionTakeTasks(fromLevel, toLevel, result3, x);
    op->applyRestrictionGive0(fromLevel, toLevel, result4, x);
    op->applyRestrictionGive(fromLevel, toLevel, result5, x);
    op->applyRestrictionGiveTasks(fromLevel, toLevel, result6, x);

    ASSERT_EQ(result1.size(), result2.size());
    ASSERT_EQ(result1.size(), result3.size());
    ASSERT_EQ(result1.size(), result4.size());
    ASSERT_EQ(result1.size(), result5.size());
    ASSERT_EQ(result1.size(), result6.size());

    for (size_t i = 0; i < result1.size(); ++i) {
        ASSERT_DOUBLE_EQ(result1[i], result2[i]);
        ASSERT_DOUBLE_EQ(result1[i], result3[i]);
        ASSERT_DOUBLE_EQ(result1[i], result4[i]);
        ASSERT_DOUBLE_EQ(result1[i], result5[i]);
        ASSERT_DOUBLE_EQ(result1[i], result6[i]);
    }
}


TEST(RestrictionTest, RestrictionAutomaticIndexing) {
    std::vector<double> fine_radii = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> fine_angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};
    double smoothing_radius = fine_radii.back() + 1.0;
    PolarGrid fine_grid(fine_radii, fine_angles, smoothing_radius);
    auto coarse_grid = coarseningGrid(fine_grid);

    Level fromLevel(0);
    fromLevel.setGrid(std::make_unique<PolarGrid>(std::move(fine_grid)));

    Level toLevel(1);
    toLevel.setGrid(std::move(coarse_grid));

    unsigned int seed = 42;
    Vector<double> x = generate_random_sample_data(fromLevel.grid(), seed);

    // Apply prolongation to both functions
    Vector<double> result1(toLevel.grid().number_of_nodes());
    Vector<double> result2(toLevel.grid().number_of_nodes());
    Vector<double> result3(toLevel.grid().number_of_nodes());
    Vector<double> result4(toLevel.grid().number_of_nodes());
    Vector<double> result5(toLevel.grid().number_of_nodes());
    Vector<double> result6(toLevel.grid().number_of_nodes());

    std::unique_ptr<Operator> op = std::make_unique<Operator>();

    op->applyRestrictionTake0(fromLevel, toLevel, result1, x);
    op->applyRestrictionTake(fromLevel, toLevel, result2, x);
    op->applyRestrictionTakeTasks(fromLevel, toLevel, result3, x);
    op->applyRestrictionGive0(fromLevel, toLevel, result4, x);
    op->applyRestrictionGive(fromLevel, toLevel, result5, x);
    op->applyRestrictionGiveTasks(fromLevel, toLevel, result6, x);

    ASSERT_EQ(result1.size(), result2.size());
    ASSERT_EQ(result1.size(), result3.size());
    ASSERT_EQ(result1.size(), result4.size());
    ASSERT_EQ(result1.size(), result5.size());
    ASSERT_EQ(result1.size(), result6.size());

    for (size_t i = 0; i < result1.size(); ++i) {
        ASSERT_DOUBLE_EQ(result1[i], result2[i]);
        ASSERT_DOUBLE_EQ(result1[i], result3[i]);
        ASSERT_DOUBLE_EQ(result1[i], result4[i]);
        ASSERT_DOUBLE_EQ(result1[i], result5[i]);
        ASSERT_DOUBLE_EQ(result1[i], result6[i]);
    }
}