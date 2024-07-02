#include <gtest/gtest.h>

#include <random>

#include "../../include/GMGPolar/gmgpolar.h"

// Function to generate sample data for vector x using random values with seed
Vector<double> generate_random_sample_data(const PolarGrid& grid, unsigned int seed) {
    Vector<double> x(grid.number_of_nodes());
    std::mt19937 gen(seed);  // Standard mersenne_twister_engine seeded with seed
    std::uniform_real_distribution<double> dist(-100.0, 100.0); 
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = dist(gen);
    }
    return x;
}


TEST(OperatorATest, applyA) {
    std::vector<double> radii = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    auto grid = std::make_unique<PolarGrid>(radii, angles);

    auto levelCache = std::make_unique<LevelCache>(*grid);

    Level level(0, std::move(grid), std::move(levelCache));

    DomainGeometry domain_geometry;
    SystemParameters system_parameters;
    bool DirBC_Interior = false;
    int maxOpenMPThreads = 1;
    int openMPTaskThreads = 1;
    
    Residual residual_operator(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);

    unsigned int seed = 42;
    Vector<double> x = generate_random_sample_data(level.grid(), seed);
    Vector<double> rhs = generate_random_sample_data(level.grid(), seed);

    Vector<double> result1(level.grid().number_of_nodes());
    Vector<double> result2(level.grid().number_of_nodes());
    Vector<double> result3(level.grid().number_of_nodes());

    residual_operator.computeResidual_V1(result1, rhs, x);
    residual_operator.computeResidual_V2(result2, rhs, x);
    residual_operator.computeResidual_V3(result3, rhs, x);

    ASSERT_EQ(result1.size(), result2.size());
    for (size_t i = 0; i < result1.size(); ++i) {
        ASSERT_NEAR(result1[i], result2[i], 1e-12);
    }

    ASSERT_EQ(result2.size(), result3.size());
    for (size_t i = 0; i < result2.size(); ++i) {
        ASSERT_NEAR(result2[i], result3[i], 1e-12);
    }

    Vector<double> result4(level.grid().number_of_nodes());
    result4 = rhs;
    residual_operator.applyATake0(result4, x, -1.0);

    ASSERT_EQ(result2.size(), result3.size());
    for (size_t i = 0; i < result1.size(); ++i) {
        MultiIndex alpha = level.grid().multiindex(i);
        if(alpha[0] == 0 && !DirBC_Interior) ASSERT_NEAR(result1[i], result4[i], 1e-8);
        else ASSERT_NEAR(result1[i], result4[i], 1e-11);
    }
}