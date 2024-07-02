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


TEST(DirectSolverTest, DirectSolverDirBC_Interior) {
    std::vector<double> radii = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    auto grid = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid);
    Level level(0, std::move(grid), std::move(levelCache));

    DomainGeometry domain_geometry;
    SystemParameters system_parameters;
    bool DirBC_Interior = true;
    int maxOpenMPThreads = 1;
    int openMPTaskThreads = 1;
    
    DirectSolver solver_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    Residual residual_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);

    const Vector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution = rhs;
    solver_op.solveInPlace(solution);

    Vector<double> residuum(level.grid().number_of_nodes());
    residual_op.computeResidual_V1(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(residuum), 0.0, 1e-11);
    ASSERT_NEAR(sqrt(l2_norm_squared(residuum)), 0.0, 1e-12);
    ASSERT_NEAR(infinity_norm(residuum), 0.0, 1e-12);
}



TEST(DirectSolverTest, DirectSolverAcrossOrigin) {
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
    
    DirectSolver solver_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    Residual residual_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);

    const Vector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution = rhs;
    solver_op.solveInPlace(solution);

    Vector<double> residuum(level.grid().number_of_nodes());
    residual_op.computeResidual_V1(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(residuum), 0.0, 1e-8);
    ASSERT_NEAR(sqrt(l2_norm_squared(residuum)), 0.0, 1e-9);
    ASSERT_NEAR(infinity_norm(residuum), 0.0, 1e-9);
}