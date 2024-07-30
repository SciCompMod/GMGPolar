#include <gtest/gtest.h>

#include "../../include/GMGPolar/gmgpolar.h"

#include "../../include/InputFunctions/domainGeometry.h"
#include "../../include/InputFunctions/systemParameters.h"
/* --------- */
/* Test Case */
/* --------- */
#include "../include/InputFunctions/DomainGeometry/czarnyGeometry.h"
#include "../include/InputFunctions/BoundaryConditions/polarR6_Boundary_CzarnyGeometry.h"
#include "../include/InputFunctions/DensityProfileCoefficients/zoniShiftedCoefficients.h"
#include "../include/InputFunctions/SourceTerms/polarR6_ZoniShifted_CzarnyGeometry.h"

#include <random>

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

TEST(SmootherTest, SmootherDirBC_Interior) {
    std::vector<double> radii = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    auto grid = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid);
    Level level(0, std::move(grid), std::move(levelCache));

    double Rmax = radii.back();
    double kappa_eps=0.3;
    double delta_e=1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    SystemParameters system_parameters(std::move(coefficients), std::move(boundary_conditions), std::move(source_term));

    bool DirBC_Interior = true;
    int maxOpenMPThreads = 1;
    int openMPTaskThreads = 1;
    
    DirectSolver solver_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    Residual residual_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    ExtrapolatedSmoother extrapolated_smoother_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);

    const Vector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    /* Set coarse values to the solution and fill the rest with random values. */
    Vector<double> smoother_solution = generate_random_sample_data(level.grid(), 42);
    for (int i_r = 0; i_r < level.grid().nr(); i_r+=2){
        for (int i_theta = 0; i_theta < level.grid().ntheta(); i_theta+=2){
           smoother_solution[level.grid().index(i_r,i_theta)] = discrete_solution[level.grid().index(i_r,i_theta)];
        }
    }

    Vector<double> temp(level.grid().number_of_nodes());
    
    for (int i = 0; i < 1000; i++){
        extrapolated_smoother_op.extrapolatedSmoothingInPlace(smoother_solution, rhs, temp);
    }

    Vector<double> error(level.grid().number_of_nodes());
    for (int i = 0; i < error.size(); i++){
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    ASSERT_NEAR(l1_norm(error), 0.0, 1e-12);
    ASSERT_NEAR(sqrt(l2_norm_squared(error)), 0.0, 1e-13);
    ASSERT_NEAR(infinity_norm(error), 0.0, 1e-13);
}

TEST(ExtrapolatedSmootherTest, ExtrapolatedSmootherAcrossOrigin) {
    std::vector<double> radii = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    auto grid = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid);
    Level level(0, std::move(grid), std::move(levelCache));

    double Rmax = radii.back();
    double kappa_eps=0.3;
    double delta_e=1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    SystemParameters system_parameters(std::move(coefficients), std::move(boundary_conditions), std::move(source_term));

    bool DirBC_Interior = false;
    int maxOpenMPThreads = 1;
    int openMPTaskThreads = 1;
    
    DirectSolver solver_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    Residual residual_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    ExtrapolatedSmoother extrapolated_smoother_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);

    const Vector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    /* Set coarse values to the solution and fill the rest with random values. */
    Vector<double> smoother_solution = generate_random_sample_data(level.grid(), 42);
    for (int i_r = 0; i_r < level.grid().nr(); i_r+=2){
        for (int i_theta = 0; i_theta < level.grid().ntheta(); i_theta+=2){
           smoother_solution[level.grid().index(i_r,i_theta)] = discrete_solution[level.grid().index(i_r,i_theta)];
        }
    }

    Vector<double> temp(level.grid().number_of_nodes());
    
    for (int i = 0; i < 1000; i++){
        extrapolated_smoother_op.extrapolatedSmoothingInPlace(smoother_solution, rhs, temp);
    }
    
    Vector<double> error(level.grid().number_of_nodes());
    for (int i = 0; i < error.size(); i++){
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    ASSERT_NEAR(l1_norm(error), 0.0, 1e-11);
    ASSERT_NEAR(sqrt(l2_norm_squared(error)), 0.0, 1e-12);
    ASSERT_NEAR(infinity_norm(error), 0.0, 1e-12);
}