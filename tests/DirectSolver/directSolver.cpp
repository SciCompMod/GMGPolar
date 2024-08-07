#include <gtest/gtest.h>

#include "../../include/GMGPolar/gmgpolar.h"

#include "../../include/InputFunctions/domainGeometry.h"
#include "../../include/InputFunctions/systemParameters.h"
/* ------ */
/* Test 1 */
/* ------ */
#include "../include/InputFunctions/DomainGeometry/circularGeometry.h"
#include "../include/InputFunctions/BoundaryConditions/cartesianR2_Boundary_CircularGeometry.h"
#include "../include/InputFunctions/DensityProfileCoefficients/sonnendruckerGyroCoefficients.h"
#include "../include/InputFunctions/SourceTerms/cartesianR2_SonnendruckerGyro_CircularGeometry.h"
/* ------ */
/* Test 2 */
/* ------ */
#include "../include/InputFunctions/DomainGeometry/shafranovGeometry.h"
#include "../include/InputFunctions/BoundaryConditions/cartesianR6_Boundary_ShafranovGeometry.h"
#include "../include/InputFunctions/DensityProfileCoefficients/zoniGyroCoefficients.h"
#include "../include/InputFunctions/SourceTerms/cartesianR6_ZoniGyro_ShafranovGeometry.h"
/* ------ */
/* Test 3 */
/* ------ */
#include "../include/InputFunctions/DomainGeometry/czarnyGeometry.h"
#include "../include/InputFunctions/BoundaryConditions/polarR6_Boundary_CzarnyGeometry.h"
#include "../include/InputFunctions/DensityProfileCoefficients/zoniShiftedCoefficients.h"
#include "../include/InputFunctions/SourceTerms/polarR6_ZoniShifted_CzarnyGeometry.h"
/* ------ */
/* Test 4 */
/* ------ */
#include "../include/InputFunctions/DomainGeometry/culhamGeometry.h"
#include "../include/InputFunctions/BoundaryConditions/refined_Boundary_CulhamGeometry.h"
#include "../include/InputFunctions/DensityProfileCoefficients/zoniShiftedGyroCoefficients.h"
#include "../include/InputFunctions/SourceTerms/refined_ZoniShiftedGyro_CulhamGeometry.h"

#include <random>

// Function to generate sample data for vector x using random values with seed
Vector<double> generate_random_sample_data_Direct_Solver(const PolarGrid& grid, unsigned int seed) {
    Vector<double> x(grid.number_of_nodes());
    std::mt19937 gen(seed);  // Standard mersenne_twister_engine seeded with seed
    std::uniform_real_distribution<double> dist(-100.0, 100.0); 
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = dist(gen);
    }
    return x;
}

/* -------- */
/* Circular */
/* -------- */

TEST(DirectSolverTest_CircularGeometry, DirectSolverDirBC_Interior_CircularGeometry) {
    std::vector<double> radii = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    auto grid = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid);
    int extrapolation = 0;
    Level level(0, std::move(grid), std::move(levelCache), extrapolation);

    double Rmax = radii.back();

    CircularGeometry domain_geometry(Rmax);

    double alpha_jump = 0.66;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<SonnendruckerGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<CartesianR2_Boundary_CircularGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CircularGeometry>(Rmax);

    SystemParameters system_parameters(std::move(coefficients), std::move(boundary_conditions), std::move(source_term));

    bool DirBC_Interior = true;
    int maxOpenMPThreads = 1;
    int openMPTaskThreads = 1;
    
    DirectSolver solver_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    Residual residual_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);

    const Vector<double> rhs = generate_random_sample_data_Direct_Solver(level.grid(), 42);
    Vector<double> solution = rhs;
    solver_op.solveInPlace(solution);

    Vector<double> residuum(level.grid().number_of_nodes());
    residual_op.computeResidual_V1(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(residuum), 0.0, 1e-11);
    ASSERT_NEAR(sqrt(l2_norm_squared(residuum)), 0.0, 1e-12);
    ASSERT_NEAR(infinity_norm(residuum), 0.0, 1e-12);
}

TEST(DirectSolverTest_CircularGeometry, DirectSolverAcrossOrigin_CircularGeometry) {
    std::vector<double> radii = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    auto grid = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid);
    int extrapolation = 0;
    Level level(0, std::move(grid), std::move(levelCache), extrapolation);

    double Rmax = radii.back();

    CircularGeometry domain_geometry(Rmax);

    double alpha_jump = 0.66;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<SonnendruckerGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<CartesianR2_Boundary_CircularGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CircularGeometry>(Rmax);

    SystemParameters system_parameters(std::move(coefficients), std::move(boundary_conditions), std::move(source_term));

    bool DirBC_Interior = false;
    int maxOpenMPThreads = 1;
    int openMPTaskThreads = 1;
    
    DirectSolver solver_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    Residual residual_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);

    const Vector<double> rhs = generate_random_sample_data_Direct_Solver(level.grid(), 42);
    Vector<double> solution = rhs;
    solver_op.solveInPlace(solution);

    Vector<double> residuum(level.grid().number_of_nodes());
    residual_op.computeResidual_V1(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(residuum), 0.0, 1e-8);
    ASSERT_NEAR(sqrt(l2_norm_squared(residuum)), 0.0, 1e-8);
    ASSERT_NEAR(infinity_norm(residuum), 0.0, 1e-8);
}

/* --------- */
/* Shafranov */
/* --------- */

TEST(DirectSolverTest_ShafranovGeometry, DirectSolverDirBC_Interior_ShafranovGeometry) {
    std::vector<double> radii = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    auto grid = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid);
    int extrapolation = 0;
    Level level(0, std::move(grid), std::move(levelCache), extrapolation);

    double Rmax = radii.back();
    double kappa_eps=0.3;
    double delta_e=0.2;

    ShafranovGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.4837;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<CartesianR6_Boundary_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR6_ZoniGyro_ShafranovGeometry>(Rmax, kappa_eps, delta_e);

    SystemParameters system_parameters(std::move(coefficients), std::move(boundary_conditions), std::move(source_term));

    bool DirBC_Interior = true;
    int maxOpenMPThreads = 1;
    int openMPTaskThreads = 1;
    
    DirectSolver solver_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    Residual residual_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);

    const Vector<double> rhs = generate_random_sample_data_Direct_Solver(level.grid(), 42);
    Vector<double> solution = rhs;
    solver_op.solveInPlace(solution);

    Vector<double> residuum(level.grid().number_of_nodes());
    residual_op.computeResidual_V1(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(residuum), 0.0, 1e-11);
    ASSERT_NEAR(sqrt(l2_norm_squared(residuum)), 0.0, 1e-12);
    ASSERT_NEAR(infinity_norm(residuum), 0.0, 1e-12);
}

TEST(DirectSolverTest_ShafranovGeometry, DirectSolverAcrossOrigin_ShafranovGeometry) {
    std::vector<double> radii = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    auto grid = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid);
    int extrapolation = 0;
    Level level(0, std::move(grid), std::move(levelCache), extrapolation);

    double Rmax = radii.back();
    double kappa_eps=0.3;
    double delta_e=0.2;

    ShafranovGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.4837;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<CartesianR6_Boundary_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR6_ZoniGyro_ShafranovGeometry>(Rmax, kappa_eps, delta_e);

    SystemParameters system_parameters(std::move(coefficients), std::move(boundary_conditions), std::move(source_term));

    bool DirBC_Interior = false;
    int maxOpenMPThreads = 1;
    int openMPTaskThreads = 1;
    
    DirectSolver solver_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    Residual residual_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);

    const Vector<double> rhs = generate_random_sample_data_Direct_Solver(level.grid(), 42);
    Vector<double> solution = rhs;
    solver_op.solveInPlace(solution);

    Vector<double> residuum(level.grid().number_of_nodes());
    residual_op.computeResidual_V1(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(residuum), 0.0, 1e-8);
    ASSERT_NEAR(sqrt(l2_norm_squared(residuum)), 0.0, 1e-8);
    ASSERT_NEAR(infinity_norm(residuum), 0.0, 1e-8);
}


/* ------ */
/* Czarny */
/* ------ */

TEST(DirectSolverTest_CzarnyGeometry, DirectSolverDirBC_Interior_CzarnyGeometry) {
    std::vector<double> radii = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    auto grid = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid);
    int extrapolation = 0;
    Level level(0, std::move(grid), std::move(levelCache), extrapolation);

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

    const Vector<double> rhs = generate_random_sample_data_Direct_Solver(level.grid(), 42);
    Vector<double> solution = rhs;
    solver_op.solveInPlace(solution);

    Vector<double> residuum(level.grid().number_of_nodes());
    residual_op.computeResidual_V1(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(residuum), 0.0, 1e-11);
    ASSERT_NEAR(sqrt(l2_norm_squared(residuum)), 0.0, 1e-12);
    ASSERT_NEAR(infinity_norm(residuum), 0.0, 1e-12);
}

TEST(DirectSolverTest_CzarnyGeometry, DirectSolverAcrossOrigin_CzarnyGeometry) {
    std::vector<double> radii = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    auto grid = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid);
    int extrapolation = 0;
    Level level(0, std::move(grid), std::move(levelCache), extrapolation);

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

    const Vector<double> rhs = generate_random_sample_data_Direct_Solver(level.grid(), 42);
    Vector<double> solution = rhs;
    solver_op.solveInPlace(solution);

    Vector<double> residuum(level.grid().number_of_nodes());
    residual_op.computeResidual_V1(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(residuum), 0.0, 1e-7);
    ASSERT_NEAR(sqrt(l2_norm_squared(residuum)), 0.0, 1e-8);
    ASSERT_NEAR(infinity_norm(residuum), 0.0, 1e-8);
}

/* ------ */
/* Culham */
/* ------ */

TEST(DirectSolverTest_CulhamGeometry, DirectSolverDirBC_Interior_CulhamGeometry) {
    std::vector<double> radii = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    auto grid = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid);
    int extrapolation = 0;
    Level level(0, std::move(grid), std::move(levelCache), extrapolation);

    double Rmax = radii.back();

    CulhamGeometry domain_geometry(Rmax);

    double alpha_jump = 0.7081;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniShiftedGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<Refined_Boundary_CulhamGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<Refined_ZoniShiftedGyro_CulhamGeometry>(Rmax);

    SystemParameters system_parameters(std::move(coefficients), std::move(boundary_conditions), std::move(source_term));

    bool DirBC_Interior = true;
    int maxOpenMPThreads = 1;
    int openMPTaskThreads = 1;
    
    DirectSolver solver_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    Residual residual_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);

    const Vector<double> rhs = generate_random_sample_data_Direct_Solver(level.grid(), 42);
    Vector<double> solution = rhs;
    solver_op.solveInPlace(solution);

    Vector<double> residuum(level.grid().number_of_nodes());
    residual_op.computeResidual_V1(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(residuum), 0.0, 1e-11);
    ASSERT_NEAR(sqrt(l2_norm_squared(residuum)), 0.0, 1e-12);
    ASSERT_NEAR(infinity_norm(residuum), 0.0, 1e-12);
}

TEST(DirectSolverTest_CulhamGeometry, DirectSolverAcrossOrigin_CulhamGeometry) {
    std::vector<double> radii = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    auto grid = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid);
    int extrapolation = 0;
    Level level(0, std::move(grid), std::move(levelCache), extrapolation);

    double Rmax = radii.back();

    CulhamGeometry domain_geometry(Rmax);

    double alpha_jump = 0.7081;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniShiftedGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<Refined_Boundary_CulhamGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<Refined_ZoniShiftedGyro_CulhamGeometry>(Rmax);

    SystemParameters system_parameters(std::move(coefficients), std::move(boundary_conditions), std::move(source_term));

    bool DirBC_Interior = false;
    int maxOpenMPThreads = 1;
    int openMPTaskThreads = 1;
    
    DirectSolver solver_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    Residual residual_op(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);

    const Vector<double> rhs = generate_random_sample_data_Direct_Solver(level.grid(), 42);
    Vector<double> solution = rhs;
    solver_op.solveInPlace(solution);

    Vector<double> residuum(level.grid().number_of_nodes());
    residual_op.computeResidual_V1(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(residuum), 0.0, 1e-7);
    ASSERT_NEAR(sqrt(l2_norm_squared(residuum)), 0.0, 1e-7);
    ASSERT_NEAR(infinity_norm(residuum), 0.0, 1e-8);
}

int main(int argc, char* argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}