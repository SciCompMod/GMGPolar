#include <gtest/gtest.h>

#include <random>
#include <vector>

#include "../../include/GMGPolar/gmgpolar.h"

#include "../../include/InputFunctions/domainGeometry.h"
#include "../../include/InputFunctions/densityProfileCoefficients.h"
#include "../../include/InputFunctions/boundaryConditions.h"
#include "../../include/InputFunctions/sourceTerm.h"
/* --------- */
/* Test Case */
/* --------- */
#include "../include/InputFunctions/DomainGeometry/czarnyGeometry.h"
#include "../include/InputFunctions/BoundaryConditions/polarR6_Boundary_CzarnyGeometry.h"
#include "../include/InputFunctions/DensityProfileCoefficients/zoniShiftedCoefficients.h"
#include "../include/InputFunctions/SourceTerms/polarR6_ZoniShifted_CzarnyGeometry.h"

#include <random>

namespace SmootherTest {
    Vector<double> generate_random_sample_data(const PolarGrid& grid, unsigned int seed) {
        Vector<double> x(grid.numberOfNodes());
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> dist(-100.0, 100.0); 
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] = dist(gen);
        }
        return x;
    }
}

using namespace SmootherTest;


TEST(SmootherTest, SequentialSmootherDirBC_Interior) {
    std::vector<double> radii = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    double Rmax = radii.back();
    double kappa_eps=0.3;
    double delta_e=1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior = true;
    int maxOpenMPThreads = 1;

    auto grid = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients);
    Level level(0, std::move(grid), std::move(levelCache), 0);

    DirectSolver solver_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);
    Residual residual_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);
    Smoother smoother_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);


    const Vector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());
    
    smoother_solution = generate_random_sample_data(level.grid(), 69);

    #pragma omp parallel for
    for (std::size_t i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations = 0;
    bool max_iterations_reached = false;
    const int max_iterations = 10000;
    const double precision = 1e-12;

    while (infinity_norm(error) > precision) {
        smoother_op.smoothingInPlace(smoother_solution, rhs, temp);

        #pragma omp parallel for
        for (std::size_t i = 0; i < error.size(); i++) {
            error[i] = discrete_solution[i] - smoother_solution[i];
        }
        iterations++;
        if (iterations >= max_iterations) {
            max_iterations_reached = true;
            std::cerr << "Max iterations reached without convergence." << std::endl;
            break;
        }
    }

    std::cout << "Convergence reached after " << iterations << " iterations." << std::endl;

    ASSERT_TRUE(!max_iterations_reached);
    ASSERT_LT(iterations, 300);
    ASSERT_NEAR(infinity_norm(error), 0.0, precision);
}


TEST(SmootherTest, ParallelSmootherDirBC_Interior) {
    std::vector<double> radii = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    double Rmax = radii.back();
    double kappa_eps=0.3;
    double delta_e=1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior = true;
    int maxOpenMPThreads = 16;

    auto grid = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients);
    Level level(0, std::move(grid), std::move(levelCache), 0);

    DirectSolver solver_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);
    Residual residual_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);
    Smoother smoother_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);


    const Vector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());
    
    smoother_solution = generate_random_sample_data(level.grid(), 69);

    #pragma omp parallel for
    for (std::size_t i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations = 0;
    bool max_iterations_reached = false;
    const int max_iterations = 10000;
    const double precision = 1e-12;

    while (infinity_norm(error) > precision) {
        smoother_op.smoothingInPlace(smoother_solution, rhs, temp);

        #pragma omp parallel for
        for (std::size_t i = 0; i < error.size(); i++) {
            error[i] = discrete_solution[i] - smoother_solution[i];
        }
        iterations++;
        if (iterations >= max_iterations) {
            max_iterations_reached = true;
            std::cerr << "Max iterations reached without convergence." << std::endl;
            break;
        }
    }

    std::cout << "Convergence reached after " << iterations << " iterations." << std::endl;

    ASSERT_TRUE(!max_iterations_reached);
    ASSERT_LT(iterations, 300);
    ASSERT_NEAR(infinity_norm(error), 0.0, precision);
}


TEST(SmootherTest, SequentialSmootherAcrossOrigin) {
    std::vector<double> radii = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    double Rmax = radii.back();
    double kappa_eps=0.3;
    double delta_e=1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior = false;
    int maxOpenMPThreads = 1;

    auto grid = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients);
    Level level(0, std::move(grid), std::move(levelCache), 0);

    DirectSolver solver_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);
    Residual residual_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);
    Smoother smoother_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);


    const Vector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());
    
    smoother_solution = generate_random_sample_data(level.grid(), 69);

    #pragma omp parallel for
    for (std::size_t i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations = 0;
    bool max_iterations_reached = false;
    const int max_iterations = 10000;
    const double precision = 1e-8;

    while (infinity_norm(error) > precision) {
        smoother_op.smoothingInPlace(smoother_solution, rhs, temp);

        #pragma omp parallel for
        for (std::size_t i = 0; i < error.size(); i++) {
            error[i] = discrete_solution[i] - smoother_solution[i];
        }
        iterations++;
        if (iterations >= max_iterations) {
            max_iterations_reached = true;
            std::cerr << "Max iterations reached without convergence." << std::endl;
            break;
        }
    }

    std::cout << "Convergence reached after " << iterations << " iterations." << std::endl;

    ASSERT_TRUE(!max_iterations_reached);
    ASSERT_LT(iterations, 600);
    ASSERT_NEAR(infinity_norm(error), 0.0, precision);
}


TEST(SmootherTest, ParallelSmootherAcrossOrigin) {
    std::vector<double> radii = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    double Rmax = radii.back();
    double kappa_eps=0.3;
    double delta_e=1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior = false;
    int maxOpenMPThreads = 16;

    auto grid = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients);
    Level level(0, std::move(grid), std::move(levelCache), 0);

    DirectSolver solver_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);
    Residual residual_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);
    Smoother smoother_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);

    const Vector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());
    
    smoother_solution = generate_random_sample_data(level.grid(), 69);

    #pragma omp parallel for
    for (std::size_t i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations = 0;
    bool max_iterations_reached = false;
    const int max_iterations = 10000;
    const double precision = 1e-8;

    while (infinity_norm(error) > 1e-8) {
        smoother_op.smoothingInPlace(smoother_solution, rhs, temp);

        #pragma omp parallel for
        for (std::size_t i = 0; i < error.size(); i++) {
            error[i] = discrete_solution[i] - smoother_solution[i];
        }
        iterations++;
        if (iterations >= max_iterations) {
            max_iterations_reached = true;
            std::cerr << "Max iterations reached without convergence." << std::endl;
            break;
        }
    }

    std::cout << "Convergence reached after " << iterations << " iterations." << std::endl;

    ASSERT_TRUE(!max_iterations_reached);
    ASSERT_LT(iterations, 600);
    ASSERT_NEAR(infinity_norm(error), 0.0, precision);
}


TEST(SmootherTest, SequentialSmootherDirBC_Interior_SmallestGrid) {
    std::vector<double> radii = {1e-5, 0.2, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/8, M_PI, M_PI+M_PI/8, M_PI+M_PI};

    double Rmax = radii.back();
    double kappa_eps=0.3;
    double delta_e=1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior = true;
    int maxOpenMPThreads = 1;

    auto grid = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients);
    Level level(0, std::move(grid), std::move(levelCache), 0);

    DirectSolver solver_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);
    Residual residual_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);
    Smoother smoother_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);


    const Vector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());
    
    smoother_solution = generate_random_sample_data(level.grid(), 69);

    #pragma omp parallel for
    for (std::size_t i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations = 0;
    bool max_iterations_reached = false;
    const int max_iterations = 10000;
    const double precision = 1e-12;

    while (infinity_norm(error) > 1e-12) {
        smoother_op.smoothingInPlace(smoother_solution, rhs, temp);

        #pragma omp parallel for
        for (std::size_t i = 0; i < error.size(); i++) {
            error[i] = discrete_solution[i] - smoother_solution[i];
        }
        iterations++;
        if (iterations >= max_iterations) {
            max_iterations_reached = true;
            std::cerr << "Max iterations reached without convergence." << std::endl;
            break;
        }
    }

    std::cout << "Convergence reached after " << iterations << " iterations." << std::endl;

    ASSERT_TRUE(!max_iterations_reached);
    ASSERT_LT(iterations, 30);
    ASSERT_NEAR(infinity_norm(error), 0.0, precision);
}


TEST(SmootherTest, ParallelSmootherDirBC_Interior_SmallestGrid) {
    std::vector<double> radii = {1e-5, 0.2, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/8, M_PI, M_PI+M_PI/8, M_PI+M_PI};

    double Rmax = radii.back();
    double kappa_eps=0.3;
    double delta_e=1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior = true;
    int maxOpenMPThreads = 16;

    auto grid = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients);
    Level level(0, std::move(grid), std::move(levelCache), 0);

    DirectSolver solver_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);
    Residual residual_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);
    Smoother smoother_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);


    const Vector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());
    
    smoother_solution = generate_random_sample_data(level.grid(), 69);

    #pragma omp parallel for
    for (std::size_t i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations = 0;
    bool max_iterations_reached = false;
    const int max_iterations = 10000;
    double precision = 1e-12;

    while (infinity_norm(error) > precision) {
        smoother_op.smoothingInPlace(smoother_solution, rhs, temp);

        #pragma omp parallel for
        for (std::size_t i = 0; i < error.size(); i++) {
            error[i] = discrete_solution[i] - smoother_solution[i];
        }
        iterations++;
        if (iterations >= max_iterations) {
            max_iterations_reached = true;
            std::cerr << "Max iterations reached without convergence." << std::endl;
            break;
        }
    }

    std::cout << "Convergence reached after " << iterations << " iterations." << std::endl;

    ASSERT_TRUE(!max_iterations_reached);
    ASSERT_LT(iterations, 30);
    ASSERT_NEAR(infinity_norm(error), 0.0, precision);
}


TEST(SmootherTest, SequentialSmootherAcrossOrigin_SmallestGrid) {
    std::vector<double> radii = {1e-5, 0.2, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/8, M_PI, M_PI+M_PI/8, M_PI+M_PI};

    double Rmax = radii.back();
    double kappa_eps=0.3;
    double delta_e=1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior = false;
    int maxOpenMPThreads = 1;

    auto grid = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients);
    Level level(0, std::move(grid), std::move(levelCache), 0);

    DirectSolver solver_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);
    Residual residual_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);
    Smoother smoother_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);


    const Vector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());
    
    smoother_solution = generate_random_sample_data(level.grid(), 69);

    #pragma omp parallel for
    for (std::size_t i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations = 0;
    bool max_iterations_reached = false;
    const int max_iterations = 10000;
    const double precision = 1e-8;

    while (infinity_norm(error) > 1e-8) {
        smoother_op.smoothingInPlace(smoother_solution, rhs, temp);

        #pragma omp parallel for
        for (std::size_t i = 0; i < error.size(); i++) {
            error[i] = discrete_solution[i] - smoother_solution[i];
        }
        iterations++;
        if (iterations >= max_iterations) {
            max_iterations_reached = true;
            std::cerr << "Max iterations reached without convergence." << std::endl;
            break;
        }
    }

    std::cout << "Convergence reached after " << iterations << " iterations." << std::endl;

    ASSERT_TRUE(!max_iterations_reached);
    ASSERT_LT(iterations, 80);
    ASSERT_NEAR(infinity_norm(error), 0.0, precision);
}


TEST(SmootherTest, ParallelSmootherAcrossOrigin_SmallestGrid) {
    std::vector<double> radii = {1e-5, 0.2, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/8, M_PI, M_PI+M_PI/8, M_PI+M_PI};

    double Rmax = radii.back();
    double kappa_eps=0.3;
    double delta_e=1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior = false;
    int maxOpenMPThreads = 16;

    auto grid = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients);
    Level level(0, std::move(grid), std::move(levelCache), 0);

    DirectSolver solver_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);
    Residual residual_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);
    Smoother smoother_op(level.grid(), level.levelCache(), domain_geometry, DirBC_Interior, maxOpenMPThreads);

    const Vector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());
    
    smoother_solution = generate_random_sample_data(level.grid(), 69);

    #pragma omp parallel for
    for (std::size_t i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations = 0;
    bool max_iterations_reached = false;
    const int max_iterations = 10000;
    const double precision = 1e-8;

    while (infinity_norm(error) > 1e-8) {
        smoother_op.smoothingInPlace(smoother_solution, rhs, temp);

        #pragma omp parallel for
        for (std::size_t i = 0; i < error.size(); i++) {
            error[i] = discrete_solution[i] - smoother_solution[i];
        }
        iterations++;
        if (iterations >= max_iterations) {
            max_iterations_reached = true;
            std::cerr << "Max iterations reached without convergence." << std::endl;
            break;
        }
    }

    std::cout << "Convergence reached after " << iterations << " iterations." << std::endl;

    ASSERT_TRUE(!max_iterations_reached);
    ASSERT_LT(iterations, 80);
    ASSERT_NEAR(infinity_norm(error), 0.0, precision);
}


int main(int argc, char* argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}