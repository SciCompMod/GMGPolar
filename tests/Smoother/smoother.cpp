#include <gtest/gtest.h>

#include <random>
#include <vector>

#include "../test_tools.h"

#include <GMGPolar/gmgpolar.h>

#include <Residual/ResidualGive/residualGive.h>
#include <DirectSolver/DirectSolverGive/directSolverGive.h>
#include <Smoother/SmootherGive/smootherGive.h>
#include <Smoother/SmootherTake/smootherTake.h>

#include <InputFunctions/domainGeometry.h>
#include <InputFunctions/densityProfileCoefficients.h>
#include <InputFunctions/boundaryConditions.h>
/* --------- */
/* Test Case */
/* --------- */
#include <InputFunctions/DomainGeometry/czarnyGeometry.h>
#include <InputFunctions/DensityProfileCoefficients/zoniShiftedCoefficients.h>

#include <random>
using namespace gmgpolar;

/* Test 1/2: */
/* Does the Take and Give Implementation match up? */

TEST(SmootherTest, smoother_DirBC_Interior)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto grid = std::make_unique<PolarGrid>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior = true;
    omp_set_num_threads(16);

    // "Take" requires cached values
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(grid), std::move(levelCache),
                                                                    ExtrapolationType::NONE, 0);

    SmootherGive smootherGive_operator(level.grid(), level.levelCache(), DirBC_Interior);
    SmootherTake smootherTake_operator(level.grid(), level.levelCache(), DirBC_Interior);

    HostVector<double> rhs   = generate_random_sample_data(level.grid(), 69);
    HostVector<double> start = generate_random_sample_data(level.grid(), 24);
    HostVector<double> temp  = generate_random_sample_data(level.grid(), 8);

    HostVector<double> solution_Give("solution_Give", start.size());
    Kokkos::deep_copy(solution_Give, start);
    smootherGive_operator.smoothing(solution_Give, rhs, temp);

    HostVector<double> solution_Take("solution_Take", start.size());
    Kokkos::deep_copy(solution_Take, start);
    smootherTake_operator.smoothing(solution_Take, rhs, temp);

    ASSERT_EQ(solution_Give.size(), solution_Take.size());
    for (uint index = 0; index < solution_Give.size(); index++) {
        int i_r, i_theta;
        level.grid().multiIndex(index, i_r, i_theta);
        if (i_r == 0 && !DirBC_Interior)
            ASSERT_NEAR(solution_Give[index], solution_Take[index], 1e-11);
        else
            ASSERT_NEAR(solution_Give[index], solution_Take[index], 1e-11);
    }
}

TEST(SmootherTest, smoother_AcrossOrigin)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto grid = std::make_unique<PolarGrid>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior = false;
    omp_set_num_threads(16);

    // "Take" requires cached values
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(grid), std::move(levelCache),
                                                                    ExtrapolationType::NONE, 0);

    SmootherGive smootherGive_operator(level.grid(), level.levelCache(), DirBC_Interior);
    SmootherTake smootherTake_operator(level.grid(), level.levelCache(), DirBC_Interior);

    HostVector<double> rhs   = generate_random_sample_data(level.grid(), 69);
    HostVector<double> start = generate_random_sample_data(level.grid(), 24);
    HostVector<double> temp  = generate_random_sample_data(level.grid(), 8);

    HostVector<double> solution_Give("solution_Give", start.size());
    Kokkos::deep_copy(solution_Give, start);
    smootherGive_operator.smoothing(solution_Give, rhs, temp);

    HostVector<double> solution_Take("solution_Take", start.size());
    Kokkos::deep_copy(solution_Take, start);
    smootherTake_operator.smoothing(solution_Take, rhs, temp);

    ASSERT_EQ(solution_Give.size(), solution_Take.size());
    for (uint index = 0; index < solution_Give.size(); index++) {
        int i_r, i_theta;
        level.grid().multiIndex(index, i_r, i_theta);
        if (i_r == 0 && !DirBC_Interior)
            ASSERT_NEAR(solution_Give[index], solution_Take[index], 1e-8);
        else
            ASSERT_NEAR(solution_Give[index], solution_Take[index], 1e-10);
    }
}

/* Test 2/2: */
/* Does the smoother converge to the directSolver solution? */

TEST(SmootherTest, SmootherDirBC_Interior)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto grid = std::make_unique<PolarGrid>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior                     = true;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(grid), std::move(levelCache),
                                                                    ExtrapolationType::NONE, 0);

    DirectSolverGive solver_op(level.grid(), level.levelCache(), DirBC_Interior);
    ResidualGive residual_op(level.grid(), level.levelCache(), DirBC_Interior);
    SmootherGive smoother_op(level.grid(), level.levelCache(), DirBC_Interior);

    HostConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    HostVector<double> discrete_solution("discrete_solution", rhs.size());
    Kokkos::deep_copy(discrete_solution, rhs);
    solver_op.solveInPlace(discrete_solution);

    HostVector<double> temp("temp", level.grid().numberOfNodes());
    HostVector<double> error("error", level.grid().numberOfNodes());
    HostVector<double> smoother_solution = generate_random_sample_data(level.grid(), 69);

    Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, error.size()),
                         [&](uint i) {
                             error(i) = discrete_solution(i) - smoother_solution(i);
                         });
    Kokkos::fence();

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-12;

    while (infinity_norm(HostConstVector<double>(error)) > precision) {
        smoother_op.smoothing(smoother_solution, rhs, temp);

        Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, error.size()),
                             [&](uint i) {
                                 error(i) = discrete_solution(i) - smoother_solution(i);
                             });
        Kokkos::fence();
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
    ASSERT_NEAR(infinity_norm(HostConstVector<double>(error)), 0.0, precision);
}

TEST(SmootherTest, SmootherAcrossOrigin)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto grid = std::make_unique<PolarGrid>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior                     = false;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(grid), std::move(levelCache),
                                                                    ExtrapolationType::NONE, 0);

    DirectSolverGive solver_op(level.grid(), level.levelCache(), DirBC_Interior);
    ResidualGive residual_op(level.grid(), level.levelCache(), DirBC_Interior);
    SmootherGive smoother_op(level.grid(), level.levelCache(), DirBC_Interior);

    HostConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    HostVector<double> discrete_solution("discrete_solution", rhs.size());
    Kokkos::deep_copy(discrete_solution, rhs);
    solver_op.solveInPlace(discrete_solution);

    HostVector<double> temp("temp", level.grid().numberOfNodes());
    HostVector<double> error("error", level.grid().numberOfNodes());
    HostVector<double> smoother_solution = generate_random_sample_data(level.grid(), 69);

    Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, error.size()),
                         [&](uint i) {
                             error(i) = discrete_solution(i) - smoother_solution(i);
                         });
    Kokkos::fence();

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-8;

    while (infinity_norm(HostConstVector<double>(error)) > precision) {
        smoother_op.smoothing(smoother_solution, rhs, temp);

        Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, error.size()),
                             [&](uint i) {
                                 error(i) = discrete_solution(i) - smoother_solution(i);
                             });
        Kokkos::fence();
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
    ASSERT_NEAR(infinity_norm(HostConstVector<double>(error)), 0.0, precision);
}

TEST(SmootherTest, SmootherDirBC_Interior_SmallestGrid)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI, M_PI + M_PI / 8, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto grid = std::make_unique<PolarGrid>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior                     = true;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(grid), std::move(levelCache),
                                                                    ExtrapolationType::NONE, 0);

    DirectSolverGive solver_op(level.grid(), level.levelCache(), DirBC_Interior);
    ResidualGive residual_op(level.grid(), level.levelCache(), DirBC_Interior);
    SmootherGive smoother_op(level.grid(), level.levelCache(), DirBC_Interior);

    HostConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    HostVector<double> discrete_solution("discrete_solution", rhs.size());
    Kokkos::deep_copy(discrete_solution, rhs);
    solver_op.solveInPlace(discrete_solution);

    HostVector<double> temp("temp", level.grid().numberOfNodes());
    HostVector<double> error("error", level.grid().numberOfNodes());
    HostVector<double> smoother_solution = generate_random_sample_data(level.grid(), 69);

    Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, error.size()),
                         [&](uint i) {
                             error(i) = discrete_solution(i) - smoother_solution(i);
                         });
    Kokkos::fence();

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    double precision            = 1e-12;

    while (infinity_norm(HostConstVector<double>(error)) > precision) {
        smoother_op.smoothing(smoother_solution, rhs, temp);

        Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, error.size()),
                             [&](uint i) {
                                 error(i) = discrete_solution(i) - smoother_solution(i);
                             });
        Kokkos::fence();
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
    ASSERT_NEAR(infinity_norm(HostConstVector<double>(error)), 0.0, precision);
}

TEST(SmootherTest, SmootherAcrossOrigin_SmallestGrid)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI, M_PI + M_PI / 8, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto grid = std::make_unique<PolarGrid>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior                     = false;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(grid), std::move(levelCache),
                                                                    ExtrapolationType::NONE, 0);

    DirectSolverGive solver_op(level.grid(), level.levelCache(), DirBC_Interior);
    ResidualGive residual_op(level.grid(), level.levelCache(), DirBC_Interior);
    SmootherGive smoother_op(level.grid(), level.levelCache(), DirBC_Interior);

    HostConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    HostVector<double> discrete_solution("discrete_solution", rhs.size());
    Kokkos::deep_copy(discrete_solution, rhs);
    solver_op.solveInPlace(discrete_solution);

    HostVector<double> temp("temp", level.grid().numberOfNodes());
    HostVector<double> error("error", level.grid().numberOfNodes());
    HostVector<double> smoother_solution = generate_random_sample_data(level.grid(), 69);

    Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, error.size()),
                         [&](uint i) {
                             error(i) = discrete_solution(i) - smoother_solution(i);
                         });
    Kokkos::fence();

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-8;

    while (infinity_norm(HostConstVector<double>(error)) > 1e-8) {
        smoother_op.smoothing(smoother_solution, rhs, temp);

        Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, error.size()),
                             [&](uint i) {
                                 error(i) = discrete_solution(i) - smoother_solution(i);
                             });
        Kokkos::fence();
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
    ASSERT_NEAR(infinity_norm(HostConstVector<double>(error)), 0.0, precision);
}

/* Using "Take" */

TEST(SmootherTest, SmootherTakeDirBC_Interior)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto grid = std::make_unique<PolarGrid>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior                     = true;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(grid), std::move(levelCache),
                                                                    ExtrapolationType::NONE, 0);

    DirectSolverGive solver_op(level.grid(), level.levelCache(), DirBC_Interior);
    ResidualGive residual_op(level.grid(), level.levelCache(), DirBC_Interior);
    SmootherTake smoother_op(level.grid(), level.levelCache(), DirBC_Interior);

    HostConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    HostVector<double> discrete_solution("discrete_solution", rhs.size());
    Kokkos::deep_copy(discrete_solution, rhs);
    solver_op.solveInPlace(discrete_solution);

    HostVector<double> temp("temp", level.grid().numberOfNodes());
    HostVector<double> error("error", level.grid().numberOfNodes());
    HostVector<double> smoother_solution = generate_random_sample_data(level.grid(), 69);

    Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, error.size()),
                         [&](uint i) {
                             error(i) = discrete_solution(i) - smoother_solution(i);
                         });
    Kokkos::fence();

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-12;

    while (infinity_norm(HostConstVector<double>(error)) > precision) {
        smoother_op.smoothing(smoother_solution, rhs, temp);

        Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, error.size()),
                             [&](uint i) {
                                 error(i) = discrete_solution(i) - smoother_solution(i);
                             });
        Kokkos::fence();
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
    ASSERT_NEAR(infinity_norm(HostConstVector<double>(error)), 0.0, precision);
}

TEST(SmootherTest, SmootherTakeAcrossOrigin)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto grid = std::make_unique<PolarGrid>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior                     = false;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(grid), std::move(levelCache),
                                                                    ExtrapolationType::NONE, 0);

    DirectSolverGive solver_op(level.grid(), level.levelCache(), DirBC_Interior);
    ResidualGive residual_op(level.grid(), level.levelCache(), DirBC_Interior);
    SmootherTake smoother_op(level.grid(), level.levelCache(), DirBC_Interior);

    HostConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    HostVector<double> discrete_solution("discrete_solution", rhs.size());
    Kokkos::deep_copy(discrete_solution, rhs);
    solver_op.solveInPlace(discrete_solution);

    HostVector<double> temp("temp", level.grid().numberOfNodes());
    HostVector<double> error("error", level.grid().numberOfNodes());
    HostVector<double> smoother_solution = generate_random_sample_data(level.grid(), 69);

    Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, error.size()),
                         [&](uint i) {
                             error(i) = discrete_solution(i) - smoother_solution(i);
                         });
    Kokkos::fence();

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-8;

    while (infinity_norm(HostConstVector<double>(error)) > precision) {
        smoother_op.smoothing(smoother_solution, rhs, temp);

        Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, error.size()),
                             [&](uint i) {
                                 error(i) = discrete_solution(i) - smoother_solution(i);
                             });
        Kokkos::fence();
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
    ASSERT_NEAR(infinity_norm(HostConstVector<double>(error)), 0.0, precision);
}

TEST(SmootherTest, SmootherTakeDirBC_Interior_SmallestGrid)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI, M_PI + M_PI / 8, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto grid = std::make_unique<PolarGrid>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior                     = true;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(grid), std::move(levelCache),
                                                                    ExtrapolationType::NONE, 0);

    DirectSolverGive solver_op(level.grid(), level.levelCache(), DirBC_Interior);
    ResidualGive residual_op(level.grid(), level.levelCache(), DirBC_Interior);
    SmootherTake smoother_op(level.grid(), level.levelCache(), DirBC_Interior);

    HostConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    HostVector<double> discrete_solution("discrete_solution", rhs.size());
    Kokkos::deep_copy(discrete_solution, rhs);
    solver_op.solveInPlace(discrete_solution);

    HostVector<double> temp("temp", level.grid().numberOfNodes());
    HostVector<double> error("error", level.grid().numberOfNodes());
    HostVector<double> smoother_solution = generate_random_sample_data(level.grid(), 69);

    Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, error.size()),
                         [&](uint i) {
                             error(i) = discrete_solution(i) - smoother_solution(i);
                         });
    Kokkos::fence();

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    double precision            = 1e-12;

    while (infinity_norm(HostConstVector<double>(error)) > precision) {
        smoother_op.smoothing(smoother_solution, rhs, temp);

        Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, error.size()),
                             [&](uint i) {
                                 error(i) = discrete_solution(i) - smoother_solution(i);
                             });
        Kokkos::fence();
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
    ASSERT_NEAR(infinity_norm(HostConstVector<double>(error)), 0.0, precision);
}

TEST(SmootherTest, SmootherTakeAcrossOrigin_SmallestGrid)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI, M_PI + M_PI / 8, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto grid = std::make_unique<PolarGrid>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior                     = false;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(grid), std::move(levelCache),
                                                                    ExtrapolationType::NONE, 0);

    DirectSolverGive solver_op(level.grid(), level.levelCache(), DirBC_Interior);
    ResidualGive residual_op(level.grid(), level.levelCache(), DirBC_Interior);
    SmootherTake smoother_op(level.grid(), level.levelCache(), DirBC_Interior);

    HostConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    HostVector<double> discrete_solution("discrete_solution", rhs.size());
    Kokkos::deep_copy(discrete_solution, rhs);
    solver_op.solveInPlace(discrete_solution);

    HostVector<double> temp("temp", level.grid().numberOfNodes());
    HostVector<double> error("error", level.grid().numberOfNodes());
    HostVector<double> smoother_solution = generate_random_sample_data(level.grid(), 69);

    Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, error.size()),
                         [&](uint i) {
                             error(i) = discrete_solution(i) - smoother_solution(i);
                         });
    Kokkos::fence();

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-8;

    while (infinity_norm(HostConstVector<double>(error)) > 1e-8) {
        smoother_op.smoothing(smoother_solution, rhs, temp);

        Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, error.size()),
                             [&](uint i) {
                                 error(i) = discrete_solution(i) - smoother_solution(i);
                             });
        Kokkos::fence();
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
    ASSERT_NEAR(infinity_norm(HostConstVector<double>(error)), 0.0, precision);
}
