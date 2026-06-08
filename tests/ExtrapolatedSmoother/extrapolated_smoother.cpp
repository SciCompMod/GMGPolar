#include <gtest/gtest.h>

#include <random>
#include <vector>

#include "../test_tools.h"

#include "../../include/GMGPolar/gmgpolar.h"

#include "../../include/Residual/ResidualGive/residualGive.h"
#include "../../include/Residual/ResidualTake/residualTake.h"
#include "../../include/DirectSolver/DirectSolverGive/directSolverGive.h"
#include "../../include/DirectSolver/DirectSolverTake/directSolverTake.h"
#include "../../include/ExtrapolatedSmoother/ExtrapolatedSmootherGive/extrapolatedSmootherGive.h"
#include "../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTake/extrapolatedSmootherTake.h"

#include "../../include/InputFunctions/domainGeometry.h"
#include "../../include/InputFunctions/densityProfileCoefficients.h"
/* --------- */
/* Test Case */
/* --------- */
#include "../include/InputFunctions/DomainGeometry/czarnyGeometry.h"
#include "../include/InputFunctions/DensityProfileCoefficients/zoniShiftedCoefficients.h"

#include <random>
using namespace gmgpolar;

/* Test 1/2: */
/* Does the Take and Give Implementation match up? */

TEST(ExtrapolatedSmootherTest, extrapolatedSmoother_DirBC_Interior)
{
    std::vector<double> radii  = {1e-2, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto tmp_grid = std::make_unique<PolarGrid<DefaultMemorySpace>>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior = true;

    // "Take" requires cached values
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *tmp_grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(tmp_grid), std::move(levelCache),
                                                                    ExtrapolationType::NONE, 0);

	const PolarGrid<DefaultMemorySpace>& grid = level.grid();
    ExtrapolatedSmootherGive smootherGive_operator(level.grid(), level.levelCache(), DirBC_Interior);
    ExtrapolatedSmootherTake smootherTake_operator(level.grid(), level.levelCache(), DirBC_Interior);

    Vector<double> rhs   = generate_random_sample_data(level.grid(), 69);
    HostVector<double> start = generate_random_sample_data(PolarGrid<Kokkos::HostSpace>(level.grid()), 24);
    Vector<double> temp  = generate_random_sample_data(level.grid(), 8);

    Vector<double> solution_Give("solution_Give", start.size());
    Kokkos::deep_copy(solution_Give, start);
    smootherGive_operator.extrapolatedSmoothing(solution_Give, rhs, temp);

    Vector<double> solution_Take("solution_Take", start.size());
    Kokkos::deep_copy(solution_Take, start);

    smootherTake_operator.extrapolatedSmoothing(solution_Take, rhs, temp);

	auto h_solution_Give = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, solution_Give);
	auto h_solution_Take = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, solution_Take);

    ASSERT_EQ(h_solution_Give.size(), h_solution_Take.size());
    for (uint index = 0; index < h_solution_Give.size(); index++) {
        int i_r, i_theta;
        level.grid().multiIndex(index, i_r, i_theta);
        if (i_r == 0 && !DirBC_Interior)
            ASSERT_NEAR(h_solution_Give[index], h_solution_Take[index], 1e-11);
        else
            ASSERT_NEAR(h_solution_Give[index], h_solution_Take[index], 1e-11);
    }
}

TEST(ExtrapolatedSmootherTest, extrapolatedSmoother_AcossOrigin)
{
    std::vector<double> radii  = {1e-2, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto tmp_grid = std::make_unique<PolarGrid<DefaultMemorySpace>>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior = false;

    // "Take" requires cached values
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *tmp_grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(tmp_grid), std::move(levelCache),
                                                                    ExtrapolationType::NONE, 0);

    ExtrapolatedSmootherGive smootherGive_operator(level.grid(), level.levelCache(), DirBC_Interior);
    ExtrapolatedSmootherTake smootherTake_operator(level.grid(), level.levelCache(), DirBC_Interior);

    Vector<double> rhs   = generate_random_sample_data(level.grid(), 69);
    HostVector<double> start = generate_random_sample_data(PolarGrid<Kokkos::HostSpace>(level.grid()), 24);
    Vector<double> temp  = generate_random_sample_data(level.grid(), 8);

    Vector<double> solution_Give("solution_Give", start.size());
    Kokkos::deep_copy(solution_Give, start);
    smootherGive_operator.extrapolatedSmoothing(solution_Give, rhs, temp);

    Vector<double> solution_Take("solution_Take", start.size());
    Kokkos::deep_copy(solution_Take, start);
    smootherTake_operator.extrapolatedSmoothing(solution_Take, rhs, temp);

	auto h_solution_Give = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, solution_Give);
	auto h_solution_Take = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, solution_Take);

    ASSERT_EQ(h_solution_Give.size(), h_solution_Take.size());
    for (uint index = 0; index < h_solution_Give.size(); index++) {
        int i_r, i_theta;
        level.grid().multiIndex(index, i_r, i_theta);
        if (i_r == 0 && !DirBC_Interior)
            ASSERT_NEAR(h_solution_Give[index], h_solution_Take[index], 1e-8);
        else
            ASSERT_NEAR(h_solution_Give[index], h_solution_Take[index], 1e-11);
    }
}

/* Test 2/2: */
/* Does the smoother converge to the DirectSolverGive solution? */

void ExtrapolatedSmootherTest_ExtrapolatedSmootherDirBC_Interior()
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto tmp_grid = std::make_unique<PolarGrid<DefaultMemorySpace>>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior = true;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *tmp_grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(tmp_grid), std::move(levelCache),
                                                                    ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

	const PolarGrid<DefaultMemorySpace>& grid = level.grid();
    DirectSolverGive solver_op(level.grid(), level.levelCache(), DirBC_Interior);
    ResidualGive residual_op(level.grid(), level.levelCache(), DirBC_Interior);
    ExtrapolatedSmootherGive extrapolated_smoother_op(level.grid(), level.levelCache(), DirBC_Interior);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    HostVector<double> h_discrete_solution("discrete_solution", rhs.size());
    Kokkos::deep_copy(h_discrete_solution, rhs);
    solver_op.solveInPlace(h_discrete_solution);
	auto discrete_solution = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_discrete_solution);

    Vector<double> temp("temp", level.grid().numberOfNodes());
    Vector<double> error("error", level.grid().numberOfNodes());
    Vector<double> smoother_solution = generate_random_sample_data(level.grid(), 69);
	        Kokkos::parallel_for(
            "Fill smoother_solution",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
                {0, 0}, // Starting point of the index space
                {(grid.nr()+1)/2, (grid.ntheta()+1)/2} // Ending point of the index space
                ),
            // Kokkos lambda function to execute for each point in the index space
            KOKKOS_LAMBDA(const int i_r_half, const int i_theta_half) {
                const int i_r = 2*i_r_half;
				const int i_theta = 2 * i_theta_half;
            smoother_solution[grid.index(i_r, i_theta)] = discrete_solution[grid.index(i_r, i_theta)];
        });

    Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, error.size()),
                         [&](uint i) {
                             error[i] = discrete_solution[i] - smoother_solution[i];
                         });
    Kokkos::fence();

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-12;

    while (infinity_norm(ConstVector<double>(error)) > precision) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

        Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, error.size()),
                             [&](uint i) {
                                 error[i] = discrete_solution[i] - smoother_solution[i];
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
    ASSERT_NEAR(infinity_norm(ConstVector<double>(error)), 0.0, precision);
}
TEST(ExtrapolatedSmootherTest, ExtrapolatedSmootherDirBC_Interior)
{ ExtrapolatedSmootherTest_ExtrapolatedSmootherDirBC_Interior();}

void ExtrapolatedSmootherTest_ExtrapolatedSmootherAcrossOrigin()
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto tmp_grid = std::make_unique<PolarGrid<DefaultMemorySpace>>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior = false;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *tmp_grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(tmp_grid), std::move(levelCache),
                                                                    ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

	const PolarGrid<DefaultMemorySpace>& grid = level.grid();
    DirectSolverGive solver_op(level.grid(), level.levelCache(), DirBC_Interior);
    ResidualGive residual_op(level.grid(), level.levelCache(), DirBC_Interior);
    ExtrapolatedSmootherGive extrapolated_smoother_op(level.grid(), level.levelCache(), DirBC_Interior);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    HostVector<double> h_discrete_solution("discrete_solution", rhs.size());
    Kokkos::deep_copy(h_discrete_solution, rhs);
    solver_op.solveInPlace(h_discrete_solution);
	auto discrete_solution = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_discrete_solution);

    Vector<double> temp("temp", level.grid().numberOfNodes());

    Vector<double> error("error", level.grid().numberOfNodes());

    Vector<double> smoother_solution = generate_random_sample_data(level.grid(), 69);
	        Kokkos::parallel_for(
            "Fill smoother_solution",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
                {0, 0}, // Starting point of the index space
                {(grid.nr()+1)/2, (grid.ntheta()+1)/2} // Ending point of the index space
                ),
            // Kokkos lambda function to execute for each point in the index space
            KOKKOS_LAMBDA(const int i_r_half, const int i_theta_half) {
                const int i_r = 2*i_r_half;
				const int i_theta = 2 * i_theta_half;
            smoother_solution[grid.index(i_r, i_theta)] = discrete_solution[grid.index(i_r, i_theta)];
    });

    Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, error.size()),
                         [&](uint i) {
                             error[i] = discrete_solution[i] - smoother_solution[i];
                         });
    Kokkos::fence();

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-8;

    while (infinity_norm(ConstVector<double>(error)) > precision) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

        Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, error.size()),
                             [&](uint i) {
                                 error[i] = discrete_solution[i] - smoother_solution[i];
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
    ASSERT_NEAR(infinity_norm(ConstVector<double>(error)), 0.0, precision);
}
TEST(ExtrapolatedSmootherTest, ExtrapolatedSmootherAcrossOrigin)
{ ExtrapolatedSmootherTest_ExtrapolatedSmootherAcrossOrigin();}

void ExtrapolatedSmootherTest_ExtrapolatedSmootherDirBC_Interior_SmallestGrid()
{
    std::vector<double> radii  = {1e-5, 0.2, 0.4, 0.45, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI, M_PI + M_PI / 8, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto tmp_grid = std::make_unique<PolarGrid<DefaultMemorySpace>>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior = true;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *tmp_grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(tmp_grid), std::move(levelCache),
                                                                    ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

	const PolarGrid<DefaultMemorySpace>& grid = level.grid();
    DirectSolverGive solver_op(level.grid(), level.levelCache(), DirBC_Interior);
    ResidualGive residual_op(level.grid(), level.levelCache(), DirBC_Interior);
    ExtrapolatedSmootherGive extrapolated_smoother_op(level.grid(), level.levelCache(), DirBC_Interior);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    HostVector<double> h_discrete_solution("discrete_solution", rhs.size());
    Kokkos::deep_copy(h_discrete_solution, rhs);
    solver_op.solveInPlace(h_discrete_solution);
	auto discrete_solution = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_discrete_solution);

    Vector<double> temp("temp", level.grid().numberOfNodes());

    Vector<double> error("error", level.grid().numberOfNodes());

    Vector<double> smoother_solution = generate_random_sample_data(level.grid(), 69);
	        Kokkos::parallel_for(
            "Fill smoother_solution",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
                {0, 0}, // Starting point of the index space
                {(grid.nr()+1)/2, (grid.ntheta()+1)/2} // Ending point of the index space
                ),
            // Kokkos lambda function to execute for each point in the index space
            KOKKOS_LAMBDA(const int i_r_half, const int i_theta_half) {
                const int i_r = 2*i_r_half;
				const int i_theta = 2 * i_theta_half;
            smoother_solution[grid.index(i_r, i_theta)] = discrete_solution[grid.index(i_r, i_theta)];
        });

    Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, error.size()),
                         [&](uint i) {
                             error[i] = discrete_solution[i] - smoother_solution[i];
                         });
    Kokkos::fence();

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    double precision            = 1e-12;

    while (infinity_norm(ConstVector<double>(error)) > precision) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

        Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, error.size()),
                             [&](uint i) {
                                 error[i] = discrete_solution[i] - smoother_solution[i];
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
    ASSERT_LT(iterations, 200);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(error)), 0.0, precision);
}
TEST(ExtrapolatedSmootherTest, ExtrapolatedSmootherDirBC_Interior_SmallestGrid)
{ ExtrapolatedSmootherTest_ExtrapolatedSmootherDirBC_Interior_SmallestGrid();}

void ExtrapolatedSmootherTest_ExtrapolatedSmootherAcrossOrigin_SmallestGrid()
{
    std::vector<double> radii  = {1e-5, 0.2, 0.4, 0.45, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI, M_PI + M_PI / 8, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto tmp_grid = std::make_unique<PolarGrid<DefaultMemorySpace>>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior = false;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *tmp_grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(tmp_grid), std::move(levelCache),
                                                                    ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

	const PolarGrid<DefaultMemorySpace>& grid = level.grid();
    DirectSolverGive solver_op(level.grid(), level.levelCache(), DirBC_Interior);
    ResidualGive residual_op(level.grid(), level.levelCache(), DirBC_Interior);
    ExtrapolatedSmootherGive extrapolated_smoother_op(level.grid(), level.levelCache(), DirBC_Interior);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    HostVector<double> h_discrete_solution("discrete_solution", rhs.size());
    Kokkos::deep_copy(h_discrete_solution, rhs);
    solver_op.solveInPlace(h_discrete_solution);
	auto discrete_solution = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_discrete_solution);

    Vector<double> temp("temp", level.grid().numberOfNodes());
    Vector<double> error("error", level.grid().numberOfNodes());
    Vector<double> smoother_solution = generate_random_sample_data(level.grid(), 69);
	        Kokkos::parallel_for(
            "Fill smoother_solution",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
                {0, 0}, // Starting point of the index space
                {(grid.nr()+1)/2, (grid.ntheta()+1)/2} // Ending point of the index space
                ),
            // Kokkos lambda function to execute for each point in the index space
            KOKKOS_LAMBDA(const int i_r_half, const int i_theta_half) {
                const int i_r = 2*i_r_half;
				const int i_theta = 2 * i_theta_half;
            smoother_solution[grid.index(i_r, i_theta)] = discrete_solution[grid.index(i_r, i_theta)];
    });

    Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, error.size()),
                         [&](uint i) {
                             error[i] = discrete_solution[i] - smoother_solution[i];
                         });
    Kokkos::fence();

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-8;

    while (infinity_norm(ConstVector<double>(error)) > 1e-8) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

        Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, error.size()),
                             [&](uint i) {
                                 error[i] = discrete_solution[i] - smoother_solution[i];
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
    ASSERT_LT(iterations, 150);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(error)), 0.0, precision);
}
TEST(ExtrapolatedSmootherTest, ExtrapolatedSmootherAcrossOrigin_SmallestGrid)
{ ExtrapolatedSmootherTest_ExtrapolatedSmootherAcrossOrigin_SmallestGrid();}

/* We now test using "Take" */

void ExtrapolatedSmootherTest_ExtrapolatedSmootherTakeDirBC_Interior()
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto tmp_grid = std::make_unique<PolarGrid<DefaultMemorySpace>>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior = true;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *tmp_grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(tmp_grid), std::move(levelCache),
                                                                    ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

	const PolarGrid<DefaultMemorySpace>& grid = level.grid();
    DirectSolverTake solver_op(level.grid(), level.levelCache(), DirBC_Interior);
    ResidualTake residual_op(level.grid(), level.levelCache(), DirBC_Interior);
    ExtrapolatedSmootherTake extrapolated_smoother_op(level.grid(), level.levelCache(), DirBC_Interior);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    HostVector<double> h_discrete_solution("discrete_solution", rhs.size());
    Kokkos::deep_copy(h_discrete_solution, rhs);
    solver_op.solveInPlace(h_discrete_solution);
	auto discrete_solution = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_discrete_solution);

    Vector<double> temp("temp", level.grid().numberOfNodes());
    Vector<double> error("error", level.grid().numberOfNodes());
    Vector<double> smoother_solution = generate_random_sample_data(level.grid(), 69);
	        Kokkos::parallel_for(
            "Fill smoother_solution",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
                {0, 0}, // Starting point of the index space
                {(grid.nr()+1)/2, (grid.ntheta()+1)/2} // Ending point of the index space
                ),
            // Kokkos lambda function to execute for each point in the index space
            KOKKOS_LAMBDA(const int i_r_half, const int i_theta_half) {
                const int i_r = 2*i_r_half;
				const int i_theta = 2 * i_theta_half;
            smoother_solution[grid.index(i_r, i_theta)] = discrete_solution[grid.index(i_r, i_theta)];
    });

    Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, error.size()),
                         [&](uint i) {
                             error[i] = discrete_solution[i] - smoother_solution[i];
                         });
    Kokkos::fence();

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-12;

    while (infinity_norm(ConstVector<double>(error)) > precision) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

        Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, error.size()),
                             [&](uint i) {
                                 error[i] = discrete_solution[i] - smoother_solution[i];
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
    ASSERT_NEAR(infinity_norm(ConstVector<double>(error)), 0.0, precision);
}
TEST(ExtrapolatedSmootherTest, ExtrapolatedSmootherTakeDirBC_Interior)
{ ExtrapolatedSmootherTest_ExtrapolatedSmootherTakeDirBC_Interior();}

void ExtrapolatedSmootherTest_ExtrapolatedSmootherTakeAcrossOrigin()
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto tmp_grid = std::make_unique<PolarGrid<DefaultMemorySpace>>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior = false;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *tmp_grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(tmp_grid), std::move(levelCache),
                                                                    ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

	const PolarGrid<DefaultMemorySpace>& grid = level.grid();
    DirectSolverTake solver_op(level.grid(), level.levelCache(), DirBC_Interior);
    ResidualTake residual_op(level.grid(), level.levelCache(), DirBC_Interior);
    ExtrapolatedSmootherTake extrapolated_smoother_op(level.grid(), level.levelCache(), DirBC_Interior);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    HostVector<double> h_discrete_solution("discrete_solution", rhs.size());
    Kokkos::deep_copy(h_discrete_solution, rhs);
    solver_op.solveInPlace(h_discrete_solution);
	auto discrete_solution = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_discrete_solution);

    Vector<double> temp("temp", level.grid().numberOfNodes());
    Vector<double> error("error", level.grid().numberOfNodes());
    Vector<double> smoother_solution = generate_random_sample_data(level.grid(), 69);
	        Kokkos::parallel_for(
            "Fill smoother_solution",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
                {0, 0}, // Starting point of the index space
                {(grid.nr()+1)/2, (grid.ntheta()+1)/2} // Ending point of the index space
                ),
            // Kokkos lambda function to execute for each point in the index space
            KOKKOS_LAMBDA(const int i_r_half, const int i_theta_half) {
                const int i_r = 2*i_r_half;
				const int i_theta = 2 * i_theta_half;
            smoother_solution[grid.index(i_r, i_theta)] = discrete_solution[grid.index(i_r, i_theta)];
    });

    Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, error.size()),
                         [&](uint i) {
                             error[i] = discrete_solution[i] - smoother_solution[i];
                         });
    Kokkos::fence();

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-8;

    while (infinity_norm(ConstVector<double>(error)) > precision) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

        Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, error.size()),
                             [&](uint i) {
                                 error[i] = discrete_solution[i] - smoother_solution[i];
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
    ASSERT_NEAR(infinity_norm(ConstVector<double>(error)), 0.0, precision);
}
TEST(ExtrapolatedSmootherTest, ExtrapolatedSmootherTakeAcrossOrigin)
{ ExtrapolatedSmootherTest_ExtrapolatedSmootherTakeAcrossOrigin();}

void ExtrapolatedSmootherTest_ExtrapolatedSmootherTakeDirBC_Interior_SmallestGrid()
{
    std::vector<double> radii  = {1e-5, 0.2, 0.4, 0.45, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI, M_PI + M_PI / 8, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto tmp_grid = std::make_unique<PolarGrid<DefaultMemorySpace>>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior = true;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *tmp_grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(tmp_grid), std::move(levelCache),
                                                                    ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

	const PolarGrid<DefaultMemorySpace>& grid = level.grid();
    DirectSolverTake solver_op(level.grid(), level.levelCache(), DirBC_Interior);
    ResidualTake residual_op(level.grid(), level.levelCache(), DirBC_Interior);
    ExtrapolatedSmootherTake extrapolated_smoother_op(level.grid(), level.levelCache(), DirBC_Interior);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    HostVector<double> h_discrete_solution("discrete_solution", rhs.size());
    Kokkos::deep_copy(h_discrete_solution, rhs);
    solver_op.solveInPlace(h_discrete_solution);
	auto discrete_solution = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_discrete_solution);

    Vector<double> temp("temp", level.grid().numberOfNodes());
    Vector<double> error("error", level.grid().numberOfNodes());
    Vector<double> smoother_solution = generate_random_sample_data(level.grid(), 69);
	        Kokkos::parallel_for(
            "Fill smoother_solution",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
                {0, 0}, // Starting point of the index space
                {(grid.nr()+1)/2, (grid.ntheta()+1)/2} // Ending point of the index space
                ),
            // Kokkos lambda function to execute for each point in the index space
            KOKKOS_LAMBDA(const int i_r_half, const int i_theta_half) {
                const int i_r = 2*i_r_half;
				const int i_theta = 2 * i_theta_half;
            smoother_solution[grid.index(i_r, i_theta)] = discrete_solution[grid.index(i_r, i_theta)];
    });

    Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, error.size()),
                         [&](uint i) {
                             error[i] = discrete_solution[i] - smoother_solution[i];
                         });
    Kokkos::fence();

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    double precision            = 1e-12;

    while (infinity_norm(ConstVector<double>(error)) > precision) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

        Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, error.size()),
                             [&](uint i) {
                                 error[i] = discrete_solution[i] - smoother_solution[i];
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
    ASSERT_LT(iterations, 200);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(error)), 0.0, precision);
}
TEST(ExtrapolatedSmootherTest, ExtrapolatedSmootherTakeDirBC_Interior_SmallestGrid)
{ ExtrapolatedSmootherTest_ExtrapolatedSmootherTakeDirBC_Interior_SmallestGrid(); }

void ExtrapolatedSmootherTest_ExtrapolatedSmootherTakeAcrossOrigin_SmallestGrid()
{
    std::vector<double> radii  = {1e-5, 0.2, 0.4, 0.45, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI, M_PI + M_PI / 8, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    auto tmp_grid = std::make_unique<PolarGrid<DefaultMemorySpace>>(radii, angles);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior = false;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *tmp_grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(tmp_grid), std::move(levelCache),
                                                                    ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

	const PolarGrid<DefaultMemorySpace>& grid = level.grid();
    DirectSolverTake solver_op(level.grid(), level.levelCache(), DirBC_Interior);
    ResidualTake residual_op(level.grid(), level.levelCache(), DirBC_Interior);
    ExtrapolatedSmootherTake extrapolated_smoother_op(level.grid(), level.levelCache(), DirBC_Interior);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    HostVector<double> h_discrete_solution("discrete_solution", rhs.size());
    Kokkos::deep_copy(h_discrete_solution, rhs);
    solver_op.solveInPlace(h_discrete_solution);
	auto discrete_solution = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_discrete_solution);

    Vector<double> temp("temp", level.grid().numberOfNodes());
    Vector<double> error("error", level.grid().numberOfNodes());
    Vector<double> smoother_solution = generate_random_sample_data(level.grid(), 69);
	        Kokkos::parallel_for(
            "Fill smoother_solution",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
                {0, 0}, // Starting point of the index space
                {(grid.nr()+1)/2, (grid.ntheta()+1)/2} // Ending point of the index space
                ),
            // Kokkos lambda function to execute for each point in the index space
            KOKKOS_LAMBDA(const int i_r_half, const int i_theta_half) {
                const int i_r = 2*i_r_half;
				const int i_theta = 2 * i_theta_half;
            smoother_solution[grid.index(i_r, i_theta)] = discrete_solution[grid.index(i_r, i_theta)];
    });

    Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, error.size()),
                         [&](uint i) {
                             error[i] = discrete_solution[i] - smoother_solution[i];
                         });
    Kokkos::fence();

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-8;

    while (infinity_norm(ConstVector<double>(error)) > 1e-8) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

        Kokkos::parallel_for("get error", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, error.size()),
                             [&](uint i) {
                                 error[i] = discrete_solution[i] - smoother_solution[i];
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
    ASSERT_LT(iterations, 150);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(error)), 0.0, precision);
}
TEST(ExtrapolatedSmootherTest, ExtrapolatedSmootherTakeAcrossOrigin_SmallestGrid)
{ ExtrapolatedSmootherTest_ExtrapolatedSmootherTakeAcrossOrigin_SmallestGrid();
}
