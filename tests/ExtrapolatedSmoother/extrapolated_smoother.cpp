#include <gtest/gtest.h>

#include <random>
#include <vector>

#include "../../include/GMGPolar/gmgpolar.h"

#include "../../include/Residual/ResidualGive/residualGive.h"
#include "../../include/Residual/ResidualTake/residualTake.h"
#include "../../include/DirectSolver/DirectSolverGiveCustomLU/directSolverGiveCustomLU.h"
#include "../../include/DirectSolver/DirectSolverTakeCustomLU/directSolverTakeCustomLU.h"
#include "../../include/ExtrapolatedSmoother/ExtrapolatedSmootherGive/extrapolatedSmootherGive.h"
#include "../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTake/extrapolatedSmootherTake.h"

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

namespace ExtrapolatedExtrapolatedSmootherTest
{
Vector<double> generate_random_sample_data(const PolarGrid& grid, unsigned int seed)
{
    Vector<double> x(grid.numberOfNodes());
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-100.0, 100.0);
    for (int i = 0; i < x.size(); ++i) {
        x[i] = dist(gen);
    }
    return x;
}
} // namespace ExtrapolatedExtrapolatedSmootherTest

using namespace ExtrapolatedExtrapolatedSmootherTest;

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

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = true;
    int maxOpenMPThreads = 16;

    // "Take" requires cached values
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    ExtrapolatedSmootherGive smootherGive_operator(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                   DirBC_Interior, maxOpenMPThreads);
    ExtrapolatedSmootherTake smootherTake_operator(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                   DirBC_Interior, maxOpenMPThreads);

    Vector<double> rhs   = generate_random_sample_data(level.grid(), 69);
    Vector<double> start = generate_random_sample_data(level.grid(), 24);
    Vector<double> temp  = generate_random_sample_data(level.grid(), 8);

    Vector<double> solution_Give = start;
    smootherGive_operator.extrapolatedSmoothing(solution_Give, rhs, temp);

    Vector<double> solution_Take = start;
    smootherTake_operator.extrapolatedSmoothing(solution_Take, rhs, temp);

    ASSERT_EQ(solution_Give.size(), solution_Take.size());
    for (int index = 0; index < solution_Give.size(); index++) {
        MultiIndex alpha = level.grid().multiIndex(index);
        if (alpha[0] == 0 && !DirBC_Interior)
            ASSERT_NEAR(solution_Give[index], solution_Take[index], 1e-11);
        else
            ASSERT_NEAR(solution_Give[index], solution_Take[index], 1e-11);
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

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = false;
    int maxOpenMPThreads = 16;

    // "Take" requires cached values
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    ExtrapolatedSmootherGive smootherGive_operator(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                   DirBC_Interior, maxOpenMPThreads);
    ExtrapolatedSmootherTake smootherTake_operator(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                   DirBC_Interior, maxOpenMPThreads);

    Vector<double> rhs   = generate_random_sample_data(level.grid(), 69);
    Vector<double> start = generate_random_sample_data(level.grid(), 24);
    Vector<double> temp  = generate_random_sample_data(level.grid(), 8);

    Vector<double> solution_Give = start;
    smootherGive_operator.extrapolatedSmoothing(solution_Give, rhs, temp);

    Vector<double> solution_Take = start;
    smootherTake_operator.extrapolatedSmoothing(solution_Take, rhs, temp);

    ASSERT_EQ(solution_Give.size(), solution_Take.size());
    for (int index = 0; index < solution_Give.size(); index++) {
        MultiIndex alpha = level.grid().multiIndex(index);
        if (alpha[0] == 0 && !DirBC_Interior)
            ASSERT_NEAR(solution_Give[index], solution_Take[index], 1e-8);
        else
            ASSERT_NEAR(solution_Give[index], solution_Take[index], 1e-11);
    }
}

/* Test 2/2: */
/* Does the smoother converge to the DirectSolverGiveCustomLU solution? */

TEST(ExtrapolatedSmootherTest, SequentialExtrapolatedSmootherDirBC_Interior)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = true;
    int maxOpenMPThreads = 1;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);
    ExtrapolatedSmootherGive extrapolated_smoother_op(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                      DirBC_Interior, maxOpenMPThreads);

    const Vector<double> rhs         = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());

    smoother_solution = generate_random_sample_data(level.grid(), 69);
    for (int i_r = 0; i_r < level.grid().nr(); i_r += 2) {
        for (int i_theta = 0; i_theta < level.grid().ntheta(); i_theta += 2) {
            smoother_solution[level.grid().index(i_r, i_theta)] = discrete_solution[level.grid().index(i_r, i_theta)];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-12;

    while (infinity_norm(error) > precision) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

#pragma omp parallel for
        for (int i = 0; i < error.size(); i++) {
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

TEST(ExtrapolatedSmootherTest, ParallelExtrapolatedSmootherDirBC_Interior)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = true;
    int maxOpenMPThreads = 16;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);
    ExtrapolatedSmootherGive extrapolated_smoother_op(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                      DirBC_Interior, maxOpenMPThreads);

    const Vector<double> rhs         = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());

    smoother_solution = generate_random_sample_data(level.grid(), 69);
    for (int i_r = 0; i_r < level.grid().nr(); i_r += 2) {
        for (int i_theta = 0; i_theta < level.grid().ntheta(); i_theta += 2) {
            smoother_solution[level.grid().index(i_r, i_theta)] = discrete_solution[level.grid().index(i_r, i_theta)];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-12;

    while (infinity_norm(error) > precision) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

#pragma omp parallel for
        for (int i = 0; i < error.size(); i++) {
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

TEST(ExtrapolatedSmootherTest, SequentialExtrapolatedSmootherAcrossOrigin)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = false;
    int maxOpenMPThreads = 1;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);
    ExtrapolatedSmootherGive extrapolated_smoother_op(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                      DirBC_Interior, maxOpenMPThreads);

    const Vector<double> rhs         = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());

    smoother_solution = generate_random_sample_data(level.grid(), 69);
    for (int i_r = 0; i_r < level.grid().nr(); i_r += 2) {
        for (int i_theta = 0; i_theta < level.grid().ntheta(); i_theta += 2) {
            smoother_solution[level.grid().index(i_r, i_theta)] = discrete_solution[level.grid().index(i_r, i_theta)];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-8;

    while (infinity_norm(error) > precision) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

#pragma omp parallel for
        for (int i = 0; i < error.size(); i++) {
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

TEST(ExtrapolatedSmootherTest, ParallelExtrapolatedSmootherAcrossOrigin)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = false;
    int maxOpenMPThreads = 16;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);
    ExtrapolatedSmootherGive extrapolated_smoother_op(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                      DirBC_Interior, maxOpenMPThreads);

    const Vector<double> rhs         = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());

    smoother_solution = generate_random_sample_data(level.grid(), 69);
    for (int i_r = 0; i_r < level.grid().nr(); i_r += 2) {
        for (int i_theta = 0; i_theta < level.grid().ntheta(); i_theta += 2) {
            smoother_solution[level.grid().index(i_r, i_theta)] = discrete_solution[level.grid().index(i_r, i_theta)];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-8;

    while (infinity_norm(error) > 1e-8) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

#pragma omp parallel for
        for (int i = 0; i < error.size(); i++) {
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

TEST(ExtrapolatedSmootherTest, SequentialExtrapolatedSmootherDirBC_Interior_SmallestGrid)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.4, 0.45, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI, M_PI + M_PI / 8, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = true;
    int maxOpenMPThreads = 1;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);
    ExtrapolatedSmootherGive extrapolated_smoother_op(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                      DirBC_Interior, maxOpenMPThreads);

    const Vector<double> rhs         = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());

    smoother_solution = generate_random_sample_data(level.grid(), 69);
    for (int i_r = 0; i_r < level.grid().nr(); i_r += 2) {
        for (int i_theta = 0; i_theta < level.grid().ntheta(); i_theta += 2) {
            smoother_solution[level.grid().index(i_r, i_theta)] = discrete_solution[level.grid().index(i_r, i_theta)];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-12;

    while (infinity_norm(error) > 1e-12) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

#pragma omp parallel for
        for (int i = 0; i < error.size(); i++) {
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
    ASSERT_LT(iterations, 200);
    ASSERT_NEAR(infinity_norm(error), 0.0, precision);
}

TEST(ExtrapolatedSmootherTest, ParallelExtrapolatedSmootherDirBC_Interior_SmallestGrid)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.4, 0.45, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI, M_PI + M_PI / 8, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = true;
    int maxOpenMPThreads = 16;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);
    ExtrapolatedSmootherGive extrapolated_smoother_op(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                      DirBC_Interior, maxOpenMPThreads);

    const Vector<double> rhs         = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());

    smoother_solution = generate_random_sample_data(level.grid(), 69);
    for (int i_r = 0; i_r < level.grid().nr(); i_r += 2) {
        for (int i_theta = 0; i_theta < level.grid().ntheta(); i_theta += 2) {
            smoother_solution[level.grid().index(i_r, i_theta)] = discrete_solution[level.grid().index(i_r, i_theta)];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    double precision            = 1e-12;

    while (infinity_norm(error) > precision) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

#pragma omp parallel for
        for (int i = 0; i < error.size(); i++) {
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
    ASSERT_LT(iterations, 200);
    ASSERT_NEAR(infinity_norm(error), 0.0, precision);
}

TEST(ExtrapolatedSmootherTest, SequentialExtrapolatedSmootherAcrossOrigin_SmallestGrid)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.4, 0.45, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI, M_PI + M_PI / 8, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = false;
    int maxOpenMPThreads = 1;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);
    ExtrapolatedSmootherGive extrapolated_smoother_op(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                      DirBC_Interior, maxOpenMPThreads);

    const Vector<double> rhs         = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());

    smoother_solution = generate_random_sample_data(level.grid(), 69);
    for (int i_r = 0; i_r < level.grid().nr(); i_r += 2) {
        for (int i_theta = 0; i_theta < level.grid().ntheta(); i_theta += 2) {
            smoother_solution[level.grid().index(i_r, i_theta)] = discrete_solution[level.grid().index(i_r, i_theta)];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-8;

    while (infinity_norm(error) > 1e-8) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

#pragma omp parallel for
        for (int i = 0; i < error.size(); i++) {
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
    ASSERT_LT(iterations, 150);
    ASSERT_NEAR(infinity_norm(error), 0.0, precision);
}

TEST(ExtrapolatedSmootherTest, ParallelExtrapolatedSmootherAcrossOrigin_SmallestGrid)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.4, 0.45, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI, M_PI + M_PI / 8, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = false;
    int maxOpenMPThreads = 16;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);
    ExtrapolatedSmootherGive extrapolated_smoother_op(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                      DirBC_Interior, maxOpenMPThreads);

    const Vector<double> rhs         = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());

    smoother_solution = generate_random_sample_data(level.grid(), 69);
    for (int i_r = 0; i_r < level.grid().nr(); i_r += 2) {
        for (int i_theta = 0; i_theta < level.grid().ntheta(); i_theta += 2) {
            smoother_solution[level.grid().index(i_r, i_theta)] = discrete_solution[level.grid().index(i_r, i_theta)];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-8;

    while (infinity_norm(error) > 1e-8) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

#pragma omp parallel for
        for (int i = 0; i < error.size(); i++) {
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
    ASSERT_LT(iterations, 150);
    ASSERT_NEAR(infinity_norm(error), 0.0, precision);
}

/* We now test using "Take" */

TEST(ExtrapolatedSmootherTest, SequentialExtrapolatedSmootherTakeDirBC_Interior)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = true;
    int maxOpenMPThreads = 1;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualTake residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);
    ExtrapolatedSmootherTake extrapolated_smoother_op(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                      DirBC_Interior, maxOpenMPThreads);

    const Vector<double> rhs         = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());

    smoother_solution = generate_random_sample_data(level.grid(), 69);
    for (int i_r = 0; i_r < level.grid().nr(); i_r += 2) {
        for (int i_theta = 0; i_theta < level.grid().ntheta(); i_theta += 2) {
            smoother_solution[level.grid().index(i_r, i_theta)] = discrete_solution[level.grid().index(i_r, i_theta)];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-12;

    while (infinity_norm(error) > precision) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

#pragma omp parallel for
        for (int i = 0; i < error.size(); i++) {
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

TEST(ExtrapolatedSmootherTest, ParallelExtrapolatedSmootherTakeDirBC_Interior)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = true;
    int maxOpenMPThreads = 16;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualTake residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);
    ExtrapolatedSmootherTake extrapolated_smoother_op(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                      DirBC_Interior, maxOpenMPThreads);

    const Vector<double> rhs         = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());

    smoother_solution = generate_random_sample_data(level.grid(), 69);
    for (int i_r = 0; i_r < level.grid().nr(); i_r += 2) {
        for (int i_theta = 0; i_theta < level.grid().ntheta(); i_theta += 2) {
            smoother_solution[level.grid().index(i_r, i_theta)] = discrete_solution[level.grid().index(i_r, i_theta)];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-12;

    while (infinity_norm(error) > precision) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

#pragma omp parallel for
        for (int i = 0; i < error.size(); i++) {
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

TEST(ExtrapolatedSmootherTest, SequentialExtrapolatedSmootherTakeAcrossOrigin)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = false;
    int maxOpenMPThreads = 1;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualTake residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);
    ExtrapolatedSmootherTake extrapolated_smoother_op(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                      DirBC_Interior, maxOpenMPThreads);

    const Vector<double> rhs         = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());

    smoother_solution = generate_random_sample_data(level.grid(), 69);
    for (int i_r = 0; i_r < level.grid().nr(); i_r += 2) {
        for (int i_theta = 0; i_theta < level.grid().ntheta(); i_theta += 2) {
            smoother_solution[level.grid().index(i_r, i_theta)] = discrete_solution[level.grid().index(i_r, i_theta)];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-8;

    while (infinity_norm(error) > precision) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

#pragma omp parallel for
        for (int i = 0; i < error.size(); i++) {
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

TEST(ExtrapolatedSmootherTest, ParallelExtrapolatedSmootherTakeAcrossOrigin)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = false;
    int maxOpenMPThreads = 16;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualTake residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);
    ExtrapolatedSmootherTake extrapolated_smoother_op(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                      DirBC_Interior, maxOpenMPThreads);

    const Vector<double> rhs         = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());

    smoother_solution = generate_random_sample_data(level.grid(), 69);
    for (int i_r = 0; i_r < level.grid().nr(); i_r += 2) {
        for (int i_theta = 0; i_theta < level.grid().ntheta(); i_theta += 2) {
            smoother_solution[level.grid().index(i_r, i_theta)] = discrete_solution[level.grid().index(i_r, i_theta)];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-8;

    while (infinity_norm(error) > 1e-8) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

#pragma omp parallel for
        for (int i = 0; i < error.size(); i++) {
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

TEST(ExtrapolatedSmootherTest, SequentialExtrapolatedSmootherTakeDirBC_Interior_SmallestGrid)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.4, 0.45, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI, M_PI + M_PI / 8, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = true;
    int maxOpenMPThreads = 1;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualTake residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);
    ExtrapolatedSmootherTake extrapolated_smoother_op(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                      DirBC_Interior, maxOpenMPThreads);

    const Vector<double> rhs         = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());

    smoother_solution = generate_random_sample_data(level.grid(), 69);
    for (int i_r = 0; i_r < level.grid().nr(); i_r += 2) {
        for (int i_theta = 0; i_theta < level.grid().ntheta(); i_theta += 2) {
            smoother_solution[level.grid().index(i_r, i_theta)] = discrete_solution[level.grid().index(i_r, i_theta)];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-12;

    while (infinity_norm(error) > 1e-12) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

#pragma omp parallel for
        for (int i = 0; i < error.size(); i++) {
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
    ASSERT_LT(iterations, 200);
    ASSERT_NEAR(infinity_norm(error), 0.0, precision);
}

TEST(ExtrapolatedSmootherTest, ParallelExtrapolatedSmootherTakeDirBC_Interior_SmallestGrid)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.4, 0.45, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI, M_PI + M_PI / 8, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = true;
    int maxOpenMPThreads = 16;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualTake residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);
    ExtrapolatedSmootherTake extrapolated_smoother_op(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                      DirBC_Interior, maxOpenMPThreads);

    const Vector<double> rhs         = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());

    smoother_solution = generate_random_sample_data(level.grid(), 69);
    for (int i_r = 0; i_r < level.grid().nr(); i_r += 2) {
        for (int i_theta = 0; i_theta < level.grid().ntheta(); i_theta += 2) {
            smoother_solution[level.grid().index(i_r, i_theta)] = discrete_solution[level.grid().index(i_r, i_theta)];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    double precision            = 1e-12;

    while (infinity_norm(error) > precision) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

#pragma omp parallel for
        for (int i = 0; i < error.size(); i++) {
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
    ASSERT_LT(iterations, 200);
    ASSERT_NEAR(infinity_norm(error), 0.0, precision);
}

TEST(ExtrapolatedSmootherTest, SequentialExtrapolatedSmootherTakeAcrossOrigin_SmallestGrid)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.4, 0.45, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI, M_PI + M_PI / 8, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = false;
    int maxOpenMPThreads = 1;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualTake residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);
    ExtrapolatedSmootherTake extrapolated_smoother_op(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                      DirBC_Interior, maxOpenMPThreads);

    const Vector<double> rhs         = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());

    smoother_solution = generate_random_sample_data(level.grid(), 69);
    for (int i_r = 0; i_r < level.grid().nr(); i_r += 2) {
        for (int i_theta = 0; i_theta < level.grid().ntheta(); i_theta += 2) {
            smoother_solution[level.grid().index(i_r, i_theta)] = discrete_solution[level.grid().index(i_r, i_theta)];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-8;

    while (infinity_norm(error) > 1e-8) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

#pragma omp parallel for
        for (int i = 0; i < error.size(); i++) {
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
    ASSERT_LT(iterations, 150);
    ASSERT_NEAR(infinity_norm(error), 0.0, precision);
}

TEST(ExtrapolatedSmootherTest, ParallelExtrapolatedSmootherTakeAcrossOrigin_SmallestGrid)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.4, 0.45, 0.9, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI, M_PI + M_PI / 8, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = false;
    int maxOpenMPThreads = 16;

    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::IMPLICIT_EXTRAPOLATION, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualTake residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);
    ExtrapolatedSmootherTake extrapolated_smoother_op(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                      DirBC_Interior, maxOpenMPThreads);

    const Vector<double> rhs         = generate_random_sample_data(level.grid(), 42);
    Vector<double> discrete_solution = rhs;
    solver_op.solveInPlace(discrete_solution);

    Vector<double> temp(level.grid().numberOfNodes());
    Vector<double> smoother_solution(level.grid().numberOfNodes());
    Vector<double> error(level.grid().numberOfNodes());

    smoother_solution = generate_random_sample_data(level.grid(), 69);
    for (int i_r = 0; i_r < level.grid().nr(); i_r += 2) {
        for (int i_theta = 0; i_theta < level.grid().ntheta(); i_theta += 2) {
            smoother_solution[level.grid().index(i_r, i_theta)] = discrete_solution[level.grid().index(i_r, i_theta)];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < error.size(); i++) {
        error[i] = discrete_solution[i] - smoother_solution[i];
    }

    int iterations              = 0;
    bool max_iterations_reached = false;
    const int max_iterations    = 10000;
    const double precision      = 1e-8;

    while (infinity_norm(error) > 1e-8) {
        extrapolated_smoother_op.extrapolatedSmoothing(smoother_solution, rhs, temp);

#pragma omp parallel for
        for (int i = 0; i < error.size(); i++) {
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
    ASSERT_LT(iterations, 150);
    ASSERT_NEAR(infinity_norm(error), 0.0, precision);
}

int main(int argc, char* argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
