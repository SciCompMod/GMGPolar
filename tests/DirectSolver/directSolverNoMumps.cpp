#include <gtest/gtest.h>

#include <random>
#include <vector>

#include "../../include/GMGPolar/gmgpolar.h"

#include "../../include/Residual/ResidualGive/residualGive.h"
#include "../../include/DirectSolver/DirectSolver-COO-MUMPS-Give/directSolverGive.h"
#include "../../include/DirectSolver/DirectSolver-COO-MUMPS-Take/directSolverTake.h"
#include "../../include/DirectSolver/DirectSolver-CSR-LU-Give/directSolverGiveCustomLU.h"
#include "../../include/DirectSolver/DirectSolver-CSR-LU-Take/directSolverTakeCustomLU.h"

#include "../../include/InputFunctions/domainGeometry.h"
#include "../../include/InputFunctions/densityProfileCoefficients.h"
#include "../../include/InputFunctions/boundaryConditions.h"
#include "../../include/InputFunctions/sourceTerm.h"
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

namespace DirectSolverTestNoMumps
{
Vector<double> generate_random_sample_data(const PolarGrid& grid, unsigned int seed)
{
    Vector<double> x("x", grid.numberOfNodes());
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-100.0, 100.0);
    for (int i = 0; i < x.size(); ++i) {
        x(i) = dist(gen);
    }
    return x;
}
} // namespace DirectSolverTestNoMumps

using namespace DirectSolverTestNoMumps;

/* Test 1/2: */
/* Does the Take and Give Implementation match up? */

TEST(DirectSolverTestNoMumps, directSolver_DirBC_Interior)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.678 * Rmax;
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
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, false);

    DirectSolverTakeCustomLU directSolverGive_operator(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                       DirBC_Interior, maxOpenMPThreads);
    DirectSolverGiveCustomLU directSolverTake_operator(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                       DirBC_Interior, maxOpenMPThreads);

    Vector<double> rhs = generate_random_sample_data(level.grid(), 69);

    Vector<double> solution_Give = rhs;
    directSolverGive_operator.solveInPlace(solution_Give);

    Vector<double> solution_Take = rhs;
    directSolverTake_operator.solveInPlace(solution_Take);

    ASSERT_EQ(solution_Give.size(), solution_Take.size());
    for (int index = 0; index < solution_Give.size(); index++) {
        MultiIndex alpha = level.grid().multiIndex(index);
        if (alpha[0] == 0 && !DirBC_Interior)
            ASSERT_NEAR(solution_Give[index], solution_Take[index], 1e-11);
        else
            ASSERT_NEAR(solution_Give[index], solution_Take[index], 1e-11);
    }
}

TEST(DirectSolverTestNoMumps, directSolver_AcrossOrigin)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.678 * Rmax;
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

    DirectSolverGiveCustomLU directSolverGive_operator(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                       DirBC_Interior, maxOpenMPThreads);
    DirectSolverTakeCustomLU directSolverTake_operator(level.grid(), level.levelCache(), domain_geometry, *coefficients,
                                                       DirBC_Interior, maxOpenMPThreads);

    Vector<double> rhs = generate_random_sample_data(level.grid(), 69);

    Vector<double> solution_Give = rhs;
    directSolverGive_operator.solveInPlace(solution_Give);

    Vector<double> solution_Take = rhs;
    directSolverTake_operator.solveInPlace(solution_Take);

    ASSERT_EQ(solution_Give.size(), solution_Take.size());
    for (int index = 0; index < solution_Give.size(); index++) {
        MultiIndex alpha = level.grid().multiIndex(index);
        if (alpha[0] == 0 && !DirBC_Interior)
            ASSERT_NEAR(solution_Give[index], solution_Take[index], 1e-8);
        else
            ASSERT_NEAR(solution_Give[index], solution_Take[index], 1e-9);
    }
}

/* Test 2/2: */
/* Are the DirectSolver and Residual are compatible with each other? */

/* -------- */
/* Circular */
/* -------- */

TEST(DirectSolverTestNoMumps_CircularGeometry, SequentialDirectSolverDirBC_Interior_CircularGeometry)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};
    double Rmax = radii.back();

    CircularGeometry domain_geometry(Rmax);

    double alpha_jump = 0.66 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<SonnendruckerGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR2_Boundary_CircularGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CircularGeometry>(Rmax);

    bool DirBC_Interior                     = true;
    int maxOpenMPThreads                    = 1;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-12);
}

TEST(DirectSolverTestNoMumps_CircularGeometry, ParallelDirectSolverDirBC_Interior_CircularGeometry)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax = radii.back();

    CircularGeometry domain_geometry(Rmax);

    double alpha_jump = 0.66 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<SonnendruckerGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR2_Boundary_CircularGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CircularGeometry>(Rmax);

    bool DirBC_Interior                     = true;
    int maxOpenMPThreads                    = 16;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-12);
}

TEST(DirectSolverTestNoMumps_CircularGeometry, SequentialDirectSolverAcrossOrigin_CircularGeometry)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax = radii.back();

    CircularGeometry domain_geometry(Rmax);

    double alpha_jump = 0.66 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<SonnendruckerGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR2_Boundary_CircularGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CircularGeometry>(Rmax);

    bool DirBC_Interior                     = false;
    int maxOpenMPThreads                    = 1;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-7);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-8);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-8);
}

TEST(DirectSolverTestNoMumps_CircularGeometry, ParallelDirectSolverAcrossOrigin_CircularGeometry)
{
    std::vector<double> radii  = {1e-5, 0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax = radii.back();

    CircularGeometry domain_geometry(Rmax);

    double alpha_jump = 0.66 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<SonnendruckerGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR2_Boundary_CircularGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CircularGeometry>(Rmax);

    bool DirBC_Interior                     = false;
    int maxOpenMPThreads                    = 16;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-7);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-8);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-8);
}

/* --------- */
/* Shafranov */
/* --------- */

TEST(DirectSolverTestNoMumps_ShafranovGeometry, DirectSolverDirBC_Interior_ShafranovGeometry)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 0.2;

    ShafranovGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump                                        = 0.4837 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR6_Boundary_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<CartesianR6_ZoniGyro_ShafranovGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior                     = true;
    int maxOpenMPThreads                    = 16;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-12);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-12);
}

TEST(DirectSolverTestNoMumps_ShafranovGeometry, DirectSolverAcrossOrigin_ShafranovGeometry)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 0.2;

    ShafranovGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump                                        = 0.4837 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR6_Boundary_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<CartesianR6_ZoniGyro_ShafranovGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior                     = false;
    int maxOpenMPThreads                    = 16;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-7);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-8);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-8);
}

/* ------ */
/* Czarny */
/* ------ */

TEST(DirectSolverTestNoMumps_CzarnyGeometry, DirectSolverDirBC_Interior_CzarnyGeometry)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.678 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior                     = true;
    int maxOpenMPThreads                    = 16;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-12);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-12);
}

TEST(DirectSolverTestNoMumps_CzarnyGeometry, DirectSolverAcrossOrigin_CzarnyGeometry)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.678 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior                     = false;
    int maxOpenMPThreads                    = 16;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-7);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-8);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-8);
}

/* ------ */
/* Culham */
/* ------ */

TEST(DirectSolverTestNoMumps_CulhamGeometry, DirectSolverDirBC_Interior_CulhamGeometry)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax = radii.back();

    CulhamGeometry domain_geometry(Rmax);

    double alpha_jump = 0.678 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<Refined_Boundary_CulhamGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<Refined_ZoniShiftedGyro_CulhamGeometry>(Rmax);

    bool DirBC_Interior                     = true;
    int maxOpenMPThreads                    = 16;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
}

TEST(DirectSolverTestNoMumps_CulhamGeometry, DirectSolverAcrossOrigin_CulhamGeometry)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax = radii.back();

    CulhamGeometry domain_geometry(Rmax);

    double alpha_jump = 0.678 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<Refined_Boundary_CulhamGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<Refined_ZoniShiftedGyro_CulhamGeometry>(Rmax);

    bool DirBC_Interior                     = false;
    int maxOpenMPThreads                    = 16;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-7);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-7);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-8);
}

/* We adjust the PolarGrid to increase the precision */

TEST(DirectSolverTestNoMumps_CircularGeometry, DirectSolverAcrossOriginHigherPrecision_CircularGeometry)
{
    std::vector<double> radii  = {1e-5,          1.441 * 1e-5,
                                  3.8833 * 1e-5, 8.7666 * 1e-5,
                                  1.8533 * 1e-4, 3.806 * 1e-4,
                                  7.713 * 1e-4,  1.55265 * 1e-3,
                                  3.1153 * 1e-3, 6.2406 * 1e-3,
                                  0.01249125,    0.0249925,
                                  0.049995,      0.1,
                                  0.2,           0.25,
                                  0.5,           0.8,
                                  0.9,           0.95,
                                  1.2,           1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax = radii.back();

    CircularGeometry domain_geometry(Rmax);

    double alpha_jump = 0.66 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<SonnendruckerGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR2_Boundary_CircularGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CircularGeometry>(Rmax);

    bool DirBC_Interior                     = false;
    int maxOpenMPThreads                    = 1;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-9);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-10);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-10);
}

TEST(DirectSolverTestNoMumps_CircularGeometry, DirectSolverAcrossOriginHigherPrecision2_CircularGeometry)
{
    std::vector<double> radii  = {0.15, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax = radii.back();

    CircularGeometry domain_geometry(Rmax);

    double alpha_jump = 0.66 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<SonnendruckerGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR2_Boundary_CircularGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CircularGeometry>(Rmax);

    bool DirBC_Interior                     = false;
    int maxOpenMPThreads                    = 1;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = false;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverGiveCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-10);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
}

/* Same test now using Take */

TEST(DirectSolverTakeCustomLUTest_CircularGeometry, SequentialDirectSolverDirBC_Interior_CircularGeometry)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};
    double Rmax = radii.back();

    CircularGeometry domain_geometry(Rmax);

    double alpha_jump = 0.66 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<SonnendruckerGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR2_Boundary_CircularGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CircularGeometry>(Rmax);

    bool DirBC_Interior                     = true;
    int maxOpenMPThreads                    = 1;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-12);
}

TEST(DirectSolverTakeCustomLUTest_CircularGeometry, ParallelDirectSolverDirBC_Interior_CircularGeometry)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax = radii.back();

    CircularGeometry domain_geometry(Rmax);

    double alpha_jump = 0.66 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<SonnendruckerGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR2_Boundary_CircularGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CircularGeometry>(Rmax);

    bool DirBC_Interior                     = true;
    int maxOpenMPThreads                    = 16;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-12);
}

TEST(DirectSolverTakeCustomLUTest_CircularGeometry, SequentialDirectSolverAcrossOrigin_CircularGeometry)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax = radii.back();

    CircularGeometry domain_geometry(Rmax);

    double alpha_jump = 0.66 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<SonnendruckerGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR2_Boundary_CircularGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CircularGeometry>(Rmax);

    bool DirBC_Interior                     = false;
    int maxOpenMPThreads                    = 1;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-7);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-8);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-8);
}

TEST(DirectSolverTakeCustomLUTest_CircularGeometry, ParallelDirectSolverAcrossOrigin_CircularGeometry)
{
    std::vector<double> radii  = {1e-5, 0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax = radii.back();

    CircularGeometry domain_geometry(Rmax);

    double alpha_jump = 0.66 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<SonnendruckerGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR2_Boundary_CircularGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CircularGeometry>(Rmax);

    bool DirBC_Interior                     = false;
    int maxOpenMPThreads                    = 16;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-7);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-8);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-8);
}

/* --------- */
/* Shafranov */
/* --------- */

TEST(DirectSolverTakeCustomLUTest_ShafranovGeometry, DirectSolverDirBC_Interior_ShafranovGeometry)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 0.2;

    ShafranovGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump                                        = 0.4837 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR6_Boundary_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<CartesianR6_ZoniGyro_ShafranovGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior                     = true;
    int maxOpenMPThreads                    = 16;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-12);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-12);
}

TEST(DirectSolverTakeCustomLUTest_ShafranovGeometry, DirectSolverAcrossOrigin_ShafranovGeometry)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 0.2;

    ShafranovGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump                                        = 0.4837 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR6_Boundary_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<CartesianR6_ZoniGyro_ShafranovGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior                     = false;
    int maxOpenMPThreads                    = 16;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-7);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-8);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-8);
}

/* ------ */
/* Czarny */
/* ------ */

TEST(DirectSolverTakeCustomLUTest_CzarnyGeometry, DirectSolverDirBC_Interior_CzarnyGeometry)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.678 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior                     = true;
    int maxOpenMPThreads                    = 16;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-12);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-12);
}

TEST(DirectSolverTakeCustomLUTest_CzarnyGeometry, DirectSolverAcrossOrigin_CzarnyGeometry)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.678 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior                     = false;
    int maxOpenMPThreads                    = 16;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-7);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-8);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-8);
}

/* ------ */
/* Culham */
/* ------ */

TEST(DirectSolverTakeCustomLUTest_CulhamGeometry, DirectSolverDirBC_Interior_CulhamGeometry)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax = radii.back();

    CulhamGeometry domain_geometry(Rmax);

    double alpha_jump = 0.678 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<Refined_Boundary_CulhamGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<Refined_ZoniShiftedGyro_CulhamGeometry>(Rmax);

    bool DirBC_Interior                     = true;
    int maxOpenMPThreads                    = 16;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
}

TEST(DirectSolverTakeCustomLUTest_CulhamGeometry, DirectSolverAcrossOrigin_CulhamGeometry)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax = radii.back();

    CulhamGeometry domain_geometry(Rmax);

    double alpha_jump = 0.678 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<Refined_Boundary_CulhamGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<Refined_ZoniShiftedGyro_CulhamGeometry>(Rmax);

    bool DirBC_Interior                     = false;
    int maxOpenMPThreads                    = 16;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-7);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-7);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-8);
}

TEST(DirectSolverTakeCustomLUTest_CircularGeometry, DirectSolverAcrossOriginHigherPrecision_CircularGeometry)
{
    std::vector<double> radii  = {1e-5,          1.441 * 1e-5,
                                  3.8833 * 1e-5, 8.7666 * 1e-5,
                                  1.8533 * 1e-4, 3.806 * 1e-4,
                                  7.713 * 1e-4,  1.55265 * 1e-3,
                                  3.1153 * 1e-3, 6.2406 * 1e-3,
                                  0.01249125,    0.0249925,
                                  0.049995,      0.1,
                                  0.2,           0.25,
                                  0.5,           0.8,
                                  0.9,           0.95,
                                  1.2,           1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax = radii.back();

    CircularGeometry domain_geometry(Rmax);

    double alpha_jump = 0.66 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<SonnendruckerGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR2_Boundary_CircularGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CircularGeometry>(Rmax);

    bool DirBC_Interior                     = false;
    int maxOpenMPThreads                    = 1;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-9);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-10);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-10);
}

TEST(DirectSolverTakeCustomLUTest_CircularGeometry, DirectSolverAcrossOriginHigherPrecision2_CircularGeometry)
{
    std::vector<double> radii  = {0.15, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax = radii.back();

    CircularGeometry domain_geometry(Rmax);

    double alpha_jump = 0.66 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<SonnendruckerGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR2_Boundary_CircularGeometry>(Rmax);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CircularGeometry>(Rmax);

    bool DirBC_Interior                     = false;
    int maxOpenMPThreads                    = 1;
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, 0);

    DirectSolverTakeCustomLU solver_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualGive residual_op(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                             maxOpenMPThreads);

    ConstVector<double> rhs = generate_random_sample_data(level.grid(), 42);
    Vector<double> solution("solution", rhs.size());
    Kokkos::deep_copy(solution, rhs);
    solver_op.solveInPlace(solution);

    Vector<double> residuum("residuum", level.grid().numberOfNodes());
    residual_op.computeResidual(residuum, rhs, solution);

    ASSERT_NEAR(l1_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
    ASSERT_NEAR(l2_norm(ConstVector<double>(residuum)), 0.0, 1e-11);
    ASSERT_NEAR(infinity_norm(ConstVector<double>(residuum)), 0.0, 1e-12);
}
