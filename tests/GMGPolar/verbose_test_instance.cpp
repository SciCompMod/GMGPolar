#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>

#include "../../include/GMGPolar/gmgpolar.h"

TEST(GMGPolar, SolvesPolarR6_WithTolerancesVerbose)
{
    /* PolarGrid */
    double R0                              = 1e-8;
    double Rmax                            = 2.0;
    int nr_exp                             = 4;
    int ntheta_exp                         = 5;
    double refinement_radius               = 0.5 * Rmax;
    int anisotropic_factor                 = 3;
    int divideBy2                          = 1;
    std::optional<double> splitting_radius = std::nullopt;
    PolarGrid grid(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, splitting_radius);

    /* PolarR6 */
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
    std::unique_ptr<DomainGeometry> domain_geometry =
        std::make_unique<ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<ExactSolution> exact_solution =
        std::make_unique<PolarR6_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniGyroCoefficients>(Rmax, refinement_radius);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniGyro_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);

    /* GMGPolar settings */
    int verbose   = 1;
    bool paraview = false;

    int maxOpenMPThreads         = 4;
    double threadReductionFactor = 1.0;
    omp_set_num_threads(maxOpenMPThreads);

    bool DirBC_Interior                                 = false;
    StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_GIVE;
    bool cacheDensityProfileCoefficients                = true;
    bool cacheDomainGeometry                            = false;

    ExtrapolationType extrapolation   = ExtrapolationType::IMPLICIT_EXTRAPOLATION;
    int maxLevels                     = -1;
    int preSmoothingSteps             = 1;
    int postSmoothingSteps            = 1;
    MultigridCycleType multigridCycle = MultigridCycleType::V_CYCLE;
    bool FMG                          = true;
    int FMG_iterations                = 3;
    MultigridCycleType FMG_cycle      = MultigridCycleType::F_CYCLE;

    int maxIterations                 = 150;
    ResidualNormType residualNormType = ResidualNormType::INFINITY_NORM;
    double absoluteTolerance          = 1e-8;
    double relativeTolerance          = 1e-8;

    // Create GMGPolar solver
    GMGPolar solver(grid, *domain_geometry.get(), *coefficients.get());
    solver.verbose(verbose);
    solver.paraview(paraview);
    solver.maxOpenMPThreads(maxOpenMPThreads);
    solver.threadReductionFactor(threadReductionFactor);
    solver.DirBC_Interior(DirBC_Interior);
    solver.stencilDistributionMethod(stencilDistributionMethod);
    solver.cacheDensityProfileCoefficients(cacheDensityProfileCoefficients);
    solver.cacheDomainGeometry(cacheDomainGeometry);
    solver.extrapolation(extrapolation);
    solver.maxLevels(maxLevels);
    solver.preSmoothingSteps(preSmoothingSteps);
    solver.postSmoothingSteps(postSmoothingSteps);
    solver.multigridCycle(multigridCycle);
    solver.FMG(FMG);
    solver.FMG_iterations(FMG_iterations);
    solver.FMG_cycle(FMG_cycle);
    solver.maxIterations(maxIterations);
    solver.residualNormType(residualNormType);
    solver.absoluteTolerance(absoluteTolerance);
    solver.relativeTolerance(relativeTolerance);

    solver.setup();
    solver.setSolution(exact_solution.get());
    solver.solve(*boundary_conditions.get(), *source_term.get());

    Vector<double>& solution     = solver.solution();
    const PolarGrid& finest_grid = solver.grid();

    EXPECT_EQ(solver.verbose(), verbose);
    EXPECT_EQ(solver.paraview(), paraview);
    EXPECT_EQ(solver.maxOpenMPThreads(), maxOpenMPThreads);
    EXPECT_DOUBLE_EQ(solver.threadReductionFactor(), threadReductionFactor);
    EXPECT_EQ(solver.DirBC_Interior(), DirBC_Interior);
    EXPECT_EQ(solver.stencilDistributionMethod(), stencilDistributionMethod);
    EXPECT_EQ(solver.cacheDensityProfileCoefficients(), cacheDensityProfileCoefficients);
    EXPECT_EQ(solver.cacheDomainGeometry(), cacheDomainGeometry);

    EXPECT_EQ(solver.extrapolation(), extrapolation);
    EXPECT_EQ(solver.maxLevels(), maxLevels);
    EXPECT_EQ(solver.preSmoothingSteps(), preSmoothingSteps);
    EXPECT_EQ(solver.postSmoothingSteps(), postSmoothingSteps);
    EXPECT_EQ(solver.multigridCycle(), multigridCycle);
    EXPECT_EQ(solver.FMG(), FMG);
    EXPECT_EQ(solver.FMG_iterations(), FMG_iterations);
    EXPECT_EQ(solver.FMG_cycle(), FMG_cycle);

    EXPECT_EQ(solver.maxIterations(), maxIterations);
    EXPECT_EQ(solver.residualNormType(), residualNormType);
    EXPECT_DOUBLE_EQ(solver.absoluteTolerance().value(), absoluteTolerance);
    EXPECT_DOUBLE_EQ(solver.relativeTolerance().value(), relativeTolerance);

    // --- Verify results ---
    int number_of_iterations                             = solver.numberOfIterations();
    std::optional<double> exact_error_weighted_euclidean = solver.exactErrorWeightedEuclidean();
    std::optional<double> exact_infinity_error           = solver.exactErrorInfinity();
    double reductionFactor                               = solver.meanResidualReductionFactor();

    ASSERT_TRUE(exact_error_weighted_euclidean.has_value());
    ASSERT_TRUE(exact_infinity_error.has_value());

    EXPECT_LE(number_of_iterations, 15);
    EXPECT_LE(exact_error_weighted_euclidean.value(), 3e-6);
    EXPECT_LE(exact_infinity_error.value(), 7e-6);
    EXPECT_GT(reductionFactor, 0.5);
    EXPECT_LT(reductionFactor, 0.7);

    // --- Verify timings (all available must be non-negative) ---
    EXPECT_GE(solver.timeSetupTotal(), 0.0);
    EXPECT_GE(solver.timeSetupCreateLevels(), 0.0);
    EXPECT_GE(solver.timeSetupRHS(), 0.0);
    EXPECT_GE(solver.timeSetupSmoother(), 0.0);
    EXPECT_GE(solver.timeSetupDirectSolver(), 0.0);

    EXPECT_GE(solver.timeSolveTotal(), 0.0);
    EXPECT_GE(solver.timeSolveInitialApproximation(), 0.0);
    EXPECT_GE(solver.timeSolveMultigridIterations(), 0.0);
    EXPECT_GE(solver.timeCheckConvergence(), 0.0);
    EXPECT_GE(solver.timeCheckExactError(), 0.0);

    EXPECT_GE(solver.timeAvgMGCTotal(), 0.0);
    EXPECT_GE(solver.timeAvgMGCPreSmoothing(), 0.0);
    EXPECT_GE(solver.timeAvgMGCPostSmoothing(), 0.0);
    EXPECT_GE(solver.timeAvgMGCResidual(), 0.0);
    EXPECT_GE(solver.timeAvgMGCDirectSolver(), 0.0);
}

TEST(GMGPolar, SolvesPolarR6_WithoutTolerances)

{
    /* PolarGrid */
    double R0                              = 1e-8;
    double Rmax                            = 2.0;
    int nr_exp                             = 4;
    int ntheta_exp                         = 5;
    double refinement_radius               = 0.5 * Rmax;
    int anisotropic_factor                 = 3;
    int divideBy2                          = 1;
    std::optional<double> splitting_radius = std::nullopt;
    PolarGrid grid(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, splitting_radius);

    /* PolarR6 */
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
    std::unique_ptr<DomainGeometry> domain_geometry =
        std::make_unique<ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<ExactSolution> exact_solution =
        std::make_unique<PolarR6_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniGyroCoefficients>(Rmax, refinement_radius);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniGyro_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);

    /* GMGPolar settings */
    int verbose   = 0;
    bool paraview = false;

    int maxOpenMPThreads         = 4;
    double threadReductionFactor = 1.0;
    omp_set_num_threads(maxOpenMPThreads);

    bool DirBC_Interior                                 = false;
    StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_GIVE;
    bool cacheDensityProfileCoefficients                = true;
    bool cacheDomainGeometry                            = false;

    ExtrapolationType extrapolation   = ExtrapolationType::IMPLICIT_EXTRAPOLATION;
    int maxLevels                     = -1;
    int preSmoothingSteps             = 1;
    int postSmoothingSteps            = 1;
    MultigridCycleType multigridCycle = MultigridCycleType::V_CYCLE;
    bool FMG                          = true;
    int FMG_iterations                = 3;
    MultigridCycleType FMG_cycle      = MultigridCycleType::F_CYCLE;

    int maxIterations                       = 150;
    ResidualNormType residualNormType       = ResidualNormType::INFINITY_NORM;
    std::optional<double> absoluteTolerance = std::nullopt;
    std::optional<double> relativeTolerance = std::nullopt;

    // Create GMGPolar solver
    GMGPolar solver(grid, *domain_geometry.get(), *coefficients.get());
    solver.verbose(verbose);
    solver.paraview(paraview);
    solver.maxOpenMPThreads(maxOpenMPThreads);
    solver.threadReductionFactor(threadReductionFactor);
    solver.DirBC_Interior(DirBC_Interior);
    solver.stencilDistributionMethod(stencilDistributionMethod);
    solver.cacheDensityProfileCoefficients(cacheDensityProfileCoefficients);
    solver.cacheDomainGeometry(cacheDomainGeometry);
    solver.extrapolation(extrapolation);
    solver.maxLevels(maxLevels);
    solver.preSmoothingSteps(preSmoothingSteps);
    solver.postSmoothingSteps(postSmoothingSteps);
    solver.multigridCycle(multigridCycle);
    solver.FMG(FMG);
    solver.FMG_iterations(FMG_iterations);
    solver.FMG_cycle(FMG_cycle);
    solver.maxIterations(maxIterations);
    solver.residualNormType(residualNormType);
    solver.absoluteTolerance(absoluteTolerance);
    solver.relativeTolerance(relativeTolerance);

    solver.setup();
    solver.setSolution(exact_solution.get());
    solver.solve(*boundary_conditions.get(), *source_term.get());

    Vector<double>& solution     = solver.solution();
    const PolarGrid& finest_grid = solver.grid();

    EXPECT_EQ(solver.verbose(), verbose);
    EXPECT_EQ(solver.paraview(), paraview);
    EXPECT_EQ(solver.maxOpenMPThreads(), maxOpenMPThreads);
    EXPECT_DOUBLE_EQ(solver.threadReductionFactor(), threadReductionFactor);
    EXPECT_EQ(solver.DirBC_Interior(), DirBC_Interior);
    EXPECT_EQ(solver.stencilDistributionMethod(), stencilDistributionMethod);
    EXPECT_EQ(solver.cacheDensityProfileCoefficients(), cacheDensityProfileCoefficients);
    EXPECT_EQ(solver.cacheDomainGeometry(), cacheDomainGeometry);

    EXPECT_EQ(solver.extrapolation(), extrapolation);
    EXPECT_EQ(solver.maxLevels(), maxLevels);
    EXPECT_EQ(solver.preSmoothingSteps(), preSmoothingSteps);
    EXPECT_EQ(solver.postSmoothingSteps(), postSmoothingSteps);
    EXPECT_EQ(solver.multigridCycle(), multigridCycle);
    EXPECT_EQ(solver.FMG(), FMG);
    EXPECT_EQ(solver.FMG_iterations(), FMG_iterations);
    EXPECT_EQ(solver.FMG_cycle(), FMG_cycle);

    EXPECT_EQ(solver.maxIterations(), maxIterations);
    EXPECT_EQ(solver.residualNormType(), residualNormType);
    EXPECT_EQ(solver.absoluteTolerance(), absoluteTolerance);
    EXPECT_EQ(solver.relativeTolerance(), relativeTolerance);

    // --- Verify results ---
    int number_of_iterations                             = solver.numberOfIterations();
    std::optional<double> exact_error_weighted_euclidean = solver.exactErrorWeightedEuclidean();
    std::optional<double> exact_infinity_error           = solver.exactErrorInfinity();
    double reductionFactor                               = solver.meanResidualReductionFactor();

    ASSERT_TRUE(exact_error_weighted_euclidean.has_value());
    ASSERT_TRUE(exact_infinity_error.has_value());

    EXPECT_EQ(number_of_iterations, 150);
    EXPECT_LE(exact_error_weighted_euclidean.value(), 3e-6);
    EXPECT_LE(exact_infinity_error.value(), 7e-6);
    EXPECT_LE(reductionFactor, 1.0);
    EXPECT_GE(reductionFactor, 0.0);

    // --- Verify timings (all available must be non-negative) ---
    EXPECT_GE(solver.timeSetupTotal(), 0.0);
    EXPECT_GE(solver.timeSetupCreateLevels(), 0.0);
    EXPECT_GE(solver.timeSetupRHS(), 0.0);
    EXPECT_GE(solver.timeSetupSmoother(), 0.0);
    EXPECT_GE(solver.timeSetupDirectSolver(), 0.0);

    EXPECT_GE(solver.timeSolveTotal(), 0.0);
    EXPECT_GE(solver.timeSolveInitialApproximation(), 0.0);
    EXPECT_GE(solver.timeSolveMultigridIterations(), 0.0);
    EXPECT_GE(solver.timeCheckConvergence(), 0.0);
    EXPECT_GE(solver.timeCheckExactError(), 0.0);

    EXPECT_GE(solver.timeAvgMGCTotal(), 0.0);
    EXPECT_GE(solver.timeAvgMGCPreSmoothing(), 0.0);
    EXPECT_GE(solver.timeAvgMGCPostSmoothing(), 0.0);
    EXPECT_GE(solver.timeAvgMGCResidual(), 0.0);
    EXPECT_GE(solver.timeAvgMGCDirectSolver(), 0.0);
}
