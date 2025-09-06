#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>

#include "../../include/GMGPolar/gmgpolar.h"

TEST(GMGPolar_IterationCount, CartesianR2_NoFMG_17x32)
{
    /* PolarGrid */
    double R0                              = 1e-8;
    double Rmax                            = 2.0;
    int nr_exp                             = 4;
    int ntheta_exp                         = 5;
    double refinement_radius               = 0.5 * Rmax;
    int anisotropic_factor                 = 0;
    int divideBy2                          = 0;
    std::optional<double> splitting_radius = std::nullopt;
    PolarGrid grid(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, splitting_radius);

    /* CartesianR2 */
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;
    std::unique_ptr<DomainGeometry> domain_geometry =
        std::make_unique<CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<ExactSolution> exact_solution =
        std::make_unique<CartesianR2_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<SonnendruckerGyroCoefficients>(Rmax, refinement_radius);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR2_Boundary_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CzarnyGeometry>(
        Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);

    /* GMGPolar settings */
    int verbose   = 0;
    bool paraview = false;

    int maxOpenMPThreads         = 4;
    double threadReductionFactor = 1.0;
    omp_set_num_threads(maxOpenMPThreads);

    bool DirBC_Interior                                 = false;
    StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_TAKE;
    bool cacheDensityProfileCoefficients                = true;
    bool cacheDomainGeometry                            = true;

    ExtrapolationType extrapolation   = ExtrapolationType::IMPLICIT_EXTRAPOLATION;
    int maxLevels                     = 7;
    int preSmoothingSteps             = 1;
    int postSmoothingSteps            = 1;
    MultigridCycleType multigridCycle = MultigridCycleType::V_CYCLE;
    bool FMG                          = false;
    int FMG_iterations                = 3;
    MultigridCycleType FMG_cycle      = MultigridCycleType::F_CYCLE;

    int maxIterations                 = 150;
    ResidualNormType residualNormType = ResidualNormType::EUCLIDEAN;
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

    int number_of_iterations                             = solver.numberOfIterations();
    std::optional<double> exact_error_weighted_euclidean = solver.exactErrorWeightedEuclidean();
    std::optional<double> exact_einfinity_error          = solver.exactErrorInfinity();
    ASSERT_TRUE(exact_error_weighted_euclidean.has_value());
    ASSERT_TRUE(exact_einfinity_error.has_value());

    EXPECT_LE(number_of_iterations, 58);
    EXPECT_LE(exact_error_weighted_euclidean.value(), 3e-2);
    EXPECT_LE(exact_einfinity_error.value(), 2e-1);
}

TEST(GMGPolar_IterationCount, CartesianR2_NoFMG_33x64)
{
    /* PolarGrid */
    double R0                              = 1e-8;
    double Rmax                            = 2.0;
    int nr_exp                             = 4;
    int ntheta_exp                         = 5;
    double refinement_radius               = 0.5 * Rmax;
    int anisotropic_factor                 = 0;
    int divideBy2                          = 1;
    std::optional<double> splitting_radius = std::nullopt;
    PolarGrid grid(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, splitting_radius);

    /* CartesianR2 */
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;
    std::unique_ptr<DomainGeometry> domain_geometry =
        std::make_unique<CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<ExactSolution> exact_solution =
        std::make_unique<CartesianR2_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<SonnendruckerGyroCoefficients>(Rmax, refinement_radius);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR2_Boundary_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CzarnyGeometry>(
        Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);

    /* GMGPolar settings */
    int verbose   = 0;
    bool paraview = false;

    int maxOpenMPThreads         = 4;
    double threadReductionFactor = 1.0;
    omp_set_num_threads(maxOpenMPThreads);

    bool DirBC_Interior                                 = false;
    StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_TAKE;
    bool cacheDensityProfileCoefficients                = true;
    bool cacheDomainGeometry                            = true;

    ExtrapolationType extrapolation   = ExtrapolationType::IMPLICIT_EXTRAPOLATION;
    int maxLevels                     = 7;
    int preSmoothingSteps             = 1;
    int postSmoothingSteps            = 1;
    MultigridCycleType multigridCycle = MultigridCycleType::V_CYCLE;
    bool FMG                          = false;
    int FMG_iterations                = 3;
    MultigridCycleType FMG_cycle      = MultigridCycleType::F_CYCLE;

    int maxIterations                 = 150;
    ResidualNormType residualNormType = ResidualNormType::EUCLIDEAN;
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

    int number_of_iterations                             = solver.numberOfIterations();
    std::optional<double> exact_error_weighted_euclidean = solver.exactErrorWeightedEuclidean();
    std::optional<double> exact_einfinity_error          = solver.exactErrorInfinity();
    ASSERT_TRUE(exact_error_weighted_euclidean.has_value());
    ASSERT_TRUE(exact_einfinity_error.has_value());

    EXPECT_LE(number_of_iterations, 53);
    EXPECT_LE(exact_error_weighted_euclidean.value(), 2e-3);
    EXPECT_LE(exact_einfinity_error.value(), 2e-2);
}

TEST(GMGPolar_IterationCount, CartesianR2_FMG_17x32)
{
    /* PolarGrid */
    double R0                              = 1e-8;
    double Rmax                            = 2.0;
    int nr_exp                             = 4;
    int ntheta_exp                         = 5;
    double refinement_radius               = 0.5 * Rmax;
    int anisotropic_factor                 = 0;
    int divideBy2                          = 0;
    std::optional<double> splitting_radius = std::nullopt;
    PolarGrid grid(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, splitting_radius);

    /* CartesianR2 */
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;
    std::unique_ptr<DomainGeometry> domain_geometry =
        std::make_unique<CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<ExactSolution> exact_solution =
        std::make_unique<CartesianR2_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<SonnendruckerGyroCoefficients>(Rmax, refinement_radius);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR2_Boundary_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CzarnyGeometry>(
        Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);

    /* GMGPolar settings */
    int verbose   = 0;
    bool paraview = false;

    int maxOpenMPThreads         = 4;
    double threadReductionFactor = 1.0;
    omp_set_num_threads(maxOpenMPThreads);

    bool DirBC_Interior                                 = false;
    StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_TAKE;
    bool cacheDensityProfileCoefficients                = true;
    bool cacheDomainGeometry                            = true;

    ExtrapolationType extrapolation   = ExtrapolationType::IMPLICIT_EXTRAPOLATION;
    int maxLevels                     = 7;
    int preSmoothingSteps             = 1;
    int postSmoothingSteps            = 1;
    MultigridCycleType multigridCycle = MultigridCycleType::V_CYCLE;
    bool FMG                          = true;
    int FMG_iterations                = 10;
    MultigridCycleType FMG_cycle      = MultigridCycleType::F_CYCLE;

    int maxIterations                 = 150;
    ResidualNormType residualNormType = ResidualNormType::EUCLIDEAN;
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

    int number_of_iterations                             = solver.numberOfIterations();
    std::optional<double> exact_error_weighted_euclidean = solver.exactErrorWeightedEuclidean();
    std::optional<double> exact_einfinity_error          = solver.exactErrorInfinity();
    ASSERT_TRUE(exact_error_weighted_euclidean.has_value());
    ASSERT_TRUE(exact_einfinity_error.has_value());

    EXPECT_LE(number_of_iterations, 58);
    EXPECT_LE(exact_error_weighted_euclidean.value(), 3e-2);
    EXPECT_LE(exact_einfinity_error.value(), 2e-1);
}

TEST(GMGPolar_IterationCount, CartesianR2_FMG_33x64)
{
    /* PolarGrid */
    double R0                              = 1e-8;
    double Rmax                            = 2.0;
    int nr_exp                             = 4;
    int ntheta_exp                         = 5;
    double refinement_radius               = 0.5 * Rmax;
    int anisotropic_factor                 = 0;
    int divideBy2                          = 1;
    std::optional<double> splitting_radius = std::nullopt;
    PolarGrid grid(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, splitting_radius);

    /* CartesianR2 */
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;
    std::unique_ptr<DomainGeometry> domain_geometry =
        std::make_unique<CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<ExactSolution> exact_solution =
        std::make_unique<CartesianR2_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<SonnendruckerGyroCoefficients>(Rmax, refinement_radius);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<CartesianR2_Boundary_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CzarnyGeometry>(
        Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);

    /* GMGPolar settings */
    int verbose   = 0;
    bool paraview = false;

    int maxOpenMPThreads         = 4;
    double threadReductionFactor = 1.0;
    omp_set_num_threads(maxOpenMPThreads);

    bool DirBC_Interior                                 = false;
    StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_TAKE;
    bool cacheDensityProfileCoefficients                = true;
    bool cacheDomainGeometry                            = true;

    ExtrapolationType extrapolation   = ExtrapolationType::IMPLICIT_EXTRAPOLATION;
    int maxLevels                     = 7;
    int preSmoothingSteps             = 1;
    int postSmoothingSteps            = 1;
    MultigridCycleType multigridCycle = MultigridCycleType::V_CYCLE;
    bool FMG                          = true;
    int FMG_iterations                = 10;
    MultigridCycleType FMG_cycle      = MultigridCycleType::F_CYCLE;

    int maxIterations                 = 150;
    ResidualNormType residualNormType = ResidualNormType::EUCLIDEAN;
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

    int number_of_iterations                             = solver.numberOfIterations();
    std::optional<double> exact_error_weighted_euclidean = solver.exactErrorWeightedEuclidean();
    std::optional<double> exact_einfinity_error          = solver.exactErrorInfinity();
    ASSERT_TRUE(exact_error_weighted_euclidean.has_value());
    ASSERT_TRUE(exact_einfinity_error.has_value());

    EXPECT_LE(number_of_iterations, 43);
    EXPECT_LE(exact_error_weighted_euclidean.value(), 2e-3);
    EXPECT_LE(exact_einfinity_error.value(), 2e-2);
}
