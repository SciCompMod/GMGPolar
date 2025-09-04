#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>

#include "../../include/GMGPolar/gmgpolar.h"

TEST(GMGPolar_StencilDistribution, Take)
{
    /* PolarGrid */
    double R0                              = 1e-8;
    double Rmax                            = 1.3;
    int nr_exp                             = 4;
    int ntheta_exp                         = -1;
    double refinement_radius               = 0.66 * Rmax;
    int anisotropic_factor                 = 3;
    int divideBy2                          = 0;
    std::optional<double> splitting_radius = std::nullopt;
    PolarGrid grid(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, splitting_radius);

    /* Czarny, PolarR6 Sonnendrucker Gyro*/
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;
    std::unique_ptr<DomainGeometry> domain_geometry =
        std::make_unique<CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<ExactSolution> exact_solution =
        std::make_unique<PolarR6_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<SonnendruckerGyroCoefficients>(Rmax, refinement_radius);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_SonnendruckerGyro_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);

    /* GMGPolar settings */
    int verbose   = 0;
    bool paraview = false;

    int maxOpenMPThreads         = 1;
    double threadReductionFactor = 1.0;
    omp_set_num_threads(maxOpenMPThreads);

    bool DirBC_Interior                                 = true;
    StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_TAKE;
    bool cacheDensityProfileCoefficients                = true;
    bool cacheDomainGeometry                            = true;

    ExtrapolationType extrapolation   = ExtrapolationType::IMPLICIT_EXTRAPOLATION;
    int maxLevels                     = 3;
    int preSmoothingSteps             = 1;
    int postSmoothingSteps            = 1;
    MultigridCycleType multigridCycle = MultigridCycleType::V_CYCLE;
    bool FMG                          = true;
    int FMG_iterations                = 3;
    MultigridCycleType FMG_cycle      = MultigridCycleType::V_CYCLE;

    int maxIterations                 = 150;
    ResidualNormType residualNormType = ResidualNormType::WEIGHTED_EUCLIDEAN;
    double absoluteTolerance          = 1e-12;
    double relativeTolerance          = 1e-10;

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

    EXPECT_LE(number_of_iterations, 36);
    EXPECT_LE(exact_error_weighted_euclidean.value(), 2e-6);
    EXPECT_LE(exact_einfinity_error.value(), 2e-5);
}

TEST(GMGPolar_StencilDistribution, GiveCache)
{
    /* PolarGrid */
    double R0                              = 1e-8;
    double Rmax                            = 1.3;
    int nr_exp                             = 4;
    int ntheta_exp                         = -1;
    double refinement_radius               = 0.66 * Rmax;
    int anisotropic_factor                 = 3;
    int divideBy2                          = 0;
    std::optional<double> splitting_radius = std::nullopt;
    PolarGrid grid(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, splitting_radius);

    /* Czarny, PolarR6 Sonnendrucker Gyro*/
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;
    std::unique_ptr<DomainGeometry> domain_geometry =
        std::make_unique<CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<ExactSolution> exact_solution =
        std::make_unique<PolarR6_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<SonnendruckerGyroCoefficients>(Rmax, refinement_radius);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_SonnendruckerGyro_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);

    /* GMGPolar settings */
    int verbose   = 0;
    bool paraview = false;

    int maxOpenMPThreads         = 1;
    double threadReductionFactor = 1.0;
    omp_set_num_threads(maxOpenMPThreads);

    bool DirBC_Interior                                 = true;
    StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_GIVE;
    bool cacheDensityProfileCoefficients                = true;
    bool cacheDomainGeometry                            = true;

    ExtrapolationType extrapolation   = ExtrapolationType::IMPLICIT_EXTRAPOLATION;
    int maxLevels                     = 3;
    int preSmoothingSteps             = 1;
    int postSmoothingSteps            = 1;
    MultigridCycleType multigridCycle = MultigridCycleType::V_CYCLE;
    bool FMG                          = true;
    int FMG_iterations                = 3;
    MultigridCycleType FMG_cycle      = MultigridCycleType::V_CYCLE;

    int maxIterations                 = 150;
    ResidualNormType residualNormType = ResidualNormType::WEIGHTED_EUCLIDEAN;
    double absoluteTolerance          = 1e-12;
    double relativeTolerance          = 1e-10;

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

    EXPECT_LE(number_of_iterations, 36);
    EXPECT_LE(exact_error_weighted_euclidean.value(), 2e-6);
    EXPECT_LE(exact_einfinity_error.value(), 2e-5);
}

TEST(GMGPolar_StencilDistribution, GiveNoCache)
{
    /* PolarGrid */
    double R0                              = 1e-8;
    double Rmax                            = 1.3;
    int nr_exp                             = 4;
    int ntheta_exp                         = -1;
    double refinement_radius               = 0.66 * Rmax;
    int anisotropic_factor                 = 3;
    int divideBy2                          = 0;
    std::optional<double> splitting_radius = std::nullopt;
    PolarGrid grid(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, splitting_radius);

    /* Czarny, PolarR6 Sonnendrucker Gyro*/
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;
    std::unique_ptr<DomainGeometry> domain_geometry =
        std::make_unique<CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<ExactSolution> exact_solution =
        std::make_unique<PolarR6_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<SonnendruckerGyroCoefficients>(Rmax, refinement_radius);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_SonnendruckerGyro_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);

    /* GMGPolar settings */
    int verbose   = 0;
    bool paraview = false;

    int maxOpenMPThreads         = 1;
    double threadReductionFactor = 1.0;
    omp_set_num_threads(maxOpenMPThreads);

    bool DirBC_Interior                                 = true;
    StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_GIVE;
    bool cacheDensityProfileCoefficients                = false;
    bool cacheDomainGeometry                            = false;

    ExtrapolationType extrapolation   = ExtrapolationType::IMPLICIT_EXTRAPOLATION;
    int maxLevels                     = 3;
    int preSmoothingSteps             = 1;
    int postSmoothingSteps            = 1;
    MultigridCycleType multigridCycle = MultigridCycleType::V_CYCLE;
    bool FMG                          = true;
    int FMG_iterations                = 3;
    MultigridCycleType FMG_cycle      = MultigridCycleType::V_CYCLE;

    int maxIterations                 = 150;
    ResidualNormType residualNormType = ResidualNormType::WEIGHTED_EUCLIDEAN;
    double absoluteTolerance          = 1e-12;
    double relativeTolerance          = 1e-10;

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

    EXPECT_LE(number_of_iterations, 36);
    EXPECT_LE(exact_error_weighted_euclidean.value(), 2e-6);
    EXPECT_LE(exact_einfinity_error.value(), 2e-5);
}

// TEST(GMGPolar_MGC_Cycles, F_Cycle)
// {
//     /* PolarGrid */
//     double R0                              = 1e-8;
//     double Rmax                            = 1.3;
//     int nr_exp                             = 4;
//     int ntheta_exp                         = -1;
//     double refinement_radius               = 0.66 * Rmax;
//     int anisotropic_factor                 = 3;
//     int divideBy2                          = 1;
//     std::optional<double> splitting_radius = std::nullopt;
//     PolarGrid grid(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, splitting_radius);

//     /* PolarR6 */

//     /* Czarny, PolarR6 Sonnendrucker Gyro*/
//     const double inverse_aspect_ratio_epsilon = 0.3;
//     const double ellipticity_e                = 1.4;
//     std::unique_ptr<DomainGeometry> domain_geometry =
//         std::make_unique<CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
//     std::unique_ptr<ExactSolution> exact_solution =
//         std::make_unique<PolarR6_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
//     std::unique_ptr<DensityProfileCoefficients> coefficients =
//         std::make_unique<SonnendruckerGyroCoefficients>(Rmax, refinement_radius);
//     std::unique_ptr<BoundaryConditions> boundary_conditions =
//         std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
//     std::unique_ptr<SourceTerm> source_term =
//         std::make_unique<PolarR6_SonnendruckerGyro_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);

//     /* GMGPolar settings */
//     int verbose   = 0;
//     bool paraview = false;

//     int maxOpenMPThreads         = 1;
//     double threadReductionFactor = 1.0;
//     omp_set_num_threads(maxOpenMPThreads);

//     bool DirBC_Interior                                 = true;
//     StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_TAKE;
//     bool cacheDensityProfileCoefficients                = true;
//     bool cacheDomainGeometry                            = true;

//     ExtrapolationType extrapolation   = ExtrapolationType::IMPLICIT_EXTRAPOLATION;
//     int maxLevels                     = 3;
//     int preSmoothingSteps             = 1;
//     int postSmoothingSteps            = 1;
//     MultigridCycleType multigridCycle = MultigridCycleType::F_CYCLE;
//     bool FMG                          = true;
//     int FMG_iterations                = 3;
//     MultigridCycleType FMG_cycle      = MultigridCycleType::F_CYCLE;

//     int maxIterations                 = 150;
//     ResidualNormType residualNormType = ResidualNormType::EUCLIDEAN;
//     double absoluteTolerance          = 1e-10;
//     double relativeTolerance          = 1e-8;

//     // Create GMGPolar solver
//     GMGPolar solver(grid, *domain_geometry.get(), *coefficients.get());
//     solver.verbose(verbose);
//     solver.paraview(paraview);
//     solver.maxOpenMPThreads(maxOpenMPThreads);
//     solver.threadReductionFactor(threadReductionFactor);
//     solver.DirBC_Interior(DirBC_Interior);
//     solver.stencilDistributionMethod(stencilDistributionMethod);
//     solver.cacheDensityProfileCoefficients(cacheDensityProfileCoefficients);
//     solver.cacheDomainGeometry(cacheDomainGeometry);
//     solver.extrapolation(extrapolation);
//     solver.maxLevels(maxLevels);
//     solver.preSmoothingSteps(preSmoothingSteps);
//     solver.postSmoothingSteps(postSmoothingSteps);
//     solver.multigridCycle(multigridCycle);
//     solver.FMG(FMG);
//     solver.FMG_iterations(FMG_iterations);
//     solver.FMG_cycle(FMG_cycle);
//     solver.maxIterations(maxIterations);
//     solver.residualNormType(residualNormType);
//     solver.absoluteTolerance(absoluteTolerance);
//     solver.relativeTolerance(relativeTolerance);

//     solver.setup();
//     solver.setSolution(exact_solution.get());
//     solver.solve(*boundary_conditions.get(), *source_term.get());

//     solver.printTimings();

//     int number_of_iterations                             = solver.numberOfIterations();
//     std::optional<double> exact_error_weighted_euclidean = solver.exactErrorWeightedEuclidean();
//     std::optional<double> exact_einfinity_error          = solver.exactErrorInfinity();
//     ASSERT_TRUE(exact_error_weighted_euclidean.has_value());
//     ASSERT_TRUE(exact_einfinity_error.has_value());

//     EXPECT_LE(number_of_iterations, 34);
//     EXPECT_LE(exact_error_weighted_euclidean.value(), 2e-7);
//     EXPECT_LE(exact_einfinity_error.value(), 2e-6);
// }

// TEST(GMGPolar_MGC_Cycles, W_Cycle)
// {
//     /* PolarGrid */
//     double R0                              = 1e-8;
//     double Rmax                            = 1.3;
//     int nr_exp                             = 4;
//     int ntheta_exp                         = -1;
//     double refinement_radius               = 0.66 * Rmax;
//     int anisotropic_factor                 = 3;
//     int divideBy2                          = 1;
//     std::optional<double> splitting_radius = std::nullopt;
//     PolarGrid grid(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, splitting_radius);

//     /* PolarR6 */

//     /* Czarny, PolarR6 Sonnendrucker Gyro*/
//     const double inverse_aspect_ratio_epsilon = 0.3;
//     const double ellipticity_e                = 1.4;
//     std::unique_ptr<DomainGeometry> domain_geometry =
//         std::make_unique<CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
//     std::unique_ptr<ExactSolution> exact_solution =
//         std::make_unique<PolarR6_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
//     std::unique_ptr<DensityProfileCoefficients> coefficients =
//         std::make_unique<SonnendruckerGyroCoefficients>(Rmax, refinement_radius);
//     std::unique_ptr<BoundaryConditions> boundary_conditions =
//         std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
//     std::unique_ptr<SourceTerm> source_term =
//         std::make_unique<PolarR6_SonnendruckerGyro_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);

//     /* GMGPolar settings */
//     int verbose   = 0;
//     bool paraview = false;

//     int maxOpenMPThreads         = 1;
//     double threadReductionFactor = 1.0;
//     omp_set_num_threads(maxOpenMPThreads);

//     bool DirBC_Interior                                 = true;
//     StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_TAKE;
//     bool cacheDensityProfileCoefficients                = true;
//     bool cacheDomainGeometry                            = true;

//     ExtrapolationType extrapolation   = ExtrapolationType::IMPLICIT_EXTRAPOLATION;
//     int maxLevels                     = 3;
//     int preSmoothingSteps             = 1;
//     int postSmoothingSteps            = 1;
//     MultigridCycleType multigridCycle = MultigridCycleType::W_CYCLE;
//     bool FMG                          = true;
//     int FMG_iterations                = 3;
//     MultigridCycleType FMG_cycle      = MultigridCycleType::W_CYCLE;

//     int maxIterations                 = 150;
//     ResidualNormType residualNormType = ResidualNormType::EUCLIDEAN;
//     double absoluteTolerance          = 1e-10;
//     double relativeTolerance          = 1e-8;

//     // Create GMGPolar solver
//     GMGPolar solver(grid, *domain_geometry.get(), *coefficients.get());
//     solver.verbose(verbose);
//     solver.paraview(paraview);
//     solver.maxOpenMPThreads(maxOpenMPThreads);
//     solver.threadReductionFactor(threadReductionFactor);
//     solver.DirBC_Interior(DirBC_Interior);
//     solver.stencilDistributionMethod(stencilDistributionMethod);
//     solver.cacheDensityProfileCoefficients(cacheDensityProfileCoefficients);
//     solver.cacheDomainGeometry(cacheDomainGeometry);
//     solver.extrapolation(extrapolation);
//     solver.maxLevels(maxLevels);
//     solver.preSmoothingSteps(preSmoothingSteps);
//     solver.postSmoothingSteps(postSmoothingSteps);
//     solver.multigridCycle(multigridCycle);
//     solver.FMG(FMG);
//     solver.FMG_iterations(FMG_iterations);
//     solver.FMG_cycle(FMG_cycle);
//     solver.maxIterations(maxIterations);
//     solver.residualNormType(residualNormType);
//     solver.absoluteTolerance(absoluteTolerance);
//     solver.relativeTolerance(relativeTolerance);

//     solver.setup();
//     solver.setSolution(exact_solution.get());
//     solver.solve(*boundary_conditions.get(), *source_term.get());

//     solver.printTimings();

//     int number_of_iterations                             = solver.numberOfIterations();
//     std::optional<double> exact_error_weighted_euclidean = solver.exactErrorWeightedEuclidean();
//     std::optional<double> exact_einfinity_error          = solver.exactErrorInfinity();
//     ASSERT_TRUE(exact_error_weighted_euclidean.has_value());
//     ASSERT_TRUE(exact_einfinity_error.has_value());

//     EXPECT_LE(number_of_iterations, 34);
//     EXPECT_LE(exact_error_weighted_euclidean.value(), 2e-7);
//     EXPECT_LE(exact_einfinity_error.value(), 2e-6);
// }

// // #include <gtest/gtest.h>

// // #include <cmath>
// // #include <cstdlib>
// // #include <ctime>
// // #include <vector>

// // #include "../../include/GMGPolar/gmgpolar.h"

// // /* Check if Take and Give yield same results */

// // class GMGPolarParamsTest : public ::testing::TestWithParam<std::tuple<StencilDistributionMethod, bool>>
// // {
// // };

// // TEST_P(GMGPolarParamsTest, PolarR6_17x32)
// // {
// //     double R0 = 1e-8, Rmax = 2.0;
// //     int nr_exp = 4, ntheta_exp = 5;
// //     double refinement_radius = 0.5 * Rmax;
// //     int anisotropic_factor = 0, divideBy2 = 1;
// //     std::optional<double> splitting_radius = std::nullopt;
// //     PolarGrid grid(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, splitting_radius);

// //     const double kappa = 0.3, delta = 0.2;
// //     auto domain_geometry     = std::make_unique<ShafranovGeometry>(Rmax, kappa, delta);
// //     auto exact_solution      = std::make_unique<PolarR6_ShafranovGeometry>(Rmax, kappa, delta);
// //     auto coefficients        = std::make_unique<ZoniGyroCoefficients>(Rmax, refinement_radius);
// //     auto boundary_conditions = std::make_unique<PolarR6_Boundary_ShafranovGeometry>(Rmax, kappa, delta);
// //     auto source_term         = std::make_unique<PolarR6_ZoniGyro_ShafranovGeometry>(Rmax, kappa, delta);

// //     int verbose                  = 1;
// //     bool paraview                = false;
// //     int maxOpenMPThreads         = 4;
// //     double threadReductionFactor = 1.0;
// //     omp_set_num_threads(maxOpenMPThreads);

// //     bool DirBC_Interior                      = false;
// //     auto [stencilDistributionMethod, useFMG] = GetParam();
// //     bool cacheDensityProfileCoefficients = true, cacheDomainGeometry = true;

// //     ExtrapolationType extrapolation = ExtrapolationType::IMPLICIT_EXTRAPOLATION;
// //     int maxLevels = 7, preSmoothingSteps = 1, postSmoothingSteps = 1;
// //     MultigridCycleType multigridCycle = MultigridCycleType::V_CYCLE;

// //     int FMG_iterations           = useFMG ? 10 : 0;
// //     MultigridCycleType FMG_cycle = useFMG ? MultigridCycleType::F_CYCLE : MultigridCycleType::V_CYCLE;

// //     int maxIterations                 = 150;
// //     ResidualNormType residualNormType = ResidualNormType::EUCLIDEAN;
// //     double absoluteTolerance = 1e-8, relativeTolerance = 1e-8;

// //     GMGPolar solver(grid, *domain_geometry, *coefficients);
// //     solver.verbose(verbose);
// //     solver.paraview(paraview);
// //     solver.maxOpenMPThreads(maxOpenMPThreads);
// //     solver.threadReductionFactor(threadReductionFactor);
// //     solver.DirBC_Interior(DirBC_Interior);
// //     solver.stencilDistributionMethod(stencilDistributionMethod);
// //     solver.cacheDensityProfileCoefficients(cacheDensityProfileCoefficients);
// //     solver.cacheDomainGeometry(cacheDomainGeometry);
// //     solver.extrapolation(extrapolation);
// //     solver.maxLevels(maxLevels);
// //     solver.preSmoothingSteps(preSmoothingSteps);
// //     solver.postSmoothingSteps(postSmoothingSteps);
// //     solver.multigridCycle(multigridCycle);
// //     solver.FMG(useFMG);
// //     solver.FMG_iterations(FMG_iterations);
// //     solver.FMG_cycle(FMG_cycle);
// //     solver.maxIterations(maxIterations);
// //     solver.residualNormType(residualNormType);
// //     solver.absoluteTolerance(absoluteTolerance);
// //     solver.relativeTolerance(relativeTolerance);

// //     solver.setup();
// //     solver.setSolution(exact_solution.get());
// //     solver.solve(*boundary_conditions, *source_term);

// //     int iters    = solver.numberOfIterations();
// //     auto err_w   = solver.exactErrorWeightedEuclidean();
// //     auto err_inf = solver.exactErrorInfinity();

// //     ASSERT_TRUE(err_w.has_value());
// //     ASSERT_TRUE(err_inf.has_value());

// //     if (useFMG) {
// //         EXPECT_LE(iters, 12);
// //     }
// //     else {
// //         EXPECT_LE(iters, 20);
// //     }
// //     EXPECT_LE(*err_w, 1e-5);
// //     EXPECT_LE(*err_inf, 1e-5);
// // }

// // INSTANTIATE_TEST_SUITE_P(StencilAndFMG, GMGPolarParamsTest,
// //                          ::testing::Values(std::make_tuple(StencilDistributionMethod::CPU_TAKE, true),
// //                                            std::make_tuple(StencilDistributionMethod::CPU_GIVE, true),
// //                                            std::make_tuple(StencilDistributionMethod::CPU_TAKE, false),
// //                                            std::make_tuple(StencilDistributionMethod::CPU_GIVE, false)));
