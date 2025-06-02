#include <gtest/gtest.h>

#include "../../include/GMGPolar/gmgpolar.h"

TEST(GeometricMultigrid, V_Cycle_Extrapolation)
{
    const double R0                           = 1e-8;
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;
    const double elongation_kappa             = 0.3;
    const double shift_delta                  = 0.2;

    /* Example 1: Polar Solution -> Higher Order 4.0 */
    const double alpha_jump = 0.4837 * Rmax;
    std::unique_ptr<DomainGeometry> domain_geometry =
        std::make_unique<ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<ExactSolution> exact_solution =
        std::make_unique<PolarR6_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniGyro_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);

    GMGPolar solver(std::move(domain_geometry), std::move(coefficients), std::move(boundary_conditions),
                    std::move(source_term));

    solver.setSolution(std::move(exact_solution));

    const int verbose   = 0;
    const bool paraview = false;

    const int maxOpenMPThreads = 16;

    const StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_GIVE;
    const bool cacheDensityProfileCoefficients                = true;
    const bool cacheDomainGeometry                            = false;

    const int nr_exp             = 3;
    const int ntheta_exp         = 5;
    const int anisotropic_factor = 2;
    int divideBy2                = 0;

    const bool DirBC_Interior = false;

    const bool FMG                     = true;
    const int FMG_iterations           = 3;
    const MultigridCycleType FMG_cycle = MultigridCycleType::V_CYCLE;

    const ExtrapolationType extrapolation   = ExtrapolationType::IMPLICIT_EXTRAPOLATION;
    const int maxLevels                     = 2;
    const int preSmoothingSteps             = 1;
    const int postSmoothingSteps            = 1;
    const MultigridCycleType multigridCycle = MultigridCycleType::V_CYCLE;

    const int maxIterations                 = 150;
    const ResidualNormType residualNormType = ResidualNormType::EUCLIDEAN;
    const double absoluteTolerance          = 1e-10;
    const double relativeTolerance          = 1e-8;

    solver.verbose(verbose);
    solver.paraview(paraview);

    solver.maxOpenMPThreads(maxOpenMPThreads);

    solver.stencilDistributionMethod(stencilDistributionMethod);
    solver.cacheDensityProfileCoefficients(cacheDensityProfileCoefficients);
    solver.cacheDomainGeometry(cacheDomainGeometry);

    solver.R0(R0);
    solver.Rmax(Rmax);
    solver.nr_exp(nr_exp);
    solver.ntheta_exp(ntheta_exp);
    solver.anisotropic_factor(anisotropic_factor);
    solver.divideBy2(divideBy2);

    solver.DirBC_Interior(DirBC_Interior);

    solver.FMG(FMG);
    solver.FMG_iterations(FMG_iterations);
    solver.FMG_cycle(FMG_cycle);

    solver.extrapolation(extrapolation);
    solver.maxLevels(maxLevels);
    solver.preSmoothingSteps(preSmoothingSteps);
    solver.postSmoothingSteps(postSmoothingSteps);
    solver.multigridCycle(multigridCycle);

    solver.maxIterations(maxIterations);
    solver.residualNormType(residualNormType);
    solver.absoluteTolerance(absoluteTolerance);
    solver.relativeTolerance(relativeTolerance);

    solver.setup();
    solver.solve();

    ASSERT_LE(solver.numberOfIterations(), 45);
}

TEST(GeometricMultigrid, V_Cycle)
{
    const double R0                           = 1e-8;
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;
    const double elongation_kappa             = 0.3;
    const double shift_delta                  = 0.2;

    /* Example 1: Polar Solution -> Higher Order 4.0 */
    const double alpha_jump = 0.4837 * Rmax;
    std::unique_ptr<DomainGeometry> domain_geometry =
        std::make_unique<ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<ExactSolution> exact_solution =
        std::make_unique<PolarR6_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniGyro_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);


    GMGPolar solver(std::move(domain_geometry), std::move(coefficients), std::move(boundary_conditions),
                    std::move(source_term));

    solver.setSolution(std::move(exact_solution));

    const int verbose   = 0;
    const bool paraview = false;

    const int maxOpenMPThreads = 16;

    const StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_TAKE;
    const bool cacheDensityProfileCoefficients                = true;
    const bool cacheDomainGeometry                            = true;

    const int nr_exp             = 3;
    const int ntheta_exp         = 5;
    const int anisotropic_factor = 2;
    int divideBy2                = 0;
    const bool DirBC_Interior = false;

    const bool FMG                     = true;
    const int FMG_iterations           = 3;
    const MultigridCycleType FMG_cycle = MultigridCycleType::V_CYCLE;

    const ExtrapolationType extrapolation   = ExtrapolationType::NONE;
    const int maxLevels                     = 2;
    const int preSmoothingSteps             = 1;
    const int postSmoothingSteps            = 1;
    const MultigridCycleType multigridCycle = MultigridCycleType::V_CYCLE;

    const int maxIterations                 = 150;
    const ResidualNormType residualNormType = ResidualNormType::WEIGHTED_EUCLIDEAN;
    const double absoluteTolerance          = 1e-12;
    const double relativeTolerance          = 1e-8;

    solver.verbose(verbose);
    solver.paraview(paraview);

    solver.maxOpenMPThreads(maxOpenMPThreads);

    solver.stencilDistributionMethod(stencilDistributionMethod);
    solver.cacheDensityProfileCoefficients(cacheDensityProfileCoefficients);
    solver.cacheDomainGeometry(cacheDomainGeometry);

    solver.R0(R0);
    solver.Rmax(Rmax);
    solver.nr_exp(nr_exp);
    solver.ntheta_exp(ntheta_exp);
    solver.anisotropic_factor(anisotropic_factor);
    solver.divideBy2(divideBy2);

    solver.DirBC_Interior(DirBC_Interior);

    solver.FMG(FMG);
    solver.FMG_iterations(FMG_iterations);
    solver.FMG_cycle(FMG_cycle);

    solver.extrapolation(extrapolation);
    solver.maxLevels(maxLevels);
    solver.preSmoothingSteps(preSmoothingSteps);
    solver.postSmoothingSteps(postSmoothingSteps);
    solver.multigridCycle(multigridCycle);

    solver.maxIterations(maxIterations);
    solver.residualNormType(residualNormType);
    solver.absoluteTolerance(absoluteTolerance);
    solver.relativeTolerance(relativeTolerance);

    solver.setup();
    solver.solve();

    ASSERT_LE(solver.numberOfIterations(), 10);
}

TEST(GeometricMultigrid, W_Cycle_Extrapolation)
{
    const double R0                           = 1e-8;
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;
    const double elongation_kappa             = 0.3;
    const double shift_delta                  = 0.2;

    /* Example 1: Polar Solution -> Higher Order 4.0 */
    const double alpha_jump = 0.4837 * Rmax;
    std::unique_ptr<DomainGeometry> domain_geometry =
        std::make_unique<ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<ExactSolution> exact_solution =
        std::make_unique<PolarR6_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniGyro_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);

    GMGPolar solver(std::move(domain_geometry), std::move(coefficients), std::move(boundary_conditions),
                    std::move(source_term));

    solver.setSolution(std::move(exact_solution));

    const int verbose   = 0;
    const bool paraview = false;

    const int maxOpenMPThreads = 16;

    const StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_GIVE;
    const bool cacheDensityProfileCoefficients                = true;
    const bool cacheDomainGeometry                            = false;

    const int nr_exp             = 3;
    const int ntheta_exp         = 5;
    const int anisotropic_factor = 2;
    int divideBy2                = 0;

    const bool DirBC_Interior = false;

    const bool FMG                     = true;
    const int FMG_iterations           = 3;
    const MultigridCycleType FMG_cycle = MultigridCycleType::W_CYCLE;

    const ExtrapolationType extrapolation   = ExtrapolationType::IMPLICIT_EXTRAPOLATION;
    const int maxLevels                     = 2;
    const int preSmoothingSteps             = 1;
    const int postSmoothingSteps            = 1;
    const MultigridCycleType multigridCycle = MultigridCycleType::W_CYCLE;

    const int maxIterations                 = 150;
    const ResidualNormType residualNormType = ResidualNormType::EUCLIDEAN;
    const double absoluteTolerance          = 1e-10;
    const double relativeTolerance          = 1e-8;

    solver.verbose(verbose);
    solver.paraview(paraview);

    solver.maxOpenMPThreads(maxOpenMPThreads);

    solver.stencilDistributionMethod(stencilDistributionMethod);
    solver.cacheDensityProfileCoefficients(cacheDensityProfileCoefficients);
    solver.cacheDomainGeometry(cacheDomainGeometry);

    solver.R0(R0);
    solver.Rmax(Rmax);
    solver.nr_exp(nr_exp);
    solver.ntheta_exp(ntheta_exp);
    solver.anisotropic_factor(anisotropic_factor);
    solver.divideBy2(divideBy2);

    solver.DirBC_Interior(DirBC_Interior);

    solver.FMG(FMG);
    solver.FMG_iterations(FMG_iterations);
    solver.FMG_cycle(FMG_cycle);

    solver.extrapolation(extrapolation);
    solver.maxLevels(maxLevels);
    solver.preSmoothingSteps(preSmoothingSteps);
    solver.postSmoothingSteps(postSmoothingSteps);
    solver.multigridCycle(multigridCycle);

    solver.maxIterations(maxIterations);
    solver.residualNormType(residualNormType);
    solver.absoluteTolerance(absoluteTolerance);
    solver.relativeTolerance(relativeTolerance);

    solver.setup();
    solver.solve();

    ASSERT_LE(solver.numberOfIterations(), 45);
}

TEST(GeometricMultigrid, W_Cycle)
{
    const double R0                           = 1e-8;
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;
    const double elongation_kappa             = 0.3;
    const double shift_delta                  = 0.2;

    /* Example 1: Polar Solution -> Higher Order 4.0 */
    const double alpha_jump = 0.4837 * Rmax;
    std::unique_ptr<DomainGeometry> domain_geometry =
        std::make_unique<ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<ExactSolution> exact_solution =
        std::make_unique<PolarR6_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniGyro_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);

    GMGPolar solver(std::move(domain_geometry), std::move(coefficients), std::move(boundary_conditions),
                    std::move(source_term));

    solver.setSolution(std::move(exact_solution));

    const int verbose   = 0;
    const bool paraview = false;

    const int maxOpenMPThreads = 16;

    const StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_GIVE;
    const bool cacheDensityProfileCoefficients                = false;
    const bool cacheDomainGeometry                            = false;

    const int nr_exp             = 3;
    const int ntheta_exp         = 5;
    const int anisotropic_factor = 2;
    int divideBy2                = 0;

    const bool DirBC_Interior = true;

    const bool FMG                     = false;
    const int FMG_iterations           = 3;
    const MultigridCycleType FMG_cycle = MultigridCycleType::W_CYCLE;

    const ExtrapolationType extrapolation   = ExtrapolationType::NONE;
    const int maxLevels                     = 2;
    const int preSmoothingSteps             = 1;
    const int postSmoothingSteps            = 1;
    const MultigridCycleType multigridCycle = MultigridCycleType::W_CYCLE;

    const int maxIterations                 = 150;
    const ResidualNormType residualNormType = ResidualNormType::INFINITY_NORM;
    const double absoluteTolerance          = 1e-12;
    const double relativeTolerance          = 1e-8;

    solver.verbose(verbose);
    solver.paraview(paraview);

    solver.maxOpenMPThreads(maxOpenMPThreads);

    solver.stencilDistributionMethod(stencilDistributionMethod);
    solver.cacheDensityProfileCoefficients(cacheDensityProfileCoefficients);
    solver.cacheDomainGeometry(cacheDomainGeometry);

    solver.R0(R0);
    solver.Rmax(Rmax);
    solver.nr_exp(nr_exp);
    solver.ntheta_exp(ntheta_exp);
    solver.anisotropic_factor(anisotropic_factor);
    solver.divideBy2(divideBy2);

    solver.DirBC_Interior(DirBC_Interior);

    solver.FMG(FMG);
    solver.FMG_iterations(FMG_iterations);
    solver.FMG_cycle(FMG_cycle);

    solver.extrapolation(extrapolation);
    solver.maxLevels(maxLevels);
    solver.preSmoothingSteps(preSmoothingSteps);
    solver.postSmoothingSteps(postSmoothingSteps);
    solver.multigridCycle(multigridCycle);

    solver.maxIterations(maxIterations);
    solver.residualNormType(residualNormType);
    solver.absoluteTolerance(absoluteTolerance);
    solver.relativeTolerance(relativeTolerance);

    solver.setup();
    solver.solve();

    ASSERT_LE(solver.numberOfIterations(), 10);
}

TEST(GeometricMultigrid, F_Cycle_Extrapolation)
{
    const double R0                           = 1e-8;
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;
    const double elongation_kappa             = 0.3;
    const double shift_delta                  = 0.2;

    /* Example 1: Polar Solution -> Higher Order 4.0 */
    const double alpha_jump = 0.4837 * Rmax;
    std::unique_ptr<DomainGeometry> domain_geometry =
        std::make_unique<ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<ExactSolution> exact_solution =
        std::make_unique<PolarR6_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniGyro_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);

    GMGPolar solver(std::move(domain_geometry), std::move(coefficients), std::move(boundary_conditions),
                    std::move(source_term));

    solver.setSolution(std::move(exact_solution));

    const int verbose   = 1;
    const bool paraview = false;

    const int maxOpenMPThreads = 16;

    const StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_GIVE;
    const bool cacheDensityProfileCoefficients                = true;
    const bool cacheDomainGeometry                            = false;

    const int nr_exp             = 3;
    const int ntheta_exp         = 5;
    const int anisotropic_factor = 2;
    int divideBy2                = 0;

    const bool DirBC_Interior = false;

    const bool FMG                     = true;
    const int FMG_iterations           = 3;
    const MultigridCycleType FMG_cycle = MultigridCycleType::F_CYCLE;

    const ExtrapolationType extrapolation   = ExtrapolationType::IMPLICIT_EXTRAPOLATION;
    const int maxLevels                     = 2;
    const int preSmoothingSteps             = 1;
    const int postSmoothingSteps            = 1;
    const MultigridCycleType multigridCycle = MultigridCycleType::F_CYCLE;

    const int maxIterations                 = 150;
    const ResidualNormType residualNormType = ResidualNormType::EUCLIDEAN;
    const double absoluteTolerance          = 1e-10;
    const double relativeTolerance          = 1e-8;

    solver.verbose(verbose);
    solver.paraview(paraview);

    solver.maxOpenMPThreads(maxOpenMPThreads);

    solver.stencilDistributionMethod(stencilDistributionMethod);
    solver.cacheDensityProfileCoefficients(cacheDensityProfileCoefficients);
    solver.cacheDomainGeometry(cacheDomainGeometry);

    solver.R0(R0);
    solver.Rmax(Rmax);
    solver.nr_exp(nr_exp);
    solver.ntheta_exp(ntheta_exp);
    solver.anisotropic_factor(anisotropic_factor);
    solver.divideBy2(divideBy2);

    solver.DirBC_Interior(DirBC_Interior);

    solver.FMG(FMG);
    solver.FMG_iterations(FMG_iterations);
    solver.FMG_cycle(FMG_cycle);

    solver.extrapolation(extrapolation);
    solver.maxLevels(maxLevels);
    solver.preSmoothingSteps(preSmoothingSteps);
    solver.postSmoothingSteps(postSmoothingSteps);
    solver.multigridCycle(multigridCycle);

    solver.maxIterations(maxIterations);
    solver.residualNormType(residualNormType);
    solver.absoluteTolerance(absoluteTolerance);
    solver.relativeTolerance(relativeTolerance);

    solver.setup();
    solver.solve();

    ASSERT_LE(solver.numberOfIterations(), 45);
}

TEST(GeometricMultigrid, F_Cycle)
{
    const double R0                           = 1e-8;
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;
    const double elongation_kappa             = 0.3;
    const double shift_delta                  = 0.2;

    /* Example 1: Polar Solution -> Higher Order 4.0 */
    const double alpha_jump = 0.4837 * Rmax;
    std::unique_ptr<DomainGeometry> domain_geometry =
        std::make_unique<ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<ExactSolution> exact_solution =
        std::make_unique<PolarR6_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniGyro_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);

    GMGPolar solver(std::move(domain_geometry), std::move(coefficients), std::move(boundary_conditions),
                    std::move(source_term));

    solver.setSolution(std::move(exact_solution));

    const int verbose   = 1;
    const bool paraview = false;

    const int maxOpenMPThreads = 16;

    const StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_GIVE;
    const bool cacheDensityProfileCoefficients                = false;
    const bool cacheDomainGeometry                            = false;

    const int nr_exp             = 3;
    const int ntheta_exp         = 5;
    const int anisotropic_factor = 2;
    int divideBy2                = 0;

    const bool DirBC_Interior = true;

    const bool FMG                     = false;
    const int FMG_iterations           = 3;
    const MultigridCycleType FMG_cycle = MultigridCycleType::F_CYCLE;

    const ExtrapolationType extrapolation   = ExtrapolationType::NONE;
    const int maxLevels                     = 2;
    const int preSmoothingSteps             = 1;
    const int postSmoothingSteps            = 1;
    const MultigridCycleType multigridCycle = MultigridCycleType::F_CYCLE;

    const int maxIterations                 = 150;
    const ResidualNormType residualNormType = ResidualNormType::INFINITY_NORM;
    const double absoluteTolerance          = 1e-12;
    const double relativeTolerance          = 1e-8;

    solver.verbose(verbose);
    solver.paraview(paraview);

    solver.maxOpenMPThreads(maxOpenMPThreads);

    solver.stencilDistributionMethod(stencilDistributionMethod);
    solver.cacheDensityProfileCoefficients(cacheDensityProfileCoefficients);
    solver.cacheDomainGeometry(cacheDomainGeometry);

    solver.R0(R0);
    solver.Rmax(Rmax);
    solver.nr_exp(nr_exp);
    solver.ntheta_exp(ntheta_exp);
    solver.anisotropic_factor(anisotropic_factor);
    solver.divideBy2(divideBy2);

    solver.DirBC_Interior(DirBC_Interior);

    solver.FMG(FMG);
    solver.FMG_iterations(FMG_iterations);
    solver.FMG_cycle(FMG_cycle);

    solver.extrapolation(extrapolation);
    solver.maxLevels(maxLevels);
    solver.preSmoothingSteps(preSmoothingSteps);
    solver.postSmoothingSteps(postSmoothingSteps);
    solver.multigridCycle(multigridCycle);

    solver.maxIterations(maxIterations);
    solver.residualNormType(residualNormType);
    solver.absoluteTolerance(absoluteTolerance);
    solver.relativeTolerance(relativeTolerance);

    solver.setup();
    solver.solve();

    ASSERT_LE(solver.numberOfIterations(), 10);
}
