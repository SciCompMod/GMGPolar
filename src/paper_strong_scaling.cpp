#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <limits>
#include <memory>
#include <numeric>

#include "../include/GMGPolar/gmgpolar.h"
#include "../include/GMGPolar/test_cases.h"

void runTest(int maxOpenMPThreads, int divideBy2, std::ofstream& outfile)
{
    const double R0                           = 1e-8;
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;
    const double elongation_kappa             = 0.3;
    const double shift_delta                  = 0.2;

    /* --------------------------------------------- */
    /* Example 1: Polar Solution -> Higher Order 4.0 */

    // const double alpha_jump = 0.4837 * Rmax;
    // std::unique_ptr<DomainGeometry> domain_geometry =
    //     std::make_unique<ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    // std::unique_ptr<ExactSolution> exact_solution =
    //     std::make_unique<PolarR6_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    // std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniGyroCoefficients>(Rmax, alpha_jump);
    // std::unique_ptr<BoundaryConditions> boundary_conditions =
    //     std::make_unique<PolarR6_Boundary_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    // std::unique_ptr<SourceTerm> source_term =
    //     std::make_unique<PolarR6_ZoniGyro_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);

    /* --------------------------------------------------------- */
    /* Example 2: Cartesian Solution (Czarny) -> Lower Order 3.5 */

    // const double alpha_jump = 0.66 * Rmax;
    // std::unique_ptr<DomainGeometry> domain_geometry = std::make_unique<CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    // std::unique_ptr<ExactSolution> exact_solution = std::make_unique<CartesianR6_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    // std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<SonnendruckerGyroCoefficients>(Rmax, alpha_jump);
    // std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<CartesianR6_Boundary_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    // std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR6_SonnendruckerGyro_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);

    /* ---------------------------------------------- */
    /* Example 3: Refined Solution -> Lower Order 3.5 */

    // const double alpha_jump = 0.9 * Rmax; // Refinement where the solution is most complex
    // std::unique_ptr<DomainGeometry> domain_geometry = std::make_unique<ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    // std::unique_ptr<ExactSolution> exact_solution = std::make_unique<Refined_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    // std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniShiftedGyroCoefficients>(Rmax, alpha_jump);
    // std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<Refined_Boundary_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    // std::unique_ptr<SourceTerm> source_term = std::make_unique<Refined_ZoniShiftedGyro_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);

    /* --------------------------------------------------------- */
    /* Example 4a: Polar Solution (Czarny with Zoni Gyro Shifted) */

    const double alpha_jump = (0.7 + 0.025 * std::log(std::sqrt(2) - 1)) * Rmax;
    std::unique_ptr<DomainGeometry> domain_geometry =
        std::make_unique<CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<ExactSolution> exact_solution =
        std::make_unique<PolarR6_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniShiftedGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShiftedGyro_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);

    /* ------------------------------------------------------------- */
    /* Example 4b: Polar Solution (Shafranov with Zoni Gyro Shifted) */

    // const double alpha_jump = (0.7 + 0.025 * std::log(std::sqrt(2) - 1)) * Rmax;
    // std::unique_ptr<DomainGeometry> domain_geometry =
    //     std::make_unique<ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    // std::unique_ptr<ExactSolution> exact_solution =
    //     std::make_unique<PolarR6_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    // std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniShiftedGyroCoefficients>(Rmax, alpha_jump);
    // std::unique_ptr<BoundaryConditions> boundary_conditions =
    //     std::make_unique<PolarR6_Boundary_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    // std::unique_ptr<SourceTerm> source_term =
    //     std::make_unique<PolarR6_ZoniShiftedGyro_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);

    /* -------------------------------------------------------------- */
    /* Example 5a: Cartesian Solution (Czarny with Zoni Gyro Shifted) */

    // const double alpha_jump = (0.7 + 0.025 * std::log(std::sqrt(2) - 1)) * Rmax;
    // std::unique_ptr<DomainGeometry> domain_geometry = std::make_unique<CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    // std::unique_ptr<ExactSolution> exact_solution = std::make_unique<CartesianR6_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    // std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniShiftedGyroCoefficients>(Rmax, alpha_jump);
    // std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<CartesianR6_Boundary_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    // std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR6_ZoniShiftedGyro_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);

    /* ----------------------------------------------------------------- */
    /* Example 5b: Cartesian Solution (Shafranov with Zoni Gyro Shifted) */

    // const double alpha_jump = (0.7 + 0.025 * std::log(std::sqrt(2) - 1)) * Rmax;
    // std::unique_ptr<DomainGeometry> domain_geometry = std::make_unique<ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    // std::unique_ptr<ExactSolution> exact_solution = std::make_unique<CartesianR6_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    // std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniShiftedGyroCoefficients>(Rmax, alpha_jump);
    // std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<CartesianR6_Boundary_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    // std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR6_ZoniShiftedGyro_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);


    std::string geometry_string = "";
    if (typeid(*domain_geometry) == typeid(CircularGeometry)) {
        geometry_string = "Circular";
    }
    else if (typeid(*domain_geometry) == typeid(ShafranovGeometry)) {
        geometry_string = "Shafranov";
    }
    else if (typeid(*domain_geometry) == typeid(CzarnyGeometry)) {
        geometry_string = "Czarny";
    }
    else if (typeid(*domain_geometry) == typeid(CulhamGeometry)) {
        geometry_string = "Culham";
    }
    else {
        geometry_string = "Unknown";
    }

    GMGPolar solver(std::move(domain_geometry), std::move(coefficients), std::move(boundary_conditions),
                    std::move(source_term));

    solver.setSolution(std::move(exact_solution));

    const int verbose   = 1;
    const bool paraview = false;

    const double threadReductionFactor = 1.0;

    const StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_GIVE;
    const bool cacheDensityProfileCoefficients                = true;
    const bool cacheDomainGeometry                            = false;

    const int nr_exp             = 4;
    const int ntheta_exp         = -1;
    const int anisotropic_factor = 3;

    const bool DirBC_Interior = false;

    const bool FMG                     = false;
    const int FMG_iterations           = 3;
    const MultigridCycleType FMG_cycle = MultigridCycleType::F_CYCLE;

    const ExtrapolationType extrapolation   = ExtrapolationType::IMPLICIT_EXTRAPOLATION;
    const int maxLevels                     = -1;
    const int preSmoothingSteps             = 1;
    const int postSmoothingSteps            = 1;
    const MultigridCycleType multigridCycle = MultigridCycleType::V_CYCLE;

    const int maxIterations                 = 300;
    const ResidualNormType residualNormType = ResidualNormType::WEIGHTED_EUCLIDEAN;
    const double absoluteTolerance          = 1e-14;
    const double relativeTolerance          = 1e-8;

    solver.verbose(verbose);
    solver.paraview(paraview);

    solver.maxOpenMPThreads(maxOpenMPThreads);
    solver.threadReductionFactor(threadReductionFactor);

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

    // Perform setup and solve
    solver.setup();
    solver.solve();

    std::string stencil_string = "";
    if (solver.stencilDistributionMethod() == StencilDistributionMethod::CPU_TAKE) {
        stencil_string = "Take";
    }
    else if (solver.stencilDistributionMethod() == StencilDistributionMethod::CPU_GIVE) {
        stencil_string = "Give";
    }

    int extrapolation_int = static_cast<int>(solver.extrapolation());

    // Write results to file
    outfile << maxOpenMPThreads << "," << divideBy2 << "," << solver.grid().nr() << "," << solver.grid().ntheta() << ","
            << geometry_string << "," << stencil_string << "," << cacheDensityProfileCoefficients << ","
            << cacheDomainGeometry << "," << FMG << "," << extrapolation_int << ","
            << solver.t_setup_total + solver.t_solve_total - solver.t_setup_rhs << ","
            << solver.t_setup_total - solver.t_setup_rhs << "," << solver.t_setup_createLevels << ","
            << solver.t_setup_smoother << "," << solver.t_setup_directSolver << "," << solver.t_solve_total << ","
            << solver.t_solve_initial_approximation << "," << solver.t_solve_multigrid_iterations << ","
            << solver.t_check_convergence << "," << solver.t_check_exact_error << "," << solver.t_avg_MGC_total << ","
            << solver.t_avg_MGC_preSmoothing << "," << solver.t_avg_MGC_postSmoothing << ","
            << solver.t_avg_MGC_residual << "," << solver.t_avg_MGC_directSolver << std::endl;
}

int main()
{
    std::ofstream outfile("strong_scaling_results.csv");
    outfile << "Threads,DivideBy2,nr,ntheta,geometry,"
            << "stencil_method,cacheDensityProfileCoefficients,cacheDomainGeometry,FMG,extrapolation_int,"
            << "TotalTime,t_setup_total,t_setup_createLevels,"
            << "t_setup_smoother,t_setup_directSolver,t_solve_total,t_solve_initial_approximation,"
            << "t_solve_multigrid_iterations,t_check_convergence,t_check_exact_error,"
            << "t_avg_MGC_total,t_avg_MGC_preSmoothing,t_avg_MGC_postSmoothing,"
            << "t_avg_MGC_residual,t_avg_MGC_directSolver\n"; // Header

    // Define the constant parameters for the problem size
    int divideBy2 = 7; // Keeping the domain division constant

    // Vary only the number of threads for strong scaling
    std::vector<int> threadCounts = {1, 2, 4, 8, 16, 32, 64};  // Thread counts to test strong scaling

    for (size_t i = 0; i < threadCounts.size(); i++) {
        runTest(threadCounts[i], divideBy2, outfile); // Keep problem size fixed, vary threads
    }

    outfile.close();
    std::cout << "Results written to strong_scaling_results.csv" << std::endl;

    return 0;
}