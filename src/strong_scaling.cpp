#include <iostream>
#include <fstream>
#include <cstdlib>

#include "../include/GMGPolar/gmgpolar.h"
#include "../include/GMGPolar/test_cases.h"
using namespace gmgpolar;

void runTest(int divideBy2, std::ostream& outfile)
{
    const double R0                           = 1e-8;
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;
    const double elongation_kappa             = 0.3;
    const double shift_delta                  = 0.2;

    /* Example 1: Polar Solution -> Higher Order 4.0 */
    // const double alpha_jump = 0.4837 * Rmax;
    // ShafranovGeometry domain_geometry(Rmax, elongation_kappa, shift_delta);
    // ZoniGyroCoefficients coefficients(Rmax, alpha_jump);
    // PolarR6_Boundary_ShafranovGeometry boundary_conditions(Rmax, elongation_kappa, shift_delta);
    // PolarR6_ShafranovGeometry exact_solution(Rmax, elongation_kappa, shift_delta);

    /* Example 2: Cartesian Solution -> Lower Order 3.5 */
    const double alpha_jump = 0.66 * Rmax;
    CzarnyGeometry domain_geometry(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    SonnendruckerGyroCoefficients coefficients(Rmax, alpha_jump);
    CartesianR2_Boundary_CzarnyGeometry boundary_conditions(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    CartesianR2_CzarnyGeometry exact_solution(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);

    /* Example 3: Refined Solution -> Lower Order 3.5 */
    // const double alpha_jump = 0.9 * Rmax; // Refinement where the solution is most complex
    // ShafranovGeometry domain_geometry(Rmax, elongation_kappa, shift_delta);
    // ZoniShiftedGyroCoefficients coefficients(Rmax, alpha_jump);
    // Refined_Boundary_ShafranovGeometry boundary_conditions(Rmax, elongation_kappa, shift_delta);
    // Refined_ShafranovGeometry exact_solution(Rmax, elongation_kappa, shift_delta);

    std::string geometry_string = "";
    if (typeid(domain_geometry) == typeid(CircularGeometry)) {
        geometry_string = "Circular";
    }
    else if (typeid(domain_geometry) == typeid(ShafranovGeometry)) {
        geometry_string = "Shafranov";
    }
    else if (typeid(domain_geometry) == typeid(CzarnyGeometry)) {
        geometry_string = "Czarny";
    }
    else if (typeid(domain_geometry) == typeid(CulhamGeometry)) {
        geometry_string = "Culham";
    }
    else {
        geometry_string = "Unknown";
    }

    const int verbose   = 1;
    const bool paraview = false;

    const StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_GIVE;
    const bool cacheDensityProfileCoefficients                = true;
    const bool cacheDomainGeometry                            = false;

    const int nr_exp             = 4;
    const int ntheta_exp         = 6;
    const int anisotropic_factor = 3;

    const bool DirBC_Interior = false;

    const bool FMG                     = true;
    const int FMG_iterations           = 3;
    const MultigridCycleType FMG_cycle = MultigridCycleType::F_CYCLE;

    const ExtrapolationType extrapolation   = ExtrapolationType::NONE;
    const int maxLevels                     = 7;
    const int preSmoothingSteps             = 1;
    const int postSmoothingSteps            = 1;
    const MultigridCycleType multigridCycle = MultigridCycleType::F_CYCLE;

    const int maxIterations                 = 10;
    const ResidualNormType residualNormType = ResidualNormType::EUCLIDEAN;
    const double absoluteTolerance          = 1e-50;
    const double relativeTolerance          = 1e-50;

    double refinement_radius               = alpha_jump;
    std::optional<double> splitting_radius = std::nullopt;
    PolarGrid<Kokkos::HostSpace> grid(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2,
                                      splitting_radius);
    GMGPolar solver(grid, domain_geometry, coefficients);

    // PolarR6_ZoniGyro_ShafranovGeometry source_term(grid, Rmax, elongation_kappa, shift_delta);
    CartesianR2_SonnendruckerGyro_CzarnyGeometry source_term(grid, Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    // Refined_ZoniShiftedGyro_ShafranovGeometry source_term(grid, Rmax, elongation_kappa, shift_delta);

    solver.verbose(verbose);
    solver.paraview(paraview);

    solver.DirBC_Interior(DirBC_Interior);
    solver.stencilDistributionMethod(stencilDistributionMethod);
    solver.cacheDensityProfileCoefficients(cacheDensityProfileCoefficients);
    solver.cacheDomainGeometry(cacheDomainGeometry);

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
    solver.setSolution(&exact_solution);
    solver.solve(boundary_conditions, source_term);

    solver.printTimings();

    std::string stencil_string = "";
    if (solver.stencilDistributionMethod() == StencilDistributionMethod::CPU_TAKE) {
        stencil_string = "Take";
    }
    else if (solver.stencilDistributionMethod() == StencilDistributionMethod::CPU_GIVE) {
        stencil_string = "Give";
    }

    int extrapolation_int = static_cast<int>(solver.extrapolation());

    // Write results to file
    outfile << Kokkos::num_threads() << "," << divideBy2 << "," << solver.grid().nr() << "," << solver.grid().ntheta()
            << "," << geometry_string << "," << stencil_string << "," << cacheDensityProfileCoefficients << ","
            << cacheDomainGeometry << "," << FMG << "," << extrapolation_int << ","
            << solver.timeSetupTotal() + solver.timeSolveTotal() << "," << solver.timeSetupTotal() << ","
            << solver.timeSetupCreateLevels() << "," << solver.timeSetupSmoother() << ","
            << solver.timeSetupDirectSolver() << "," << solver.timeSolveTotal() << ","
            << solver.timeSolveInitialApproximation() << "," << solver.timeSolveMultigridIterations() << ","
            << solver.timeCheckConvergence() << "," << solver.timeCheckExactError() << "," << solver.timeAvgMGCTotal()
            << "," << solver.timeAvgMGCPreSmoothing() << "," << solver.timeAvgMGCPostSmoothing() << ","
            << solver.timeAvgMGCResidual() << "," << solver.timeAvgMGCDirectSolver() << std::endl;
}

int main(int argc, char* argv[])
{
    Kokkos::ScopeGuard kokkos_scope(argc, argv);

    if (argc < 2) {
        std::cerr << "Usage: strong_scaling_run <output_csv>\n";
        return 1;
    }

    std::ofstream outfile(argv[1], std::ios::app);
    runTest(7, outfile);
    return 0;
}
