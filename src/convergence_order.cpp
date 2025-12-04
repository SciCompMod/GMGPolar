#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include "../include/GMGPolar/gmgpolar.h"
#include "../include/GMGPolar/test_cases.h"

int main(int argc, char* argv[])
{
    omp_set_num_threads(omp_get_max_threads());

    const int verbose   = 0;
    const bool paraview = false;

    const int maxOpenMPThreads         = 16;
    const double threadReductionFactor = 1.0;

    const StencilDistributionMethod stencilDistributionMethod = StencilDistributionMethod::CPU_TAKE;
    const bool cacheDensityProfileCoefficients                = true;
    const bool cacheDomainGeometry                            = true;

    const int nr_exp             = 4;
    const int ntheta_exp         = 6;
    const int anisotropic_factor = 3;
    int divideBy2                = 0;
    const int MAX_DIVIDE_BY_2    = 6;

    const bool DirBC_Interior = false;

    const bool FMG                     = true;
    const int FMG_iterations           = 3;
    const MultigridCycleType FMG_cycle = MultigridCycleType::F_CYCLE;

    const ExtrapolationType extrapolation   = ExtrapolationType::IMPLICIT_EXTRAPOLATION;
    const int maxLevels                     = 6;
    const int preSmoothingSteps             = 1;
    const int postSmoothingSteps            = 1;
    const MultigridCycleType multigridCycle = MultigridCycleType::F_CYCLE;

    const int maxIterations                 = 150;
    const ResidualNormType residualNormType = ResidualNormType::EUCLIDEAN;
    const double absoluteTolerance          = 1e-10;
    const double relativeTolerance          = 1e-8;

    const double R0                           = 1e-8;
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;
    const double elongation_kappa             = 0.3;
    const double shift_delta                  = 0.2;

    /* Example 1: Polar Solution -> Higher Order 4.0 */
    const double alpha_jump = 0.4837 * Rmax;
    ShafranovGeometry domain_geometry(Rmax, elongation_kappa, shift_delta);
    ZoniGyroCoefficients coefficients(Rmax, alpha_jump);
    PolarR6_Boundary_ShafranovGeometry boundary_conditions(Rmax, elongation_kappa, shift_delta);
    PolarR6_ZoniGyro_ShafranovGeometry source_term(Rmax, elongation_kappa, shift_delta);
    PolarR6_ShafranovGeometry exact_solution(Rmax, elongation_kappa, shift_delta);

    /* Example 2: Cartesian Solution -> Lower Order 3.5 */
    // const double alpha_jump = 0.66 * Rmax;
    // CzarnyGeometry domain_geometry(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    // SonnendruckerGyroCoefficients coefficients(Rmax, alpha_jump);
    // CartesianR2_Boundary_CzarnyGeometry boundary_conditions(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    // CartesianR2_SonnendruckerGyro_CzarnyGeometry source_term(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    // CartesianR2_CzarnyGeometry exact_solution(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);

    /* Example 3: Refined Solution -> Lower Order 3.5 */
    // const double alpha_jump = 0.9 * Rmax; // Refinement where the solution is most complex
    // ShafranovGeometry domain_geometry(Rmax, elongation_kappa, shift_delta);
    // ZoniShiftedGyroCoefficients coefficients(Rmax, alpha_jump);
    // Refined_Boundary_ShafranovGeometry boundary_conditions(Rmax, elongation_kappa, shift_delta);
    // Refined_ZoniShiftedGyro_ShafranovGeometry source_term(Rmax, elongation_kappa, shift_delta);
    // Refined_ShafranovGeometry exact_solution(Rmax, elongation_kappa, shift_delta);

    std::vector<int> table_nr(MAX_DIVIDE_BY_2);
    std::vector<int> table_ntheta(MAX_DIVIDE_BY_2);
    std::vector<int> table_dofs(MAX_DIVIDE_BY_2);
    std::vector<int> table_iterations(MAX_DIVIDE_BY_2);
    std::vector<double> table_reduction_factor(MAX_DIVIDE_BY_2);
    std::vector<double> table_exact_error_weighted_euclidean(MAX_DIVIDE_BY_2);
    std::vector<double> table_exact_error_weighted_euclidean_order(MAX_DIVIDE_BY_2);
    std::vector<double> table_exact_error_infinity(MAX_DIVIDE_BY_2);
    std::vector<double> table_exact_error_infinity_order(MAX_DIVIDE_BY_2);
    std::vector<double> table_total_time(MAX_DIVIDE_BY_2);

    for (divideBy2 = 0; divideBy2 < MAX_DIVIDE_BY_2; divideBy2++) {
        double refinement_radius               = alpha_jump;
        std::optional<double> splitting_radius = std::nullopt;
        PolarGrid grid(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2,
                       splitting_radius);
        GMGPolar solver(grid, domain_geometry, coefficients);

        solver.verbose(verbose);
        solver.paraview(paraview);

        solver.maxOpenMPThreads(maxOpenMPThreads);
        solver.threadReductionFactor(threadReductionFactor);

        omp_set_num_threads(maxOpenMPThreads); // Global OpenMP thread limit

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

        solver.setup();
        solver.setSolution(&exact_solution);
        solver.solve(boundary_conditions, source_term);

        table_nr[divideBy2]                             = solver.grid().nr();
        table_ntheta[divideBy2]                         = solver.grid().ntheta();
        table_dofs[divideBy2]                           = solver.grid().numberOfNodes();
        table_iterations[divideBy2]                     = solver.numberOfIterations();
        table_reduction_factor[divideBy2]               = solver.meanResidualReductionFactor();
        table_exact_error_weighted_euclidean[divideBy2] = solver.exactErrorWeightedEuclidean().value();
        table_exact_error_infinity[divideBy2]           = solver.exactErrorInfinity().value();
        table_total_time[divideBy2]                     = solver.timeSetupTotal() + solver.timeSolveTotal();
    }

    table_exact_error_weighted_euclidean_order[0] = std::numeric_limits<double>::max();

    table_exact_error_infinity_order[0] = std::numeric_limits<double>::max();

    for (int i = 1; i < MAX_DIVIDE_BY_2; i++) {
        table_exact_error_weighted_euclidean_order[i] =
            log(table_exact_error_weighted_euclidean[i - 1] / table_exact_error_weighted_euclidean[i]) /
            log(sqrt(static_cast<double>(table_dofs[i]) / static_cast<double>(table_dofs[i - 1])));

        table_exact_error_infinity_order[i] =
            log(table_exact_error_infinity[i - 1] / table_exact_error_infinity[i]) /
            log(sqrt(static_cast<double>(table_dofs[i]) / static_cast<double>(table_dofs[i - 1])));
    }

    // Set up table formatting
    std::cout << std::fixed << std::setprecision(2);

    // Print the header
    std::cout << std::setw(12) << "nr x nθ" << std::setw(9) << "its" << std::setw(7) << "ρ" << std::setw(25)
              << "||err||_2/sqrt(m)" << std::setw(6) << "ord" << std::setw(22) << "||err||_∞" << std::setw(6) << "ord"
              << std::setw(15) << "time[s]" << std::endl;

    // Print the table rows with data
    for (int i = 0; i < MAX_DIVIDE_BY_2; i++) {
        // Start the row
        std::cout << std::setw(6) << table_nr[i] << " x " << table_ntheta[i] << std::setw(7.5) << table_iterations[i]
                  << std::setw(9) << table_reduction_factor[i];

        // Print ||err||_2/sqrt(m) in scientific notation
        std::cout << std::scientific << std::setprecision(2);
        std::cout << std::setw(20) << table_exact_error_weighted_euclidean[i];

        // Print ord for ||err||_2/sqrt(m)
        std::cout << std::fixed << std::setprecision(2); // Go back to fixed notation
        std::cout << std::setw(10);
        if (table_exact_error_weighted_euclidean_order[i] != std::numeric_limits<double>::max())
            std::cout << table_exact_error_weighted_euclidean_order[i];
        else
            std::cout << "-";

        // Print ||err||_∞ in scientific notation
        std::cout << std::scientific << std::setprecision(2);
        std::cout << std::setw(18) << table_exact_error_infinity[i];

        // Print ord for ||err||_∞
        std::cout << std::fixed << std::setprecision(2); // Back to fixed
        std::cout << std::setw(8);
        if (table_exact_error_infinity_order[i] != std::numeric_limits<double>::max())
            std::cout << table_exact_error_infinity_order[i];
        else
            std::cout << "-";

        // Print time[s]
        std::cout << std::setw(12) << table_total_time[i] << std::endl;
    }

    return 0;
}
