#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <chrono>
#include <limits>

#include "../include/GMGPolar/gmgpolar.h"
#include "../include/GMGPolar/test_cases.h"

int main(int argc, char* argv[]){
    #ifdef NDEBUG
        std::cout << "Build Type: Release\n"<< std::endl;
    #else
        std::cout << "Build Type: Debug\n"<< std::endl;
    #endif

    const double R0 = 1e-8; const double Rmax = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3; const double ellipticity_e = 1.4;
    const double elongation_kappa = 0.3; const double shift_delta = 0.2;
    const double alpha_jump = 0.66;

    /* Example 1: Cartesian Solution -> Lower Order */
    std::unique_ptr<DomainGeometry> domain_geometry = std::make_unique<CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<ExactSolution> exact_solution = std::make_unique<CartesianR2_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<SonnendruckerGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<CartesianR2_Boundary_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);

    /* Example 2: Polar Solution -> Higher Order */
    // std::unique_ptr<DomainGeometry> domain_geometry = std::make_unique<ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    // std::unique_ptr<ExactSolution> exact_solution = std::make_unique<PolarR6_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    // std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniGyroCoefficients>(Rmax, alpha_jump);
    // std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<PolarR6_Boundary_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);
    // std::unique_ptr<SourceTerm> source_term = std::make_unique<PolarR6_ZoniGyro_ShafranovGeometry>(Rmax, elongation_kappa, shift_delta);

    GMGPolar solver(std::move(domain_geometry), std::move(coefficients), std::move(boundary_conditions), std::move(source_term));

    solver.setSolution(std::move(exact_solution));

    const int verbose = 0;
    const bool paraview = false;

    const int maxOpenMPThreads = 32;
    const double threadReductionFactor = 1.0;

    const ImplementationType implementationType = ImplementationType::CPU_GIVE;
    const bool cacheDensityProfileCoefficients = true;
    const bool cacheDomainGeometry = false;

    const int nr_exp = 4;
    const int ntheta_exp = 6;
    const int anisotropic_factor = 3;
    int divideBy2 = 0; const int MAX_DIVIDE_BY_2 = 6;

    const bool DirBC_Interior = false;

    const bool FMG = true;
    const int FMG_iterations = 3;
    const MultigridCycleType FMG_cycle = MultigridCycleType::F_CYCLE;

    const ExtrapolationType extrapolation = ExtrapolationType::IMPLICIT_EXTRAPOLATION;
    const int maxLevels = 6;
    const int preSmoothingSteps = 1; 
    const int postSmoothingSteps = 1;
    const MultigridCycleType multigridCycle = MultigridCycleType::F_CYCLE;

    const int maxIterations = 150;
    const ResidualNormType residualNormType = ResidualNormType::EUCLIDEAN;
    const double absoluteTolerance = 1e-10;
    const double relativeTolerance = 1e-8;

    solver.verbose(verbose);
    solver.paraview(paraview);

    solver.maxOpenMPThreads(maxOpenMPThreads);
    solver.threadReductionFactor(threadReductionFactor);

    solver.implementationType(implementationType);
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

    for (divideBy2 = 0; divideBy2 < MAX_DIVIDE_BY_2; divideBy2++)
    {
        solver.divideBy2(divideBy2);

        auto start = std::chrono::high_resolution_clock::now();

        solver.setup();
        solver.solve();

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed_seconds = end - start;
        double total_time_seconds = elapsed_seconds.count(); 

        table_nr[divideBy2] = solver.grid().nr();
        table_ntheta[divideBy2] = solver.grid().ntheta();
        table_dofs[divideBy2] = solver.grid().numberOfNodes();
        table_iterations[divideBy2] = solver.numberOfIterations();
        table_reduction_factor[divideBy2] = solver.meanResidualReductionFactor();
        table_exact_error_weighted_euclidean[divideBy2] = solver.exactErrorWeightedEuclidean().value();
        table_exact_error_infinity[divideBy2] = solver.exactErrorInfinity().value();
        table_total_time[divideBy2] = total_time_seconds;
    }

    table_exact_error_weighted_euclidean_order[0] = std::numeric_limits<double>::max();;
    table_exact_error_infinity_order[0] = std::numeric_limits<double>::max();;

    for (int i = 1; i < MAX_DIVIDE_BY_2; i++)
    {
        table_exact_error_weighted_euclidean_order[i] = 
            log(table_exact_error_weighted_euclidean[i-1] / table_exact_error_weighted_euclidean[i]) /
            log(sqrt(static_cast<double>(table_dofs[i]) / static_cast<double>(table_dofs[i-1])));

        table_exact_error_infinity_order[i] = 
            log(table_exact_error_infinity[i-1] / table_exact_error_infinity[i]) /
            log(sqrt(static_cast<double>(table_dofs[i]) / static_cast<double>(table_dofs[i-1])));
    }
    

    // Set up table formatting
    std::cout << std::fixed << std::setprecision(2);

    // Print the header
    std::cout << std::setw(12) << "nr x nθ"
              << std::setw(8)  << "its"
              << std::setw(8)  << "ρ̂"
              << std::setw(25) << "||err||_2/sqrt(m)"
              << std::setw(6)  << "ord"
              << std::setw(20) << "||err||_∞"
              << std::setw(6)  << "ord"
              << std::setw(15) << "time[s]" << std::endl;

    // Print the table rows with data
    for (int i = 0; i < MAX_DIVIDE_BY_2; i++) {
        // Start the row
        std::cout << std::setw(6)  << table_nr[i] << " x " << table_ntheta[i]
                  << std::setw(7)  << table_iterations[i]
                  << std::setw(8)  << table_reduction_factor[i];

        // Print ||err||_2/sqrt(m) in scientific notation
        std::cout << std::scientific << std::setprecision(3);
        std::cout << std::setw(20) << table_exact_error_weighted_euclidean[i];

        // Print ord for ||err||_2/sqrt(m)
        std::cout << std::fixed << std::setprecision(1);  // Go back to fixed notation
        std::cout << std::setw(10);
        if (table_exact_error_weighted_euclidean_order[i] != std::numeric_limits<double>::max())
            std::cout << table_exact_error_weighted_euclidean_order[i];
        else
            std::cout << "-";

        // Print ||err||_∞ in scientific notation
        std::cout << std::scientific << std::setprecision(3);
        std::cout << std::setw(18) << table_exact_error_infinity[i];

        // Print ord for ||err||_∞
        std::cout << std::fixed << std::setprecision(1);  // Back to fixed
        std::cout << std::setw(8);
        if (table_exact_error_infinity_order[i] != std::numeric_limits<double>::max())
            std::cout << table_exact_error_infinity_order[i];
        else
            std::cout << "-";

        // Print time[s]
        std::cout << std::setw(12) << table_total_time[i]
                  << std::endl;
    }
    std::cout<<"\n";

    return 0;



    
}