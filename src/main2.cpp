#include <iostream>
#include <chrono>
#include <memory>
#include <random>

#include "../include/InputFunctions/domainGeometry.h"
#include "../include/InputFunctions/DomainGeometry/circularGeometry.h"
#include "../include/InputFunctions/DomainGeometry/shafranovGeometry.h"
#include "../include/InputFunctions/DomainGeometry/czarnyGeometry.h"
#include "../include/InputFunctions/DomainGeometry/culhamGeometry.h"

#include "../include/InputFunctions/exactSolution.h"
#include "../include/InputFunctions/densityProfileCoefficients.h"
#include "../include/InputFunctions/boundaryConditions.h"
#include "../include/InputFunctions/sourceTerm.h"

#include "../include/InputFunctions/ExactSolution/cartesianR2_CzarnyGeometry.h"

#include "../include/InputFunctions/BoundaryConditions/cartesianR2_Boundary_CzarnyGeometry.h"
#include "../include/InputFunctions/DensityProfileCoefficients/sonnendruckerGyroCoefficients.h"
#include "../include/InputFunctions/SourceTerms/cartesianR2_SonnendruckerGyro_CzarnyGeometry.h"

#include "../include/GMGPolar/gmgpolar.h"

int main(int argc, char* argv[]){
    #ifdef NDEBUG
        std::cout << "Build Type: Release\n"<< std::endl;
    #else
        std::cout << "Build Type: Debug\n"<< std::endl;
    #endif

    double R0 = 1e-5;
    double Rmax = 1.3;
    double inverse_aspect_ratio_epsilon = 0.3;
    double ellipticity_e = 1.4;
    double alpha_jump = 0.66;

    std::unique_ptr<DomainGeometry> domain_geometry = std::make_unique<CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<ExactSolution> exact_solution = std::make_unique<CartesianR2_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);

    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<SonnendruckerGyroCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<CartesianR2_Boundary_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CzarnyGeometry>(Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);

    // Initialize the solver with domain geometry and system parameters
    GMGPolar solver(std::move(domain_geometry), std::move(coefficients), std::move(boundary_conditions), std::move(source_term));

    // Set the exact solution if available
    solver.setSolution(std::move(exact_solution));

    int verbose = 0;
    int maxOpenMPThreads = 8;
    double threadReductionFactor = 1.0;

    int nr_exp = 4;
    int ntheta_exp = 4;
    int anisotropic_factor = 0;
    int divideBy2 = 0;

    bool DirBC_Interior = true;

    int extrapolation = 1;
    int maxLevels = 6;
    int preSmoothingSteps = 1;
    int postSmoothingSteps = 1;
    MultigridCycleType multigridCycle = MultigridCycleType::V_CYCLE;

    int maxIterations = 150;
    ResidualNormType residualNormType = ResidualNormType::EUCLIDEAN;
    double absoluteTolerance = 1e-12;
    double relativeTolerance = 1e-5;

    /* --------------- */
    /* Grid Parameters */
    solver.R0(R0);
    solver.Rmax(Rmax);
    solver.nr_exp(nr_exp);
    solver.ntheta_exp(ntheta_exp);
    solver.anisotropic_factor(anisotropic_factor);
    solver.divideBy2(divideBy2);
    
    solver.write_grid_file(false);
    solver.load_grid_file(false);
    solver.file_grid_radii("");
    solver.file_grid_angles("");

    /* ------------------- */
    /* Geometry Parameters */
    solver.DirBC_Interior(DirBC_Interior);

    /* -------------------- */
    /* Multigrid Parameters */
    solver.extrapolation(extrapolation);

    solver.maxLevels(maxLevels);
    solver.multigridCycle(multigridCycle);
    
    solver.preSmoothingSteps(preSmoothingSteps);
    solver.postSmoothingSteps(postSmoothingSteps);

    solver.maxIterations(maxIterations);
    solver.residualNormType(residualNormType);
    solver.absoluteTolerance(absoluteTolerance);
    solver.relativeTolerance(relativeTolerance);

    /* ------------------ */
    /* Control Parameters */
    solver.verbose(verbose);
    solver.maxOpenMPThreads(maxOpenMPThreads);
    solver.threadReductionFactor(threadReductionFactor);

    int benchmark_size = 8;

    std::vector<double> setup_times(benchmark_size);
    std::vector<double> solve_times(benchmark_size);
    std::vector<double> sqrt_n(benchmark_size);

    for (int divideBy2 = 0; divideBy2 < benchmark_size; divideBy2++)
    {
        solver.divideBy2(divideBy2);

        auto start_setup = std::chrono::high_resolution_clock::now();
        solver.setup();
        auto end_setup = std::chrono::high_resolution_clock::now();
        
        auto setup_duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_setup - start_setup).count();
        setup_times[divideBy2] = setup_duration;

        auto start_solve = std::chrono::high_resolution_clock::now();
        solver.solve();
        auto end_solve = std::chrono::high_resolution_clock::now();
        
        auto solve_duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_solve - start_solve).count();
        solve_times[divideBy2] = solve_duration;

        // Recover solution
        Vector<double>& solution = solver.solution();
        const PolarGrid& grid = solver.grid();

        int num_nodes = grid.numberOfNodes(); // Assuming this method exists and returns an integer
        sqrt_n[divideBy2] = std::sqrt(num_nodes);
    }

    std::ofstream file("timings.txt");
    if (file.is_open()) {
        file << "Benchmark Size: " << benchmark_size << "\n";
        file << "Setup Times (seconds):\n";
        for (const auto& time : setup_times) {
            file << time << "\n";
        }
        file << "Solve Times (seconds):\n";
        for (const auto& time : solve_times) {
            file << time << "\n";
        }
        file << "Square Root of Number of Nodes:\n";
        for (const auto& sqrt_nodes : sqrt_n) {
            file << sqrt_nodes << "\n";
        }
        file.close();
        std::cout << "Timings and square root of number of nodes have been written to timings.txt\n";
    } else {
        std::cerr << "Unable to open file";
    }
    
    return 0;
}