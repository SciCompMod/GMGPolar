#include <iostream>
#include <chrono>
#include <memory>
#include <random>

// #include "../include/InputFunctions/domainGeometry.h"
// #include "../include/InputFunctions/DomainGeometry/circularGeometry.h"
// #include "../include/InputFunctions/DomainGeometry/shafranovGeometry.h"
// #include "../include/InputFunctions/DomainGeometry/czarnyGeometry.h"
// #include "../include/InputFunctions/DomainGeometry/culhamGeometry.h"

// #include "../include/InputFunctions/exactSolution.h"
// #include "../include/InputFunctions/systemParameters.h"

// #include "../include/InputFunctions/ExactSolution/cartesianR2_CircularGeometry.h"

// #include "../include/InputFunctions/BoundaryConditions/cartesianR2_Boundary_CircularGeometry.h"
// #include "../include/InputFunctions/DensityProfileCoefficients/sonnendruckerGyroCoefficients.h"
// #include "../include/InputFunctions/SourceTerms/cartesianR2_SonnendruckerGyro_CircularGeometry.h"

#include "../include/GMGPolar/gmgpolar.h"

int main(int argc, char* argv[]){
    #ifdef NDEBUG
        std::cout << "Build Type: Release\n"<<std::endl;
    #else
        std::cout << "Build Type: Debug\n"<<std::endl;
    #endif

    std::cout << std::scientific  << std::setprecision(10);

    GMGPolar solver;

    // Configure solver parameters from command-line arguments
    solver.setParameters(argc, argv);
    
    // Setup and solve the problem
    solver.setup();
    solver.solve();

    // Recover solution
    Vector<double>& solution = solver.solution();
    const PolarGrid& grid = solver.grid();

    // Print Timings
    solver.printTimings();





    // ShafranovGeometry domain_geometry;
    
    // std::unique_ptr<ExactSolution> exact_solution = std::make_unique<CartesianR2_CircularGeometry>(Rmax);

    // std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<SonnendruckerGyroCoefficients>(Rmax, alpha_jump);
    // std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<CartesianR2_Boundary_CircularGeometry>(Rmax);
    // std::unique_ptr<SourceTerm> source_term = std::make_unique<CartesianR2_SonnendruckerGyro_CircularGeometry>(Rmax);

    // SystemParameters system_parameters(std::move(coefficients), std::move(boundary_conditions), std::move(source_term));

    // // Optional exact solution
    // ExactSolution exact_solution; 

    // // Initialize the solver with domain geometry and system parameters
    // GMGPolar solver(domain_geometry, system_parameters);


    // // Set the exact solution if available
    // solver.setSolution(exact_solution);

    // // Configure solver parameters from command-line arguments
    // solver.setParameters(argc, argv);
    
    // // Time the setup section
    // auto start = std::chrono::high_resolution_clock::now();
    
    // // Setup and solve the problem
    // solver.setup();
    // auto endSetup = std::chrono::high_resolution_clock::now();
    
    // solver.solve();
    // auto endSolve = std::chrono::high_resolution_clock::now();
    
    // // Calculate durations
    // auto setupDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endSetup - start);
    // auto solveDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endSolve - endSetup);
    // auto totalDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endSolve - start);

    // // Print durations
    // std::cout << "Setup time: " << setupDuration.count() << " ms" << std::endl;
    // std::cout << "Solve time: " << solveDuration.count() << " ms" << std::endl;
    // std::cout << "Total time: " << totalDuration.count() << " ms" << std::endl;

    return 0;
}