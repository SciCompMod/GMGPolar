#include <iostream>
#include <chrono>

#include "../include/InputFunctions/domainGeometry.h"
#include "../include/InputFunctions/systemParameters.h"
#include "../include/InputFunctions/exactSolution.h"

#include "../include/GMGPolar/gmgpolar.h"

#include "../include/LinearAlgebra/symmetricTridiagonalSolver.h"

int main(int argc, char* argv[]){
    #ifdef NDEBUG
        std::cout << "Build Type: Release\n";
    #else
        std::cout << "Build Type: Debug\n";
    #endif

    // Define the problem to be solved
    DomainGeometry domain_geometry;
    SystemParameters system_parameters;

    // Optional exact solution
    ExactSolution exact_solution; 

    // Initialize the solver with domain geometry and system parameters
    GMGPolar solver(domain_geometry, system_parameters);

    // Set the exact solution if available
    solver.setSolution(exact_solution);

    // Configure solver parameters from command-line arguments
    solver.setParameters(argc, argv);

    // // Setup and solve the problem
    // solver.setup();
    // solver.solve();


    // Time the setup section
    auto start = std::chrono::high_resolution_clock::now();
    
    // Setup and solve the problem
    solver.setup();
    auto endSetup = std::chrono::high_resolution_clock::now();
    
    solver.solve();
    auto endSolve = std::chrono::high_resolution_clock::now();
    
    // Calculate durations
    auto setupDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endSetup - start);
    auto solveDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endSolve - endSetup);
    auto totalDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endSolve - start);

    // Print durations
    std::cout << "Setup time: " << setupDuration.count() << " ms" << std::endl;
    std::cout << "Solve time: " << solveDuration.count() << " ms" << std::endl;
    std::cout << "Total time: " << totalDuration.count() << " ms" << std::endl;

    return 0;
}