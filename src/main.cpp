#include <iostream>
#include <chrono>
#include <memory>
#include <random>

#include "../include/GMGPolar/gmgpolar.h"

int main(int argc, char* argv[]){
    #ifdef NDEBUG
        std::cout << "Build Type: Release\n"<<std::endl;
    #else
        std::cout << "Build Type: Debug\n"<<std::endl;
    #endif

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

    return 0;
}
