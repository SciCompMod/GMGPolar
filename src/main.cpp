#include <iostream>

#include "../include/GMGPolar/gmgpolar.h"

int main(int argc, char* argv[])
{
// Display Build Type
#ifdef NDEBUG
    std::cout << "Build Type: Release" << std::endl;
#else
    std::cout << "Build Type: Debug" << std::endl;
#endif

// Check Likwid Status
#ifdef GMGPOLAR_USE_LIKWID
    std::cout << "Likwid: ON" << std::endl;
#else
    std::cout << "Likwid: OFF" << std::endl;
#endif

// Check MUMPS Status
#ifdef GMGPOLAR_USE_MUMPS
    std::cout << "MUMPS: ON\n" << std::endl;
#else
    std::cout << "MUMPS: OFF\n" << std::endl;
#endif

    // Initialize LIKWID markers if enabled
    LIKWID_INIT();

    // Initialize solver and set parameters from command-line arguments
    GMGPolar solver;
    solver.setParameters(argc, argv);
    // Run Solver Setup with optional LIKWID markers
    solver.setup();
    // Execute Solve Phase with optional LIKWID markers
    solver.solve();

    // Finalize LIKWID markers if enabled
    LIKWID_CLOSE();

    // Retrieve and print solution and timings
    Vector<double>& solution = solver.solution();
    const PolarGrid& grid    = solver.grid();

    solver.printTimings();

    return 0;
}
