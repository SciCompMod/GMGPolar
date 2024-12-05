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
    std::cout << "Likwid: ON\n" << std::endl;
#else
    std::cout << "Likwid: OFF\n" << std::endl;
#endif

    // Initialize solver and set parameters from command-line arguments
    GMGPolar solver;
    solver.setParameters(argc, argv);

    // Initialize LIKWID markers if enabled
    LIKWID_INIT();
    LIKWID_REGISTER("Setup");
    LIKWID_REGISTER("Solve");

    // Run Solver Setup with optional LIKWID markers
    LIKWID_START("Setup");
    solver.setup();
    LIKWID_STOP("Setup");

    // Execute Solve Phase with optional LIKWID markers
    LIKWID_START("Solve");
    solver.solve();
    LIKWID_STOP("Solve");

    // Finalize LIKWID markers if enabled
    LIKWID_CLOSE();

    // Retrieve and print solution and timings
    Vector<double>& solution = solver.solution();
    const PolarGrid& grid    = solver.grid();

    solver.printTimings();

    return 0;
}
