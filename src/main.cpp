#include <iostream>

#include "../include/GMGPolar/gmgpolar.h"

int main(int argc, char* argv[])
{
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
