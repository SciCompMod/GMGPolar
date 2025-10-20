#include <iostream>

#include "../include/ConfigParser/config_parser.h"
#include "../include/GMGPolar/gmgpolar.h"

int main(int argc, char* argv[])
{
    Kokkos::ScopeGuard kokkos_scope(argc, argv);
    // Initialize LIKWID markers if enabled
    LIKWID_INIT();

    // Parse command-line arguments to extract problem configuration
    ConfigParser parser;
    parser.parse(argc, argv);

    // Create GMGPolar solver
    GMGPolar solver(parser.grid(), parser.domainGeometry(), parser.densityProfileCoefficients());

    // --- General solver output and visualization settings --- //
    solver.verbose(parser.verbose()); // Enable/disable verbose output
    solver.paraview(parser.paraview()); // Enable/disable ParaView output

    // --- Parallelization and threading settings --- //
    solver.maxOpenMPThreads(parser.maxOpenMPThreads()); // Maximum OpenMP threads to use
    solver.threadReductionFactor(parser.threadReductionFactor()); // Reduce threads on coarser grids

    omp_set_num_threads(parser.maxOpenMPThreads()); // Global OpenMP thread limit

    // --- Numerical method setup --- //
    solver.DirBC_Interior(parser.DirBC_Interior()); // Interior boundary conditions: Dirichlet, Across-the-origin,
    solver.stencilDistributionMethod(parser.stencilDistributionMethod()); // Stencil distribution strategy: Take, Give
    solver.cacheDensityProfileCoefficients(
        parser.cacheDensityProfileCoefficients()); // Cache density profile coefficients: alpha, beta
    solver.cacheDomainGeometry(parser.cacheDomainGeometry()); // Cache domain geometry data: arr, att, art, detDF

    // --- Multigrid settings --- //
    solver.extrapolation(parser.extrapolation()); // Enable/disable extrapolation
    solver.maxLevels(parser.maxLevels()); // Max multigrid levels (-1 = use deepest possible)
    solver.preSmoothingSteps(parser.preSmoothingSteps()); // Smoothing before coarse-grid correction
    solver.postSmoothingSteps(parser.postSmoothingSteps()); // Smoothing after coarse-grid correction
    solver.multigridCycle(parser.multigridCycle()); // Multigrid cycle type
    solver.FMG(parser.FMG()); // Full Multigrid mode on/off
    solver.FMG_iterations(parser.FMG_iterations()); // FMG iteration count
    solver.FMG_cycle(parser.FMG_cycle()); // FMG cycle type

    // --- Iterative solver controls --- //
    solver.maxIterations(parser.maxIterations()); // Max number of iterations
    solver.residualNormType(parser.residualNormType()); // Residual norm type (L2, weighted-L2, L∞)
    solver.absoluteTolerance(parser.absoluteTolerance()); // Absolute residual tolerance
    solver.relativeTolerance(parser.relativeTolerance()); // Relative residual tolerance

    // --- Finalize solver setup --- //
    solver.setup(); // (allocates internal data, prepares operators, etc.)

    // --- Provide optional exact solution --- //
    solver.setSolution(&parser.exactSolution());
    // --- Solve Phase --- //
    solver.solve(parser.boundaryConditions(), parser.sourceTerm());

    // --- Retrieve solution and associated grid --- //
    Vector<double> solution = solver.solution();
    const PolarGrid& grid   = solver.grid();

    // Finalize LIKWID performance markers
    LIKWID_CLOSE();

    // Print timing statistics for each solver phase
    solver.printTimings();

    return 0;
}
