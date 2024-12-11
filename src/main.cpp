#include <iostream>

#include "../include/GMGPolar/gmgpolar.h"

#include "../include/DirectSolver/directSolver.h"

int main(int argc, char* argv[]) {
    // Display Build Type
    #ifdef NDEBUG
        std::cout << "Build Type: Release" << std::endl;
    #else
        std::cout << "Build Type: Debug" << std::endl;
    #endif

    auto domain_geometry = std::make_unique<DomainGeometry>();
    auto exact_solution = std::make_unique<const ExactSolution>();
    auto coefficients = std::make_unique<const DensityProfileCoefficients>();
    auto boundary_conditions = std::make_unique<const BoundaryConditions>();
    auto source_term = std::make_unique<const SourceTerm>();

    GMGPolar solver(std::move(domain_geometry), std::move(coefficients), std::move(boundary_conditions), std::move(source_term));
    solver.setSolution(std::move(exact_solution));

    solver.setParameters(argc, argv);

    solver.setup();

    solver.solve();

    // solver.printTimings();

    return 0;
}
