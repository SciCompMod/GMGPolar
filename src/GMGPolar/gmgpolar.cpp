#include "../../include/GMGPolar/gmgpolar.h"

// #include <iomanip>

GMGPolar::GMGPolar(std::unique_ptr<const DomainGeometry> domain_geometry, 
                   std::unique_ptr<const DensityProfileCoefficients> density_profile_coefficients,
                   std::unique_ptr<const BoundaryConditions> boundary_conditions,
                   std::unique_ptr<const SourceTerm> source_term) :

    domain_geometry_(std::move(domain_geometry)), 
    density_profile_coefficients_(std::move(density_profile_coefficients)),
    boundary_conditions_(std::move(boundary_conditions)), 
    source_term_(std::move(source_term)),
    parser_()
{
    resetTimings();
    initializeGrid(); initializeMultigrid(); initializeGeneral();
    parseGrid(); parseMultigrid(); parseGeneral();
}


void GMGPolar::setParameters(int argc, char* argv[]) {
    if(argc != 0){
        try {
            parser_.parse_check(argc, argv);
        } catch (const cmdline::cmdline_error& parse_error) {
            std::cerr << "Error: " << parse_error.what() << std::endl;
            std::cerr << "Usage: " << parser_.usage() << std::endl;
        }
    }
    parseGrid(); parseMultigrid(); parseGeneral();
}


void GMGPolar::setSolution(std::unique_ptr<const ExactSolution> exact_solution) {
    exact_solution_ = std::move(exact_solution);
}

void GMGPolar::printTimings() const {
    std::cout << "\n------------------"<< std::endl;
    std::cout << "Timing Information" << std::endl;
    std::cout << "------------------"<< std::endl;
    std::cout << "Setup Time: " << t_setup_total-t_setup_rhs << " seconds" << std::endl;
    std::cout << "    Create Levels: " << t_setup_createLevels << " seconds" << std::endl;
    std::cout << "    Smoother: " << t_setup_smoother << " seconds" << std::endl;
    std::cout << "    Direct Solver: " << t_setup_directSolver << " seconds" << std::endl;
    std::cout << "    (Build rhs: " << t_setup_rhs << " seconds)" << std::endl;
    std::cout << "\nSolve Time: " << t_solve_total << " seconds" << std::endl;
    std::cout << "    Initial Approximation: " << t_solve_initial_approximation << " seconds" << std::endl;
    std::cout << "    Multigrid Iteration: " << t_solve_multigrid_iterations << " seconds" << std::endl;
    std::cout << "    Check Convergence: " << t_check_convergence << " seconds" << std::endl;
    std::cout << "    (Check Exact Error: " << t_check_exact_error << " seconds)" << std::endl;
    std::cout << "\nAverage Multigrid Iteration: " << t_avg_MGC_total << " seconds" << std::endl;
    std::cout << "    Pre Smoothing: " << t_avg_MGC_preSmoothing << " seconds" << std::endl;
    std::cout << "    Post Smoothing: " << t_avg_MGC_postSmoothing << " seconds" << std::endl;
    std::cout << "    Residual: " << t_avg_MGC_residual << " seconds" << std::endl;
    std::cout << "    Direct Solve: " << t_avg_MGC_directSolver << " seconds" << std::endl;
    std::cout << "    Other Computations: " << std::max(t_avg_MGC_total - t_avg_MGC_preSmoothing - t_avg_MGC_postSmoothing - t_avg_MGC_residual - t_avg_MGC_directSolver, 0.0) << " seconds" << std::endl;
    std::cout <<"\n"<< std::endl;
}