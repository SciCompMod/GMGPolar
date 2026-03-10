#include "../../include/GMGPolar/gmgpolar.h"

// =============================================================================
//   Convergence Functions
// =============================================================================

bool IGMGPolar::converged(double residual_norm, double relative_residual_norm)
{
    if (relative_tolerance_.has_value()) {
        if (!(relative_residual_norm > relative_tolerance_.value())) {
            return true;
        }
    }
    if (absolute_tolerance_.has_value()) {
        if (!(residual_norm > absolute_tolerance_.value())) {
            return true;
        }
    }
    return false;
}

// =============================================================================
//   Output and Logging Functions
// =============================================================================

void IGMGPolar::printIterationHeader(const ExactSolution* exact_solution)
{
    if (verbose_ <= 0)
        return;

    const int table_spacing = 4;
    std::cout << "------------------------------\n";
    std::cout << "---- Multigrid Iterations ----\n";
    std::cout << "------------------------------\n";
    std::cout << std::left;
    std::cout << std::setw(3 + table_spacing) << "it";
    if (absolute_tolerance_.has_value() || relative_tolerance_.has_value()) {
        std::cout << std::setw(9 + table_spacing) << "||r_k||";
        if (!FMG_)
            std::cout << std::setw(15 + table_spacing) << "||r_k||/||r_0||";
        else
            std::cout << std::setw(15 + table_spacing) << "||r_k||/||rhs||";
    }
    if (exact_solution != nullptr) {
        std::cout << std::setw(12 + table_spacing) << "||u-u_k||_l2";
        std::cout << std::setw(13 + table_spacing) << "||u-u_k||_inf";
    }
    std::cout << "\n";
    std::cout << std::right << std::setfill(' ');
}

void IGMGPolar::printIterationInfo(int iteration, double current_residual_norm, double current_relative_residual_norm,
                                   const ExactSolution* exact_solution)
{
    if (verbose_ <= 0)
        return;

    const int table_spacing = 4;
    std::cout << std::left << std::scientific << std::setprecision(2);
    std::cout << std::setw(3 + table_spacing) << iteration;
    if (absolute_tolerance_.has_value() || relative_tolerance_.has_value()) {
        std::cout << std::setw(9 + table_spacing + 2) << current_residual_norm;
        std::cout << std::setw(15 + table_spacing) << current_relative_residual_norm;
    }
    if (exact_solution != nullptr) {
        std::cout << std::setw(12 + table_spacing) << exact_errors_.back().first;
        std::cout << std::setw(13 + table_spacing) << exact_errors_.back().second;
    }
    std::cout << "\n";
    std::cout << std::right << std::defaultfloat << std::setprecision(6) << std::setfill(' ');
}
