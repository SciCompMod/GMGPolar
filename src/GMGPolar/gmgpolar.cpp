#include "../../include/GMGPolar/gmgpolar.h"

/* ---------------------------------------------------------------------- */
/* Constructor & Initialization                                           */
/* ---------------------------------------------------------------------- */
IGMGPolar::IGMGPolar(const PolarGrid& grid, const DensityProfileCoefficients& density_profile_coefficients)
    : grid_(grid)
    , density_profile_coefficients_(density_profile_coefficients)
    , exact_solution_(nullptr)
    // General solver output and visualization settings
    , verbose_(0)
    , paraview_(false)
    // Parallelization and threading settings
    , max_omp_threads_(omp_get_max_threads())
    , thread_reduction_factor_(1.0)
    // Numerical method setup
    , DirBC_Interior_(true)
    , stencil_distribution_method_(StencilDistributionMethod::CPU_GIVE)
    , cache_density_profile_coefficients_(true)
    , cache_domain_geometry_(false)
    // Multigrid settings
    , extrapolation_(ExtrapolationType::IMPLICIT_EXTRAPOLATION)
    , max_levels_(-1)
    , pre_smoothing_steps_(1)
    , post_smoothing_steps_(1)
    , multigrid_cycle_(MultigridCycleType::V_CYCLE)
    // FMG settings
    , FMG_(false)
    , FMG_iterations_(3)
    , FMG_cycle_(MultigridCycleType::F_CYCLE)
    // Convergence settings
    , max_iterations_(300)
    , residual_norm_type_(ResidualNormType::WEIGHTED_EUCLIDEAN)
    , absolute_tolerance_(1e-8)
    , relative_tolerance_(1e-8)
    // Level management and internal solver data
    , number_of_levels_(0)
    , interpolation_(nullptr)
    , full_grid_smoothing_(false)
    , number_of_iterations_(0)
    , mean_residual_reduction_factor_(1.0)
{
    resetAllTimings();
    LIKWID_REGISTER("Setup");
    LIKWID_REGISTER("Solve");
}

void IGMGPolar::setSolution(const ExactSolution* exact_solution)
{
    exact_solution_ = exact_solution;
}

/* ---------------------------------------------------------------------- */
/* General output & visualization                                         */
/* ---------------------------------------------------------------------- */
int IGMGPolar::verbose() const
{
    return verbose_;
}
void IGMGPolar::verbose(int verbose)
{
    verbose_ = verbose;
}

bool IGMGPolar::paraview() const
{
    return paraview_;
}
void IGMGPolar::paraview(bool paraview)
{
    paraview_ = paraview;
}

/* ---------------------------------------------------------------------- */
/* Parallelization & threading                                            */
/* ---------------------------------------------------------------------- */
int IGMGPolar::maxOpenMPThreads() const
{
    return max_omp_threads_;
}
void IGMGPolar::maxOpenMPThreads(int max_omp_threads)
{
    max_omp_threads_ = max_omp_threads;
}

double IGMGPolar::threadReductionFactor() const
{
    return thread_reduction_factor_;
}
void IGMGPolar::threadReductionFactor(double thread_reduction_factor)
{
    thread_reduction_factor_ = thread_reduction_factor;
}

/* ---------------------------------------------------------------------- */
/* Numerical method options                                               */
/* ---------------------------------------------------------------------- */
bool IGMGPolar::DirBC_Interior() const
{
    return DirBC_Interior_;
}
void IGMGPolar::DirBC_Interior(bool DirBC_Interior)
{
    DirBC_Interior_ = DirBC_Interior;
}

StencilDistributionMethod IGMGPolar::stencilDistributionMethod() const
{
    return stencil_distribution_method_;
}
void IGMGPolar::stencilDistributionMethod(StencilDistributionMethod stencil_distribution_method)
{
    stencil_distribution_method_ = stencil_distribution_method;
}

bool IGMGPolar::cacheDensityProfileCoefficients() const
{
    return cache_density_profile_coefficients_;
}
void IGMGPolar::cacheDensityProfileCoefficients(bool cache_density_profile_coefficients)
{
    cache_density_profile_coefficients_ = cache_density_profile_coefficients;
}

bool IGMGPolar::cacheDomainGeometry() const
{
    return cache_domain_geometry_;
}
void IGMGPolar::cacheDomainGeometry(bool cache_domain_geometry)
{
    cache_domain_geometry_ = cache_domain_geometry;
}

/* ---------------------------------------------------------------------- */
/* Multigrid controls                                                     */
/* ---------------------------------------------------------------------- */
ExtrapolationType IGMGPolar::extrapolation() const
{
    return extrapolation_;
}
void IGMGPolar::extrapolation(ExtrapolationType extrapolation)
{
    extrapolation_ = extrapolation;
}

int IGMGPolar::maxLevels() const
{
    return max_levels_;
}
void IGMGPolar::maxLevels(int max_levels)
{
    max_levels_ = max_levels;
}

MultigridCycleType IGMGPolar::multigridCycle() const
{
    return multigrid_cycle_;
}
void IGMGPolar::multigridCycle(MultigridCycleType multigrid_cycle)
{
    multigrid_cycle_ = multigrid_cycle;
}

int IGMGPolar::preSmoothingSteps() const
{
    return pre_smoothing_steps_;
}
void IGMGPolar::preSmoothingSteps(int pre_smoothing_steps)
{
    pre_smoothing_steps_ = pre_smoothing_steps;
}

int IGMGPolar::postSmoothingSteps() const
{
    return post_smoothing_steps_;
}
void IGMGPolar::postSmoothingSteps(int post_smoothing_steps)
{
    post_smoothing_steps_ = post_smoothing_steps;
}

bool IGMGPolar::FMG() const
{
    return FMG_;
}
void IGMGPolar::FMG(bool FMG)
{
    FMG_ = FMG;
}

int IGMGPolar::FMG_iterations() const
{
    return FMG_iterations_;
}
void IGMGPolar::FMG_iterations(int FMG_iterations)
{
    FMG_iterations_ = FMG_iterations;
}

MultigridCycleType IGMGPolar::FMG_cycle() const
{
    return FMG_cycle_;
}
void IGMGPolar::FMG_cycle(MultigridCycleType FMG_cycle)
{
    FMG_cycle_ = FMG_cycle;
}

/* ---------------------------------------------------------------------- */
/* Iterative solver termination                                           */
/* ---------------------------------------------------------------------- */

int IGMGPolar::maxIterations() const
{
    return max_iterations_;
}
void IGMGPolar::maxIterations(int maxIterations)
{
    max_iterations_ = maxIterations;
}

ResidualNormType IGMGPolar::residualNormType() const
{
    return residual_norm_type_;
}
void IGMGPolar::residualNormType(ResidualNormType residualNormType)
{
    residual_norm_type_ = residualNormType;
}

std::optional<double> IGMGPolar::absoluteTolerance() const
{
    return absolute_tolerance_;
}
void IGMGPolar::absoluteTolerance(std::optional<double> tol)
{
    absolute_tolerance_ = tol.has_value() && tol.value() >= 0.0 ? tol : std::nullopt;
}

std::optional<double> IGMGPolar::relativeTolerance() const
{
    return relative_tolerance_;
}
void IGMGPolar::relativeTolerance(std::optional<double> tol)
{
    relative_tolerance_ = tol.has_value() && tol.value() >= 0.0 ? tol : std::nullopt;
}

/* ---------------------------------------------------------------------- */
/* Solution & Grid Access                                                 */
/* ---------------------------------------------------------------------- */
Vector<double> IGMGPolar::solution()
{
    int level_depth = 0;
    return levels_[level_depth].solution();
}
ConstVector<double> IGMGPolar::solution() const
{
    int level_depth = 0;
    return levels_[level_depth].solution();
}

const PolarGrid& IGMGPolar::grid() const
{
    return grid_;
}

/* ---------------------------------------------------------------------- */
/* Setup timings                                                          */
/* ---------------------------------------------------------------------- */
double IGMGPolar::timeSetupTotal() const
{
    return t_setup_total_;
}
double IGMGPolar::timeSetupCreateLevels() const
{
    return t_setup_createLevels_;
}
double IGMGPolar::timeSetupRHS() const
{
    return t_setup_rhs_;
}
double IGMGPolar::timeSetupSmoother() const
{
    return t_setup_smoother_;
}
double IGMGPolar::timeSetupDirectSolver() const
{
    return t_setup_directSolver_;
}

/* ---------------------------------------------------------------------- */
/* Solve timings                                                          */
/* ---------------------------------------------------------------------- */
double IGMGPolar::timeSolveTotal() const
{
    return t_solve_total_;
}
double IGMGPolar::timeSolveInitialApproximation() const
{
    return t_solve_initial_approximation_;
}
double IGMGPolar::timeSolveMultigridIterations() const
{
    return t_solve_multigrid_iterations_;
}
double IGMGPolar::timeCheckConvergence() const
{
    return t_check_convergence_;
}
double IGMGPolar::timeCheckExactError() const
{
    return t_check_exact_error_;
}

/* ---------------------------------------------------------------------- */
/* Average Multigrid Cycle timings                                        */
/* ---------------------------------------------------------------------- */
double IGMGPolar::timeAvgMGCTotal() const
{
    return t_avg_MGC_total_;
}
double IGMGPolar::timeAvgMGCPreSmoothing() const
{
    return t_avg_MGC_preSmoothing_;
}
double IGMGPolar::timeAvgMGCPostSmoothing() const
{
    return t_avg_MGC_postSmoothing_;
}
double IGMGPolar::timeAvgMGCResidual() const
{
    return t_avg_MGC_residual_;
}
double IGMGPolar::timeAvgMGCDirectSolver() const
{
    return t_avg_MGC_directSolver_;
}

/* ---------------------------------------------------------------------- */
/* Reset timings                                                          */
/* ---------------------------------------------------------------------- */

void IGMGPolar::resetAllTimings()
{
    resetSetupPhaseTimings();
    resetSolvePhaseTimings();
    resetAvgMultigridCycleTimings();
}

void IGMGPolar::resetSetupPhaseTimings()
{
    t_setup_total_        = 0.0;
    t_setup_createLevels_ = 0.0;
    t_setup_rhs_          = 0.0;
    t_setup_smoother_     = 0.0;
    t_setup_directSolver_ = 0.0;
}

void IGMGPolar::resetSolvePhaseTimings()
{
    t_solve_total_                 = 0.0;
    t_solve_initial_approximation_ = 0.0;
    t_solve_multigrid_iterations_  = 0.0;
    t_check_convergence_           = 0.0;
    t_check_exact_error_           = 0.0;
}

void IGMGPolar::resetAvgMultigridCycleTimings()
{
    t_avg_MGC_total_         = 0.0;
    t_avg_MGC_preSmoothing_  = 0.0;
    t_avg_MGC_postSmoothing_ = 0.0;
    t_avg_MGC_residual_      = 0.0;
    t_avg_MGC_directSolver_  = 0.0;
}

/* ---------------------------------------------------------------------- */
/* Diagnostics & statistics                                               */
/* ---------------------------------------------------------------------- */
// Print timing breakdown for setup, smoothing, coarse solve, etc.
void IGMGPolar::printTimings() const
{
    // t_setup_rhs_ is neither included in t_setup_total_ and t_solve_total_.
    std::cout << "\n------------------" << std::endl;
    std::cout << "Timing Information" << std::endl;
    std::cout << "------------------" << std::endl;
    std::cout << "Setup Time: " << t_setup_total_ << " seconds" << std::endl;
    std::cout << "    Create Levels: " << t_setup_createLevels_ << " seconds" << std::endl;
    std::cout << "    Smoother: " << t_setup_smoother_ << " seconds" << std::endl;
    std::cout << "    Direct Solver: " << t_setup_directSolver_ << " seconds" << std::endl;
    std::cout << "    (Build rhs: " << t_setup_rhs_ << " seconds)" << std::endl;
    std::cout << "\nSolve Time: " << t_solve_total_ << " seconds" << std::endl;
    std::cout << "    Initial Approximation: " << t_solve_initial_approximation_ << " seconds" << std::endl;
    std::cout << "    Multigrid Iteration: " << t_solve_multigrid_iterations_ << " seconds" << std::endl;
    std::cout << "    Check Convergence: " << t_check_convergence_ << " seconds" << std::endl;
    std::cout << "    (Check Exact Error: " << t_check_exact_error_ << " seconds)" << std::endl;
    std::cout << "\nAverage Multigrid Iteration: " << t_avg_MGC_total_ << " seconds" << std::endl;
    std::cout << "    PreSmoothing: " << t_avg_MGC_preSmoothing_ << " seconds" << std::endl;
    std::cout << "    PostSmoothing: " << t_avg_MGC_postSmoothing_ << " seconds" << std::endl;
    std::cout << "    Residual: " << t_avg_MGC_residual_ << " seconds" << std::endl;
    std::cout << "    DirectSolve: " << t_avg_MGC_directSolver_ << " seconds" << std::endl;
    std::cout << "    Other Computations: "
              << std::max(t_avg_MGC_total_ - t_avg_MGC_preSmoothing_ - t_avg_MGC_postSmoothing_ - t_avg_MGC_residual_ -
                              t_avg_MGC_directSolver_,
                          0.0)
              << " seconds" << std::endl;
}

// Number of iterations taken by last solve.
int IGMGPolar::numberOfIterations() const
{
    return number_of_iterations_;
}

// Mean residual reduction factor per iteration.
double IGMGPolar::meanResidualReductionFactor() const
{
    return mean_residual_reduction_factor_;
}

// Error norms (only available if exact solution was set).
std::optional<double> IGMGPolar::exactErrorWeightedEuclidean() const
{
    if (exact_solution_) {
        return exact_errors_.back().first;
    }
    return std::nullopt;
}
std::optional<double> IGMGPolar::exactErrorInfinity() const
{
    if (exact_solution_) {
        return exact_errors_.back().second;
    }
    return std::nullopt;
}
