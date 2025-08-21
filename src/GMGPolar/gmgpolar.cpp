#include "../../include/GMGPolar/gmgpolar.h"

#include <iomanip>

GMGPolar::GMGPolar(const PolarGrid& grid, const DomainGeometry& domain_geometry,
                   const DensityProfileCoefficients& density_profile_coefficients)
    : grid_(grid)
    , domain_geometry_(domain_geometry)
    , density_profile_coefficients_(density_profile_coefficients)
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

/* ---------------------------------------------------------------------- */
/* General output & visualization                                         */
/* ---------------------------------------------------------------------- */
int GMGPolar::verbose() const
{
    return verbose_;
}
void GMGPolar::verbose(int verbose)
{
    verbose_ = verbose;
}

bool GMGPolar::paraview() const
{
    return paraview_;
}
void GMGPolar::paraview(bool paraview)
{
    paraview_ = paraview;
}

/* ---------------------------------------------------------------------- */
/* Parallelization & threading                                            */
/* ---------------------------------------------------------------------- */
int GMGPolar::maxOpenMPThreads() const
{
    return max_omp_threads_;
}
void GMGPolar::maxOpenMPThreads(int max_omp_threads_in)
{
    max_omp_threads_ = max_omp_threads_in;
}

double GMGPolar::threadReductionFactor() const
{
    return thread_reduction_factor_;
}
void GMGPolar::threadReductionFactor(double factor)
{
    thread_reduction_factor_ = factor;
}

/* ---------------------------------------------------------------------- */
/* Numerical method options                                               */
/* ---------------------------------------------------------------------- */
bool GMGPolar::DirBC_Interior() const
{
    return DirBC_Interior_;
}
void GMGPolar::DirBC_Interior(bool val)
{
    DirBC_Interior_ = val;
}

StencilDistributionMethod GMGPolar::stencilDistributionMethod() const
{
    return stencil_distribution_method_;
}
void GMGPolar::stencilDistributionMethod(StencilDistributionMethod method)
{
    stencil_distribution_method_ = method;
}

bool GMGPolar::cacheDensityProfileCoefficients() const
{
    return cache_density_profile_coefficients_;
}
void GMGPolar::cacheDensityProfileCoefficients(bool cache)
{
    cache_density_profile_coefficients_ = cache;
}

bool GMGPolar::cacheDomainGeometry() const
{
    return cache_domain_geometry_;
}
void GMGPolar::cacheDomainGeometry(bool cache)
{
    cache_domain_geometry_ = cache;
}

/* ---------------------------------------------------------------------- */
/* Multigrid controls                                                     */
/* ---------------------------------------------------------------------- */
ExtrapolationType GMGPolar::extrapolation() const
{
    return extrapolation_;
}
void GMGPolar::extrapolation(ExtrapolationType type)
{
    extrapolation_ = type;
}

int GMGPolar::maxLevels() const
{
    return max_levels_;
}
void GMGPolar::maxLevels(int levels)
{
    max_levels_ = levels;
}

MultigridCycleType GMGPolar::multigridCycle() const
{
    return multigrid_cycle_;
}
void GMGPolar::multigridCycle(MultigridCycleType type)
{
    multigrid_cycle_ = type;
}

int GMGPolar::preSmoothingSteps() const
{
    return pre_smoothing_steps_;
}
void GMGPolar::preSmoothingSteps(int steps)
{
    pre_smoothing_steps_ = steps;
}

int GMGPolar::postSmoothingSteps() const
{
    return post_smoothing_steps_;
}
void GMGPolar::postSmoothingSteps(int steps)
{
    post_smoothing_steps_ = steps;
}

bool GMGPolar::FMG() const
{
    return FMG_;
}
void GMGPolar::FMG(bool enable)
{
    FMG_ = enable;
}

int GMGPolar::FMG_iterations() const
{
    return FMG_iterations_;
}
void GMGPolar::FMG_iterations(int iterations)
{
    FMG_iterations_ = iterations;
}

MultigridCycleType GMGPolar::FMG_cycle() const
{
    return FMG_cycle_;
}
void GMGPolar::FMG_cycle(MultigridCycleType type)
{
    FMG_cycle_ = type;
}

/* ---------------------------------------------------------------------- */
/* Iterative solver termination                                           */
/* ---------------------------------------------------------------------- */
int GMGPolar::maxIterations() const
{
    return max_iterations_;
}
void GMGPolar::maxIterations(int max_iterations_in)
{
    max_iterations_ = max_iterations_in;
}

ResidualNormType GMGPolar::residualNormType() const
{
    return residual_norm_type_;
}
void GMGPolar::residualNormType(ResidualNormType type)
{
    residual_norm_type_ = type;
}

std::optional<double> GMGPolar::absoluteTolerance() const
{
    return absolute_tolerance_;
}
void GMGPolar::absoluteTolerance(std::optional<double> tol)
{
    absolute_tolerance_ = tol.has_value() && tol.value() >= 0.0 ? tol : std::nullopt;
}

std::optional<double> GMGPolar::relativeTolerance() const
{
    return relative_tolerance_;
}
void GMGPolar::relativeTolerance(std::optional<double> tol)
{
    relative_tolerance_ = tol.has_value() && tol.value() >= 0.0 ? tol : std::nullopt;
}

/* ---------------------------------------------------------------------- */
/* Solution & Grid Access                                                 */
/* ---------------------------------------------------------------------- */
Vector<double>& GMGPolar::solution()
{
    int level_depth = 0;
    return levels_[level_depth].solution();
}
const Vector<double>& GMGPolar::solution() const
{
    int level_depth = 0;
    return levels_[level_depth].solution();
}

const PolarGrid& GMGPolar::grid() const
{
    return grid_;
}

/* ---------------------------------------------------------------------- */
/* Diagnostics & statistics                                               */
/* ---------------------------------------------------------------------- */
void GMGPolar::printTimings() const
{
    std::cout << "\n------------------------------" << std::endl;
    std::cout << "----- Timing Information -----" << std::endl;
    std::cout << "------------------------------" << std::endl;
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

int GMGPolar::numberOfIterations() const
{
    return number_of_iterations_;
}
double GMGPolar::meanResidualReductionFactor() const
{
    return mean_residual_reduction_factor_;
}

std::optional<double> GMGPolar::exactErrorWeightedEuclidean() const
{
    if (exact_errors_.empty())
        return std::nullopt;
    return exact_errors_.back().first;
}

std::optional<double> GMGPolar::exactErrorInfinity() const
{
    if (exact_errors_.empty())
        return std::nullopt;
    return exact_errors_.back().second;
}

// -------------------------
// Setup timings
// -------------------------
double GMGPolar::timeSetupTotal() const
{
    return t_setup_total_;
}
double GMGPolar::timeSetupCreateLevels() const
{
    return t_setup_createLevels_;
}
double GMGPolar::timeSetupRHS() const
{
    return t_setup_rhs_;
}
double GMGPolar::timeSetupSmoother() const
{
    return t_setup_smoother_;
}
double GMGPolar::timeSetupDirectSolver() const
{
    return t_setup_directSolver_;
}

// -------------------------
// Solve timings
// -------------------------
double GMGPolar::timeSolveTotal() const
{
    return t_solve_total_;
}
double GMGPolar::timeSolveInitialApproximation() const
{
    return t_solve_initial_approximation_;
}
double GMGPolar::timeSolveMultigridIterations() const
{
    return t_solve_multigrid_iterations_;
}
double GMGPolar::timeCheckConvergence() const
{
    return t_check_convergence_;
}
double GMGPolar::timeCheckExactError() const
{
    return t_check_exact_error_;
}

// -------------------------
// Average Multigrid Cycle timings
// -------------------------
double GMGPolar::timeAvgMGCTotal() const
{
    return t_avg_MGC_total_;
}
double GMGPolar::timeAvgMGCPreSmoothing() const
{
    return t_avg_MGC_preSmoothing_;
}
double GMGPolar::timeAvgMGCPostSmoothing() const
{
    return t_avg_MGC_postSmoothing_;
}
double GMGPolar::timeAvgMGCResidual() const
{
    return t_avg_MGC_residual_;
}
double GMGPolar::timeAvgMGCDirectSolver() const
{
    return t_avg_MGC_directSolver_;
}