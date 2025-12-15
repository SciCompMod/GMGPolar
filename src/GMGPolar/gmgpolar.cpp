
/* ---------------------------------------------------------------------- */
/* Constructor & Initialization                                           */
/* ---------------------------------------------------------------------- */
template<DomainGeometryConcept DomainGeometry>
GMGPolar<DomainGeometry>::GMGPolar(const PolarGrid& grid, const DomainGeometry& domain_geometry,
                   const DensityProfileCoefficients& density_profile_coefficients)
    : grid_(grid)
    , domain_geometry_(domain_geometry)
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

template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::setSolution(const ExactSolution* exact_solution)
{
    exact_solution_ = exact_solution;
}

/* ---------------------------------------------------------------------- */
/* General output & visualization                                         */
/* ---------------------------------------------------------------------- */
template<DomainGeometryConcept DomainGeometry>
int GMGPolar<DomainGeometry>::verbose() const
{
    return verbose_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::verbose(int verbose)
{
    verbose_ = verbose;
}

template<DomainGeometryConcept DomainGeometry>
bool GMGPolar<DomainGeometry>::paraview() const
{
    return paraview_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::paraview(bool paraview)
{
    paraview_ = paraview;
}

/* ---------------------------------------------------------------------- */
/* Parallelization & threading                                            */
/* ---------------------------------------------------------------------- */
template<DomainGeometryConcept DomainGeometry>
int GMGPolar<DomainGeometry>::maxOpenMPThreads() const
{
    return max_omp_threads_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::maxOpenMPThreads(int max_omp_threads)
{
    max_omp_threads_ = max_omp_threads;
}

template<DomainGeometryConcept DomainGeometry>
double GMGPolar<DomainGeometry>::threadReductionFactor() const
{
    return thread_reduction_factor_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::threadReductionFactor(double thread_reduction_factor)
{
    thread_reduction_factor_ = thread_reduction_factor;
}

/* ---------------------------------------------------------------------- */
/* Numerical method options                                               */
/* ---------------------------------------------------------------------- */
template<DomainGeometryConcept DomainGeometry>
bool GMGPolar<DomainGeometry>::DirBC_Interior() const
{
    return DirBC_Interior_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::DirBC_Interior(bool DirBC_Interior)
{
    DirBC_Interior_ = DirBC_Interior;
}

template<DomainGeometryConcept DomainGeometry>
StencilDistributionMethod GMGPolar<DomainGeometry>::stencilDistributionMethod() const
{
    return stencil_distribution_method_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::stencilDistributionMethod(StencilDistributionMethod stencil_distribution_method)
{
    stencil_distribution_method_ = stencil_distribution_method;
}

template<DomainGeometryConcept DomainGeometry>
bool GMGPolar<DomainGeometry>::cacheDensityProfileCoefficients() const
{
    return cache_density_profile_coefficients_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::cacheDensityProfileCoefficients(bool cache_density_profile_coefficients)
{
    cache_density_profile_coefficients_ = cache_density_profile_coefficients;
}

template<DomainGeometryConcept DomainGeometry>
bool GMGPolar<DomainGeometry>::cacheDomainGeometry() const
{
    return cache_domain_geometry_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::cacheDomainGeometry(bool cache_domain_geometry)
{
    cache_domain_geometry_ = cache_domain_geometry;
}

/* ---------------------------------------------------------------------- */
/* Multigrid controls                                                     */
/* ---------------------------------------------------------------------- */
template<DomainGeometryConcept DomainGeometry>
ExtrapolationType GMGPolar<DomainGeometry>::extrapolation() const
{
    return extrapolation_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::extrapolation(ExtrapolationType extrapolation)
{
    extrapolation_ = extrapolation;
}

template<DomainGeometryConcept DomainGeometry>
int GMGPolar<DomainGeometry>::maxLevels() const
{
    return max_levels_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::maxLevels(int max_levels)
{
    max_levels_ = max_levels;
}

template<DomainGeometryConcept DomainGeometry>
MultigridCycleType GMGPolar<DomainGeometry>::multigridCycle() const
{
    return multigrid_cycle_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::multigridCycle(MultigridCycleType multigrid_cycle)
{
    multigrid_cycle_ = multigrid_cycle;
}

template<DomainGeometryConcept DomainGeometry>
int GMGPolar<DomainGeometry>::preSmoothingSteps() const
{
    return pre_smoothing_steps_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::preSmoothingSteps(int pre_smoothing_steps)
{
    pre_smoothing_steps_ = pre_smoothing_steps;
}

template<DomainGeometryConcept DomainGeometry>
int GMGPolar<DomainGeometry>::postSmoothingSteps() const
{
    return post_smoothing_steps_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::postSmoothingSteps(int post_smoothing_steps)
{
    post_smoothing_steps_ = post_smoothing_steps;
}

template<DomainGeometryConcept DomainGeometry>
bool GMGPolar<DomainGeometry>::FMG() const
{
    return FMG_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::FMG(bool FMG)
{
    FMG_ = FMG;
}

template<DomainGeometryConcept DomainGeometry>
int GMGPolar<DomainGeometry>::FMG_iterations() const
{
    return FMG_iterations_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::FMG_iterations(int FMG_iterations)
{
    FMG_iterations_ = FMG_iterations;
}

template<DomainGeometryConcept DomainGeometry>
MultigridCycleType GMGPolar<DomainGeometry>::FMG_cycle() const
{
    return FMG_cycle_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::FMG_cycle(MultigridCycleType FMG_cycle)
{
    FMG_cycle_ = FMG_cycle;
}

/* ---------------------------------------------------------------------- */
/* Iterative solver termination                                           */
/* ---------------------------------------------------------------------- */

template<DomainGeometryConcept DomainGeometry>
int GMGPolar<DomainGeometry>::maxIterations() const
{
    return max_iterations_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::maxIterations(int maxIterations)
{
    max_iterations_ = maxIterations;
}

template<DomainGeometryConcept DomainGeometry>
ResidualNormType GMGPolar<DomainGeometry>::residualNormType() const
{
    return residual_norm_type_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::residualNormType(ResidualNormType residualNormType)
{
    residual_norm_type_ = residualNormType;
}

template<DomainGeometryConcept DomainGeometry>
std::optional<double> GMGPolar<DomainGeometry>::absoluteTolerance() const
{
    return absolute_tolerance_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::absoluteTolerance(std::optional<double> tol)
{
    absolute_tolerance_ = tol.has_value() && tol.value() >= 0.0 ? tol : std::nullopt;
}

template<DomainGeometryConcept DomainGeometry>
std::optional<double> GMGPolar<DomainGeometry>::relativeTolerance() const
{
    return relative_tolerance_;
}
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::relativeTolerance(std::optional<double> tol)
{
    relative_tolerance_ = tol.has_value() && tol.value() >= 0.0 ? tol : std::nullopt;
}

/* ---------------------------------------------------------------------- */
/* Solution & Grid Access                                                 */
/* ---------------------------------------------------------------------- */
template<DomainGeometryConcept DomainGeometry>
Vector<double> GMGPolar<DomainGeometry>::solution()
{
    int level_depth = 0;
    return levels_[level_depth].solution();
}
template<DomainGeometryConcept DomainGeometry>
ConstVector<double> GMGPolar<DomainGeometry>::solution() const
{
    int level_depth = 0;
    return levels_[level_depth].solution();
}

template<DomainGeometryConcept DomainGeometry>
const PolarGrid& GMGPolar<DomainGeometry>::grid() const
{
    return grid_;
}

/* ---------------------------------------------------------------------- */
/* Setup timings                                                          */
/* ---------------------------------------------------------------------- */
template<DomainGeometryConcept DomainGeometry>
double GMGPolar<DomainGeometry>::timeSetupTotal() const
{
    return t_setup_total_;
}
template<DomainGeometryConcept DomainGeometry>
double GMGPolar<DomainGeometry>::timeSetupCreateLevels() const
{
    return t_setup_createLevels_;
}
template<DomainGeometryConcept DomainGeometry>
double GMGPolar<DomainGeometry>::timeSetupRHS() const
{
    return t_setup_rhs_;
}
template<DomainGeometryConcept DomainGeometry>
double GMGPolar<DomainGeometry>::timeSetupSmoother() const
{
    return t_setup_smoother_;
}
template<DomainGeometryConcept DomainGeometry>
double GMGPolar<DomainGeometry>::timeSetupDirectSolver() const
{
    return t_setup_directSolver_;
}

/* ---------------------------------------------------------------------- */
/* Solve timings                                                          */
/* ---------------------------------------------------------------------- */
template<DomainGeometryConcept DomainGeometry>
double GMGPolar<DomainGeometry>::timeSolveTotal() const
{
    return t_solve_total_;
}
template<DomainGeometryConcept DomainGeometry>
double GMGPolar<DomainGeometry>::timeSolveInitialApproximation() const
{
    return t_solve_initial_approximation_;
}
template<DomainGeometryConcept DomainGeometry>
double GMGPolar<DomainGeometry>::timeSolveMultigridIterations() const
{
    return t_solve_multigrid_iterations_;
}
template<DomainGeometryConcept DomainGeometry>
double GMGPolar<DomainGeometry>::timeCheckConvergence() const
{
    return t_check_convergence_;
}
template<DomainGeometryConcept DomainGeometry>
double GMGPolar<DomainGeometry>::timeCheckExactError() const
{
    return t_check_exact_error_;
}

/* ---------------------------------------------------------------------- */
/* Average Multigrid Cycle timings                                        */
/* ---------------------------------------------------------------------- */
template<DomainGeometryConcept DomainGeometry>
double GMGPolar<DomainGeometry>::timeAvgMGCTotal() const
{
    return t_avg_MGC_total_;
}
template<DomainGeometryConcept DomainGeometry>
double GMGPolar<DomainGeometry>::timeAvgMGCPreSmoothing() const
{
    return t_avg_MGC_preSmoothing_;
}
template<DomainGeometryConcept DomainGeometry>
double GMGPolar<DomainGeometry>::timeAvgMGCPostSmoothing() const
{
    return t_avg_MGC_postSmoothing_;
}
template<DomainGeometryConcept DomainGeometry>
double GMGPolar<DomainGeometry>::timeAvgMGCResidual() const
{
    return t_avg_MGC_residual_;
}
template<DomainGeometryConcept DomainGeometry>
double GMGPolar<DomainGeometry>::timeAvgMGCDirectSolver() const
{
    return t_avg_MGC_directSolver_;
}

/* ---------------------------------------------------------------------- */
/* Reset timings                                                          */
/* ---------------------------------------------------------------------- */

template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::resetAllTimings()
{
    resetSetupPhaseTimings();
    resetSolvePhaseTimings();
    resetAvgMultigridCycleTimings();
}

template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::resetSetupPhaseTimings()
{
    t_setup_total_        = 0.0;
    t_setup_createLevels_ = 0.0;
    t_setup_rhs_          = 0.0;
    t_setup_smoother_     = 0.0;
    t_setup_directSolver_ = 0.0;
}

template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::resetSolvePhaseTimings()
{
    t_solve_total_                 = 0.0;
    t_solve_initial_approximation_ = 0.0;
    t_solve_multigrid_iterations_  = 0.0;
    t_check_convergence_           = 0.0;
    t_check_exact_error_           = 0.0;
}

template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::resetAvgMultigridCycleTimings()
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
template<DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::printTimings() const
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
template<DomainGeometryConcept DomainGeometry>
int GMGPolar<DomainGeometry>::numberOfIterations() const
{
    return number_of_iterations_;
}

// Mean residual reduction factor per iteration.
template<DomainGeometryConcept DomainGeometry>
double GMGPolar<DomainGeometry>::meanResidualReductionFactor() const
{
    return mean_residual_reduction_factor_;
}

// Error norms (only available if exact solution was set).
template<DomainGeometryConcept DomainGeometry>
std::optional<double> GMGPolar<DomainGeometry>::exactErrorWeightedEuclidean() const
{
    if (exact_solution_) {
        return exact_errors_.back().first;
    }
    return std::nullopt;
}
template<DomainGeometryConcept DomainGeometry>
std::optional<double> GMGPolar<DomainGeometry>::exactErrorInfinity() const
{
    if (exact_solution_) {
        return exact_errors_.back().second;
    }
    return std::nullopt;
}
