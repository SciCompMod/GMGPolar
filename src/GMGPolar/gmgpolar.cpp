#include "../../include/GMGPolar/gmgpolar.h"

#include <iomanip>

GMGPolar::GMGPolar()
    : parser_()
{
    // Initialize LIKWID markers if enabled
    LIKWID_REGISTER("Setup");
    LIKWID_REGISTER("Solve");
    resetAllTimings();
    initializeGrid();
    initializeGeometry();
    initializeMultigrid();
    initializeGeneral();

    setParameters(0, nullptr);
}

GMGPolar::GMGPolar(std::unique_ptr<const DomainGeometry> domain_geometry,
                   std::unique_ptr<const DensityProfileCoefficients> density_profile_coefficients,
                   std::unique_ptr<const BoundaryConditions> boundary_conditions,
                   std::unique_ptr<const SourceTerm> source_term)
    :

    domain_geometry_(std::move(domain_geometry))
    , density_profile_coefficients_(std::move(density_profile_coefficients))
    , boundary_conditions_(std::move(boundary_conditions))
    , source_term_(std::move(source_term))
    , parser_()
{
    LIKWID_REGISTER("Setup");
    LIKWID_REGISTER("Solve");
    resetAllTimings();
    initializeGrid();
    initializeGeometry();
    initializeMultigrid();
    initializeGeneral();

    parseGrid(); /* Removed: parseGeometry(); */
    parseMultigrid();
    parseGeneral();
}

void GMGPolar::setParameters(int argc, char* argv[])
{
    if (argc != 0) {
        try {
            parser_.parse_check(argc, argv);
        }
        catch (const cmdline::cmdline_error& parse_error) {
            std::cerr << "Error: " << parse_error.what() << std::endl;
            std::cerr << "Usage: " << parser_.usage() << std::endl;
        }
    }

    parseGrid();
    parseGeometry();
    parseMultigrid();
    parseGeneral();
}

void GMGPolar::setSolution(std::unique_ptr<const ExactSolution> exact_solution)
{
    exact_solution_ = std::move(exact_solution);
}

/* ----------------- */
/* GMGPolar Solution */
Vector<double>& GMGPolar::solution()
{
    return levels_[0].solution();
}
const Vector<double>& GMGPolar::solution() const
{
    return levels_[0].solution();
}

const PolarGrid& GMGPolar::grid() const
{
    return levels_[0].grid();
}

/* Solve Properties */
int GMGPolar::numberOfIterations() const
{
    return number_of_iterations_;
}
double GMGPolar::meanResidualReductionFactor() const
{
    return mean_residual_reduction_factor_;
}
// Only when exact solution provided
std::optional<double> GMGPolar::exactErrorWeightedEuclidean() const
{
    if (exact_solution_) {
        return exact_errors_.back().first;
    }
    return std::nullopt;
}
std::optional<double> GMGPolar::exactErrorInfinity() const
{
    if (exact_solution_) {
        return exact_errors_.back().second;
    }
    return std::nullopt;
}

void GMGPolar::printTimings() const
{
    std::cout << "\n------------------" << std::endl;
    std::cout << "Timing Information" << std::endl;
    std::cout << "------------------" << std::endl;
    std::cout << "Setup Time: " << t_setup_total_ - t_setup_rhs_ << " seconds" << std::endl;
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
    std::cout << "\n" << std::endl;
}

/* --------------- */
/* Grid Parameters */
double GMGPolar::R0() const
{
    return R0_;
}

void GMGPolar::R0(double R0)
{
    R0_ = R0;
}

double GMGPolar::Rmax() const
{
    return Rmax_;
}

void GMGPolar::Rmax(double Rmax)
{
    Rmax_ = Rmax;
}

int GMGPolar::nr_exp() const
{
    return nr_exp_;
}

void GMGPolar::nr_exp(int nr_exp)
{
    nr_exp_ = nr_exp;
}

int GMGPolar::ntheta_exp() const
{
    return ntheta_exp_;
}

void GMGPolar::ntheta_exp(int ntheta_exp)
{
    ntheta_exp_ = ntheta_exp;
}

int GMGPolar::anisotropic_factor() const
{
    return anisotropic_factor_;
}

void GMGPolar::anisotropic_factor(int anisotropic_factor)
{
    anisotropic_factor_ = anisotropic_factor;
}

int GMGPolar::divideBy2() const
{
    return divideBy2_;
}

void GMGPolar::divideBy2(int divideBy2)
{
    divideBy2_ = divideBy2;
}

/* ------------------- */
/* Geometry Parameters */
bool GMGPolar::DirBC_Interior() const
{
    return DirBC_Interior_;
}

void GMGPolar::DirBC_Interior(bool DirBC_Interior)
{
    DirBC_Interior_ = DirBC_Interior;
}

/* -------------------- */
/* Multigrid Parameters */

bool GMGPolar::FMG() const
{
    return FMG_;
}

void GMGPolar::FMG(bool FMG)
{
    FMG_ = FMG;
}

int GMGPolar::FMG_iterations() const
{
    return FMG_iterations_;
}

void GMGPolar::FMG_iterations(int FMG_iterations)
{
    FMG_iterations_ = FMG_iterations;
}

MultigridCycleType GMGPolar::FMG_cycle() const
{
    return FMG_cycle_;
}

void GMGPolar::FMG_cycle(MultigridCycleType FMG_cycle)
{
    FMG_cycle_ = FMG_cycle;
}

ExtrapolationType GMGPolar::extrapolation() const
{
    return extrapolation_;
}

void GMGPolar::extrapolation(ExtrapolationType extrapolation)
{
    extrapolation_ = extrapolation;
}

int GMGPolar::maxLevels() const
{
    return max_levels_;
}

void GMGPolar::maxLevels(int max_levels)
{
    max_levels_ = max_levels;
}

MultigridCycleType GMGPolar::multigridCycle() const
{
    return multigrid_cycle_;
}

void GMGPolar::multigridCycle(MultigridCycleType multigrid_cycle)
{
    multigrid_cycle_ = multigrid_cycle;
}

int GMGPolar::preSmoothingSteps() const
{
    return pre_smoothing_steps_;
}

void GMGPolar::preSmoothingSteps(int pre_smoothing_steps)
{
    pre_smoothing_steps_ = pre_smoothing_steps;
}

int GMGPolar::postSmoothingSteps() const
{
    return post_smoothing_steps_;
}

void GMGPolar::postSmoothingSteps(int post_smoothing_steps)
{
    post_smoothing_steps_ = post_smoothing_steps;
}

int GMGPolar::maxIterations() const
{
    return max_iterations_;
}

void GMGPolar::maxIterations(int maxIterations)
{
    max_iterations_ = maxIterations;
}

ResidualNormType GMGPolar::residualNormType() const
{
    return residual_norm_type_;
}
void GMGPolar::residualNormType(ResidualNormType residualNormType)
{
    residual_norm_type_ = residualNormType;
}

double GMGPolar::absoluteTolerance() const
{
    if (absolute_tolerance_.has_value()) {
        return absolute_tolerance_.value();
    }
    else {
        return -1.0;
    }
}

void GMGPolar::absoluteTolerance(double absolute_tolerance)
{
    if (absolute_tolerance > 0) {
        absolute_tolerance_ = absolute_tolerance;
    }
    else {
        absolute_tolerance_ = std::nullopt;
    }
}

double GMGPolar::relativeTolerance() const
{
    if (relative_tolerance_.has_value()) {
        return relative_tolerance_.value();
    }
    else {
        return -1.0;
    }
}

void GMGPolar::relativeTolerance(double relative_tolerance)
{
    if (relative_tolerance > 0) {
        relative_tolerance_ = relative_tolerance;
    }
    else {
        relative_tolerance_ = std::nullopt;
    }
}

/* ------------------ */
/* Control Parameters */
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

int GMGPolar::maxOpenMPThreads() const
{
    return max_omp_threads_;
}

void GMGPolar::maxOpenMPThreads(int max_omp_threads)
{
    max_omp_threads_ = max_omp_threads;
}

double GMGPolar::threadReductionFactor() const
{
    return thread_reduction_factor_;
}

void GMGPolar::threadReductionFactor(double thread_reduction_factor)
{
    thread_reduction_factor_ = thread_reduction_factor;
}

StencilDistributionMethod GMGPolar::stencilDistributionMethod() const
{
    return stencil_distribution_method_;
}

void GMGPolar::stencilDistributionMethod(StencilDistributionMethod stencil_distribution_method)
{
    stencil_distribution_method_ = stencil_distribution_method;
}

bool GMGPolar::cacheDensityProfileCoefficients() const
{
    return cache_density_profile_coefficients_;
}

void GMGPolar::cacheDensityProfileCoefficients(bool cache_density_profile_coefficients)
{
    cache_density_profile_coefficients_ = cache_density_profile_coefficients;
}

bool GMGPolar::cacheDomainGeometry() const
{
    return cache_domain_geometry_;
}

void GMGPolar::cacheDomainGeometry(bool cache_domain_geometry)
{
    cache_domain_geometry_ = cache_domain_geometry;
}

/* ---------------------------------------------------------------------- */
/* Setup timings                                                          */
/* ---------------------------------------------------------------------- */
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

/* ---------------------------------------------------------------------- */
/* Solve timings                                                          */
/* ---------------------------------------------------------------------- */
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

/* ---------------------------------------------------------------------- */
/* Average Multigrid Cycle timings                                        */
/* ---------------------------------------------------------------------- */
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

/* ---------------------------------------------------------------------- */
/* Reset timings                                                          */
/* ---------------------------------------------------------------------- */

void GMGPolar::resetAllTimings()
{
    resetSetupPhaseTimings();
    resetSolvePhaseTimings();
    resetAvgMultigridCycleTimings();
}

void GMGPolar::resetSetupPhaseTimings()
{
    t_setup_total_        = 0.0;
    t_setup_createLevels_ = 0.0;
    t_setup_rhs_          = 0.0;
    t_setup_smoother_     = 0.0;
    t_setup_directSolver_ = 0.0;
}

void GMGPolar::resetSolvePhaseTimings()
{
    t_solve_total_                 = 0.0;
    t_solve_initial_approximation_ = 0.0;
    t_solve_multigrid_iterations_  = 0.0;
    t_check_convergence_           = 0.0;
    t_check_exact_error_           = 0.0;
}

void GMGPolar::resetAvgMultigridCycleTimings()
{
    t_avg_MGC_total_         = 0.0;
    t_avg_MGC_preSmoothing_  = 0.0;
    t_avg_MGC_postSmoothing_ = 0.0;
    t_avg_MGC_residual_      = 0.0;
    t_avg_MGC_directSolver_  = 0.0;
}