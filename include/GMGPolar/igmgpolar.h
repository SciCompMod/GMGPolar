#pragma once

#include <chrono>
#include <filesystem>
#include <iostream>
#include <memory>
#include <omp.h>
#include <optional>
#include <utility>

#include "../InputFunctions/domainGeometry.h"

template <concepts::DomainGeometry DomainGeometry>
class LevelCache;

template <concepts::DomainGeometry DomainGeometry>
class Level;

#include "../InputFunctions/boundaryConditions.h"
#include "../InputFunctions/densityProfileCoefficients.h"
#include "../InputFunctions/exactSolution.h"
#include "../InputFunctions/sourceTerm.h"
#include "../Interpolation/interpolation.h"
#include "../Level/level.h"
#include "../LinearAlgebra/Vector/vector.h"
#include "../LinearAlgebra/Vector/vector_operations.h"
#include "../PolarGrid/polargrid.h"
#include "../Definitions/global_definitions.h"
#include "test_cases.h"

class IGMGPolar
{
public:
    /* ---------------------------------------------------------------------- */
    /* Constructor & Initialization                                           */
    /* ---------------------------------------------------------------------- */
    // Construct a polar PDE multigrid solver for the Poisson-like equation:
    // - \nabla \cdot (\alpha \nabla u) + \beta u = rhs_f  in \Omega,
    // with Dirichlet boundary condition        u = u_D    on \partial \Omega.
    // Parameters:
    // - grid: Cartesian mesh discretizing the computational domain.
    IGMGPolar(const PolarGrid& grid);

    /* ---------------------------------------------------------------------- */
    /* General output & visualization                                         */
    /* ---------------------------------------------------------------------- */
    // Verbose level (0 = silent, >0 = increasingly detailed logs).
    int verbose() const;
    void verbose(int verbose);

    // Enable/disable ParaView output for visualization.
    bool paraview() const;
    void paraview(bool paraview);

    /* ---------------------------------------------------------------------- */
    /* Parallelization                                                        */
    /* ---------------------------------------------------------------------- */
    // Maximum number of OpenMP threads to use.
    int maxOpenMPThreads() const;
    void maxOpenMPThreads(int max_omp_threads);

    /* ---------------------------------------------------------------------- */
    /* Numerical method options                                               */
    /* ---------------------------------------------------------------------- */
    // Treatment of the interior boundary at the origin:
    //     true  -> use Dirichlet boundary condition
    //     false -> use discretization across the origin
    bool DirBC_Interior() const;
    void DirBC_Interior(bool DirBC_Interior);

    // Strategy for distributing the stencil (Take, Give).
    StencilDistributionMethod stencilDistributionMethod() const;
    void stencilDistributionMethod(StencilDistributionMethod stencil_distribution_method);

    // Cache density profile coefficients (alpha, beta).
    bool cacheDensityProfileCoefficients() const;
    void cacheDensityProfileCoefficients(bool cache_density_profile_coefficients);

    // Cache domain geometry data (arr, att, art, detDF).
    bool cacheDomainGeometry() const;
    void cacheDomainGeometry(bool cache_domain_geometry);

    /* ---------------------------------------------------------------------- */
    /* Multigrid controls                                                     */
    /* ---------------------------------------------------------------------- */
    // Implicit extrapolation to increase the order of convergence
    ExtrapolationType extrapolation() const;
    void extrapolation(ExtrapolationType extrapolation);

    // Maximum multigrid levels (-1 = use deepest possible).
    int maxLevels() const;
    void maxLevels(int max_levels);

    // Multigrid cycle type (V-cycle, W-cycle, F-cycle).
    MultigridCycleType multigridCycle() const;
    void multigridCycle(MultigridCycleType multigrid_cycle);

    // Pre-/post-smoothing steps per level.
    int preSmoothingSteps() const;
    void preSmoothingSteps(int pre_smoothing_steps);
    int postSmoothingSteps() const;
    void postSmoothingSteps(int post_smoothing_steps);

    // Full Multigrid (FMG) control.
    bool FMG() const;
    void FMG(bool FMG);
    int FMG_iterations() const;
    void FMG_iterations(int FMG_iterations);
    MultigridCycleType FMG_cycle() const;
    void FMG_cycle(MultigridCycleType FMG_cycle);

    /* ---------------------------------------------------------------------- */
    /* Iterative solver termination                                           */
    /* ---------------------------------------------------------------------- */
    // Maximum number of iterations for the solver.
    int maxIterations() const;
    void maxIterations(int max_iterations);

    // Type of residual norm used to check convergence.
    ResidualNormType residualNormType() const;
    void residualNormType(ResidualNormType residual_norm_type);

    // Absolute residual tolerance for convergence.
    std::optional<double> absoluteTolerance() const;
    void absoluteTolerance(std::optional<double> tol);

    // Relative residual tolerance (relative to initial residual).
    std::optional<double> relativeTolerance() const;
    void relativeTolerance(std::optional<double> tol);

    // If an exact solution is provided, the solver will compute the exact error at each iteration.
    void setSolution(const ExactSolution* exact_solution);

    /* ---------------------------------------------------------------------- */
    /* Solution & Grid Access                                                 */
    /* ---------------------------------------------------------------------- */
    // Return a reference to the computed solution vector.
    Vector<double> solution();
    ConstVector<double> solution() const;

    // Return the underlying cartesian mesh used for discretization.
    const PolarGrid& grid() const;

    /* ---------------------------------------------------------------------- */
    /* Diagnostics & statistics                                               */
    /* ---------------------------------------------------------------------- */
    // Print timing breakdown for setup, smoothing, coarse solve, etc.
    void printTimings() const;

    // Number of iterations taken by last solve.
    int numberOfIterations() const;

    // Mean residual reduction factor per iteration.
    double meanResidualReductionFactor() const;

    // Error norms (only available if exact solution was set).
    std::optional<double> exactErrorWeightedEuclidean() const;
    std::optional<double> exactErrorInfinity() const;

    /* ---------------------------------------------------------------------- */
    /* Timing Statistics                                                      */
    /* ---------------------------------------------------------------------- */
    double timeSetupTotal() const;
    double timeSetupCreateLevels() const;
    double timeSetupRHS() const;
    double timeSetupSmoother() const;
    double timeSetupDirectSolver() const;

    double timeSolveTotal() const;
    double timeSolveInitialApproximation() const;
    double timeSolveMultigridIterations() const;
    double timeCheckConvergence() const;
    double timeCheckExactError() const;

    double timeAvgMGCTotal() const;
    double timeAvgMGCPreSmoothing() const;
    double timeAvgMGCPostSmoothing() const;
    double timeAvgMGCResidual() const;
    double timeAvgMGCDirectSolver() const;

protected:
    /* ------------------------------------ */
    /* Grid Configuration & Input Functions */
    /* ------------------------------------ */
    const PolarGrid& grid_;
    const ExactSolution* exact_solution_; // Optional exact solution for validation

    /* ------------------ */
    /* Control Parameters */
    /* ------------------ */
    // General solver output and visualization settings
    int verbose_;
    bool paraview_;
    // Parallelization and threading settings
    int max_omp_threads_;
    // Numerical method setup
    bool DirBC_Interior_;
    StencilDistributionMethod stencil_distribution_method_;
    bool cache_density_profile_coefficients_;
    bool cache_domain_geometry_;
    // Multigrid settings
    ExtrapolationType extrapolation_;
    int max_levels_;
    int pre_smoothing_steps_;
    int post_smoothing_steps_;
    MultigridCycleType multigrid_cycle_;
    // FMG settings
    bool FMG_;
    int FMG_iterations_;
    MultigridCycleType FMG_cycle_;
    // Convergence settings
    int max_iterations_;
    ResidualNormType residual_norm_type_;
    std::optional<double> absolute_tolerance_;
    std::optional<double> relative_tolerance_;

    /* ---------------- */
    /* Multigrid levels */
    int number_of_levels_;

    std::vector<std::pair<double, double>> exact_errors_;

    /* ---------------------- */
    /* Interpolation operator */
    std::unique_ptr<Interpolation> interpolation_;

    /* ------------------------------------------------------------------------- */
    /* Chooses if full grid smoothing is active on level 0 for extrapolation > 0 */
    bool full_grid_smoothing_ = false;

    /* -------------------- */
    /* Convergence criteria */
    int number_of_iterations_;
    std::vector<double> residual_norms_;
    double mean_residual_reduction_factor_;
    bool converged(double current_residual_norm, double first_residual_norm);

    /* ------------------------------------------------------------------------- */
    /* Compute the extrapolated residual: res_ex = 4/3 res_fine - 1/3 res_coarse */
    void extrapolatedResidual(int current_level, Vector<double> residual, ConstVector<double> residual_next_level);

    /* --------------- */
    /* Setup Functions */
    int chooseNumberOfLevels(const PolarGrid& finest_grid);

    /* ----------------- */
    /* Print information */
    void printSettings(const PolarGrid& finest_grid, const PolarGrid& coarsest_grid) const;
    void printIterationHeader(const ExactSolution* exact_solution);
    void printIterationInfo(int iteration, double current_residual_norm, double current_relative_residual_norm,
                            const ExactSolution* exact_solution);

    /* ------------------------------ */
    /* Timing statistics for GMGPolar */
    void resetAllTimings();

    void resetSetupPhaseTimings();
    double t_setup_total_;
    double t_setup_createLevels_;
    double t_setup_rhs_;
    double t_setup_smoother_;
    double t_setup_directSolver_;

    void resetSolvePhaseTimings();
    double t_solve_total_;
    double t_solve_initial_approximation_;
    double t_solve_multigrid_iterations_;
    double t_check_convergence_;
    double t_check_exact_error_;

    void resetAvgMultigridCycleTimings();
    double t_avg_MGC_total_;
    double t_avg_MGC_preSmoothing_;
    double t_avg_MGC_postSmoothing_;
    double t_avg_MGC_residual_;
    double t_avg_MGC_directSolver_;
};
