#pragma once

#include <chrono>
#include <filesystem>
#include <iostream>
#include <memory>
#include <omp.h>
#include <optional>
#include <utility>

#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/boundaryConditions.h"
#include "../InputFunctions/exactSolution.h"
#include "../InputFunctions/sourceTerm.h"
#include "../Interpolation/interpolation.h"
#include "../LinearAlgebra/Vector/vector.h"
#include "../LinearAlgebra/Vector/vector_operations.h"
#include "../PolarGrid/polargrid.h"
#include "../Definitions/global_definitions.h"

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

    // Preconditioned Conjugate Gradient (PCG) control.
    bool PCG() const;
    void PCG(bool PCG);
    bool PCG_FMG() const;
    void PCG_FMG(bool PCG_FMG);
    int PCG_FMG_iterations() const;
    void PCG_FMG_iterations(int PCG_FMG_iterations);
    MultigridCycleType PCG_FMG_cycle() const;
    void PCG_FMG_cycle(MultigridCycleType PCG_FMG_cycle);
    int PCG_MG_iterations() const;
    void PCG_MG_iterations(int PCG_MG_iterations);
    MultigridCycleType PCG_MG_cycle() const;
    void PCG_MG_cycle(MultigridCycleType PCG_MG_cycle);

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
    double timeConjugateGradient() const;

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
    // PCG settings
    bool PCG_;
    bool PCG_FMG_;
    int PCG_FMG_iterations_;
    MultigridCycleType PCG_FMG_cycle_;
    int PCG_MG_iterations_;
    MultigridCycleType PCG_MG_cycle_;
    // Convergence settings
    int max_iterations_;
    ResidualNormType residual_norm_type_;
    std::optional<double> absolute_tolerance_;
    std::optional<double> relative_tolerance_;

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
    double t_conjugate_gradient_;

    void resetAvgMultigridCycleTimings();
    double t_avg_MGC_total_;
    double t_avg_MGC_preSmoothing_;
    double t_avg_MGC_postSmoothing_;
    double t_avg_MGC_residual_;
    double t_avg_MGC_directSolver_;
};
