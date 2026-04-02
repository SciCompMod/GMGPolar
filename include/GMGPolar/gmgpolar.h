#pragma once

#include <filesystem>
#include <iostream>
#include <omp.h>
#include <utility>

#include "../InputFunctions/densityProfileCoefficients.h"
#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/sourceTerm.h"
#include "../Level/level.h"

#include "igmgpolar.h"

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
class GMGPolar : public IGMGPolar
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
    // - domain_geometry: Mapping from the reference domain to the physical domain \Omega.
    // - density_profile_coefficients: Coefficients \alpha and \beta defining the PDE.
    GMGPolar(const PolarGrid& grid, const DomainGeometry& domain_geometry,
             const DensityProfileCoefficients& density_profile_coefficients)
        : IGMGPolar(grid)
        , domain_geometry_(domain_geometry)
        , density_profile_coefficients_(density_profile_coefficients)
        , exact_solution_(nullptr)
        // Level management and internal solver data
        , number_of_levels_(0)
        , interpolation_(nullptr)
        , full_grid_smoothing_(false)
        , number_of_iterations_(0)
        , mean_residual_reduction_factor_(1.0)
    {
    }

    // If an exact solution is provided, the solver will compute the exact error at each iteration.
    void setSolution(const ExactSolution* exact_solution);

    /* ---------------------------------------------------------------------- */
    /* Setup & Solve                                                          */
    /* ---------------------------------------------------------------------- */
    // Finalize solver setup (allocate data, build operators, etc.).
    void setup();

    // Solve system with given boundary conditions and source term.
    // Multiple solves with different inputs are supported.
    template <concepts::BoundaryConditions BoundaryConditions, concepts::SourceTerm SourceTerm>
    void solve(const BoundaryConditions& boundary_conditions, const SourceTerm& source_term);

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

private:
    /* ------------------------------------ */
    /* Grid Configuration & Input Functions */
    /* ------------------------------------ */
    const DomainGeometry& domain_geometry_;
    const DensityProfileCoefficients& density_profile_coefficients_;
    const ExactSolution* exact_solution_; // Optional exact solution for validation

    /* ---------------- */
    /* Multigrid levels */
    int number_of_levels_;
    std::vector<Level<DomainGeometry, DensityProfileCoefficients>> levels_;

    /* ---------------------- */
    /* Interpolation operator */
    std::unique_ptr<Interpolation> interpolation_;

    /* ------------------------------------------------------------------------- */
    /* Chooses if full grid smoothing is active on level 0 for extrapolation > 0 */
    bool full_grid_smoothing_ = false;

    /* -------------------------------------------------- */
    /* Vectors for PCG (Preconditioned Conjugate Gradient)
    * https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method
    *
    * Dedicated vectors:
    *   x  (solution)            -> pcg_solution_
    *   p  (search direction)    -> pcg_search_direction_
    *
    * Reused vectors (to avoid extra allocations):
    *   r    (residual)                       -> level.rhs()
    *   z    (preconditioned residual)        -> level.solution()
    *   A*p  (matrix applied to search dir.)  -> level.residual()
    */
    AllocatableVector<double> pcg_solution_; // x (solution)
    AllocatableVector<double> pcg_search_direction_; // p (search direction)

    /* -------------------- */
    /* Convergence criteria */
    int number_of_iterations_;
    std::vector<double> residual_norms_;
    double mean_residual_reduction_factor_;
    bool converged(double current_residual_norm, double first_residual_norm);

    /* ---------------------------------------------------- */
    /* Compute exact error if an exact solution is provided */
    // The results are stored as a pair: (weighted L2 error, infinity error).
    std::vector<std::pair<double, double>> exact_errors_;
    std::pair<double, double> computeExactError(Level<DomainGeometry, DensityProfileCoefficients>& level,
                                                ConstVector<double> solution, Vector<double> error,
                                                const ExactSolution& exact_solution);

    /* ------------------------------------------------------------------------- */
    /* Compute the extrapolated residual: res_ex = 4/3 res_fine - 1/3 res_coarse */
    void extrapolatedResidual(int current_level, Vector<double> residual, ConstVector<double> residual_next_level);

    /* --------------- */
    /* Setup Functions */
    int chooseNumberOfLevels(const PolarGrid& finest_grid);
    template <concepts::BoundaryConditions BoundaryConditions, concepts::SourceTerm SourceTerm>
    void build_rhs_f(const Level<DomainGeometry, DensityProfileCoefficients>& level, Vector<double> rhs_f,
                     const BoundaryConditions& boundary_conditions, const SourceTerm& source_term);
    void discretize_rhs_f(const Level<DomainGeometry, DensityProfileCoefficients>& level, Vector<double> rhs_f);
    bool checkUniformRefinement(const PolarGrid& grid, double tolerance) const;

    /* --------------- */
    /* Solve Functions */
    void fullMultigridApproximation(MultigridCycleType FMG_cycle, int FMG_iterations);
    void solveMultigrid(double& initial_residual_norm, double& current_residual_norm,
                        double& current_relative_residual_norm);
    void solvePCG(double& initial_residual_norm, double& current_residual_norm, double& current_relative_residual_norm);
    double residualNorm(const ResidualNormType& norm_type,
                        const Level<DomainGeometry, DensityProfileCoefficients>& level,
                        ConstVector<double> residual) const;
    void evaluateExactError(Level<DomainGeometry, DensityProfileCoefficients>& level,
                            const ExactSolution& exact_solution);
    void updateResidualNorms(Level<DomainGeometry, DensityProfileCoefficients>& level, int iteration,
                             double& initial_residual_norm, double& current_residual_norm,
                             double& current_relative_residual_norm);
    void initRhsHierarchy(Vector<double> rhs);
    void applyMultigridIterations(Level<DomainGeometry, DensityProfileCoefficients>& level, MultigridCycleType cycle,
                                  int iterations);
    void applyExtrapolatedMultigridIterations(Level<DomainGeometry, DensityProfileCoefficients>& level,
                                              MultigridCycleType cycle, int iterations);

    /* ----------------- */
    /* Print information */
    void printSettings(const PolarGrid& finest_grid, const PolarGrid& coarsest_grid) const;
    void printIterationHeader(const ExactSolution* exact_solution);
    void printIterationInfo(int iteration, double current_residual_norm, double current_relative_residual_norm,
                            const ExactSolution* exact_solution);

    /* ------------------- */
    /* Multigrid Functions */
    void multigrid_V_Cycle(int level_depth, Vector<double> solution, ConstVector<double> rhs, Vector<double> residual);
    void multigrid_W_Cycle(int level_depth, Vector<double> solution, ConstVector<double> rhs, Vector<double> residual);
    void multigrid_F_Cycle(int level_depth, Vector<double> solution, ConstVector<double> rhs, Vector<double> residual);
    void extrapolated_multigrid_V_Cycle(int level_depth, Vector<double> solution, ConstVector<double> rhs,
                                        Vector<double> residual);
    void extrapolated_multigrid_W_Cycle(int level_depth, Vector<double> solution, ConstVector<double> rhs,
                                        Vector<double> residual);
    void extrapolated_multigrid_F_Cycle(int level_depth, Vector<double> solution, ConstVector<double> rhs,
                                        Vector<double> residual);

    /* ----------------------- */
    /* Interpolation functions */
    void prolongation(int current_level, Vector<double> result, ConstVector<double> x) const;
    void restriction(int current_level, Vector<double> result, ConstVector<double> x) const;
    void injection(int current_level, Vector<double> result, ConstVector<double> x) const;
    void extrapolatedProlongation(int current_level, Vector<double> result, ConstVector<double> x) const;
    void extrapolatedRestriction(int current_level, Vector<double> result, ConstVector<double> x) const;
    void FMGInterpolation(int current_level, Vector<double> result, ConstVector<double> x) const;

    /* ------------- */
    /* Visualization */
    void writeToVTK(const std::filesystem::path& file_path, const PolarGrid& grid);
    void writeToVTK(const std::filesystem::path& file_path,
                    const Level<DomainGeometry, DensityProfileCoefficients>& level, ConstVector<double> grid_function);
};

#include "utils.h"
#include "setup.h"
#include "solver.h"
#include "MultigridMethods/extrapolated_multigrid_F_Cycle.h"
#include "MultigridMethods/extrapolated_multigrid_W_Cycle.h"
#include "MultigridMethods/multigrid_V_Cycle.h"
#include "MultigridMethods/extrapolated_multigrid_V_Cycle.h"
#include "MultigridMethods/multigrid_F_Cycle.h"
#include "MultigridMethods/multigrid_W_Cycle.h"
