#pragma once

#include <filesystem>
#include <iostream>
#include <omp.h>
#include <utility>

#include "../InputFunctions/densityProfileCoefficients.h"
#include "../InputFunctions/domainGeometry.h"

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
    {
    }

    /* ---------------------------------------------------------------------- */
    /* Setup & Solve                                                          */
    /* ---------------------------------------------------------------------- */
    // Finalize solver setup (allocate data, build operators, etc.).
    void setup();

    // Solve system with given boundary conditions and source term.
    // Multiple solves with different inputs are supported.
    template <concepts::BoundaryConditions BoundaryConditions>
    void solve(const BoundaryConditions& boundary_conditions, const SourceTerm& source_term);

    // Return the computed solution vector (hides IGMGPolar::solution()).
    Vector<double> solution();
    ConstVector<double> solution() const;

private:
    /* ------------------------------------ */
    /* Grid Configuration & Input Functions */
    /* ------------------------------------ */
    const DomainGeometry& domain_geometry_;
    const DensityProfileCoefficients& density_profile_coefficients_;

    /* ---------------- */
    /* Multigrid levels */
    std::vector<Level<DomainGeometry>> levels_;

    /* ---------------------------------------------------- */
    /* Compute exact error if an exact solution is provided */
    // The results are stored as a pair: (weighted L2 error, infinity error).
    std::pair<double, double> computeExactError(Level<DomainGeometry>& level, ConstVector<double> solution, Vector<double> error,
                                                const ExactSolution& exact_solution);

    /* --------------- */
    /* Setup Functions */
    void discretize_rhs_f(const Level<DomainGeometry>& level, Vector<double> rhs_f);
    template <concepts::BoundaryConditions BoundaryConditions>
    void build_rhs_f(const Level<DomainGeometry>& level, Vector<double> rhs_f, const BoundaryConditions& boundary_conditions,
                     const SourceTerm& source_term);

    /* --------------- */
    /* Solve Functions */
    void initializeSolution();

    double residualNorm(const ResidualNormType& norm_type, const Level<DomainGeometry>& level, ConstVector<double> residual) const;
    void evaluateExactError(Level<DomainGeometry>& level, const ExactSolution& exact_solution);
    void updateResidualNorms(Level<DomainGeometry>& level, int iteration, double& initial_residual_norm, double& current_residual_norm,
                             double& current_relative_residual_norm);
    void extrapolatedResidual(int current_level, Vector<double> residual, ConstVector<double> residual_next_level);
    void fullMultigridApproximation(MultigridCycleType FMG_cycle, int FMG_iterations);
    void initRhsHierarchy(Vector<double> rhs);
    void solveMultigrid(double& initial_residual_norm, double& current_residual_norm,
                        double& current_relative_residual_norm);
    void solvePCG(double& initial_residual_norm, double& current_residual_norm,
                  double& current_relative_residual_norm);
    void applyMultigridIterations(Level<DomainGeometry>& level, MultigridCycleType cycle, int iterations);
    void applyExtrapolatedMultigridIterations(Level<DomainGeometry>& level, MultigridCycleType cycle, int iterations);

    /* ------------------- */
    /* Multigrid Functions */
    void multigrid_V_Cycle(int level_depth, Vector<double> solution, Vector<double> rhs, Vector<double> residual);
    void multigrid_W_Cycle(int level_depth, Vector<double> solution, Vector<double> rhs, Vector<double> residual);
    void multigrid_F_Cycle(int level_depth, Vector<double> solution, Vector<double> rhs, Vector<double> residual);
    void extrapolated_multigrid_V_Cycle(int level_depth, Vector<double> solution, Vector<double> rhs,
                                        Vector<double> residual);
    void extrapolated_multigrid_W_Cycle(int level_depth, Vector<double> solution, Vector<double> rhs,
                                        Vector<double> residual);
    void extrapolated_multigrid_F_Cycle(int level_depth, Vector<double> solution, Vector<double> rhs,
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
    void writeToVTK(const std::filesystem::path& file_path, const Level<DomainGeometry>& level,
                    ConstVector<double> grid_function);
};

#include "build_rhs_f.h"
#include "setup.h"
#include "solver.inl"
#include "solver.h"
#include "writeToVTK.h"
