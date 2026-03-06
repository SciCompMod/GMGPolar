#pragma once

#include <memory>
#include <optional>

#include "../GMGPolar/gmgpolar.h"
#include "../ConfigParser/test_selection.h"

// GMGPolarDispatcher is an application-layer helper for command-line tools.
// It wraps a GMGPolar<DG,DC> solver (created via double std::visit over DG x DC
// variants) and exposes a variant-based solve() so callers don't need to know
// the concrete geometry and coefficient types.
//
// Library users who know their types at compile time should use GMGPolar<DG,DC>
// directly; GMGPolarDispatcher is intended only for driver programs that receive
// geometry/problem selections at runtime (e.g. main.cpp, scaling benchmarks).
class GMGPolarDispatcher
{
public:
    GMGPolarDispatcher(const PolarGrid& grid, const DomainGeometryVariant& domain_geometry,
                       const DensityProfileCoefficientsVariant& density_profile_coefficients)
    {
        solver_ = std::visit(
            [&grid](auto const& dg, auto const& dc) -> std::unique_ptr<IGMGPolar> {
                using DG = std::decay_t<decltype(dg)>;
                using DC = std::decay_t<decltype(dc)>;
                return std::make_unique<GMGPolar<DG, DC>>(grid, dg, dc);
            },
            domain_geometry, density_profile_coefficients);
    }

    /* ---------------------------------------------------------------------- */
    /* General output & visualization                                         */
    /* ---------------------------------------------------------------------- */
    int verbose() const { return solver_->verbose(); }
    void verbose(int v) { solver_->verbose(v); }

    bool paraview() const { return solver_->paraview(); }
    void paraview(bool p) { solver_->paraview(p); }

    /* ---------------------------------------------------------------------- */
    /* Parallelization                                                        */
    /* ---------------------------------------------------------------------- */
    int maxOpenMPThreads() const { return solver_->maxOpenMPThreads(); }
    void maxOpenMPThreads(int n) { solver_->maxOpenMPThreads(n); }

    /* ---------------------------------------------------------------------- */
    /* Numerical method options                                               */
    /* ---------------------------------------------------------------------- */
    bool DirBC_Interior() const { return solver_->DirBC_Interior(); }
    void DirBC_Interior(bool b) { solver_->DirBC_Interior(b); }

    StencilDistributionMethod stencilDistributionMethod() const
    {
        return solver_->stencilDistributionMethod();
    }
    void stencilDistributionMethod(StencilDistributionMethod m) { solver_->stencilDistributionMethod(m); }

    bool cacheDensityProfileCoefficients() const { return solver_->cacheDensityProfileCoefficients(); }
    void cacheDensityProfileCoefficients(bool b) { solver_->cacheDensityProfileCoefficients(b); }

    bool cacheDomainGeometry() const { return solver_->cacheDomainGeometry(); }
    void cacheDomainGeometry(bool b) { solver_->cacheDomainGeometry(b); }

    /* ---------------------------------------------------------------------- */
    /* Multigrid controls                                                     */
    /* ---------------------------------------------------------------------- */
    ExtrapolationType extrapolation() const { return solver_->extrapolation(); }
    void extrapolation(ExtrapolationType e) { solver_->extrapolation(e); }

    int maxLevels() const { return solver_->maxLevels(); }
    void maxLevels(int n) { solver_->maxLevels(n); }

    MultigridCycleType multigridCycle() const { return solver_->multigridCycle(); }
    void multigridCycle(MultigridCycleType c) { solver_->multigridCycle(c); }

    int preSmoothingSteps() const { return solver_->preSmoothingSteps(); }
    void preSmoothingSteps(int n) { solver_->preSmoothingSteps(n); }

    int postSmoothingSteps() const { return solver_->postSmoothingSteps(); }
    void postSmoothingSteps(int n) { solver_->postSmoothingSteps(n); }

    bool FMG() const { return solver_->FMG(); }
    void FMG(bool b) { solver_->FMG(b); }

    int FMG_iterations() const { return solver_->FMG_iterations(); }
    void FMG_iterations(int n) { solver_->FMG_iterations(n); }

    MultigridCycleType FMG_cycle() const { return solver_->FMG_cycle(); }
    void FMG_cycle(MultigridCycleType c) { solver_->FMG_cycle(c); }

    /* ---------------------------------------------------------------------- */
    /* Iterative solver termination                                           */
    /* ---------------------------------------------------------------------- */
    int maxIterations() const { return solver_->maxIterations(); }
    void maxIterations(int n) { solver_->maxIterations(n); }

    ResidualNormType residualNormType() const { return solver_->residualNormType(); }
    void residualNormType(ResidualNormType t) { solver_->residualNormType(t); }

    std::optional<double> absoluteTolerance() const { return solver_->absoluteTolerance(); }
    void absoluteTolerance(std::optional<double> tol) { solver_->absoluteTolerance(tol); }

    std::optional<double> relativeTolerance() const { return solver_->relativeTolerance(); }
    void relativeTolerance(std::optional<double> tol) { solver_->relativeTolerance(tol); }

    /* ---------------------------------------------------------------------- */
    /* Setup & Solve                                                          */
    /* ---------------------------------------------------------------------- */
    void setup() { solver_->setup(); }

    void setSolution(const ExactSolution* exact_solution) { solver_->setSolution(exact_solution); }

    // Solve with variant-based boundary conditions (variant types stay in the app layer).
    void solve(const BoundaryConditionsVariant& boundary_conditions, const SourceTerm& source_term)
    {
        std::visit([&](auto const& bc) { solver_->solve(bc, source_term); }, boundary_conditions);
    }

    /* ---------------------------------------------------------------------- */
    /* Solution & Grid Access                                                 */
    /* ---------------------------------------------------------------------- */
    Vector<double> solution() { return solver_->solution(); }
    ConstVector<double> solution() const { return solver_->solution(); }
    const PolarGrid& grid() const { return solver_->grid(); }

    /* ---------------------------------------------------------------------- */
    /* Diagnostics & statistics                                               */
    /* ---------------------------------------------------------------------- */
    void printTimings() const { solver_->printTimings(); }
    int numberOfIterations() const { return solver_->numberOfIterations(); }
    double meanResidualReductionFactor() const { return solver_->meanResidualReductionFactor(); }
    std::optional<double> exactErrorWeightedEuclidean() const { return solver_->exactErrorWeightedEuclidean(); }
    std::optional<double> exactErrorInfinity() const { return solver_->exactErrorInfinity(); }

private:
    std::unique_ptr<IGMGPolar> solver_;
};
