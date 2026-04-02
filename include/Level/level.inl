#pragma once
#include "../../include/Residual/ResidualGive/residualGive.h"
#include "../../include/Residual/ResidualTake/residualTake.h"

#include "../../include/DirectSolver/DirectSolver-COO-MUMPS-Give/directSolverGive.h"
#include "../../include/DirectSolver/DirectSolver-COO-MUMPS-Take/directSolverTake.h"
#include "../../include/DirectSolver/DirectSolver-CSR-LU-Give/directSolverGiveCustomLU.h"
#include "../../include/DirectSolver/DirectSolver-CSR-LU-Take/directSolverTakeCustomLU.h"

#include "../../include/Smoother/SmootherGive/smootherGive.h"
#include "../../include/Smoother/SmootherTake/smootherTake.h"

#include "../../include/ExtrapolatedSmoother/ExtrapolatedSmootherGive/extrapolatedSmootherGive.h"
#include "../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTake/extrapolatedSmootherTake.h"

// ----------- //
// Constructor //
template <concepts::DomainGeometry DomainGeometry>
Level<DomainGeometry>::Level(const int level_depth, std::unique_ptr<const PolarGrid> grid,
                             std::unique_ptr<const LevelCache<DomainGeometry>> level_cache,
                             const ExtrapolationType extrapolation, const bool FMG, const bool PCG_FMG)
    : level_depth_(level_depth)
    , grid_(std::move(grid))
    , level_cache_(std::move(level_cache))
    , rhs_("rhs", (FMG || PCG_FMG || level_depth == 0 || (level_depth == 1 && extrapolation != ExtrapolationType::NONE))
                      ? grid_->numberOfNodes()
                      : 0)
    , solution_("solution", grid_->numberOfNodes())
    , residual_("residual", grid_->numberOfNodes())
    , error_correction_("err_correction", (level_depth > 0) ? grid_->numberOfNodes() : 0)
{
}

// ---------------- //
// Getter Functions //
template <concepts::DomainGeometry DomainGeometry>
int Level<DomainGeometry>::level_depth() const
{
    return level_depth_;
}

template <concepts::DomainGeometry DomainGeometry>
const PolarGrid& Level<DomainGeometry>::grid() const
{
    return *grid_;
}

template <concepts::DomainGeometry DomainGeometry>
const LevelCache<DomainGeometry>& Level<DomainGeometry>::levelCache() const
{
    return *level_cache_;
}

template <concepts::DomainGeometry DomainGeometry>
Vector<double> Level<DomainGeometry>::rhs()
{
    return rhs_;
}

template <concepts::DomainGeometry DomainGeometry>
ConstVector<double> Level<DomainGeometry>::rhs() const
{
    return rhs_;
}

template <concepts::DomainGeometry DomainGeometry>
Vector<double> Level<DomainGeometry>::solution()
{
    return solution_;
}

template <concepts::DomainGeometry DomainGeometry>
ConstVector<double> Level<DomainGeometry>::solution() const
{
    return solution_;
}

template <concepts::DomainGeometry DomainGeometry>
Vector<double> Level<DomainGeometry>::residual()
{
    return residual_;
}

template <concepts::DomainGeometry DomainGeometry>
ConstVector<double> Level<DomainGeometry>::residual() const
{
    return residual_;
}

template <concepts::DomainGeometry DomainGeometry>
Vector<double> Level<DomainGeometry>::error_correction()
{
    return error_correction_;
}

template <concepts::DomainGeometry DomainGeometry>
ConstVector<double> Level<DomainGeometry>::error_correction() const
{
    return error_correction_;
}

// -------------- //
// Apply Residual //
template <concepts::DomainGeometry DomainGeometry>
void Level<DomainGeometry>::initializeResidual(const DensityProfileCoefficients& density_profile_coefficients,
                                               const bool DirBC_Interior, const int num_omp_threads,
                                               const StencilDistributionMethod stencil_distribution_method)
{
    if (stencil_distribution_method == StencilDistributionMethod::CPU_TAKE) {
        op_residual_ = std::make_unique<ResidualTake<DomainGeometry>>(
            *grid_, *level_cache_, density_profile_coefficients, DirBC_Interior, num_omp_threads);
    }
    else if (stencil_distribution_method == StencilDistributionMethod::CPU_GIVE) {
        op_residual_ = std::make_unique<ResidualGive<DomainGeometry>>(
            *grid_, *level_cache_, density_profile_coefficients, DirBC_Interior, num_omp_threads);
    }
    if (!op_residual_)
        throw std::runtime_error("Failed to initialize Residual.");
}

template <concepts::DomainGeometry DomainGeometry>
void Level<DomainGeometry>::computeResidual(Vector<double> result, ConstVector<double> rhs, ConstVector<double> x) const
{
    if (!op_residual_)
        throw std::runtime_error("Residual not initialized.");
    op_residual_->computeResidual(result, rhs, x);
}
template <concepts::DomainGeometry DomainGeometry>
void Level<DomainGeometry>::applySystemOperator(Vector<double> result, ConstVector<double> x) const
{
    if (!op_residual_)
        throw std::runtime_error("Residual not initialized.");
    op_residual_->applySystemOperator(result, x);
}

// ------------------- //
// Solve coarse System //
template <concepts::DomainGeometry DomainGeometry>
void Level<DomainGeometry>::initializeDirectSolver(const DensityProfileCoefficients& density_profile_coefficients,
                                                   const bool DirBC_Interior, const int num_omp_threads,
                                                   const StencilDistributionMethod stencil_distribution_method)
{
#ifdef GMGPOLAR_USE_MUMPS
    if (stencil_distribution_method == StencilDistributionMethod::CPU_TAKE) {
        op_directSolver_ = std::make_unique<DirectSolver_COO_MUMPS_Take<DomainGeometry>>(
            *grid_, *level_cache_, density_profile_coefficients, DirBC_Interior, num_omp_threads);
    }
    else if (stencil_distribution_method == StencilDistributionMethod::CPU_GIVE) {
        op_directSolver_ = std::make_unique<DirectSolver_COO_MUMPS_Give<DomainGeometry>>(
            *grid_, *level_cache_, density_profile_coefficients, DirBC_Interior, num_omp_threads);
    }
#else
    if (stencil_distribution_method == StencilDistributionMethod::CPU_TAKE) {
        op_directSolver_ = std::make_unique<DirectSolver_CSR_LU_Take<DomainGeometry>>(
            *grid_, *level_cache_, density_profile_coefficients, DirBC_Interior, num_omp_threads);
    }
    else if (stencil_distribution_method == StencilDistributionMethod::CPU_GIVE) {
        op_directSolver_ = std::make_unique<DirectSolver_CSR_LU_Give<DomainGeometry>>(
            *grid_, *level_cache_, density_profile_coefficients, DirBC_Interior, num_omp_threads);
    }
#endif
    if (!op_directSolver_)
        throw std::runtime_error("Failed to initialize Direct Solver.");
}

template <concepts::DomainGeometry DomainGeometry>
void Level<DomainGeometry>::directSolveInPlace(Vector<double> x) const
{
    if (!op_directSolver_)
        throw std::runtime_error("Coarse Solver not initialized.");
    op_directSolver_->solveInPlace(x);
}

// --------------- //
// Apply Smoothing //
template <concepts::DomainGeometry DomainGeometry>
void Level<DomainGeometry>::initializeSmoothing(const DensityProfileCoefficients& density_profile_coefficients,
                                                const bool DirBC_Interior, const int num_omp_threads,
                                                const StencilDistributionMethod stencil_distribution_method)
{
    if (stencil_distribution_method == StencilDistributionMethod::CPU_TAKE) {
        op_smoother_ = std::make_unique<SmootherTake<DomainGeometry>>(
            *grid_, *level_cache_, density_profile_coefficients, DirBC_Interior, num_omp_threads);
    }
    else if (stencil_distribution_method == StencilDistributionMethod::CPU_GIVE) {
        op_smoother_ = std::make_unique<SmootherGive<DomainGeometry>>(
            *grid_, *level_cache_, density_profile_coefficients, DirBC_Interior, num_omp_threads);
    }
    if (!op_smoother_)
        throw std::runtime_error("Failed to initialize Smoother.");
}

template <concepts::DomainGeometry DomainGeometry>
void Level<DomainGeometry>::smoothing(Vector<double> x, ConstVector<double> rhs, Vector<double> temp) const
{
    if (!op_smoother_)
        throw std::runtime_error("Smoother not initialized.");
    op_smoother_->smoothing(x, rhs, temp);
}

// ---------------------------- //
// Apply Extrapolated Smoothing //
template <concepts::DomainGeometry DomainGeometry>
void Level<DomainGeometry>::initializeExtrapolatedSmoothing(
    const DensityProfileCoefficients& density_profile_coefficients, const bool DirBC_Interior,
    const int num_omp_threads, const StencilDistributionMethod stencil_distribution_method)
{
    if (stencil_distribution_method == StencilDistributionMethod::CPU_TAKE) {
        op_extrapolated_smoother_ = std::make_unique<ExtrapolatedSmootherTake<DomainGeometry>>(
            *grid_, *level_cache_, density_profile_coefficients, DirBC_Interior, num_omp_threads);
    }
    else if (stencil_distribution_method == StencilDistributionMethod::CPU_GIVE) {
        op_extrapolated_smoother_ = std::make_unique<ExtrapolatedSmootherGive<DomainGeometry>>(
            *grid_, *level_cache_, density_profile_coefficients, DirBC_Interior, num_omp_threads);
    }
    if (!op_extrapolated_smoother_)
        throw std::runtime_error("Failed to initialize Extrapolated Smoother.");
}

template <concepts::DomainGeometry DomainGeometry>
void Level<DomainGeometry>::extrapolatedSmoothing(Vector<double> x, ConstVector<double> rhs, Vector<double> temp) const
{
    if (!op_extrapolated_smoother_)
        throw std::runtime_error("Extrapolated Smoother not initialized.");
    op_extrapolated_smoother_->extrapolatedSmoothing(x, rhs, temp);
}
