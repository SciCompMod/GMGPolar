#include "../../include/Level/level.h"

#include "../../include/Residual/ResidualGive/residualGive.h"
#include "../../include/Residual/ResidualTake/residualTake.h"

#include "../../include/DirectSolver/DirectSolverGive/directSolverGive.h"
#include "../../include/DirectSolver/DirectSolverTake/directSolverTake.h"

#include "../../include/Smoother/SmootherGive/smootherGive.h"
#include "../../include/Smoother/SmootherTake/smootherTake.h"

#include "../../include/ExtrapolatedSmoother/ExtrapolatedSmootherGive/extrapolatedSmootherGive.h"
#include "../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTake/extrapolatedSmootherTake.h"

// ----------- //
// Constructor //
Level::Level(const int level,
             std::unique_ptr<const PolarGrid> grid,
             std::unique_ptr<const LevelCache> level_cache,
             const ExtrapolationType extrapolation,
             const bool FMG)
    : level_(level)
    , grid_(std::move(grid))
    , level_cache_(std::move(level_cache))
    , rhs_((FMG || level == 0 || (level == 1 && extrapolation != ExtrapolationType::NONE)) ? grid_->numberOfNodes() : 0)
    , solution_(grid_->numberOfNodes())
    , residual_(grid_->numberOfNodes())
    , error_correction_((level > 0) ? grid_->numberOfNodes() : 0)
{
}

// ---------------- //
// Getter Functions //
int Level::level() const
{
    return level_;
}

const PolarGrid& Level::grid() const
{
    return *grid_;
}

const LevelCache& Level::levelCache() const
{
    return *level_cache_;
}

Vector<double>& Level::rhs()
{
    return rhs_;
}
const Vector<double>& Level::rhs() const
{
    return rhs_;
}
Vector<double>& Level::solution()
{
    return solution_;
}
const Vector<double>& Level::solution() const
{
    return solution_;
}
Vector<double>& Level::residual()
{
    return residual_;
}
const Vector<double>& Level::residual() const
{
    return residual_;
}
Vector<double>& Level::error_correction()
{
    return error_correction_;
}
const Vector<double>& Level::error_correction() const
{
    return error_correction_;
}

// -------------- //
// Apply Residual //
void Level::initializeResidual(const DomainGeometry& domain_geometry,
                               const DensityProfileCoefficients& density_profile_coefficients,
                               const bool DirBC_Interior,
                               const int num_omp_threads,
                               const StencilDistributionMethod stencil_distribution_method)
{
    if (stencil_distribution_method == StencilDistributionMethod::CPU_TAKE)
    {
        op_residual_ =
            std::make_unique<ResidualTake>(*grid_, *level_cache_, domain_geometry, density_profile_coefficients, DirBC_Interior, num_omp_threads);
    }
    else if (stencil_distribution_method == StencilDistributionMethod::CPU_GIVE)
    {
        op_residual_ =
            std::make_unique<ResidualGive>(*grid_, *level_cache_, domain_geometry, density_profile_coefficients, DirBC_Interior, num_omp_threads);
    }
    if (!op_residual_) throw std::runtime_error("Failed to initialize Residual.");
}
void Level::computeResidual(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const
{
    if (!op_residual_) throw std::runtime_error("Residual not initialized.");
    op_residual_->computeResidual(result, rhs, x);
}

// ------------------- //
// Solve coarse System //
void Level::initializeDirectSolver(const DomainGeometry& domain_geometry,
                                   const DensityProfileCoefficients& density_profile_coefficients,
                                   const bool DirBC_Interior,
                                   const int num_omp_threads,
                                   const StencilDistributionMethod stencil_distribution_method)
{
    if (stencil_distribution_method == StencilDistributionMethod::CPU_TAKE)
    {
        op_directSolver_ =
            std::make_unique<DirectSolverTake>(*grid_, *level_cache_, domain_geometry, density_profile_coefficients, DirBC_Interior, num_omp_threads);
    }
    else if (stencil_distribution_method == StencilDistributionMethod::CPU_GIVE)
    {
        op_directSolver_ =
            std::make_unique<DirectSolverGive>(*grid_, *level_cache_, domain_geometry, density_profile_coefficients, DirBC_Interior, num_omp_threads);
    }
    if (!op_directSolver_) throw std::runtime_error("Failed to initialize Direct Solver.");
}
void Level::directSolveInPlace(Vector<double>& x) const
{
    if (!op_directSolver_) throw std::runtime_error("Coarse Solver not initialized.");
    op_directSolver_->solveInPlace(x);
}

// --------------- //
// Apply Smoothing //
void Level::initializeSmoothing(const DomainGeometry& domain_geometry,
                                const DensityProfileCoefficients& density_profile_coefficients,
                                const bool DirBC_Interior,
                                const int num_omp_threads,
                                const StencilDistributionMethod stencil_distribution_method)
{
    if (stencil_distribution_method == StencilDistributionMethod::CPU_TAKE)
    {
        op_smoother_ =
            std::make_unique<SmootherTake>(*grid_, *level_cache_, domain_geometry, density_profile_coefficients, DirBC_Interior, num_omp_threads);
    }
    else if (stencil_distribution_method == StencilDistributionMethod::CPU_GIVE)
    {
        op_smoother_ =
            std::make_unique<SmootherGive>(*grid_, *level_cache_, domain_geometry, density_profile_coefficients, DirBC_Interior, num_omp_threads);
    }
    if (!op_smoother_) throw std::runtime_error("Failed to initialize Smoother.");
}
void Level::smoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp) const
{
    if (!op_smoother_) throw std::runtime_error("Smoother not initialized.");
    op_smoother_->smoothingInPlace(x, rhs, temp);
}

// ---------------------------- //
// Apply Extrapolated Smoothing //
void Level::initializeExtrapolatedSmoothing(const DomainGeometry& domain_geometry,
                                            const DensityProfileCoefficients& density_profile_coefficients,
                                            const bool DirBC_Interior,
                                            const int num_omp_threads,
                                            const StencilDistributionMethod stencil_distribution_method)
{
    if (stencil_distribution_method == StencilDistributionMethod::CPU_TAKE)
    {
        op_extrapolated_smoother_ = std::make_unique<ExtrapolatedSmootherTake>(*grid_, *level_cache_, domain_geometry, density_profile_coefficients,
                                                                               DirBC_Interior, num_omp_threads);
    }
    else if (stencil_distribution_method == StencilDistributionMethod::CPU_GIVE)
    {
        op_extrapolated_smoother_ = std::make_unique<ExtrapolatedSmootherGive>(*grid_, *level_cache_, domain_geometry, density_profile_coefficients,
                                                                               DirBC_Interior, num_omp_threads);
    }
    if (!op_extrapolated_smoother_) throw std::runtime_error("Failed to initialize Extrapolated Smoother.");
}
void Level::extrapolatedSmoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp) const
{
    if (!op_extrapolated_smoother_) throw std::runtime_error("Extrapolated Smoother not initialized.");
    op_extrapolated_smoother_->extrapolatedSmoothingInPlace(x, rhs, temp);
}
