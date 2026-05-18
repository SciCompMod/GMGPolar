#pragma once

namespace gmgpolar
{

// Required to prevent circular dependencies.
template <class LevelCacheType>
class DirectSolver;
template <class LevelCacheType>
class Residual;
template <class LevelCacheType>
class Smoother;
template <class LevelCacheType>
class ExtrapolatedSmoother;

} // namespace gmgpolar

#include <memory>
#include <omp.h>
#include <vector>

#include <Kokkos_Core.hpp>

#include "../PolarGrid/polargrid.h"

#include "../InputFunctions/boundaryConditions.h"
#include "../InputFunctions/densityProfileCoefficients.h"
#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/sourceTerm.h"

#include "../DirectSolver/directSolver.h"
#include "../ExtrapolatedSmoother/extrapolatedSmoother.h"
#include "../Residual/residual.h"
#include "../Smoother/smoother.h"

#include "../Definitions/geometry_helper.h"

#include "../../include/Residual/ResidualGive/residualGive.h"
#include "../../include/Residual/ResidualTake/residualTake.h"

#include "../../include/DirectSolver/DirectSolverGive/directSolverGive.h"
#include "../../include/DirectSolver/DirectSolverTake/directSolverTake.h"

#include "../../include/Smoother/SmootherGive/smootherGive.h"
#include "../../include/Smoother/SmootherTake/smootherTake.h"

#include "../../include/ExtrapolatedSmoother/ExtrapolatedSmootherGive/extrapolatedSmootherGive.h"
#include "../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTake/extrapolatedSmootherTake.h"

namespace gmgpolar
{

// The `Level` class represents a single level of a multigrid method.
// In multigrid solvers, the computational domain is divided into different levels, where each level corresponds to a grid with a different resolution.
// The `Level` class manages the specific data structures and operations needed to solve a problem at a given level, including residual computation, direct solving, and smoothing.
// It holds information for the solution, residuals, right-hand side, and error corrections used in the multigrid method.

// The `LevelCache` class is responsible for caching auxiliary data required for solving a problem at a specific level of a multigrid method.
// It stores essential data such as profile coefficients (e.g., `alpha`, `beta`)
// that are frequently used in the solution process. Additionally, depending on the stencil distribution strategy, it can store transformation
// coefficients (`arr`, `att`, `art`, `detDF`) related to the domain geometry. These coefficients are critical for efficient matrix-free stencil operations
// and contribute to the accuracy and performance of the multigrid solver.

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
class LevelCache;

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
class Level
{
public:
    using LevelCacheType = LevelCache<DomainGeometry, DensityProfileCoefficients>;

public:
    // ----------- //
    // Constructor //
    explicit Level(const int level_depth, std::unique_ptr<const PolarGrid> grid,
                   std::unique_ptr<const LevelCacheType> level_cache, const ExtrapolationType extrapolation,
                   const bool FMG, const bool PCG_FMG = false);

    // ---------------- //
    // Getter Functions //
    int level_depth() const;
    const PolarGrid& grid() const;
    const LevelCacheType& levelCache() const;

    HostVector<double> rhs();
    HostConstVector<double> rhs() const;
    HostVector<double> solution();
    HostConstVector<double> solution() const;
    HostVector<double> residual();
    HostConstVector<double> residual() const;
    HostVector<double> error_correction();
    HostConstVector<double> error_correction() const;

    // -------------- //
    // Apply Residual //
    void initializeResidual(const bool DirBC_Interior, const StencilDistributionMethod stencil_distribution_method);
    void computeResidual(HostVector<double> result, HostConstVector<double> rhs, HostConstVector<double> x) const;
    void applySystemOperator(HostVector<double> result, HostConstVector<double> x) const;

    // ------------------- //
    // Solve coarse System //
    void initializeDirectSolver(const bool DirBC_Interior, const StencilDistributionMethod stencil_distribution_method);
    // Note: The rhs (right-hand side) vector gets overwritten by the solution.
    void directSolveInPlace(HostVector<double> x) const;

    // --------------- //
    // Apply Smoothing //
    void initializeSmoothing(const bool DirBC_Interior, const StencilDistributionMethod stencil_distribution_method);
    void smoothing(HostVector<double> x, HostConstVector<double> rhs, HostVector<double> temp) const;

    // ---------------------------- //
    // Apply Extrapolated Smoothing //
    void initializeExtrapolatedSmoothing(const bool DirBC_Interior,
                                         const StencilDistributionMethod stencil_distribution_method);
    void extrapolatedSmoothing(HostVector<double> x, HostConstVector<double> rhs, HostVector<double> temp) const;

private:
    const int level_depth_;
    std::unique_ptr<const PolarGrid> grid_;
    std::unique_ptr<const LevelCacheType> level_cache_;

    std::unique_ptr<DirectSolver<LevelCacheType>> op_directSolver_;
    std::unique_ptr<Residual<LevelCacheType>> op_residual_;
    std::unique_ptr<Smoother<LevelCacheType>> op_smoother_;
    std::unique_ptr<ExtrapolatedSmoother<LevelCacheType>> op_extrapolated_smoother_;

    HostVector<double> rhs_;
    HostVector<double> solution_;
    HostVector<double> residual_;
    HostVector<double> error_correction_;
};

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
class LevelCache
{
public:
    explicit LevelCache(const PolarGrid& grid, const DensityProfileCoefficients& density_profile_coefficients,
                        const DomainGeometry& domain_geometry, const bool cache_density_profile_coefficients,
                        const bool cache_domain_geometry);

    const DomainGeometry& domainGeometry() const;
    const DensityProfileCoefficients& densityProfileCoefficients() const;

    bool cacheDensityProfileCoefficients() const;
    HostConstVector<double> coeff_alpha() const;
    HostConstVector<double> coeff_beta() const;

    bool cacheDomainGeometry() const;
    HostConstVector<double> arr() const;
    HostConstVector<double> att() const;
    HostConstVector<double> art() const;
    HostConstVector<double> detDF() const;

    KOKKOS_INLINE_FUNCTION void obtainValues(const int i_r, const int i_theta, const int global_index, double r,
                                             double theta, double& coeff_beta, double& arr, double& att, double& art,
                                             double& detDF) const
    {
        coeff_beta = cache_density_profile_coefficients_ ? coeff_beta_[global_index]
                                                         : density_profile_coefficients_.beta(r, theta);

        if (cache_domain_geometry_) {
            arr   = arr_[global_index];
            att   = att_[global_index];
            art   = art_[global_index];
            detDF = detDF_[global_index];
        }
        else {
            double coeff_alpha = cache_density_profile_coefficients_ ? coeff_alpha_[global_index]
                                                                     : density_profile_coefficients_.alpha(r, theta);

            compute_jacobian_elements(domain_geometry_, r, theta, coeff_alpha, arr, att, art, detDF);
        }
    }

private:
    const DomainGeometry& domain_geometry_;
    const DensityProfileCoefficients& density_profile_coefficients_;

    bool cache_density_profile_coefficients_; // cache alpha(r, theta), beta(r, theta)
    HostVector<double> coeff_alpha_;
    HostVector<double> coeff_beta_;

    bool cache_domain_geometry_; // cache arr, att, art, detDF
    HostVector<double> arr_;
    HostVector<double> att_;
    HostVector<double> art_;
    HostVector<double> detDF_;
};

#include "levelCache.inl"
#include "level.inl"

} // namespace gmgpolar
