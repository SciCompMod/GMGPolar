#pragma once

// Required to prevent circular dependencies.
class DirectSolver;
class Residual;
class Smoother;
class ExtrapolatedSmoother;

#include <memory>
#include <omp.h>
#include <vector>

#include "../PolarGrid/polargrid.h"

#include "../InputFunctions/boundaryConditions.h"
#include "../InputFunctions/densityProfileCoefficients.h"
#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/sourceTerm.h"

#include "../DirectSolver/directSolver.h"
#include "../ExtrapolatedSmoother/extrapolatedSmoother.h"
#include "../Residual/residual.h"
#include "../Smoother/smoother.h"

#include "../common/geometry_helper.h"

// The `Level` class represents a single level of a multigrid method.
// In multigrid solvers, the computational domain is divided into different levels, where each level corresponds to a grid with a different resolution.
// The `Level` class manages the specific data structures and operations needed to solve a problem at a given level, including residual computation, direct solving, and smoothing.
// It holds information for the solution, residuals, right-hand side, and error corrections used in the multigrid method.

// The `LevelCache` class is responsible for caching auxiliary data required for solving a problem at a specific level of a multigrid method.
// It stores essential data such as profile coefficients (e.g., `alpha`, `beta`)
// that are frequently used in the solution process. Additionally, depending on the stencil distribution strategy, it can store transformation
// coefficients (`arr`, `att`, `art`, `detDF`) related to the domain geometry. These coefficients are critical for efficient matrix-free stencil operations
// and contribute to the accuracy and performance of the multigrid solver.

class LevelCache;

class Level
{
public:
    // ----------- //
    // Constructor //
    explicit Level(const int level_depth, std::unique_ptr<const PolarGrid> grid,
                   std::unique_ptr<const LevelCache> level_cache, const ExtrapolationType extrapolation,
                   const bool FMG);

    // ---------------- //
    // Getter Functions //
    int level_depth() const;
    const PolarGrid& grid() const;
    const LevelCache& levelCache() const;

    Vector<double> rhs();
    ConstVector<double> rhs() const;
    Vector<double> solution();
    ConstVector<double> solution() const;
    Vector<double> residual();
    ConstVector<double> residual() const;
    Vector<double> error_correction();
    ConstVector<double> error_correction() const;

    // -------------- //
    // Apply Residual //
    void initializeResidual(const DomainGeometry& domain_geometry,
                            const DensityProfileCoefficients& density_profile_coefficients, const bool DirBC_Interior,
                            const int num_omp_threads, const StencilDistributionMethod stencil_distribution_method);
    void computeResidual(Vector<double> result, ConstVector<double> rhs, ConstVector<double> x) const;

    // ------------------- //
    // Solve coarse System //
    void initializeDirectSolver(const DomainGeometry& domain_geometry,
                                const DensityProfileCoefficients& density_profile_coefficients,
                                const bool DirBC_Interior, const int num_omp_threads,
                                const StencilDistributionMethod stencil_distribution_method);
    // Note: The rhs (right-hand side) vector gets overwritten by the solution.
    void directSolveInPlace(Vector<double> x) const;

    // --------------- //
    // Apply Smoothing //
    void initializeSmoothing(const DomainGeometry& domain_geometry,
                             const DensityProfileCoefficients& density_profile_coefficients, const bool DirBC_Interior,
                             const int num_omp_threads, const StencilDistributionMethod stencil_distribution_method);
    void smoothing(Vector<double> x, ConstVector<double> rhs, Vector<double> temp) const;

    // ---------------------------- //
    // Apply Extrapolated Smoothing //
    void initializeExtrapolatedSmoothing(const DomainGeometry& domain_geometry,
                                         const DensityProfileCoefficients& density_profile_coefficients,
                                         const bool DirBC_Interior, const int num_omp_threads,
                                         const StencilDistributionMethod stencil_distribution_method);
    void extrapolatedSmoothing(Vector<double> x, ConstVector<double> rhs, Vector<double> temp) const;

private:
    const int level_depth_;
    std::unique_ptr<const PolarGrid> grid_;
    std::unique_ptr<const LevelCache> level_cache_;

    std::unique_ptr<DirectSolver> op_directSolver_;
    std::unique_ptr<Residual> op_residual_;
    std::unique_ptr<Smoother> op_smoother_;
    std::unique_ptr<ExtrapolatedSmoother> op_extrapolated_smoother_;

    Vector<double> rhs_;
    Vector<double> solution_;
    Vector<double> residual_;
    Vector<double> error_correction_;
};

class LevelCache
{
public:
    explicit LevelCache(const PolarGrid& grid, const DensityProfileCoefficients& density_profile_coefficients,
                        const DomainGeometry& domain_geometry, const bool cache_density_profile_coefficients,
                        const bool cache_domain_geometry);
    explicit LevelCache(const Level& previous_level, const PolarGrid& current_grid);

    const DomainGeometry& domainGeometry() const;
    const DensityProfileCoefficients& densityProfileCoefficients() const;

    bool cacheDensityProfileCoefficients() const;
    ConstVector<double> coeff_alpha() const;
    ConstVector<double> coeff_beta() const;

    bool cacheDomainGeometry() const;
    ConstVector<double> arr() const;
    ConstVector<double> att() const;
    ConstVector<double> art() const;
    ConstVector<double> detDF() const;

    inline void obtainValues(const int i_r, const int i_theta, const int global_index, const double& r,
                             const double& theta, double& coeff_beta, double& arr, double& att, double& art,
                             double& detDF) const
    {
        if (cache_density_profile_coefficients_)
            coeff_beta = coeff_beta_[global_index];
        else
            coeff_beta = density_profile_coefficients_.beta(r, theta);

        double coeff_alpha;
        if (!cache_domain_geometry_) {
            if (cache_density_profile_coefficients_)
                coeff_alpha = coeff_alpha_[global_index];
            else
                coeff_alpha = density_profile_coefficients_.alpha(r, theta);
        }

        if (cache_domain_geometry_) {
            arr   = arr_[global_index];
            att   = att_[global_index];
            art   = art_[global_index];
            detDF = detDF_[global_index];
        }
        else {
            compute_jacobian_elements(domain_geometry_, r, theta, coeff_alpha, arr, att, art, detDF);
        }
    }

private:
    const DomainGeometry& domain_geometry_;
    const DensityProfileCoefficients& density_profile_coefficients_;

    bool cache_density_profile_coefficients_; // cache alpha(r_i), beta(r_i)
    Vector<double> coeff_alpha_;
    Vector<double> coeff_beta_;

    bool cache_domain_geometry_; // cache arr, att, art, detDF
    Vector<double> arr_;
    Vector<double> att_;
    Vector<double> art_;
    Vector<double> detDF_;
};
