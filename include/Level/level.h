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

// The `Level` class represents a single level of a multigrid method. 
// In multigrid solvers, the computational domain is divided into different levels, where each level corresponds to a grid with a different resolution. 
// The `Level` class manages the specific data structures and operations needed to solve a problem at a given level, including residual computation, direct solving, and smoothing. 
// It holds information for the solution, residuals, right-hand side, and error corrections used in the multigrid method.

// The `LevelCache` class is responsible for caching auxiliary data required for solving a problem at a specific level of a multigrid method. 
// It stores essential data such as trigonometric values (e.g., `sin_theta` and `cos_theta`) and profile coefficients (e.g., `alpha`, `beta`) 
// that are frequently used in the solution process. Additionally, depending on the stencil distribution strategy, it can store transformation 
// coefficients (`arr`, `att`, `art`) related to the domain geometry. These coefficients are critical for efficient matrix-free stencil operations 
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

    Vector<double>& rhs();
    const Vector<double>& rhs() const;
    Vector<double>& solution();
    const Vector<double>& solution() const;
    Vector<double>& residual();
    const Vector<double>& residual() const;
    Vector<double>& error_correction();
    const Vector<double>& error_correction() const;

    // -------------- //
    // Apply Residual //
    void initializeResidual(const DomainGeometry& domain_geometry,
                            const DensityProfileCoefficients& density_profile_coefficients, const bool DirBC_Interior,
                            const int num_omp_threads, const StencilDistributionMethod stencil_distribution_method);
    void computeResidual(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const;

    // ------------------- //
    // Solve coarse System //
    void initializeDirectSolver(const DomainGeometry& domain_geometry,
                                const DensityProfileCoefficients& density_profile_coefficients,
                                const bool DirBC_Interior, const int num_omp_threads,
                                const StencilDistributionMethod stencil_distribution_method);
    // Note: The rhs (right-hand side) vector gets overwritten during the solution process.
    void directSolveInPlace(Vector<double>& x) const;

    // --------------- //
    // Apply Smoothing //
    void initializeSmoothing(const DomainGeometry& domain_geometry,
                             const DensityProfileCoefficients& density_profile_coefficients, const bool DirBC_Interior,
                             const int num_omp_threads, const StencilDistributionMethod stencil_distribution_method);
    void smoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp) const;

    // ---------------------------- //
    // Apply Extrapolated Smoothing //
    void initializeExtrapolatedSmoothing(const DomainGeometry& domain_geometry,
                                         const DensityProfileCoefficients& density_profile_coefficients,
                                         const bool DirBC_Interior, const int num_omp_threads,
                                         const StencilDistributionMethod stencil_distribution_method);
    void extrapolatedSmoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp) const;

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

    const std::vector<double>& sin_theta() const;
    const std::vector<double>& cos_theta() const;

    bool cacheDensityProfileCoefficients() const;
    const std::vector<double>& coeff_alpha() const;
    const std::vector<double>& coeff_beta() const;

    bool cacheDomainGeometry() const;
    const Vector<double>& arr() const;
    const Vector<double>& att() const;
    const Vector<double>& art() const;
    const Vector<double>& detDF() const;

private:
    std::vector<double> sin_theta_;
    std::vector<double> cos_theta_;

    bool cache_density_profile_coefficients_;
    std::vector<double> coeff_alpha_;
    std::vector<double> coeff_beta_;

    bool cache_domain_geometry_;
    Vector<double> arr_;
    Vector<double> att_;
    Vector<double> art_;
    Vector<double> detDF_;
};
