#pragma once

class DirectSolver;
class Residual;
class Smoother;
class ExtrapolatedSmoother;

#include <memory>
#include <omp.h>
#include <vector>

#include "../PolarGrid/polargrid.h"

#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/densityProfileCoefficients.h"
#include "../InputFunctions/boundaryConditions.h"
#include "../InputFunctions/sourceTerm.h"

#include "../DirectSolver/directSolver.h"
#include "../Residual/residual.h"
#include "../Smoother/smoother.h"
#include "../ExtrapolatedSmoother/extrapolated_smoother.h"

class LevelCache;

class Level {
public:
    // ----------- //
    // Constructor //
    explicit Level(int level, std::unique_ptr<const PolarGrid> grid, std::unique_ptr<const LevelCache> level_cache, int extrapolation);

    // ---------------- //
    // Getter Functions //
    int level() const;
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
    void initializeResidual(const DomainGeometry& domain_geometry, const bool DirBC_Interior, const int num_omp_threads);
    void computeResidual(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const;

    // ------------------- //
    // Solve coarse System //
    void initializeDirectSolver(const DomainGeometry& domain_geometry, const bool DirBC_Interior, const int num_omp_threads);
    void directSolveInPlace(Vector<double>& x) const;

    // --------------- //
    // Apply Smoothing //
    void initializeSmoothing(const DomainGeometry& domain_geometry, const bool DirBC_Interior, const int num_omp_threads);
    void smoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp) const;

    // ---------------------------- //
    // Apply Extrapolated Smoothing //
    void initializeExtrapolatedSmoothing(const DomainGeometry& domain_geometry, const bool DirBC_Interior, const int num_omp_threads);
    void extrapolatedSmoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp) const;

private:
    const int level_;
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



class LevelCache {
public:
    explicit LevelCache(const PolarGrid& grid, const DensityProfileCoefficients& density_profile_coefficients);
    explicit LevelCache(const Level& previous_level, const PolarGrid& current_grid);

    const std::vector<double>& sin_theta() const;
    const std::vector<double>& cos_theta() const;

    const std::vector<double>& coeff_alpha() const;
    const std::vector<double>& coeff_beta() const;

private:
    std::vector<double> sin_theta_;
    std::vector<double> cos_theta_;

    std::vector<double> coeff_alpha_;
    std::vector<double> coeff_beta_;
};