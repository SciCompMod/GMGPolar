#pragma once

class DirectSolver;
class Residual;
class Smoother;
class ExtrapolatedSmoother;

#include "../PolarGrid/polargrid.h"

#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/systemParameters.h"

#include "../DirectSolver/directSolver.h"
#include "../Residual/residual.h"
#include "../Smoother/smoother.h"
#include "../ExtrapolatedSmoother/extrapolated_smoother.h"

#include <memory>
#include <omp.h>
#include <vector>

class LevelCache;

class Level {
public:
    // ----------- //
    // Constructor //

    explicit Level(int level, std::unique_ptr<const PolarGrid> grid, std::unique_ptr<LevelCache> level_cache, int extrapolation);

    // ---------------- //
    // Getter Functions //
    int level() const;
    const PolarGrid& grid() const;
    const LevelCache& levelCache() const;
    LevelCache& levelCache();

    Vector<double>& rhs();
    const Vector<double>& rhs() const;
    Vector<double>& solution();
    const Vector<double>& solution() const;
    Vector<double>& residual();
    const Vector<double>& residual() const;
    Vector<double>& rhs_error();
    const Vector<double>& rhs_error() const;

    // -------------- //
    // Apply Residual //
    void initializeResidual(
        const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, 
        const int maxOpenMPThreads, const int openMPTaskThreads);
    void computeResidual(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const;

    // ------------------- //
    // Solve coarse System //
    void initializeDirectSolver(
        const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior,
        const int maxOpenMPThreads, const int openMPTaskThreads);
    void directSolveInPlace(Vector<double>& x) const;

    // --------------- //
    // Apply Smoothing //
    void initializeSmoothing(
        const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior,
        const int maxOpenMPThreads, const int openMPTaskThreads);
    void smoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp) const;

    // ---------------------------- //
    // Apply Extrapolated Smoothing //
    void initializeExtrapolatedSmoothing(
        const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior,
        const int maxOpenMPThreads, const int openMPTaskThreads);
    void extrapolatedSmoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp) const;

private:
    const int level_;
    std::unique_ptr<const PolarGrid> grid_;
    std::unique_ptr<LevelCache> level_cache_;

    std::unique_ptr<DirectSolver> op_directSolver_;
    std::unique_ptr<Residual> op_residual_;
    std::unique_ptr<Smoother> op_smoother_;
    std::unique_ptr<ExtrapolatedSmoother> op_extrapolated_smoother_;

    Vector<double> rhs_;
    Vector<double> solution_;
    Vector<double> residual_;
    Vector<double> rhs_error_;
};



class LevelCache {
public:
    explicit LevelCache(const PolarGrid& grid);
    explicit LevelCache(const Level& previous_level, const PolarGrid& current_grid);

    const std::vector<double>& sin_theta() const;
    const std::vector<double>& cos_theta() const;
private:
    std::vector<double> sin_theta_;
    std::vector<double> cos_theta_;
};