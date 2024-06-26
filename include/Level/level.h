#pragma once

class CoarseSolver;
class Residual;
class Smoother;

#include "../PolarGrid/polargrid.h"

#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/systemParameters.h"

#include "../CoarseSolver/coarseSolver.h"
#include "../Residual/residual.h"
#include "../Smoother/smoother.h"

#include <memory>
#include <omp.h>
#include <vector>

class LevelCache;

class Level {
public:
    // ----------- //
    // Constructor //
    explicit Level(const int level, std::unique_ptr<const PolarGrid> grid, std::unique_ptr<LevelCache> level_cache);

    // ---------------- //
    // Getter Functions //
    int level() const;
    const PolarGrid& grid() const;
    const LevelCache& levelCache() const;
    LevelCache& levelCache();

    Vector<double>& solution();
    Vector<double>& residual();

    // -------------- //
    // Apply Residual //
    void initializeResidual(
        const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, 
        const int maxOpenMPThreads, const int openMPTaskThreads);
    void computeResidual(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const;

    // ------------------- //
    // Solve coarse System //
    void initializeCoarseSolver(
        const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior,
        const int maxOpenMPThreads, const int openMPTaskThreads);
    void coarseSolveInPlace(Vector<double>& x) const;

    // --------------- //
    // Apply Smoothing //
    void initializeSmoothing(
        const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior,
        const int maxOpenMPThreads, const int openMPTaskThreads);
    void smoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp) const;


    // void initializeExtrapolatedSmoothing(
    //     const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior,
    //     const int maxOpenMPThreads, const int openMPTaskThreads);
    // void extrapolatedSmoothing(Vector<double>& result, const Vector<double>& x);

private:
    const int level_;
    std::unique_ptr<const PolarGrid> grid_;
    std::unique_ptr<LevelCache> level_cache_;

    std::unique_ptr<CoarseSolver> op_coarseSolver_;
    std::unique_ptr<Residual> op_residual_;
    std::unique_ptr<Smoother> op_smoother_;

    // std::unique_ptr<const ExtrapolatedSmoother> 

    Vector<double> solution_;
    Vector<double> residual_;
};



class LevelCache {
public:
    explicit LevelCache(const PolarGrid& grid, const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior);
    explicit LevelCache(const Level& previous_level, const PolarGrid& current_grid, const bool DirBC_Interior);

    const std::vector<double>& sin_theta() const;
    const std::vector<double>& cos_theta() const;

    const Vector<double>& rhs() const;
    Vector<double>& rhs();
private:
    std::vector<double> sin_theta_;
    std::vector<double> cos_theta_;

    /* Vector<double> and is used when the size matches grid.numberOfNodes(). */
    Vector<double> rhs_;
};