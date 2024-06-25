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
    explicit Level(const int level, std::unique_ptr<const PolarGrid> grid, std::unique_ptr<const LevelCache> level_cache);

    // ---------------- //
    // Getter Functions //
    int level() const;
    const PolarGrid& grid() const;
    const LevelCache& levelCache() const;

    // -------------- //
    // Apply Residual //
    void initializeResidual(
        const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, 
        const int maxOpenMPThreads, const int openMPTaskThreads);
    void computeResidual(Vector<double>& result, const Vector<double>& x);

    // ------------------- //
    // Solve coarse System //
    void initializeCoarseSolver(
        const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior,
        const int maxOpenMPThreads, const int openMPTaskThreads);
    void coarseSolveInPlace(Vector<double>& x);

    // --------------- //
    // Apply Smoothing //
    void initializeSmoothing(
        const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior,
        const int maxOpenMPThreads, const int openMPTaskThreads);
    void smoothingInPlace(Vector<double>& x, Vector<double>& temp_rhs);


    // void initializeExtrapolatedSmoothing(
    //     const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior,
    //     const int maxOpenMPThreads, const int openMPTaskThreads);
    // void extrapolatedSmoothing(Vector<double>& result, const Vector<double>& x);

private:
    const int level_;
    std::unique_ptr<const PolarGrid> grid_;
    std::unique_ptr<const LevelCache> level_cache_;

    std::unique_ptr<CoarseSolver> coarseSolver_;
    std::unique_ptr<Residual> residual_;
    std::unique_ptr<Smoother> smoother_;

    // std::unique_ptr<const ExtrapolatedSmoother> 
};



class LevelCache {
public:
    explicit LevelCache(const PolarGrid& grid, const SystemParameters& system_parameters, const bool DirBC_Interior);
    explicit LevelCache(const Level& previous_level, const PolarGrid& current_grid, const bool DirBC_Interior);

    const std::vector<double>& sin_theta() const;
    const std::vector<double>& cos_theta() const;

    const Vector<double>& rhs_f() const;
    const std::vector<double>& u_D() const;
    const std::vector<double>& u_D_Interior() const;
private:
    std::vector<double> sin_theta_;
    std::vector<double> cos_theta_;

    /* Vector<double> and is used when the size matches grid.numberOfNodes(). */
    Vector<double> rhs_f_;
    std::vector<double> u_D_;
    std::vector<double> u_D_Interior_;

    /* Note that we don't store the discretization of rhs_f. */
    // Vector<double> discretization_rhs_f_;
};