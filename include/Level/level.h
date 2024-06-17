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
    explicit Level(const int level, std::unique_ptr<const PolarGrid> grid);

    const PolarGrid& grid() const;
    int level() const;

    void initializeResidual(
        const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, 
        const int maxOpenMPThreads, const int openMPTaskThreads);
    void computeResidual(Vector<double>& result, const Vector<double>& x);

    void initializeCoarseSolver(
        const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior,
        const int maxOpenMPThreads, const int openMPTaskThreads);
    void coarseSolveInPlace(Vector<double>& x);

    void initializeSmoothing(
        const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior,
        const int maxOpenMPThreads, const int openMPTaskThreads);
    void smoothingInPlace(Vector<double>& x, Vector<double>& temp_rhs);

    // void initializeExtrapolatedSmoothing(
    //     const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior,
    //     const int maxOpenMPThreads, const int openMPTaskThreads);
    // void extrapolatedSmoothing(Vector<double>& result, const Vector<double>& x);

    const LevelCache& levelData() const;

private:
    const int level_;
    std::unique_ptr<const PolarGrid> grid_;
    std::unique_ptr<const LevelCache> level_data_;

    std::unique_ptr<CoarseSolver> coarseSolver_;
    std::unique_ptr<Residual> residual_;
    std::unique_ptr<Smoother> smoother_;

    // std::unique_ptr<const ExtrapolatedSmoother> 
};


class LevelCache {
public:
    explicit LevelCache(const int level, const PolarGrid& grid);

    const std::vector<double>& sin_theta() const;
    const std::vector<double>& cos_theta() const;
private:
    std::vector<double> sin_theta_;
    std::vector<double> cos_theta_;
};