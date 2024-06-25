#include "../../include/Level/level.h"

// ----------- //
// Constructor //
Level::Level(const int level, std::unique_ptr<const PolarGrid> grid, std::unique_ptr<const LevelCache> level_cache) :
    level_(level), 
    grid_(std::move(grid)), 
    level_cache_(std::move(level_cache))
{}


// ---------------- //
// Getter Functions //
int Level::level() const {
    return level_;
}

const PolarGrid& Level::grid() const {
    return *grid_;
}

const LevelCache& Level::levelCache() const {
    return *level_cache_;
}


// -------------- //
// Apply Residual //
void Level::initializeResidual(
    const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior,
    const int maxOpenMPThreads, const int openMPTaskThreads)
{
    residual_ = std::make_unique<Residual>(*grid_, *level_cache_, domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    if (!residual_) throw std::runtime_error("Failed to initialize Residual.");
}
void Level::computeResidual(Vector<double>& result, const Vector<double>& x){
    if (!residual_) throw std::runtime_error("Residual not initialized.");
    residual_ -> computeResidual_V1(result, x);
}

// ------------------- //
// Solve coarse System //
void Level::initializeCoarseSolver(
    const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, 
    const int maxOpenMPThreads, const int openMPTaskThreads)
{
    coarseSolver_ = std::make_unique<CoarseSolver>(*grid_, *level_cache_, domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    if (!coarseSolver_) throw std::runtime_error("Failed to initialize Coarse Solver.");
}
void Level::coarseSolveInPlace(Vector<double>& x){
    if (!coarseSolver_) throw std::runtime_error("Coarse Solver not initialized.");
    coarseSolver_->solveInPlace(x);
}

// --------------- //
// Apply Smoothing //
void Level::initializeSmoothing(
    const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, 
    const int maxOpenMPThreads, const int openMPTaskThreads)
{
    smoother_ = std::make_unique<Smoother>(*grid_, *level_cache_, domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    if (!smoother_) throw std::runtime_error("Failed to initialize Smoother.");
}
void Level::smoothingInPlace(Vector<double>& x, Vector<double>& temp_rhs){
    if (!smoother_) throw std::runtime_error("Smoother not initialized.");
    smoother_ -> smoothingInPlace(x, temp_rhs);
}

