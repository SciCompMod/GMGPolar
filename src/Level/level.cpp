#include "../../include/Level/level.h"

// ----------- //
// Constructor //
Level::Level(const int level, std::unique_ptr<const PolarGrid> grid, std::unique_ptr<LevelCache> level_cache) :
    level_(level), 
    grid_(std::move(grid)), 
    level_cache_(std::move(level_cache)),
    solution_(grid_->number_of_nodes()),
    residual_(grid_->number_of_nodes())
{}


// ---------------- //
// Getter Functions //
int Level::level() const {
    return level_;
}

const PolarGrid& Level::grid() const {
    return *grid_;
}

const LevelCache& Level::levelCache() const{
    return *level_cache_;
}

LevelCache& Level::levelCache() {
    return *level_cache_;
}

Vector<double>& Level::solution() {
    return solution_;
}
Vector<double>& Level::residual() {
    return residual_;
}

// -------------- //
// Apply Residual //
void Level::initializeResidual(
    const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior,
    const int maxOpenMPThreads, const int openMPTaskThreads)
{
    op_residual_ = std::make_unique<Residual>(*grid_, *level_cache_, domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    if (!op_residual_) throw std::runtime_error("Failed to initialize Residual.");
}
void Level::computeResidual(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const {
    if (!op_residual_) throw std::runtime_error("Residual not initialized.");
    op_residual_ -> computeResidual_V1(result, rhs, x);
}

// ------------------- //
// Solve coarse System //
void Level::initializeCoarseSolver(
    const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, 
    const int maxOpenMPThreads, const int openMPTaskThreads)
{
    op_coarseSolver_ = std::make_unique<CoarseSolver>(*grid_, *level_cache_, domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    if (!op_coarseSolver_) throw std::runtime_error("Failed to initialize Coarse Solver.");
}
void Level::coarseSolveInPlace(Vector<double>& x) const {
    if (!op_coarseSolver_) throw std::runtime_error("Coarse Solver not initialized.");
    op_coarseSolver_->solveInPlace(x);
}

// --------------- //
// Apply Smoothing //
void Level::initializeSmoothing(
    const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, 
    const int maxOpenMPThreads, const int openMPTaskThreads)
{
    op_smoother_ = std::make_unique<Smoother>(*grid_, *level_cache_, domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    if (!op_smoother_) throw std::runtime_error("Failed to initialize Smoother.");
}
void Level::smoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp) const {
    if (!op_smoother_) throw std::runtime_error("Smoother not initialized.");
    op_smoother_ -> smoothingInPlace(x, rhs, temp);
}

