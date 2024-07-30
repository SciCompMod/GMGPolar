#include "../../include/Level/level.h"

// ----------- //
// Constructor //
Level::Level(int level, std::unique_ptr<const PolarGrid> grid, std::unique_ptr<LevelCache> level_cache, int extrapolation) :
    level_(level),
    grid_(std::move(grid)),
    level_cache_(std::move(level_cache)),
    rhs_((level == 0 || level == 1 && extrapolation > 0) ? grid_->number_of_nodes() : 0),
    solution_(grid_->number_of_nodes()),
    residual_(grid_->number_of_nodes()),
    rhs_error_((level > 0) ? grid_->number_of_nodes() : 0)
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

Vector<double>& Level::rhs() {
    return rhs_;
}
const Vector<double>& Level::rhs() const {
    return rhs_;
}
Vector<double>& Level::solution() {
    return solution_;
}
const Vector<double>& Level::solution() const {
    return solution_;
}
Vector<double>& Level::residual() {
    return residual_;
}
const Vector<double>& Level::residual() const {
    return residual_;
}
Vector<double>& Level::rhs_error() {
    return rhs_error_;
}
const Vector<double>& Level::rhs_error() const {
    return rhs_error_;
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
void Level::initializeDirectSolver(
    const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, 
    const int maxOpenMPThreads, const int openMPTaskThreads)
{
    op_directSolver_ = std::make_unique<DirectSolver>(*grid_, *level_cache_, domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    if (!op_directSolver_) throw std::runtime_error("Failed to initialize Direct Solver.");
}
void Level::directSolveInPlace(Vector<double>& x) const {
    if (!op_directSolver_) throw std::runtime_error("Coarse Solver not initialized.");
    op_directSolver_->solveInPlace(x);
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

// ---------------------------- //
// Apply Extrapolated Smoothing //
void Level::initializeExtrapolatedSmoothing(
    const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, 
    const int maxOpenMPThreads, const int openMPTaskThreads)
{
    op_extrapolated_smoother_ = std::make_unique<ExtrapolatedSmoother>(*grid_, *level_cache_, domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    if (!op_extrapolated_smoother_) throw std::runtime_error("Failed to initialize Extrapolated Smoother.");
}
void Level::extrapolatedSmoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp) const {
    if (!op_extrapolated_smoother_) throw std::runtime_error("Extrapolated Smoother not initialized.");
    op_extrapolated_smoother_ -> extrapolatedSmoothingInPlace(x, rhs, temp);
}
