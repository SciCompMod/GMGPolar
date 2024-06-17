#include "../../include/Level/level.h"

Level::Level(const int level, std::unique_ptr<const PolarGrid> grid) :
    level_(level), 
    grid_(std::move(grid)), 
    level_data_(std::make_unique<const LevelCache>(level_, *grid_)) 
{}

int Level::level() const {
    return level_;
}

const PolarGrid& Level::grid() const {
    return *grid_;
}


LevelCache::LevelCache(const int level, const PolarGrid& grid) : 
    sin_theta_(grid.ntheta()), 
    cos_theta_(grid.ntheta()) 
{
    #pragma omp parallel for
    for (int i = 0; i < grid.ntheta(); i++) {
        sin_theta_[i] = sin(grid.theta(i));
        cos_theta_[i] = cos(grid.theta(i));
    }
}

const std::vector<double>& LevelCache::sin_theta() const {
    return sin_theta_;
}
const std::vector<double>& LevelCache::cos_theta() const {
    return cos_theta_;
}

void Level::initializeResidual(const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, const int maxOpenMPThreads, const int openMPTaskThreads) {
    residual_ = std::make_unique<Residual>(*grid_, *level_data_, domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    if (!residual_) throw std::runtime_error("Failed to initialize Residual.");
}
void Level::computeResidual(Vector<double>& result, const Vector<double>& x){
    if (!residual_) throw std::runtime_error("Residual not initialized.");
    residual_ -> computeResidual_V1(result, x);
}


void Level::initializeCoarseSolver(const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, const int maxOpenMPThreads, const int openMPTaskThreads){
    coarseSolver_ = std::make_unique<CoarseSolver>(*grid_, *level_data_, domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    if (!coarseSolver_) throw std::runtime_error("Failed to initialize Coarse Solver.");
}
void Level::coarseSolveInPlace(Vector<double>& x){
    if (!coarseSolver_) throw std::runtime_error("Coarse Solver not initialized.");
    coarseSolver_->solveInPlace(x);
}


void Level::initializeSmoothing(const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, const int maxOpenMPThreads, const int openMPTaskThreads){
    smoother_ = std::make_unique<Smoother>(*grid_, *level_data_, domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);
    if (!smoother_) throw std::runtime_error("Failed to initialize Smoother.");
}
void Level::smoothingInPlace(Vector<double>& x, Vector<double>& temp_rhs){
    if (!smoother_) throw std::runtime_error("Smoother not initialized.");
    smoother_ -> smoothing(x, temp_rhs);
}

// void Level::initializeExtrapolatedSmoothing(const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, const int maxOpenMPThreads, const int openMPTaskThreads){}
// void Level::extrapolatedSmoothing(Vector<double>& result, const Vector<double>& x){}


const LevelCache& Level::levelData() const {
    return *level_data_;
}