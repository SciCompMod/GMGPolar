#include "../../include/Level/level.h"

#include "../../include/Residual/ResidualTakeCPU/residual.h"
#include "../../include/Residual/ResidualTakeGPU/residual.h"

#include "../../include/DirectSolver/directSolver.h"

#include "../../include/Smoother/SmootherTakeCPU/smoother.h"
#include "../../include/Smoother/SmootherTakeGPU/smoother.h"

#include "../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTakeCPU/extrapolatedSmoother.h"
#include "../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTakeGPU/extrapolatedSmoother.h"

Level::Level(const int level, const ProcessingType processing_type, std::unique_ptr<const PolarGrid> grid,
             std::unique_ptr<const LevelCache> level_cache, const ExtrapolationType extrapolation, const bool FMG)
    : level_(level)
    , processing_type_(processing_type)
    , grid_(std::move(grid))
    , level_cache_(std::move(level_cache))
    , extrapolation_(extrapolation)
    , FMG_(FMG)
    , device_grid_(nullptr)
{
    const auto num_nodes = grid_->numberOfNodes();
    const auto rhs_size =
        (FMG_ || level_ == 0 || (level_ == 1 && extrapolation_ != ExtrapolationType::NONE)) ? num_nodes : 0;
    const auto error_size = (level_ > 0) ? num_nodes : 0;

    if (processing_type == ProcessingType::CPU || processing_type == ProcessingType::CPU_HYBRID) {
        rhs_              = Vector<double>(rhs_size);
        solution_         = Vector<double>(num_nodes);
        residual_         = Vector<double>(num_nodes);
        error_correction_ = Vector<double>(error_size);
    }
    else {
        rhs_              = Vector<double>(0);
        solution_         = Vector<double>(0);
        residual_         = Vector<double>(0);
        error_correction_ = Vector<double>(0);
    }

    if (processing_type == ProcessingType::GPU || processing_type == ProcessingType::CPU_HYBRID) {
        cudaMalloc(&device_grid_, sizeof(PolarGrid));
        cudaMemcpy(device_grid_, grid_.get(), sizeof(PolarGrid), cudaMemcpyHostToDevice);
        gpu_rhs_              = GPU_Vector<double>(rhs_size);
        gpu_solution_         = GPU_Vector<double>(num_nodes);
        gpu_residual_         = GPU_Vector<double>(num_nodes);
        gpu_error_correction_ = GPU_Vector<double>(error_size);
    }
    else {
        gpu_rhs_              = GPU_Vector<double>(0);
        gpu_solution_         = GPU_Vector<double>(0);
        gpu_residual_         = GPU_Vector<double>(0);
        gpu_error_correction_ = GPU_Vector<double>(0);
    }
}

Level::Level(Level&& other) noexcept
    : level_(other.level_)
    , processing_type_(other.processing_type_)
    , grid_(std::move(other.grid_))
    , device_grid_(other.device_grid_)
    , level_cache_(std::move(other.level_cache_))
    , extrapolation_(other.extrapolation_)
    , FMG_(other.FMG_)
    , rhs_(std::move(other.rhs_))
    , solution_(std::move(other.solution_))
    , residual_(std::move(other.residual_))
    , error_correction_(std::move(other.error_correction_))
    , gpu_rhs_(std::move(other.gpu_rhs_))
    , gpu_solution_(std::move(other.gpu_solution_))
    , gpu_residual_(std::move(other.gpu_residual_))
    , gpu_error_correction_(std::move(other.gpu_error_correction_))
{
    other.device_grid_ = nullptr;
}

// Move Assignment Operator
Level& Level::operator=(Level&& other) noexcept
{
    if (this != &other) {
        // Free existing GPU resources if necessary
        if (device_grid_) {
            cudaFree(device_grid_);
        }

        // Move data members
        level_         = other.level_;
        grid_          = std::move(other.grid_);
        device_grid_   = other.device_grid_;
        level_cache_   = std::move(other.level_cache_);
        extrapolation_ = other.extrapolation_;
        FMG_           = other.FMG_;

        rhs_                  = std::move(other.rhs_);
        solution_             = std::move(other.solution_);
        residual_             = std::move(other.residual_);
        error_correction_     = std::move(other.error_correction_);
        gpu_rhs_              = std::move(other.gpu_rhs_);
        gpu_solution_         = std::move(other.gpu_solution_);
        gpu_residual_         = std::move(other.gpu_residual_);
        gpu_error_correction_ = std::move(other.gpu_error_correction_);

        // Leave the source object in a valid state
        other.device_grid_ = nullptr;
    }
    return *this;
}
// Destructor
Level::~Level()
{
if (device_grid_) {
    cudaPointerAttributes attributes;
    cudaError_t err = cudaPointerGetAttributes(&attributes, device_grid_);
    if (err != cudaSuccess) {
        std::cerr << "Invalid device pointer before free: " << cudaGetErrorString(err) << std::endl;
    }
    cudaFree(device_grid_);
    device_grid_ = nullptr;
}
}

// ---------------- //
// Getter Functions //
int Level::level() const
{
    return level_;
}

const PolarGrid& Level::grid() const
{
    return *grid_;
}

ProcessingType Level::processingType() const
{
    return processing_type_;
}

PolarGrid* Level::device_grid() const
{
    return device_grid_;
}

const LevelCache& Level::levelCache() const
{
    return *level_cache_;
}

Vector<double>& Level::rhs()
{
    return rhs_;
}
const Vector<double>& Level::rhs() const
{
    return rhs_;
}
Vector<double>& Level::solution()
{
    return solution_;
}
const Vector<double>& Level::solution() const
{
    return solution_;
}
Vector<double>& Level::residual()
{
    return residual_;
}
const Vector<double>& Level::residual() const
{
    return residual_;
}
Vector<double>& Level::error_correction()
{
    return error_correction_;
}
const Vector<double>& Level::error_correction() const
{
    return error_correction_;
}

GPU_Vector<double>& Level::GPU_rhs()
{
    return gpu_rhs_;
}
const GPU_Vector<double>& Level::GPU_rhs() const
{
    return gpu_rhs_;
}
GPU_Vector<double>& Level::GPU_solution()
{
    return gpu_solution_;
}
const GPU_Vector<double>& Level::GPU_solution() const
{
    return gpu_solution_;
}
GPU_Vector<double>& Level::GPU_residual()
{
    return gpu_residual_;
}
const GPU_Vector<double>& Level::GPU_residual() const
{
    return gpu_residual_;
}
GPU_Vector<double>& Level::GPU_error_correction()
{
    return gpu_error_correction_;
}
const GPU_Vector<double>& Level::GPU_error_correction() const
{
    return gpu_error_correction_;
}


// -------------- //
// Apply Residual //
void Level::initializeResidual(const DomainGeometry& domain_geometry,
                               const DensityProfileCoefficients& density_profile_coefficients,
                               const bool DirBC_Interior)
{
    if(processing_type_ == ProcessingType::GPU){
        op_residual_GPU_ = std::make_unique<ResidualTakeGPU>(*this, domain_geometry, density_profile_coefficients, DirBC_Interior);
        if (!op_residual_GPU_) throw std::runtime_error("Failed to initialize GPU Residual.");
    }
    else{
        op_residual_CPU_ = std::make_unique<ResidualTakeCPU>(*this, domain_geometry, density_profile_coefficients, DirBC_Interior);
        if (!op_residual_CPU_) throw std::runtime_error("Failed to initialize CPU Residual.");
    }
}

void Level::computeResidual(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const
{
    if (!op_residual_CPU_) throw std::runtime_error("CPU Residual not initialized.");
    op_residual_CPU_->computeResidual(result, rhs, x);
}
void Level::computeResidual(GPU_Vector<double>& result, const GPU_Vector<double>& rhs, const GPU_Vector<double>& x) const
{
    if (!op_residual_GPU_) throw std::runtime_error("GPU Residual not initialized.");
    op_residual_GPU_->computeResidual(result, rhs, x);
}


// ------------------- //
// Solve coarse System //
void Level::initializeDirectSolver(const DomainGeometry& domain_geometry,
                                   const DensityProfileCoefficients& density_profile_coefficients,
                                   const bool DirBC_Interior)
{

    op_directSolver_ = std::make_unique<DirectSolver>(*this, domain_geometry, density_profile_coefficients, DirBC_Interior);
    if (!op_directSolver_) throw std::runtime_error("Failed to initialize Direct Solver.");
}
void Level::directSolveInPlace(Vector<double>& x) const
{
    if (!op_directSolver_) throw std::runtime_error("Coarse Solver not initialized.");
    op_directSolver_->solveInPlace(x);
}


// --------------- //
// Apply Smoothing //
void Level::initializeSmoothing(const DomainGeometry& domain_geometry,
                                const DensityProfileCoefficients& density_profile_coefficients,
                                const bool DirBC_Interior)
{
    if(processing_type_ == ProcessingType::GPU){
        op_smoother_GPU_ = std::make_unique<SmootherTakeGPU>(*this, domain_geometry, density_profile_coefficients, DirBC_Interior);
        if (!op_smoother_GPU_) throw std::runtime_error("Failed to initialize GPU Smoother.");
    }
    else{
        op_smoother_CPU_ = std::make_unique<SmootherTakeCPU>(*this, domain_geometry, density_profile_coefficients, DirBC_Interior);
        if (!op_smoother_CPU_) throw std::runtime_error("Failed to initialize CPU Smoother.");
    }
}

void Level::smoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp) const
{
    if (!op_smoother_CPU_) throw std::runtime_error("CPU Smoother not initialized.");
    op_smoother_CPU_->smoothingInPlace(x, rhs, x);
}
void Level::smoothingInPlace(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp) const
{
    if (!op_smoother_GPU_) throw std::runtime_error("GPU Smoother not initialized.");
    op_smoother_GPU_->smoothingInPlace(x, rhs, temp);
}

// ---------------------------- //
// Apply Extrapolated Smoothing //
void Level::initializeExtrapolatedSmoothing(const DomainGeometry& domain_geometry,
                                const DensityProfileCoefficients& density_profile_coefficients,
                                const bool DirBC_Interior)
{
    if(processing_type_ == ProcessingType::GPU){
        op_extrapolated_smoother_GPU_ = std::make_unique<ExtrapolatedSmootherTakeGPU>(*this, domain_geometry, density_profile_coefficients, DirBC_Interior);
        if (!op_extrapolated_smoother_GPU_) throw std::runtime_error("Failed to initialize GPU Extrapolated Smoother.");
    }
    else{
        op_extrapolated_smoother_CPU_ = std::make_unique<ExtrapolatedSmootherTakeCPU>(*this, domain_geometry, density_profile_coefficients, DirBC_Interior);
        if (!op_extrapolated_smoother_CPU_) throw std::runtime_error("Failed to initialize CPU Extrapolated Smoother.");
    }
}

void Level::extrapolatedSmoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp) const
{
    if (!op_extrapolated_smoother_CPU_) throw std::runtime_error("CPU Extrapolated Smoother not initialized.");
    op_extrapolated_smoother_CPU_->extrapolatedSmoothingInPlace(x, rhs, x);
}
void Level::extrapolatedSmoothingInPlace(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp) const
{
    if (!op_extrapolated_smoother_GPU_) throw std::runtime_error("GPU Extrapolated Smoother not initialized.");
    op_extrapolated_smoother_GPU_->extrapolatedSmoothingInPlace(x, rhs, temp);
}

