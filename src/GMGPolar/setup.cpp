#include "../../include/GMGPolar/gmgpolar.h"

#include <chrono>
#include <iostream>

// clang-format off
void GMGPolar::setup() {
    auto start_setup = std::chrono::high_resolution_clock::now();

    resetTimings();

    auto start_setup_createLevels = std::chrono::high_resolution_clock::now();

    // -------------------------------- //
    // Create the finest mesh (level 0) //
    // -------------------------------- //
    auto finest_grid = std::make_unique<PolarGrid>(createFinestGrid()); /* Implementation below */
    if(verbose_ > 0) {
        std::cout << "System of size (nr x ntheta) = (" << finest_grid->nr() << " x " << finest_grid->ntheta() << ")\n";
        std::cout << "on the coordinates (r x theta): (" << R0_ << ", " << Rmax_ << ") x (" << 0 << ", " << 2 * M_PI << ")\n";

        std::cout << "Anisotropy factor: " << anisotropic_factor_ << std::endl;
        std::cout << "Dirichlet boundary (interior): " << DirBC_Interior_ << std::endl; 
    }
    if(paraview_) writeToVTK("output_finest_grid", *finest_grid);

    // ----------------------------------------- //
    // Definining the number of multigrid levels //
    // ----------------------------------------- //

    number_of_levels_ = chooseNumberOfLevels(*finest_grid); /* Implementation below */
    assert(number_of_levels_ >= 2);

    int actual_gpu_levels;
    if(gpu_levels_ < 0 || gpu_levels_ >= number_of_levels_){
        actual_gpu_levels = number_of_levels_ - 1;
    } else{
        actual_gpu_levels = gpu_levels_;
    }

    if(verbose_ > 0) std::cout <<"Number of levels: "<<number_of_levels_<<"\n";
    if(verbose_ > 0) std::cout <<"Number of GPU levels: "<<actual_gpu_levels<<"\n";

    // ---------------------------------------------------------- //
    // Building PolarGrid and LevelCache for all multigrid levels //
    // ---------------------------------------------------------- //

    levels_.clear(); levels_.reserve(number_of_levels_);

    /* Example: number_of_levels_ = 5, actual_gpu_levels = 2. */
    /* processing_type = [GPU, GPU; CPU_HYBRID, CPU, CPU]; */
    std::vector<ProcessingType> processing_type(number_of_levels_);
    for (int level = 0; level < number_of_levels_; level++){
        if (level < actual_gpu_levels) {
            processing_type[level] = ProcessingType::GPU;
        } else if (level == actual_gpu_levels) {
            processing_type[level] = ProcessingType::CPU_HYBRID;
        } else {
            processing_type[level] = ProcessingType::CPU;
        }
    }
    /* If no GPU levels exist, CPU_Hybrid is set to CPU. */
    if(processing_type[0] != ProcessingType::GPU){
        for (int level = 0; level < number_of_levels_; level++)
            processing_type[level] = ProcessingType::CPU;
    }

    int finest_level = 0;
    auto finest_levelCache = std::make_unique<LevelCache>(processing_type[finest_level], *finest_grid, *density_profile_coefficients_, *domain_geometry_);
    levels_.emplace_back(finest_level, processing_type[finest_level], std::move(finest_grid), std::move(finest_levelCache), extrapolation_, FMG_);

    for(int current_level = 1; current_level < number_of_levels_; current_level++) {
        auto current_grid = std::make_unique<PolarGrid>(coarseningGrid(levels_[current_level-1].grid()));
        auto current_levelCache = std::make_unique<LevelCache>(processing_type[current_level], *current_grid, *density_profile_coefficients_, *domain_geometry_);
        levels_.emplace_back(current_level, processing_type[current_level], std::move(current_grid), std::move(current_levelCache), extrapolation_, FMG_);
    }

    auto end_setup_createLevels = std::chrono::high_resolution_clock::now();
    t_setup_createLevels += std::chrono::duration<double>(end_setup_createLevels - start_setup_createLevels).count();

    if(paraview_) writeToVTK("output_coarsest_grid", levels_.back().grid());

    if(verbose_ > 0) std::cout <<"Maxmimum number of threads: "<<max_omp_threads_<<"\n";
    
    interpolation_ = std::make_unique<Interpolation>(DirBC_Interior_);

    auto start_setup_rhs = std::chrono::high_resolution_clock::now();

    // ------------------------------------- //
    // Build rhs_f on Level 0 (finest Level) //
    // ------------------------------------- //
    if(levels_[0].processingType() == ProcessingType::GPU){
        build_rhs_f(levels_[0], levels_[0].GPU_rhs());
    } else{
        build_rhs_f(levels_[0], levels_[0].rhs());
    }

    /* ---------------- */
    /* Discretize rhs_f */
    /* ---------------- */
    int initial_rhs_f_levels = FMG_ ? number_of_levels_ : (extrapolation_ == ExtrapolationType::NONE ? 1 : 2);
    /* Loop through the levels, injecting and discretizing rhs */
    for (int level_idx = 0; level_idx < initial_rhs_f_levels; level_idx++) 
    {
        Level& current_level = levels_[level_idx];
        // Inject rhs if there is a next level
        if (level_idx + 1 < initial_rhs_f_levels) {
            Level& next_level = levels_[level_idx + 1];

            /* Injection */
            if(current_level.processingType() == ProcessingType::GPU){
                assert(next_level.processingType() != ProcessingType::CPU);
                injection(level_idx, next_level.GPU_rhs(), current_level.GPU_rhs());
            } else{
                assert(next_level.processingType() == ProcessingType::CPU);
                injection(level_idx, next_level.rhs(), current_level.rhs());
            }
            if(next_level.processingType() == ProcessingType::CPU_HYBRID){
                copyDeviceToHost(next_level.GPU_rhs(), next_level.rhs());
            }

        }
        /* Discretize the rhs for the current level */
        if(current_level.processingType() == ProcessingType::GPU){
            discretize_rhs_f(current_level, current_level.GPU_rhs());
        } else{
            discretize_rhs_f(current_level, current_level.rhs());
        }
    } 

    auto end_setup_rhs = std::chrono::high_resolution_clock::now();
    t_setup_rhs += std::chrono::duration<double>(end_setup_rhs - start_setup_rhs).count();

    // -------------------------------------------------------
    // Initializing various operators based on the level index
    for (int current_level = 0; current_level < number_of_levels_; current_level++){
        // ---------------------- //
        // Level 0 (finest Level) //
        // ---------------------- //
        if(current_level == 0){
            auto start_setup_smoother = std::chrono::high_resolution_clock::now();
            switch(extrapolation_) {
                case ExtrapolationType::NONE:
                    full_grid_smoothing_ = true;
                    levels_[current_level].initializeSmoothing(*domain_geometry_, *density_profile_coefficients_, DirBC_Interior_);
                    break;
                case ExtrapolationType::IMPLICIT_EXTRAPOLATION:
                    full_grid_smoothing_ = false;
                    levels_[current_level].initializeExtrapolatedSmoothing(*domain_geometry_, *density_profile_coefficients_, DirBC_Interior_);
                    break;
                case ExtrapolationType::IMPLICIT_FULL_GRID_SMOOTHING:
                    full_grid_smoothing_ = true;
                    levels_[current_level].initializeSmoothing(*domain_geometry_, *density_profile_coefficients_, DirBC_Interior_);
                    break;
                case ExtrapolationType::COMBINED:
                    full_grid_smoothing_ = true;
                    levels_[current_level].initializeSmoothing(*domain_geometry_, *density_profile_coefficients_, DirBC_Interior_);
                    levels_[current_level].initializeExtrapolatedSmoothing(*domain_geometry_, *density_profile_coefficients_, DirBC_Interior_);
                    break;
                default:
                    full_grid_smoothing_ = false;
                    levels_[current_level].initializeSmoothing(*domain_geometry_, *density_profile_coefficients_, DirBC_Interior_);
                    levels_[current_level].initializeExtrapolatedSmoothing(*domain_geometry_, *density_profile_coefficients_, DirBC_Interior_);
                    break;
            }
            auto end_setup_smoother = std::chrono::high_resolution_clock::now();
            t_setup_smoother += std::chrono::duration<double>(end_setup_smoother - start_setup_smoother).count();
            levels_[current_level].initializeResidual(*domain_geometry_, *density_profile_coefficients_, DirBC_Interior_);
        }
        // -------------------------- //
        // Level n-1 (coarsest Level) //
        // -------------------------- //
        else if(current_level == number_of_levels_ - 1){
            assert(levels_[current_level].processingType() != ProcessingType::GPU);
            auto start_setup_directSolver = std::chrono::high_resolution_clock::now();
            levels_[current_level].initializeDirectSolver(*domain_geometry_, *density_profile_coefficients_, DirBC_Interior_);
            auto end_setup_directSolver = std::chrono::high_resolution_clock::now();
            t_setup_directSolver += std::chrono::duration<double>(end_setup_directSolver - start_setup_directSolver).count();
            levels_[current_level].initializeResidual(*domain_geometry_, *density_profile_coefficients_, DirBC_Interior_);
        }
        // ------------------- //
        // Intermediate levels //
        // ------------------- //
        else{
            auto start_setup_smoother = std::chrono::high_resolution_clock::now();
            levels_[current_level].initializeSmoothing(*domain_geometry_, *density_profile_coefficients_, DirBC_Interior_);
            auto end_setup_smoother = std::chrono::high_resolution_clock::now();
            t_setup_smoother += std::chrono::duration<double>(end_setup_smoother - start_setup_smoother).count();
            levels_[current_level].initializeResidual(*domain_geometry_, *density_profile_coefficients_, DirBC_Interior_);
        }
    }

    auto end_setup = std::chrono::high_resolution_clock::now();
    t_setup_total += std::chrono::duration<double>(end_setup - start_setup).count();
}



PolarGrid GMGPolar::createFinestGrid() {
    const double& refinement_radius = density_profile_coefficients_->getAlphaJump(); /* Radius of anisotropic grid refinement */
    std::optional<double> splitting_radius = std::nullopt; /* (Automatic) line splitting radius for the smoother */
    PolarGrid finest_grid(R0_, Rmax_, nr_exp_, ntheta_exp_, refinement_radius, anisotropic_factor_, divideBy2_, splitting_radius);
    return finest_grid;
}

int GMGPolar::chooseNumberOfLevels(const PolarGrid& finestGrid) {
    const int minRadialNodes = 5;
    const int minAngularDivisions = 4;

    // Minimum level for Multigrid
    const int multigridMinLevel = 2;

    // Calculate radial maximum level
    int radialNodes = finestGrid.nr();
    int radialMaxLevel = 1;
    while ((radialNodes + 1) / 2 >= minRadialNodes && (radialNodes + 1) % 2 == 0) {
        radialNodes = (radialNodes + 1) / 2;
        radialMaxLevel++;
    }

    // Calculate angular maximum level
    int angularDivisions = finestGrid.ntheta();
    int angularMaxLevel = 1;
    while (angularDivisions / 2 >= minAngularDivisions && angularDivisions % 2 == 0 && (angularDivisions/2) % 2 == 0) {
        angularDivisions = angularDivisions / 2;
        angularMaxLevel++;
    }

    /* Currently unused: Number of levels which guarantee linear scalability */
    const int linear_complexity_levels = 
        std::ceil( (2.0 * std::log(static_cast<double>(finestGrid.numberOfNodes())) - std::log(3.0)) / (3.0 * std::log(4.0)));

    // Determine the number of levels as the minimum of radial maximum level, angular maximum level, 
    // and the maximum levels specified.
    int levels = std::min(radialMaxLevel, angularMaxLevel);
    if(max_levels_ > 0) levels = std::min(max_levels_, levels);

    // Check if levels is less than Multigrid minimum level and throw an error
    if (levels < multigridMinLevel) {
        throw std::runtime_error("Number of possible levels is less than Multigrid minimum level");
    }

    return levels;
}

void GMGPolar::resetTimings(){
    t_setup_total = 0.0;
    t_setup_createLevels = 0.0;
    t_setup_rhs = 0.0;
    t_setup_smoother = 0.0;
    t_setup_directSolver = 0.0;

    t_solve_total = 0.0;
    t_solve_initial_approximation = 0.0;
    t_solve_multigrid_iterations = 0.0;
    t_check_convergence = 0.0;
    t_check_exact_error = 0.0;

    t_avg_MGC_total = 0.0;
    t_avg_MGC_preSmoothing = 0.0;
    t_avg_MGC_postSmoothing = 0.0;
    t_avg_MGC_residual = 0.0;
    t_avg_MGC_directSolver = 0.0;
}