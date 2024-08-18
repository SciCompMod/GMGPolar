#include "../../include/GMGPolar/gmgpolar.h"

void GMGPolar::setup() {
    resetTimings();
    auto start_setup = std::chrono::high_resolution_clock::now();

    auto start_setup_createLevels = std::chrono::high_resolution_clock::now();

    // --------------------------------
    // Create the finest mesh (level 0)
    auto finest_grid = std::make_unique<PolarGrid>(createFinestGrid()); /* Implementation below */
    std::cout << "System of size (nr x ntheta) = (" << finest_grid->nr() << " x " << finest_grid->ntheta() << ")\n";
    std::cout << "on the coordinates (r x theta): (" << R0_ << ", " << Rmax_ << ") x (" << 0 << ", " << 2 * M_PI << ")\n";

    // ----------------------------------------------------------
    // Building PolarGrid and LevelCache for all multigrid levels
    number_of_levels_ = chooseNumberOfLevels(*finest_grid); /* Implementation below */
    levels_.clear(); levels_.reserve(number_of_levels_);

    int current_level = 0;
    auto finest_levelCache = std::make_unique<LevelCache>(*finest_grid, *density_profile_coefficients_);
    levels_.emplace_back(current_level, std::move(finest_grid), std::move(finest_levelCache), extrapolation_);

    for(current_level = 1; current_level < number_of_levels_; current_level++) {
        auto current_grid = std::make_unique<PolarGrid>(coarseningGrid(levels_[current_level-1].grid()));
        auto current_levelCache = std::make_unique<LevelCache>(levels_[current_level-1], *current_grid);
        levels_.emplace_back(current_level, std::move(current_grid), std::move(current_levelCache), extrapolation_);
    }

    auto end_setup_createLevels = std::chrono::high_resolution_clock::now();
    t_setup_createLevels += std::chrono::duration<double>(end_setup_createLevels - start_setup_createLevels).count();

    if(write_grid_file_) {
        const int precision = 18;
        levels_.back().grid().writeToFile("_radii_coarse.txt", "_angles_coarse.txt", precision);
    }

    // -----------------------------------------------------------
    // Initializing the optimal number of threads for OpenMP tasks 
    threads_per_level_.resize(number_of_levels_, max_omp_threads_);

    for (int current_level = 0; current_level < number_of_levels_; current_level++){
        threads_per_level_[current_level] = std::max(1, std::min(max_omp_threads_, 
            static_cast<int>(std::floor(max_omp_threads_ * std::pow(thread_reduction_factor_, current_level)))
        ));
    }
    
    interpolation_ = std::make_unique<Interpolation>(threads_per_level_);

    auto start_setup_rhs = std::chrono::high_resolution_clock::now();

    // ------------------------------------- //
    // Build rhs_f on Level 0 (finest Level) //
    build_rhs_f(levels_[0], levels_[0].rhs());

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
            if(extrapolation_ > 0){
                levels_[current_level].initializeExtrapolatedSmoothing(*domain_geometry_, DirBC_Interior_, threads_per_level_[current_level]);
                if(extrapolation_ > 1){
                    levels_[current_level].initializeSmoothing(*domain_geometry_, DirBC_Interior_, threads_per_level_[current_level]);
                }
            } else{
                levels_[current_level].initializeSmoothing(*domain_geometry_, DirBC_Interior_, threads_per_level_[current_level]);
            }
            auto end_setup_smoother = std::chrono::high_resolution_clock::now();
            t_setup_smoother += std::chrono::duration<double>(end_setup_smoother - start_setup_smoother).count();
            levels_[current_level].initializeResidual(*domain_geometry_, DirBC_Interior_, threads_per_level_[current_level]);
        }
        // -------------------------- //
        // Level n-1 (coarsest Level) //
        // -------------------------- //
        else if(current_level == number_of_levels_ - 1){
            auto start_setup_directSolver = std::chrono::high_resolution_clock::now();
            levels_[current_level].initializeDirectSolver(*domain_geometry_, DirBC_Interior_, threads_per_level_[current_level]);
            auto end_setup_directSolver = std::chrono::high_resolution_clock::now();
            t_setup_directSolver += std::chrono::duration<double>(end_setup_directSolver - start_setup_directSolver).count();
            levels_[current_level].initializeResidual(*domain_geometry_, DirBC_Interior_, threads_per_level_[current_level]);
        }
        // ------------------- //
        // Intermediate levels //
        // ------------------- //
        else{
            auto start_setup_smoother = std::chrono::high_resolution_clock::now();
            levels_[current_level].initializeSmoothing(*domain_geometry_, DirBC_Interior_, threads_per_level_[current_level]);
            auto end_setup_smoother = std::chrono::high_resolution_clock::now();
            t_setup_smoother += std::chrono::duration<double>(end_setup_smoother - start_setup_smoother).count();
            levels_[current_level].initializeResidual(*domain_geometry_, DirBC_Interior_, threads_per_level_[current_level]);
        }
    }
    auto end_setup = std::chrono::high_resolution_clock::now();
    t_setup_total += std::chrono::duration<double>(end_setup - start_setup).count();
}


PolarGrid GMGPolar::createFinestGrid() {
    PolarGrid finest_grid;

    if(load_grid_file_) {
        assert(!file_grid_radii_.empty() && !file_grid_angles_.empty());
        finest_grid = PolarGrid(file_grid_radii_, file_grid_angles_);
    } else {
        const double& refinement_radius = density_profile_coefficients_->getAlphaJump(); /* Radius of anisotropic grid refinement */
        std::optional<double> splitting_radius = std::nullopt; /* (Automatic) line splitting radius for the smoother */
        finest_grid = PolarGrid(R0_, Rmax_, nr_exp_, ntheta_exp_, refinement_radius, anisotropic_factor_, divideBy2_, splitting_radius);
    }
    if(write_grid_file_) {
        const int precision = 18;
        finest_grid.writeToFile(file_grid_radii_, file_grid_angles_, precision);
    }
    return finest_grid;
}


int GMGPolar::chooseNumberOfLevels(const PolarGrid& finestGrid) {
    const int minRadialNodes = 5;
    const int minAngularDivisions = 4;

    // Minimum level for Multigrid
    const int multigridMinLevel = extrapolation_ ? 2 : 2;

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
    t_solve_multigrid_iterations = 0.0;
    t_check_convergence = 0.0;
    t_check_exact_error = 0.0;

    t_avg_MGC_total = 0.0;
    t_avg_MGC_preSmoothing = 0.0;
    t_avg_MGC_postSmoothing = 0.0;
    t_avg_MGC_residual = 0.0;
    t_avg_MGC_directSolver = 0.0;
}