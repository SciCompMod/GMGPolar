#include "../../include/GMGPolar/gmgpolar.h"

void GMGPolar::setup() {
    auto start_setup = std::chrono::high_resolution_clock::now();

    auto start_setup_createLevels = std::chrono::high_resolution_clock::now();

    // --------------------------------
    // Create the finest mesh (level 0)
    auto finest_grid = std::make_unique<PolarGrid>(createFinestGrid()); /* Implementation below */
    std::cout << "System of size (nr x ntheta) = (" << finest_grid->nr() << " x " << finest_grid->ntheta() << ")\n";
    std::cout << "on the coordinates (r x theta): (" << R0_ << ", " << Rmax_ << ") x (" << 0 << ", " << 2 * M_PI << ")\n";

    // ----------------------------------------------------------
    // Building PolarGrid and LevelCache for all multigrid levels
    numberOflevels_ = chooseNumberOfLevels(*finest_grid); /* Implementation below */
    levels_.reserve(numberOflevels_);

    int current_level = 0;
    auto finest_levelCache = std::make_unique<LevelCache>(*finest_grid);
    levels_.emplace_back(current_level, std::move(finest_grid), std::move(finest_levelCache));

    for(current_level = 1; current_level < numberOflevels_; current_level++) {
        auto current_grid = std::make_unique<PolarGrid>(coarseningGrid(levels_[current_level-1].grid()));
        auto current_levelCache = std::make_unique<LevelCache>(levels_[current_level-1], *current_grid);
        levels_.emplace_back(current_level, std::move(current_grid), std::move(current_levelCache));
    }

    auto end_setup_createLevels = std::chrono::high_resolution_clock::now();
    t_setup_createLevels += std::chrono::duration<double>(end_setup_createLevels - start_setup_createLevels).count();

    if(write_grid_file_) {
        const int precision = 18;
        levels_.back().grid().writeToFile("_radii_coarse.txt", "_angles_coarse.txt", precision);
    }

    // -----------------------------------------------------------
    // Initializing the optimal number of threads for OpenMP tasks 
    taskingThreads_.resize(numberOflevels_, maxOpenMPThreads_);
    if(finestLevelThreads_ > 0){
        for (int current_level = 0; current_level < numberOflevels_; current_level++){
            taskingThreads_[current_level] = std::max(1, std::min(maxOpenMPThreads_, 
                static_cast<int>(std::floor(finestLevelThreads_ * std::pow(threadReductionFactor_, current_level)))
            ));
        }
    }

    interpolation_ = std::make_unique<Interpolation>(maxOpenMPThreads_, taskingThreads_);

    auto start_setup_rhs = std::chrono::high_resolution_clock::now();

    // ------------------------------------- //
    // Build rhs_f on Level 0 (finest Level) //
    build_rhs_f(levels_[0], levels_[0].rhs());

    auto end_setup_rhs = std::chrono::high_resolution_clock::now();
    t_setup_rhs += std::chrono::duration<double>(end_setup_rhs - start_setup_rhs).count();

    // -------------------------------------------------------
    // Initializing various operators based on the level index
    for (int current_level = 0; current_level < numberOflevels_; current_level++){
        // ---------------------- //
        // Level 0 (finest Level) //
        // ---------------------- //
        if(current_level == 0){
            if(extrapolation_){
                auto start_setup_smoother = std::chrono::high_resolution_clock::now();
                levels_[current_level].initializeExtrapolatedSmoothing(*domain_geometry_, *system_parameters_, DirBC_Interior_, 
                    maxOpenMPThreads_, taskingThreads_[current_level]
                );

                levels_[current_level].initializeSmoothing(*domain_geometry_, *system_parameters_, DirBC_Interior_, 
                    maxOpenMPThreads_, taskingThreads_[current_level]
                );

                auto end_setup_smoother = std::chrono::high_resolution_clock::now();
                t_setup_smoother += std::chrono::duration<double>(end_setup_smoother - start_setup_smoother).count();
            } else{
                auto start_setup_smoother = std::chrono::high_resolution_clock::now();
                levels_[current_level].initializeSmoothing(*domain_geometry_, *system_parameters_, DirBC_Interior_, 
                    maxOpenMPThreads_, taskingThreads_[current_level]
                );
                auto end_setup_smoother = std::chrono::high_resolution_clock::now();
                t_setup_smoother += std::chrono::duration<double>(end_setup_smoother - start_setup_smoother).count();
            }
            levels_[current_level].initializeResidual(*domain_geometry_, *system_parameters_, DirBC_Interior_, 
                maxOpenMPThreads_, taskingThreads_[current_level]
            );
        }
        // -------------------------- //
        // Level n-1 (coarsest Level) //
        // -------------------------- //
        else if(current_level == numberOflevels_ - 1){
            auto start_setup_directSolver = std::chrono::high_resolution_clock::now();
            levels_[current_level].initializeDirectSolver(*domain_geometry_, *system_parameters_, DirBC_Interior_, 
                maxOpenMPThreads_, taskingThreads_[current_level]
            );
            auto end_setup_directSolver = std::chrono::high_resolution_clock::now();
            t_setup_directSolver += std::chrono::duration<double>(end_setup_directSolver - start_setup_directSolver).count();
        }
        // ------------------- //
        // Intermediate levels //
        // ------------------- //
        else{
            auto start_setup_smoother = std::chrono::high_resolution_clock::now();
            levels_[current_level].initializeSmoothing(*domain_geometry_, *system_parameters_, DirBC_Interior_, 
                maxOpenMPThreads_, taskingThreads_[current_level]
            );
            auto end_setup_smoother = std::chrono::high_resolution_clock::now();
            t_setup_smoother += std::chrono::duration<double>(end_setup_smoother - start_setup_smoother).count();

            levels_[current_level].initializeResidual(*domain_geometry_, *system_parameters_, DirBC_Interior_, 
                maxOpenMPThreads_, taskingThreads_[current_level]
            );
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
        const double& refinement_radius = system_parameters_->getAlphaJump(); /* Radius of anisotropic grid refinement */
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
    const int minAngularDivisions = 8;

    // Minimum level for Multigrid
    const int multigridMinLevel = extrapolation_ ? 3 : 2;

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

    /* Currently unused */
    const int linear_complexity_levels = 
        std::ceil( (2.0 * std::log(static_cast<double>(finestGrid.number_of_nodes())) - std::log(3.0)) / (3.0 * std::log(4.0)));

    // Determine the number of levels as the minimum of radial maximum level, angular maximum level, 
    // and the maximum levels specified.
    int levels = std::min(radialMaxLevel, angularMaxLevel);
    if(maxLevels_ > 0) levels = std::min(maxLevels_, levels);

    // Check if levels is less than Multigrid minimum level and throw an error
    if (levels < multigridMinLevel) {
        throw std::runtime_error("Number of possible levels is less than Multigrid minimum level");
    }

    return levels;
}

