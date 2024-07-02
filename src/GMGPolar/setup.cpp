#include "../../include/GMGPolar/gmgpolar.h"

void GMGPolar::setup() {
    // -------------------------------
    // Create the finest mesh (level 0)
    auto finest_grid = std::make_unique<PolarGrid>(createFinestGrid());
    std::cout << "System of size (nr x ntheta) = (" << finest_grid->nr() << " x " << finest_grid->ntheta() << ")\n";
    std::cout << "on the coordinates (r x theta): (" << R0 << ", " << Rmax << ") x (" << 0 << ", " << 2 * M_PI << ")\n";

    // ----------------------------------------------------------
    // Building PolarGrid and LevelCache for all multigrid levels
    numberOflevels_ = chooseNumberOfLevels(*finest_grid);
    levels_.reserve(numberOflevels_);

    int current_level = 0;
    auto finest_levelCache = std::make_unique<LevelCache>(*finest_grid);
    levels_.emplace_back(current_level, std::move(finest_grid), std::move(finest_levelCache));

    for(current_level = 1; current_level < numberOflevels_; current_level++) {
        auto current_grid = std::make_unique<PolarGrid>(coarseningGrid(levels_[current_level-1].grid()));
        auto current_levelCache = std::make_unique<LevelCache>(levels_[current_level-1], *current_grid);
        levels_.emplace_back(current_level, std::move(current_grid), std::move(current_levelCache));
    }

    if(write_grid_file) {
        const int precision = 18;
        levels_.back().grid().writeToFile("_radii_coarse.txt", "_angles_coarse.txt", precision);
    }

    // -----------------------------------------------------------
    // Initializing the optimal number of threads for OpenMP tasks 
    taskingThreads_.resize(numberOflevels_, maxOpenMPThreads);
    if(finestLevelThreads > 0){
        for (int current_level = 0; current_level < numberOflevels_; current_level++){
            taskingThreads_[current_level] = std::max(1, std::min(maxOpenMPThreads, 
                static_cast<int>(std::floor(finestLevelThreads * std::pow(threadReductionFactor, current_level)))
            ));
        }
    }

    interpolation_ = std::make_unique<Interpolation>(maxOpenMPThreads, taskingThreads_);

    // -------------------------------------------------------
    // Initializing various operators based on the level index
    for (int current_level = 0; current_level < numberOflevels_; current_level++){
        // ---------------------- //
        // Level 0 (finest Level) //
        // ---------------------- //
        if(current_level == 0){
            if(extrapolation){
                // levels_[current_level].initializeExtrapolatedSmoothing(domain_geometry_, system_parameters_, DirBC_Interior, numThreadsUsed);
            } else{
                levels_[current_level].initializeSmoothing(domain_geometry_, system_parameters_, DirBC_Interior, 
                    maxOpenMPThreads, taskingThreads_[current_level]
                );
            }
            levels_[current_level].initializeResidual(domain_geometry_, system_parameters_, DirBC_Interior, 
                maxOpenMPThreads, taskingThreads_[current_level]
            );
        }
        // -------------------------- //
        // Level n-1 (coarsest Level) //
        // -------------------------- //
        else if(current_level == numberOflevels_ - 1){
            levels_[current_level].initializeDirectSolver(domain_geometry_, system_parameters_, DirBC_Interior, 
                maxOpenMPThreads, taskingThreads_[current_level]
            );
        }
        // ------------------- //
        // Intermediate levels //
        // ------------------- //
        else{
            levels_[current_level].initializeSmoothing(domain_geometry_, system_parameters_, DirBC_Interior, 
                maxOpenMPThreads, taskingThreads_[current_level]
            );
            levels_[current_level].initializeResidual(domain_geometry_, system_parameters_, DirBC_Interior, 
                maxOpenMPThreads, taskingThreads_[current_level]
            );
        }
    }
}

PolarGrid GMGPolar::createFinestGrid() {
    PolarGrid finest_grid;
    if(load_grid_file) {
        assert(!file_grid_r.empty() && !file_grid_theta.empty());
        finest_grid = PolarGrid(file_grid_r, file_grid_theta);
    } else {
        const double& refinement_radius = system_parameters_.getAlphaJump();
        std::optional<double> splitting_radius = std::nullopt; // Automatic smoother splitting radius
        finest_grid = PolarGrid(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, splitting_radius);
    }
    if(write_grid_file) {
        const int precision = 18;
        finest_grid.writeToFile(file_grid_r, file_grid_theta, precision);
    }
    return finest_grid;
}

int GMGPolar::chooseNumberOfLevels(const PolarGrid& finestGrid) {
    const int minRadialNodes = 4;
    const int minAngularDivisions = 4;

    // Minimum level for Multigrid
    const int multigridMinLevel = extrapolation ? 3 : 2;

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

    const int linear_complexity_levels = 
        std::ceil( (2.0 * std::log(static_cast<double>(finestGrid.number_of_nodes())) - std::log(3.0)) / (3.0 * std::log(4.0)));

    // Determine the number of levels as the minimum of radial maximum level, angular maximum level, 
    // and the maximum levels specified.
    int levels = std::min(radialMaxLevel, angularMaxLevel);
    if(maxLevels > 0) levels = std::min(maxLevels, levels);
    if(maxLevels < 0) levels = std::min(linear_complexity_levels, levels);

    // Check if levels is less than Multigrid minimum level and throw an error
    if (levels < multigridMinLevel) {
        throw std::runtime_error("Number of possible levels is less than Multigrid minimum level");
    }

    return levels;
}

