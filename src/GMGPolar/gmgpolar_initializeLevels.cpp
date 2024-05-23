#include "../../include/GMGPolar/gmgpolar.h"

PolarGrid GMGPolar::createFinestGrid() {
    PolarGrid finest_grid;
    if(load_grid_file) {
        assert(!file_grid_r.empty() && !file_grid_theta.empty());
        finest_grid = PolarGrid(file_grid_r, file_grid_theta);
    } else {
        finest_grid = PolarGrid(R0, Rmax, nr_exp, ntheta_exp, r_jump_, anisotropic_factor, alpha, divideBy2);
    }
    if(write_grid_file) {
        const int precision = 18;
        finest_grid.writeToFile("_radii_divisions.txt", "_theta_divisions.txt", precision);
    }
    return finest_grid;
}


int GMGPolar::numberOfLevels(const PolarGrid& finestGrid) {
    // Minimum number of radial nodes and angular divisions
    const int minRadialNodes = 3;
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

    // Determine the number of levels as the minimum of radial maximum level, angular maximum level, 
    // and the maximum levels specified.
    int levels = std::min(radialMaxLevel, angularMaxLevel);
    if(maxLevels > 0) levels = std::min(maxLevels, levels);

    // Check if levels is less than Multigrid minimum level and throw an error
    if (levels < multigridMinLevel) {
        throw std::runtime_error("Number of possible levels is less than Multigrid minimum level");
    }

    return levels;
}

void GMGPolar::initializeLevels(const int levels, PolarGrid& finest_grid) {
    const auto exactFunctions = selectExactFunctionsClass();

    levels_.reserve(levels);

    levels_.emplace_back(0, std::make_unique<PolarGrid>(std::move(finest_grid)), exactFunctions);
    levels_[0].setOperator(*this);

    for(int current_level = 1; current_level < levels; current_level++) {
        levels_.emplace_back(current_level, coarseningGrid(levels_[current_level-1].grid()), exactFunctions);
        levels_[current_level].setOperator(*this);
    }

    if(write_grid_file) {
        const int precision = 18;
        levels_.back().grid().writeToFile("_radii_divisions_coarse.txt", "_theta_divisions_coarse.txt", precision);
    }
}