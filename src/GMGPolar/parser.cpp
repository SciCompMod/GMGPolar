#include "../../include/GMGPolar/gmgpolar.h"

/* Defines if user input is required */
enum {
    OPTIONAL = 0,
    REQUIRED = 1
};

void GMGPolar::parseGrid() {
    R0 = parser_.get<double>("R0");
    Rmax = parser_.get<double>("Rmax");
    nr_exp = parser_.get<int>("nr_exp");
    ntheta_exp = parser_.get<int>("ntheta_exp");
    anisotropic_factor = parser_.get<int>("anisotropic_factor");
    divideBy2 = parser_.get<int>("divideBy2");
    write_grid_file = parser_.get<int>("write_grid_file") != 0;
    load_grid_file = parser_.get<int>("load_grid_file") != 0;
    file_grid_r = parser_.get<std::string>("file_grid_r");
    file_grid_theta = parser_.get<std::string>("file_grid_theta");
}

void GMGPolar::parseGeometry() {
    DirBC_Interior = parser_.get<int>("DirBC_Interior") != 0;
}

void GMGPolar::parseMultigrid() {
    extrapolation = parser_.get<int>("extrapolation");
    maxLevels = parser_.get<int>("maxLevels");
    v1 = parser_.get<int>("v1");
    v2 = parser_.get<int>("v2");
    cycle = parser_.get<int>("cycle");
}

void GMGPolar::parseGeneral() {
    maxOpenMPThreads = parser_.get<int>("maxOpenMPThreads");
    finestLevelThreads = parser_.get<int>("finestLevelThreads");
    threadReductionFactor = parser_.get<double>("threadReductionFactor");
    omp_set_num_threads(maxOpenMPThreads);
}

void GMGPolar::initializeGrid() {
    parser_.add<double>(
        "R0", 'r', "Interior radius of the disk", 
        OPTIONAL, 1e-5
    );
    parser_.add<double>(
        "Rmax", 'R', "Exterior radius of the disk", 
        OPTIONAL, 1.3
    );
    parser_.add<int>(
        "nr_exp", 'n', 
        "Number of nodes (exponents) in the radial direction.", 
        OPTIONAL, 4
    );
    parser_.add<int>(
        "ntheta_exp", '\0', 
        "Number of nodes (exponents) in the angular direction.", 
        OPTIONAL, -1
    );
    parser_.add<int>(
        "anisotropic_factor", '\0', 
        "Defines anisotropic discretization in r-direction.", 
        OPTIONAL, 0
    );
    parser_.add<int>(
        "divideBy2", '\0', 
        "Refines the grid globally divideBy2 times.", 
        OPTIONAL, 0
    );
    parser_.add<int>(
        "write_grid_file", '\0', "Enable writing the finest PolarGrid to a file.",
        OPTIONAL, 0, cmdline::oneof(0,1)
    );
    parser_.add<int>(
        "load_grid_file", '\0', "Enable loading the finest PolarGrid from a file.",
        OPTIONAL, 0, cmdline::oneof(0,1)
    );
    parser_.add<std::string>(
        "file_grid_r", '\0', "Path to the file containing radii values for grid divisions in the r-direction.",
        OPTIONAL, ""
    );
    parser_.add<std::string>(
        "file_grid_theta", '\0', "Path to the file containing theta values for grid divisions in the theta-direction.",
        OPTIONAL, ""
    );
}

void GMGPolar::initializeGeometry() {
    parser_.add<int>(
        "DirBC_Interior", '\0', "Defines the boundary condition on the interior circle. Across-origin(0), Dirichlet-boundary(1).",
        OPTIONAL, 0, cmdline::oneof(0,1)
    );
}


void GMGPolar::initializeMultigrid() {
    parser_.add<int>(
        "extrapolation", 'e', 
        "Specifies if extrapolation is used used.", 
        OPTIONAL, 0, cmdline::oneof(0,1)
    );
    parser_.add<int>(
        "maxLevels", 'l', 
        "Defines the maximum number of levels used in the multigrid scheme.", 
        OPTIONAL, -1
    );
    parser_.add<int>(
        "v1", '\0', 
        "Defines the number of pre-smoothing steps.", 
        OPTIONAL, 1
    );
    parser_.add<int>(
        "v2", '\0', 
        "Defines the number of post-smoothing steps.", 
        OPTIONAL, 1
    );
    // Type of Multigrid Cycle
    // - 1: V-cycle (default setting)
    // - 2: W-cycle
    parser_.add<int>(
        "cycle", '\0', 
        "Type of Multigrid Cycle.", 
        OPTIONAL, 1
    );
}

void GMGPolar::initializeGeneral() {
    parser_.add<int>(
        "maxOpenMPThreads", '\0', 
        "Defines the maximum number of OpenMP threads used.", 
        REQUIRED, 1
    );
    parser_.add<int>(
        "finestLevelThreads", '\0', 
        "Optimal number of OpenMP threads for the finest level, can exceed maxOpenMPThreads.", 
        OPTIONAL, -1
    );
    parser_.add<double>(
        "threadReductionFactor", '\0', 
        "Thread reduction factor to coarser levels.", 
        OPTIONAL, 1.0
    );
}
