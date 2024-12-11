#include "../../include/GMGPolar/gmgpolar.h"

/* Specifies whether user input is required */
enum
{
    OPTIONAL = 0,
    REQUIRED = 1
};

void GMGPolar::parseGrid() {
    R0_ = parser_.get<double>("R0");
    Rmax_ = parser_.get<double>("Rmax");
    nr_exp_ = parser_.get<int>("nr_exp");
    ntheta_exp_ = parser_.get<int>("ntheta_exp");
    anisotropic_factor_ = parser_.get<int>("anisotropic_factor");
    divideBy2_ = parser_.get<int>("divideBy2");
    DirBC_Interior_ = parser_.get<int>("DirBC_Interior") != 0;
}

void GMGPolar::parseMultigrid() {
    FMG_ = parser_.get<int>("FMG") != 0;
    FMG_iterations_ = parser_.get<int>("FMG_iterations");
    const int FMG_cycleValue = parser_.get<int>("FMG_cycle");
    if (FMG_cycleValue == static_cast<int>(MultigridCycleType::V_CYCLE) ||
        FMG_cycleValue == static_cast<int>(MultigridCycleType::W_CYCLE) ||
        FMG_cycleValue == static_cast<int>(MultigridCycleType::F_CYCLE)) {
        FMG_cycle_ = static_cast<MultigridCycleType>(FMG_cycleValue);
    } else {
        throw std::runtime_error("Invalid extrapolation value.\n");
    }
    const int extrapolationValue = parser_.get<int>("extrapolation");
    if (extrapolationValue == static_cast<int>(ExtrapolationType::NONE) ||
        extrapolationValue == static_cast<int>(ExtrapolationType::IMPLICIT_EXTRAPOLATION) ||
        extrapolationValue == static_cast<int>(ExtrapolationType::IMPLICIT_FULL_GRID_SMOOTHING) ||
        extrapolationValue == static_cast<int>(ExtrapolationType::COMBINED)) {
        extrapolation_ = static_cast<ExtrapolationType>(extrapolationValue);
    } else {
        throw std::runtime_error("Invalid extrapolation value.\n");
    }
    max_levels_ = parser_.get<int>("maxLevels");
    gpu_levels_ = parser_.get<int>("gpuLevels");
    pre_smoothing_steps_ = parser_.get<int>("preSmoothingSteps");
    post_smoothing_steps_ = parser_.get<int>("postSmoothingSteps");

    const int cycleValue = parser_.get<int>("multigridCycle");
    if (cycleValue == static_cast<int>(MultigridCycleType::V_CYCLE) ||
        cycleValue == static_cast<int>(MultigridCycleType::W_CYCLE) ||
        cycleValue == static_cast<int>(MultigridCycleType::F_CYCLE)) {
        multigrid_cycle_ = static_cast<MultigridCycleType>(cycleValue);
    } else {
        throw std::runtime_error("Invalid multigrid cycle value.\n");
    }

    max_iterations_ = parser_.get<int>("maxIterations");

    const int normValue = parser_.get<int>("residualNormType");
    if (normValue == static_cast<int>(ResidualNormType::EUCLIDEAN) ||
        normValue == static_cast<int>(ResidualNormType::WEIGHTED_EUCLIDEAN) ||
        normValue == static_cast<int>(ResidualNormType::INFINITY_NORM)) {
        residual_norm_type_ = static_cast<ResidualNormType>(normValue);
    } else {
        throw std::runtime_error("Invalid residual norm type.\n");
    }

    double absTol = parser_.get<double>("absoluteTolerance");
    if (absTol < 0) {
        absolute_tolerance_ = std::nullopt;
    } else {
        absolute_tolerance_ = absTol;
    }
    double relTol = parser_.get<double>("relativeTolerance");
    if (relTol < 0) {
        relative_tolerance_ = std::nullopt;
    } else {
        relative_tolerance_ = relTol;
    }
}

void GMGPolar::parseGeneral() {
    verbose_ = parser_.get<int>("verbose");
    paraview_ = parser_.get<int>("paraview") != 0;
    max_omp_threads_ = parser_.get<int>("maxOpenMPThreads"); omp_set_num_threads(max_omp_threads_);
}


void GMGPolar::initializeGrid() {
    parser_.add<double>(
        "R0", 'r', "Interior radius of the disk.", 
        OPTIONAL, 1e-5
    );
    parser_.add<double>(
        "Rmax", 'R', "Exterior radius of the disk.", 
        OPTIONAL, 1.3
    );
    parser_.add<int>(
        "nr_exp", 'n', 
        "Number of nodes (exponents) in the radial direction.", 
        OPTIONAL, 5
    );
    parser_.add<int>(
        "ntheta_exp", '\0', 
        "Number of nodes (exponents) in the angular direction.", 
        OPTIONAL, -1
    );
    parser_.add<int>(
        "anisotropic_factor", '\0', 
        "Defines anisotropic discretization in the radial direction.",
        OPTIONAL, 0
    );
    parser_.add<int>(
        "divideBy2", '\0', 
        "Refines the grid globally `divideBy2` times.",
        OPTIONAL, 0
    );
    parser_.add<int>(
        "DirBC_Interior", '\0', "Defines the boundary condition on the interior circle: Across-origin (0), Dirichlet boundary (1).",
        OPTIONAL, 0, cmdline::oneof(0,1)
    );
}

void GMGPolar::initializeMultigrid() {
    parser_.add<int>(
        "FMG", '\0', 
        "Specifies whether the initial approximation is obtained by nested iteration.", 
        OPTIONAL, 0, cmdline::oneof(0,1)
    );
    parser_.add<int>(
        "FMG_iterations", '\0', 
        "Specifies the number of FMG iterations.", 
        OPTIONAL, 2
    );
    parser_.add<int>(
        "FMG_cycle", '\0', 
        "Specifies the type of FMG Cycle: V-cycle (0), W-cycle (1), F-cycle (2).", 
        OPTIONAL, 0, cmdline::oneof(0,1,2)
    );
    parser_.add<int>(
        "extrapolation", 'e', 
        "Specifies the type of extrapolation: No extrapolation (0), Implicit extrapolation (1), Implicit extrapolation with full grid smoothing (2), Combination of both methods (3).",
        OPTIONAL, 0, cmdline::oneof(0,1,2,3)
    );
    parser_.add<int>(
        "maxLevels", 'l', 
        "Defines the maximum number of levels in the multigrid scheme.", 
        OPTIONAL, -1
    );
    parser_.add<int>(
        "gpuLevels", '\0', 
        "Defines the number of GPU multigrid levels.", 
        OPTIONAL, -1
    );
    parser_.add<int>(
        "preSmoothingSteps", '\0', 
        "Number of pre-smoothing steps.", 
        OPTIONAL, 1
    );
    parser_.add<int>(
        "postSmoothingSteps", '\0', 
        "Number of post-smoothing steps.", 
        OPTIONAL, 1
    );
    parser_.add<int>(
        "multigridCycle", '\0', 
        "Type of Multigrid Cycle: V-cycle (0), W-cycle (1), F-cycle (2).", 
        OPTIONAL, 0, cmdline::oneof(0, 1, 2)
    );

    parser_.add<int>(
        "residualNormType", '\0', 
        "Type of Residual Norm: Euclidean (0), Weighted Euclidean (1), Infinity (2).", 
        OPTIONAL, 0, cmdline::oneof(0, 1, 2)
    );
    parser_.add<int>(
        "maxIterations", '\0', 
        "Maximum number of Multigrid iterations.", 
        OPTIONAL, 150
    );
    parser_.add<double>(
        "absoluteTolerance", '\0', 
        "Convergence achieved when absolute tolerance is reached.", 
        OPTIONAL, 1e-8
    );
    parser_.add<double>(
        "relativeTolerance", '\0', 
        "Convergence achieved when relative tolerance is reached.", 
        OPTIONAL, 1e-8
    );
}

void GMGPolar::initializeGeneral() {
    parser_.add<int>(
        "verbose", '\0', 
        "Controls the verbosity of the output. Higher values produce more detailed diagnostic information.", 
        OPTIONAL, 1
    );
    parser_.add<int>(
        "paraview", '\0', 
        "Specifies whether to generate Paraview output files.", 
        OPTIONAL, 0
    );
    parser_.add<int>(
        "maxOpenMPThreads", '\0', 
        "Defines the maximum number of OpenMP threads used.", 
        OPTIONAL, 1
    );
}