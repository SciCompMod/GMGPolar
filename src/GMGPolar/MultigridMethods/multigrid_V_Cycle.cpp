#include "../../../include/GMGPolar/gmgpolar.h"

void GMGPolar::multigrid_V_Cycle(const int level_depth, Vector<double>& solution, Vector<double>& rhs, Vector<double>& residual) {
    assert(0 <= level_depth && level_depth < numberOflevels_-1);

    auto start_MGC = std::chrono::high_resolution_clock::now();

    Level& level = levels_[level_depth];
    Level& next_level = levels_[level_depth+1];

    auto start_MGC_preSmoothing = std::chrono::high_resolution_clock::now();

    /* ------------ */
    /* Presmoothing */
    for (int i = 0; i < preSmoothingSteps_; i++){
        level.smoothingInPlace(solution, rhs, residual);
    }

    auto end_MGC_preSmoothing = std::chrono::high_resolution_clock::now();
    t_avg_MGC_preSmoothing += std::chrono::duration<double>(end_MGC_preSmoothing - start_MGC_preSmoothing).count();

    /* ---------------------- */
    /* Coarse grid correction */
    /* ---------------------- */

    auto start_MGC_residual = std::chrono::high_resolution_clock::now();

    /* Compute the residual */
    level.computeResidual(residual, rhs, solution);

    auto end_MGC_residual = std::chrono::high_resolution_clock::now();
    t_avg_MGC_residual += std::chrono::duration<double>(end_MGC_residual - start_MGC_residual).count();

    /* -------------------------- */
    /* Solve A * error = residual */
    if(level_depth+1 == numberOflevels_-1){
        /* --------------------- */
        /* Using a direct solver */
        /* --------------------- */

        /* Step 1: Restrict the residual */
        restrictToLowerLevel(level_depth, next_level.residual(), residual);

        /* Step 2: Solve for the error in place */
        auto start_MGC_directSolver = std::chrono::high_resolution_clock::now(); 

        next_level.directSolveInPlace(next_level.residual());

        auto end_MGC_directSolver = std::chrono::high_resolution_clock::now();
        t_avg_MGC_directSolver += std::chrono::duration<double>(end_MGC_directSolver - start_MGC_directSolver).count();
    } else{
        /* ------------------------------------------ */
        /* By recursively calling the multigrid cycle */
        /* ------------------------------------------ */

        /* Step 1: Restrict the residual. */
        restrictToLowerLevel(level_depth, next_level.rhs_error(), residual);

        /* Step 2: Set starting error to zero. */
        assign(next_level.residual(), 0.0);

        /* Step 3: Solve for the error by recursively calling the multigrid cycle. */
        multigrid_V_Cycle(level_depth+1, next_level.residual(), next_level.rhs_error(), next_level.solution());
    }

    /* Interpolate the correction */
    prolongateToUpperLevel(level_depth+1, residual, next_level.residual());

    /* Compute the corrected approximation: u = u + error */
    add(solution, residual);

    auto start_MGC_postSmoothing = std::chrono::high_resolution_clock::now();

    /* ------------- */
    /* Postsmoothing */
    for (int i = 0; i < postSmoothingSteps_; i++){
        level.smoothingInPlace(solution, rhs, residual);
    }

    auto end_MGC_postSmoothing = std::chrono::high_resolution_clock::now();
    t_avg_MGC_postSmoothing += std::chrono::duration<double>(end_MGC_postSmoothing - start_MGC_postSmoothing).count();

    auto end_MGC = std::chrono::high_resolution_clock::now();
    t_avg_MGC_total += std::chrono::duration<double>(end_MGC - start_MGC).count();
}