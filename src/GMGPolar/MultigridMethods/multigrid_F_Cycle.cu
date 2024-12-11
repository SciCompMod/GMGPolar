#include "../../../include/GMGPolar/gmgpolar.h"

#include "../../../include/LinearAlgebra/Vector/gpu_vector_operations.h"

void GMGPolar::multigrid_F_Cycle(const int level_depth, Vector<double>& solution, Vector<double>& rhs,
                                 Vector<double>& residual)
{
    assert(0 <= level_depth && level_depth < number_of_levels_ - 1);

    auto start_MGC = std::chrono::high_resolution_clock::now();

    Level& level      = levels_[level_depth];
    Level& next_level = levels_[level_depth + 1];

    assert(level.processingType() != ProcessingType::GPU);
    assert(next_level.processingType() != ProcessingType::GPU);

    auto start_MGC_preSmoothing = std::chrono::high_resolution_clock::now();

    /* ------------ */
    /* Presmoothing */
    for (int i = 0; i < pre_smoothing_steps_; i++) {
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
    if (level_depth + 1 == number_of_levels_ - 1) {
        /* --------------------- */
        /* Using a direct solver */
        /* --------------------- */

        /* Step 1: Restrict the residual */
        restriction(level_depth, next_level.residual(), residual);

        /* Step 2: Solve for the error in place */
        auto start_MGC_directSolver = std::chrono::high_resolution_clock::now();

        next_level.directSolveInPlace(next_level.residual());

        auto end_MGC_directSolver = std::chrono::high_resolution_clock::now();
        t_avg_MGC_directSolver += std::chrono::duration<double>(end_MGC_directSolver - start_MGC_directSolver).count();
    }
    else {
        /* ------------------------------------------ */
        /* By recursively calling the multigrid cycle */
        /* ------------------------------------------ */

        /* Step 1: Restrict the residual. */
        restriction(level_depth, next_level.error_correction(), residual);

        /* Step 2: Set starting error to zero. */
        assign(next_level.residual(), 0.0);

        /* Step 3: Solve for the error by recursively calling the multigrid cycle. */
        multigrid_F_Cycle(level_depth + 1, next_level.residual(), next_level.error_correction(), next_level.solution());
        multigrid_V_Cycle(level_depth + 1, next_level.residual(), next_level.error_correction(), next_level.solution());
    }

    /* Interpolate the correction */
    prolongation(level_depth + 1, residual, next_level.residual());

    /* Compute the corrected approximation: u = u + error */
    add(solution, residual);

    auto start_MGC_postSmoothing = std::chrono::high_resolution_clock::now();

    /* ------------- */
    /* Postsmoothing */
    for (int i = 0; i < post_smoothing_steps_; i++) {
        level.smoothingInPlace(solution, rhs, residual);
    }

    auto end_MGC_postSmoothing = std::chrono::high_resolution_clock::now();
    t_avg_MGC_postSmoothing += std::chrono::duration<double>(end_MGC_postSmoothing - start_MGC_postSmoothing).count();

    auto end_MGC = std::chrono::high_resolution_clock::now();
    t_avg_MGC_total += std::chrono::duration<double>(end_MGC - start_MGC).count();
}



void GMGPolar::multigrid_F_Cycle(const int level_depth, GPU_Vector<double>& solution, GPU_Vector<double>& rhs,
                                 GPU_Vector<double>& residual)
{
    assert(0 <= level_depth && level_depth < number_of_levels_ - 1);

    auto start_MGC = std::chrono::high_resolution_clock::now();

    Level& level      = levels_[level_depth];
    Level& next_level = levels_[level_depth + 1];

    assert(level.processingType() == ProcessingType::GPU);
    assert(next_level.processingType() != ProcessingType::CPU);

    auto start_MGC_preSmoothing = std::chrono::high_resolution_clock::now();

    /* ------------ */
    /* Presmoothing */
    for (int i = 0; i < pre_smoothing_steps_; i++) {
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
    if (level_depth + 1 == number_of_levels_ - 1) {
        /* --------------------- */
        /* Using a direct solver */
        /* --------------------- */

        /* Step 1: Restrict the residual */
        restriction(level_depth, next_level.GPU_residual(), residual);
        if(next_level.processingType() == ProcessingType::CPU_HYBRID){
            copyDeviceToHost(next_level.GPU_residual(), next_level.residual());
        } 

        /* Step 2: Solve for the error in place */
        auto start_MGC_directSolver = std::chrono::high_resolution_clock::now();

        next_level.directSolveInPlace(next_level.residual());

        auto end_MGC_directSolver = std::chrono::high_resolution_clock::now();
        t_avg_MGC_directSolver += std::chrono::duration<double>(end_MGC_directSolver - start_MGC_directSolver).count();
    }
    else {
        /* ------------------------------------------ */
        /* By recursively calling the multigrid cycle */
        /* ------------------------------------------ */

        /* Step 1: Restrict the residual. */
        restriction(level_depth, next_level.GPU_error_correction(), residual);
        if(next_level.processingType() == ProcessingType::CPU_HYBRID){
            copyDeviceToHost(next_level.GPU_error_correction(), next_level.error_correction());
        } 

        /* Step 2: Set starting error to zero. */
        if(next_level.processingType() == ProcessingType::GPU){
            assign(next_level.GPU_residual(), 0.0);
        } else {
            assign(next_level.residual(), 0.0);
        }

        /* Step 3: Solve for the error by recursively calling the multigrid cycle. */
        if(next_level.processingType() == ProcessingType::GPU){
            multigrid_F_Cycle(level_depth + 1, next_level.GPU_residual(), next_level.GPU_error_correction(), next_level.GPU_solution());
            multigrid_V_Cycle(level_depth + 1, next_level.GPU_residual(), next_level.GPU_error_correction(), next_level.GPU_solution());
        } else {
            multigrid_F_Cycle(level_depth + 1, next_level.residual(), next_level.error_correction(), next_level.solution());
            multigrid_V_Cycle(level_depth + 1, next_level.residual(), next_level.error_correction(), next_level.solution());
        }
    }

    /* Interpolate the correction */
    if(next_level.processingType() == ProcessingType::CPU_HYBRID){
        copyHostToDevice(next_level.residual(), next_level.GPU_residual());
    }
    prolongation(level_depth + 1, residual, next_level.GPU_residual());
  
    /* Compute the corrected approximation: u = u + error */
    add(solution, residual);

    auto start_MGC_postSmoothing = std::chrono::high_resolution_clock::now();

    /* ------------- */
    /* Postsmoothing */
    for (int i = 0; i < post_smoothing_steps_; i++) {
        level.smoothingInPlace(solution, rhs, residual);
    }

    auto end_MGC_postSmoothing = std::chrono::high_resolution_clock::now();
    t_avg_MGC_postSmoothing += std::chrono::duration<double>(end_MGC_postSmoothing - start_MGC_postSmoothing).count();

    auto end_MGC = std::chrono::high_resolution_clock::now();
    t_avg_MGC_total += std::chrono::duration<double>(end_MGC - start_MGC).count();
}