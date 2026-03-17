#pragma once

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::multigrid_F_Cycle(int level_depth, Vector<double> solution,
                                                                             ConstVector<double> rhs,
                                                                             Vector<double> residual)
{
    assert(0 <= level_depth && level_depth < number_of_levels_);

    std::chrono::high_resolution_clock::time_point start_MGC;
    if (level_depth == 0) {
        start_MGC = std::chrono::high_resolution_clock::now();
    }

    /* ------------------------ */
    /* Solve A * solution = rhs */
    /* ------------------------ */
    if (level_depth == number_of_levels_ - 1) {
        /* ------------------------------ */
        /* Coarsest level: solve directly */
        /* ------------------------------ */
        Level<DomainGeometry>& level = levels_[level_depth];

        /* Step 1: Copy rhs in solution */
        Kokkos::deep_copy(solution, rhs);

        /* Step 2: Solve for the solution in place */
        auto start_MGC_directSolver = std::chrono::high_resolution_clock::now();

        level.directSolveInPlace(solution);

        auto end_MGC_directSolver = std::chrono::high_resolution_clock::now();
        t_avg_MGC_directSolver_ += std::chrono::duration<double>(end_MGC_directSolver - start_MGC_directSolver).count();
    }
    else {
        /* ------------------------------------------------------- */
        /* Multigrid V-cycle on current level:                     */
        /* presmoothing, coarse-grid correction, and postsmoothing */
        /* ------------------------------------------------------- */
        Level<DomainGeometry>& level      = levels_[level_depth];
        Level<DomainGeometry>& next_level = levels_[level_depth + 1];

        auto start_MGC_preSmoothing = std::chrono::high_resolution_clock::now();

        /* ------------ */
        /* Presmoothing */
        for (int i = 0; i < pre_smoothing_steps_; i++) {
            level.smoothing(solution, rhs, residual);
        }

        auto end_MGC_preSmoothing = std::chrono::high_resolution_clock::now();
        t_avg_MGC_preSmoothing_ += std::chrono::duration<double>(end_MGC_preSmoothing - start_MGC_preSmoothing).count();

        /* ----------------------------- */
        /* Compute fine-grid residual    */
        /* residual = rhs - A * solution */
        /* ----------------------------- */
        auto start_MGC_residual = std::chrono::high_resolution_clock::now();

        level.computeResidual(residual, rhs, solution);

        auto end_MGC_residual = std::chrono::high_resolution_clock::now();
        t_avg_MGC_residual_ += std::chrono::duration<double>(end_MGC_residual - start_MGC_residual).count();

        /* -------------------------- */
        /* Solve A * error = residual */
        /* -------------------------- */
        /* By recursively calling the multigrid cycle */

        /* Step 1: Restrict the residual. */
        restriction(level_depth, next_level.error_correction(), residual);

        /* Step 2: Set starting error to zero. */
        assign(next_level.residual(), 0.0);

        /* Step 3: Solve for the error by recursively calling the multigrid cycle. */
        multigrid_F_Cycle(level_depth + 1,
                          next_level.residual(), // error (solution)
                          next_level.error_correction(), // coarse residual (rhs)
                          next_level.solution()); // workspace

        // Don't do a second recursion on the coarsest level since the DirectSolver is exact.
        if (level_depth + 1 != number_of_levels - 1) {
            multigrid_V_Cycle(level_depth + 1,
                              next_level.residual(), // error (solution)
                              next_level.error_correction(), // coarse residual (rhs)
                              next_level.solution()); // workspace
        }

        /* Interpolate the correction */
        prolongation(level_depth + 1, residual, next_level.residual());

        /* Compute the corrected approximation: u = u + error */
        add(solution, ConstVector<double>(residual));

        auto start_MGC_postSmoothing = std::chrono::high_resolution_clock::now();

        /* ------------- */
        /* Postsmoothing */
        for (int i = 0; i < post_smoothing_steps_; i++) {
            level.smoothing(solution, rhs, residual);
        }

        auto end_MGC_postSmoothing = std::chrono::high_resolution_clock::now();
        t_avg_MGC_postSmoothing_ +=
            std::chrono::duration<double>(end_MGC_postSmoothing - start_MGC_postSmoothing).count();
    }

    if (level_depth == 0) {
        auto end_MGC = std::chrono::high_resolution_clock::now();
        t_avg_MGC_total_ += std::chrono::duration<double>(end_MGC - start_MGC).count();
    }
}
