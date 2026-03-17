#pragma once

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::extrapolated_multigrid_F_Cycle(int level_depth,
                                                                                          Vector<double> solution,
                                                                                          ConstVector<double> rhs,
                                                                                          Vector<double> residual)
{
    assert(0 <= level_depth && level_depth < number_of_levels_);

    std::chrono::high_resolution_clock::time_point start_MGC;
    if (level_depth == 0) {
        start_MGC = std::chrono::high_resolution_clock::now();
    }

    if (level_depth == number_of_levels_ - 1) {
        /* ---------------------------------------------------- */
        /* Coarsest level: solve A * x = rhs using DirectSolver */
        /* ---------------------------------------------------- */
        // Note: A_ex * x = rhs is only solvable when level_depth >= 2.

        Level<DomainGeometry>& coarsest_level = levels_[level_depth];

        /* Step 1: Copy rhs in solution */
        Kokkos::deep_copy(solution, rhs);

        /* Step 2: Solve for the solution in place */
        auto start_MGC_directSolver = std::chrono::high_resolution_clock::now();

        coarsest_level.directSolveInPlace(solution);

        auto end_MGC_directSolver = std::chrono::high_resolution_clock::now();
        t_avg_MGC_directSolver_ += std::chrono::duration<double>(end_MGC_directSolver - start_MGC_directSolver).count();
    }
    else {
        /* ------------------------------ */
        /* Extrapolated multigrid F-cycle */
        /* ------------------------------ */
        Level<DomainGeometry>& level      = levels_[level_depth];
        Level<DomainGeometry>& next_level = levels_[level_depth + 1];

        /* ------------ */
        /* Presmoothing */
        /* ------------ */
        auto start_MGC_preSmoothing = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < pre_smoothing_steps_; i++) {
            if (level.level_depth() == 0 && !full_grid_smoothing_) {
                level.extrapolatedSmoothing(solution, rhs, residual);
            }
            else {
                level.smoothing(solution, rhs, residual);
            }
        }

        auto end_MGC_preSmoothing = std::chrono::high_resolution_clock::now();
        t_avg_MGC_preSmoothing_ += std::chrono::duration<double>(end_MGC_preSmoothing - start_MGC_preSmoothing).count();

        /* -------------------------------------------- */
        /* Compute the restricted extrapolated residual */
        /* -------------------------------------------- */
        auto start_MGC_residual = std::chrono::high_resolution_clock::now();

        // P_ex^T (f_l - A_l*u_l)
        level.computeResidual(residual, rhs, solution);
        extrapolatedRestriction(level.level_depth(), next_level.residual(), residual);

        // f_{l-1} - A_{l-1}* Inject(u_l)
        injection(level.level_depth(), next_level.solution(), solution);
        next_level.computeResidual(next_level.error_correction(), next_level.rhs(), next_level.solution());

        // res_ex = 4/3 * P_ex^T (f_l - A_l*u_l) - 1/3 * (f_{l-1} - A_{l-1}* Inject(u_l))
        linear_combination(next_level.residual(), 4.0 / 3.0, ConstVector<double>(next_level.error_correction()),
                           -1.0 / 3.0);

        auto end_MGC_residual = std::chrono::high_resolution_clock::now();
        t_avg_MGC_residual_ += std::chrono::duration<double>(end_MGC_residual - start_MGC_residual).count();

        /* -------------------------------------------------- */
        /* Solve A * error = restricted extrapolated residual */
        /* -------------------------------------------------- */
        // Note: We deliberately use the non-extrapolated multigrid cycle here.

        /* Set starting error to zero. */
        assign(next_level.error_correction(), 0.0);

        /* Solve for the error by recursively calling the multigrid cycle. */
        multigrid_F_Cycle(next_level.level_depth(), next_level.error_correction(), next_level.residual(),
                          next_level.solution());

        /* Don't do a second recursion on the coarsest level since the DirectSolver is exact. */
        if (next_level.level_depth() != number_of_levels_ - 1) {
            multigrid_V_Cycle(next_level.level_depth(), next_level.error_correction(), next_level.residual(),
                              next_level.solution());
        }

        /* -------------------------- */
        /* Interpolate the correction */
        /* -------------------------- */
        // Use 'residual' instead of 'level.error_correction()' as a temporary buffer.
        // Note: 'level.error_correction()' has size 0 at level depth = 0.
        extrapolatedProlongation(next_level.level_depth(), residual, next_level.error_correction());

        /* ----------------------------------- */
        /* Compute the corrected approximation */
        /* ----------------------------------- */
        add(solution, ConstVector<double>(residual));

        /* ------------- */
        /* Postsmoothing */
        /* ------------- */
        auto start_MGC_postSmoothing = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < post_smoothing_steps_; i++) {
            if (level.level_depth() == 0 && !full_grid_smoothing_) {
                level.extrapolatedSmoothing(solution, rhs, residual);
            }
            else {
                level.smoothing(solution, rhs, residual);
            }
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
