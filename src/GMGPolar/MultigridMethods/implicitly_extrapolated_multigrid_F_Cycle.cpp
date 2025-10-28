#include "../../../include/GMGPolar/gmgpolar.h"

void GMGPolar::implicitlyExtrapolatedMultigrid_F_Cycle(const int level_depth, Vector<double> solution,
                                                       Vector<double> const rhs, Vector<double> residual)
{
    assert(0 <= level_depth && level_depth < number_of_levels_ - 1);

    std::chrono::high_resolution_clock::time_point start_MGC;
    if (level_depth == 0) {
        start_MGC = std::chrono::high_resolution_clock::now();
    }

    Level& level      = levels_[level_depth];
    Level& next_level = levels_[level_depth + 1];

    auto start_MGC_preSmoothing = std::chrono::high_resolution_clock::now();

    /* ------------ */
    /* Presmoothing */
    for (int i = 0; i < pre_smoothing_steps_; i++) {
        if (level_depth == 0 && !full_grid_smoothing_) {
            level.extrapolatedSmoothing(solution, rhs, residual);
        }
        else {
            level.smoothing(solution, rhs, residual);
        }
    }

    auto end_MGC_preSmoothing = std::chrono::high_resolution_clock::now();
    t_avg_MGC_preSmoothing_ += std::chrono::duration<double>(end_MGC_preSmoothing - start_MGC_preSmoothing).count();

    /* ---------------------- */
    /* Coarse grid correction */
    /* ---------------------- */

    /* -------------------------- */
    /* Solve A * error = residual */
    if (level_depth + 1 == number_of_levels_ - 1) {
        /* --------------------- */
        /* Using a direct solver */
        /* --------------------- */

        /* Step 1: Compute extrapolated residual */
        auto start_MGC_residual = std::chrono::high_resolution_clock::now();

        // P_ex^T (f_l - A_l*u_l)
        level.computeResidual(residual, rhs, solution);
        extrapolatedRestriction(level_depth, next_level.residual(), residual);

        // f_{l-1} - A_{l-1}* Inject(u_l)
        injection(level_depth, next_level.solution(), solution);
        next_level.computeResidual(next_level.error_correction(), next_level.rhs(), next_level.solution());

        // res_ex = 4/3 * P_ex^T (f_l - A_l*u_l) - 1/3 * (f_{l-1} - A_{l-1}* Inject(u_l))
        linear_combination(next_level.residual(), 4.0 / 3.0, ConstVector<double>(next_level.error_correction()),
                           -1.0 / 3.0);

        auto end_MGC_residual = std::chrono::high_resolution_clock::now();
        t_avg_MGC_residual_ += std::chrono::duration<double>(end_MGC_residual - start_MGC_residual).count();

        /* Step 2: Solve for the error in place */
        auto start_MGC_directSolver = std::chrono::high_resolution_clock::now();

        next_level.directSolveInPlace(next_level.residual());

        auto end_MGC_directSolver = std::chrono::high_resolution_clock::now();
        t_avg_MGC_directSolver_ += std::chrono::duration<double>(end_MGC_directSolver - start_MGC_directSolver).count();
    }
    else {
        /* ------------------------------------------ */
        /* By recursively calling the multigrid cycle */
        /* ------------------------------------------ */

        /* Step 1: Compute extrapolated residual. */
        auto start_MGC_residual = std::chrono::high_resolution_clock::now();

        // P_ex^T (f_l - A_l*u_l)
        level.computeResidual(residual, rhs, solution);
        extrapolatedRestriction(level_depth, next_level.error_correction(), residual);

        // f_{l-1} - A_{l-1}* Inject(u_l)
        injection(level_depth, next_level.solution(), solution);
        next_level.computeResidual(next_level.residual(), next_level.rhs(), next_level.solution());

        // res_ex = 4/3 * P_ex^T (f_l - A_l*u_l) - 1/3 * (f_{l-1} - A_{l-1}* Inject(u_l))
        linear_combination(next_level.error_correction(), 4.0 / 3.0, ConstVector<double>(next_level.residual()),
                           -1.0 / 3.0);

        auto end_MGC_residual = std::chrono::high_resolution_clock::now();
        t_avg_MGC_residual_ += std::chrono::duration<double>(end_MGC_residual - start_MGC_residual).count();

        /* Step 2: Set starting error to zero. */
        assign(next_level.residual(), 0.0);

        /* Step 3: Solve for the error by recursively calling the multigrid cycle. */
        multigrid_F_Cycle(level_depth + 1, next_level.residual(), next_level.error_correction(), next_level.solution());
        multigrid_V_Cycle(level_depth + 1, next_level.residual(), next_level.error_correction(), next_level.solution());
    }

    /* Interpolate the correction */
    extrapolatedProlongation(level_depth + 1, residual, next_level.residual());

    /* Compute the corrected approximation: u = u + error */
    add(solution, ConstVector<double>(residual));

    auto start_MGC_postSmoothing = std::chrono::high_resolution_clock::now();

    /* ------------- */
    /* Postsmoothing */
    for (int i = 0; i < post_smoothing_steps_; i++) {
        if (level_depth == 0 && !full_grid_smoothing_) {
            level.extrapolatedSmoothing(solution, rhs, residual);
        }
        else {
            level.smoothing(solution, rhs, residual);
        }
    }

    auto end_MGC_postSmoothing = std::chrono::high_resolution_clock::now();
    t_avg_MGC_postSmoothing_ += std::chrono::duration<double>(end_MGC_postSmoothing - start_MGC_postSmoothing).count();

    if (level_depth == 0) {
        auto end_MGC = std::chrono::high_resolution_clock::now();
        t_avg_MGC_total_ += std::chrono::duration<double>(end_MGC - start_MGC).count();
    }
}