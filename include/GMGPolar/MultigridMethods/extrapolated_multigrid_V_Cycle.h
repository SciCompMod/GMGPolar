#pragma once

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::extrapolated_multigrid_V_Cycle(int level_depth,
                                                                                          HostVector<double> h_solution,
                                                                                          HostConstVector<double> h_rhs,
                                                                                          HostVector<double> h_residual)
{
    assert(level_depth == 0);
    assert(number_of_levels_ >= 2);

    std::chrono::high_resolution_clock::time_point start_MGC;
    if (level_depth == 0) {
        start_MGC = std::chrono::high_resolution_clock::now();
    }

    /* ------------------------------ */
    /* Extrapolated multigrid V-cycle */
    /* ------------------------------ */
    Level<DomainGeometry, DensityProfileCoefficients>& level      = levels_[level_depth];
    Level<DomainGeometry, DensityProfileCoefficients>& next_level = levels_[level_depth + 1];

    /* ------------ */
    /* Presmoothing */
    /* ------------ */
    auto start_MGC_preSmoothing = std::chrono::high_resolution_clock::now();

    auto residual = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_residual);
    auto solution = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_solution);
    auto rhs      = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_rhs); // const

    for (int i = 0; i < pre_smoothing_steps_; i++) {
        if (level.level_depth() == 0 && !full_grid_smoothing_) {
            level.extrapolatedSmoothing(solution, rhs, residual);
        }
        else {
            level.smoothing(solution, rhs, residual);
        }
    }
	Kokkos::deep_copy(h_residual, residual);
	Kokkos::deep_copy(h_solution, solution);

    auto end_MGC_preSmoothing = std::chrono::high_resolution_clock::now();
    t_avg_MGC_preSmoothing_ += std::chrono::duration<double>(end_MGC_preSmoothing - start_MGC_preSmoothing).count();

    /* -------------------------------------------- */
    /* Compute the restricted extrapolated h_residual */
    /* -------------------------------------------- */
    auto start_MGC_residual = std::chrono::high_resolution_clock::now();

    // P_ex^T (f_l - A_l*u_l)
    level.computeResidual(h_residual, h_rhs, h_solution);
    auto next_level_residual = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), next_level.residual());
    Kokkos::deep_copy(residual, h_residual);
    extrapolatedRestriction(level.level_depth(), next_level_residual, residual);
    Kokkos::deep_copy(next_level.residual(), next_level_residual);

    // f_{l-1} - A_{l-1}* Inject(u_l)
    Kokkos::deep_copy(solution, h_solution);
    auto next_level_solution = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), next_level.solution());
    injection(level.level_depth(), next_level_solution, solution);
    Kokkos::deep_copy(next_level.solution(), next_level_solution);
    next_level.computeResidual(next_level.error_correction(), next_level.rhs(), next_level.solution());

    // res_ex = 4/3 * P_ex^T (f_l - A_l*u_l) - 1/3 * (f_{l-1} - A_{l-1}* Inject(u_l))
    linear_combination(next_level.residual(), 4.0 / 3.0, HostConstVector<double>(next_level.error_correction()),
                       -1.0 / 3.0);

    auto end_MGC_residual = std::chrono::high_resolution_clock::now();
    t_avg_MGC_residual_ += std::chrono::duration<double>(end_MGC_residual - start_MGC_residual).count();

    /* -------------------------------------------------- */
    /* Solve A * error = restricted extrapolated h_residual */
    /* -------------------------------------------------- */
    // Note: We deliberately use the non-extrapolated multigrid cycle here.

    /* Set starting error to zero. */
    assign(next_level.error_correction(), 0.0);

    /* Solve for the error by recursively calling the multigrid cycle. */
    multigrid_V_Cycle(next_level.level_depth(), next_level.error_correction(), next_level.residual(),
                      next_level.solution());

    /* -------------------------- */
    /* Interpolate the correction */
    /* -------------------------- */
    // Use 'h_residual' instead of 'level.error_correction()' as a temporary buffer.
    // Note: 'level.error_correction()' has size 0 at level depth = 0.
    extrapolatedProlongation(next_level.level_depth(), h_residual, next_level.error_correction());

    /* ----------------------------------- */
    /* Compute the corrected approximation */
    /* ----------------------------------- */
    add(h_solution, HostConstVector<double>(h_residual));

    /* ------------- */
    /* Postsmoothing */
    /* ------------- */
    auto start_MGC_postSmoothing = std::chrono::high_resolution_clock::now();

	Kokkos::deep_copy(residual, h_residual);
	Kokkos::deep_copy(solution, h_solution);
    for (int i = 0; i < post_smoothing_steps_; i++) {
        if (level.level_depth() == 0 && !full_grid_smoothing_) {
            level.extrapolatedSmoothing(solution, rhs, residual);
        }
        else {
            level.smoothing(solution, rhs, residual);
        }
    }
	Kokkos::deep_copy(h_residual, residual);
	Kokkos::deep_copy(h_solution, solution);

    auto end_MGC_postSmoothing = std::chrono::high_resolution_clock::now();
    t_avg_MGC_postSmoothing_ += std::chrono::duration<double>(end_MGC_postSmoothing - start_MGC_postSmoothing).count();

    if (level_depth == 0) {
        auto end_MGC = std::chrono::high_resolution_clock::now();
        t_avg_MGC_total_ += std::chrono::duration<double>(end_MGC - start_MGC).count();
    }
}
