#pragma once

namespace residual_give
{

template <class LevelCacheType>
static inline void node_apply_a_give(int i_r, int i_theta, const PolarGrid& grid, const LevelCacheType& level_cache,
                                     bool DirBC_Interior, Vector<double>& result, ConstVector<double>& x)
{
    /* ---------------------------------------- */
    /* Compute or retrieve stencil coefficients */
    /* ---------------------------------------- */
    const int center   = grid.index(i_r, i_theta);
    const double r     = grid.radius(i_r);
    const double theta = grid.theta(i_theta);

    double coeff_beta, arr, att, art, detDF;
    level_cache.obtainValues(i_r, i_theta, center, r, theta, coeff_beta, arr, att, art, detDF);

    /* -------------------- */
    /* Node in the interior */
    /* -------------------- */
    if (i_r > 0 && i_r < grid.nr() - 1) {
        const double h1 = grid.radialSpacing(i_r - 1);
        const double h2 = grid.radialSpacing(i_r);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff1 = 0.5 * (k1 + k2) / h1;
        const double coeff2 = 0.5 * (k1 + k2) / h2;
        const double coeff3 = 0.5 * (h1 + h2) / k1;
        const double coeff4 = 0.5 * (h1 + h2) / k2;
        const double coeff5 = 0.25 * (h1 + h2) * (k1 + k2);

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int left   = grid.index(i_r - 1, i_theta);
        const int right  = grid.index(i_r + 1, i_theta);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int top    = grid.index(i_r, i_theta_P1);

        /* Fill result(i,j) */
        result[center] += (coeff5 * coeff_beta * std::fabs(detDF) * x[center] /* beta_{i,j} */
                           - coeff1 * arr * x[left] /* Left */
                           - coeff2 * arr * x[right] /* Right */
                           - coeff3 * att * x[bottom] /* Bottom */
                           - coeff4 * att * x[top] /* Top */
                           + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) *
                                 x[center]); /* Center: (Left, Right, Bottom, Top) */

        /* Fill result(i-1,j) */
        if (i_r > 1 || (i_r == 1 && !DirBC_Interior)) { /* Don't give to the inner dirichlet boundary! */
            result[left] += (-coeff1 * arr * x[center] /* Right */
                             + coeff1 * arr * x[left] /* Center: (Right) */
                             - 0.25 * art * x[top] /* Top Right */
                             + 0.25 * art * x[bottom]); /* Bottom Right */
        }
        /* Fill result(i+1,j) */
        if (i_r < grid.nr() - 2) { /* Don't give to the outer dirichlet boundary! */
            result[right] += (-coeff2 * arr * x[center] /* Left */
                              + coeff2 * arr * x[right] /* Center: (Left) */
                              + 0.25 * art * x[top] /* Top Left */
                              - 0.25 * art * x[bottom]); /* Bottom Left */
        }
        /* Fill result(i,j-1) */
        result[bottom] += (-coeff3 * att * x[center] /* Top */
                           + coeff3 * att * x[bottom] /* Center: (Top) */
                           - 0.25 * art * x[right] /* Top Right */
                           + 0.25 * art * x[left]); /* Top Left */
        /* Fill result(i,j+1) */
        result[top] += (-coeff4 * att * x[center] /* Bottom */
                        + coeff4 * att * x[top] /* Center: (Bottom) */
                        + 0.25 * art * x[right] /* Bottom Right */
                        - 0.25 * art * x[left]); /* Bottom Left */
    }
    /* -------------------------- */
    /* Node on the inner boundary */
    /* -------------------------- */
    else if (i_r == 0) {
        /* ------------------------------------------------ */
        /* Case 1: Dirichlet boundary on the inner boundary */
        /* ------------------------------------------------ */
        if (DirBC_Interior) {
            /* Fill result(i,j) */
            result[center] += x[center];

            /* Give value to the interior nodes! */
            double h2 = grid.radialSpacing(i_r);
            double k1 = grid.angularSpacing(i_theta - 1);
            double k2 = grid.angularSpacing(i_theta);

            double coeff2 = 0.5 * (k1 + k2) / h2;

            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

            const int right  = grid.index(i_r + 1, i_theta);
            const int bottom = grid.index(i_r, i_theta_M1);
            const int top    = grid.index(i_r, i_theta_P1);

            /* Fill result(i+1,j) */
            result[right] += (-coeff2 * arr * x[center] /* Left */
                              + coeff2 * arr * x[right] /* Center: (Left) */
                              + 0.25 * art * x[top] /* Top Left */
                              - 0.25 * art * x[bottom]); /* Bottom Left */
        }
        else {
            /* ------------------------------------------------------------- */
            /* Case 2: Across origin discretization on the interior boundary */
            /* ------------------------------------------------------------- */
            /* h1 gets replaced with 2 * R0. */
            /* (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()/2)). */
            /* Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil. */
            const double h1 = 2.0 * grid.radius(0);
            const double h2 = grid.radialSpacing(i_r);
            const double k1 = grid.angularSpacing(i_theta - 1);
            const double k2 = grid.angularSpacing(i_theta);

            const double coeff1 = 0.5 * (k1 + k2) / h1;
            const double coeff2 = 0.5 * (k1 + k2) / h2;
            const double coeff3 = 0.5 * (h1 + h2) / k1;
            const double coeff4 = 0.5 * (h1 + h2) / k2;
            const double coeff5 = 0.25 * (h1 + h2) * (k1 + k2);

            const int i_theta_M1     = grid.wrapThetaIndex(i_theta - 1);
            const int i_theta_P1     = grid.wrapThetaIndex(i_theta + 1);
            const int i_theta_Across = grid.wrapThetaIndex(i_theta + grid.ntheta() / 2);

            const int left   = grid.index(i_r, i_theta_Across);
            const int right  = grid.index(i_r + 1, i_theta);
            const int bottom = grid.index(i_r, i_theta_M1);
            const int top    = grid.index(i_r, i_theta_P1);

            /* Fill result(i,j) */
            result[center] += (coeff5 * coeff_beta * std::fabs(detDF) * x[center] /* beta_{i,j} */
                               - coeff1 * arr * x[left] /* Left */
                               - coeff2 * arr * x[right] /* Right */
                               - coeff3 * att * x[bottom] /* Bottom */
                               - coeff4 * att * x[top] /* Top */
                               + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) *
                                     x[center]); /* Center: (Left, Right, Bottom, Top) */
            /* Fill result(i-1,j) */
            /* From view the view of the across origin node, the directions are roatated by 180 degrees in the stencil! */
            result[left] += (-coeff1 * arr * x[center] /* Right -> Left */
                             + coeff1 * arr * x[left]); /* Center: (Right) -> Center: (Left)*/
            /* + 0.25 * art * x[top]; // Top Right -> Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
            /* - 0.25 * art * x[bottom]; // Bottom Right -> Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
            /* Fill result(i+1,j) */
            result[right] += (-coeff2 * arr * x[center] /* Left */
                              + coeff2 * arr * x[right] /* Center: (Left) */
                              + 0.25 * art * x[top] /* Top Left */
                              - 0.25 * art * x[bottom]); /* Bottom Left */
            /* Fill result(i,j-1) */
            result[bottom] += (-coeff3 * att * x[center] /* Top */
                               + coeff3 * att * x[bottom] /* Center: (Top) */
                               - 0.25 * art * x[right]); /* Top Right */
            /* + 0.25 * art * x[left]; // Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
            /* Fill result(i,j+1) */
            result[top] += (-coeff4 * att * x[center] /* Bottom */
                            + coeff4 * att * x[top] /* Center: (Bottom) */
                            + 0.25 * art * x[right]); /* Bottom Right */
            /* - 0.25 * art * x[left]; // Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
        }
    }
    /* ----------------------------- */
    /* Node on to the outer boundary */
    /* ----------------------------- */
    else if (i_r == grid.nr() - 1) {
        /* Fill result of (i,j) */
        result[center] += x[center];

        /* Give value to the interior nodes! */
        double h1 = grid.radialSpacing(i_r - 1);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int left   = grid.index(i_r - 1, i_theta);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int top    = grid.index(i_r, i_theta_P1);

        /* Fill result(i-1,j) */
        result[left] += (-coeff1 * arr * x[center] /* Right */
                         + coeff1 * arr * x[left] /* Center: (Right) */
                         - 0.25 * art * x[top] /* Top Right */
                         + 0.25 * art * x[bottom]); /* Bottom Right */
    }
}

} // namespace residual_give

template <class LevelCacheType>
void ResidualGive<LevelCacheType>::applySystemOperator(Vector<double> result, ConstVector<double> x) const
{
    assert(result.size() == x.size());

    const PolarGrid& grid             = Residual<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = Residual<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = Residual<LevelCacheType>::DirBC_Interior_;

    using residual_give::node_apply_a_give;

    assign(result, 0.0);

    /* ---------------- */
    /* Circular section */
    /* ---------------- */
    // We parallelize over i_r (step 3) to avoid data race conditions between adjacent circles.
    // The i_theta loop is sequential inside the kernel.
    const int num_circle_tasks = grid.numberSmootherCircles();

    for (int start_circle = 0; start_circle < 3; ++start_circle) {
        const int num_circular_tasks = (num_circle_tasks - start_circle + 2) / 3;
        Kokkos::parallel_for(
            "ResidualGive: ApplyA (Circular)", Kokkos::RangePolicy<>(0, num_circular_tasks),
            KOKKOS_LAMBDA(const int circle_task) {
                const int i_r = start_circle + circle_task * 3;
                for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                    node_apply_a_give(i_r, i_theta, grid, level_cache, DirBC_Interior, result, x);
                }
            });
        Kokkos::fence();
    }

    /* -------------- */
    /* Radial section */
    /* -------------- */
    // We parallelize over i_theta (step 3) to avoid data race conditions between adjacent radial lines.
    // The i_r loop is sequential inside the kernel.
    // Due to periodicity in the angular direction, handle up to 2 additional
    // radial lines (i_theta = 0 and 1) before the parallel passes.
    const int additional_radial_tasks = grid.ntheta() % 3;
    const int num_radial_tasks        = grid.ntheta() - additional_radial_tasks;

    for (int i_theta = 0; i_theta < additional_radial_tasks; i_theta++) {
        Kokkos::parallel_for(
            "ResidualGive: ApplyA (Radial, additional)", Kokkos::RangePolicy<>(0, 1), KOKKOS_LAMBDA(const int) {
                for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                    node_apply_a_give(i_r, i_theta, grid, level_cache, DirBC_Interior, result, x);
                }
            });
        Kokkos::fence();
    }

    for (int start_radial = 0; start_radial < 3; ++start_radial) {
        const int num_radial_batches = (num_radial_tasks - start_radial + 2) / 3;
        Kokkos::parallel_for(
            "ResidualGive: ApplyA (Radial)", Kokkos::RangePolicy<>(0, num_radial_batches),
            KOKKOS_LAMBDA(const int radial_task) {
                const int i_theta = start_radial + radial_task * 3;
                for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                    node_apply_a_give(i_r, i_theta, grid, level_cache, DirBC_Interior, result, x);
                }
            });
        Kokkos::fence();
    }
}
// clang-format on