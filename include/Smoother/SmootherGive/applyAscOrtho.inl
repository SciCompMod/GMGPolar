#pragma once

namespace smoother_give
{

template <class LevelCacheType>
static KOKKOS_INLINE_FUNCTION void nodeApplyAscOrthoCircleGiveInside(int i_r, int i_theta, const PolarGrid& grid,
                                                                     const LevelCacheType& level_cache,
                                                                     bool DirBC_Interior, ConstVector<double>& x,
                                                                     ConstVector<double>& rhs, Vector<double>& result)
{
    assert(i_r >= 0 && i_r < grid.numberSmootherCircles());

    /* ---------------------------------------- */
    /* Compute or retrieve stencil coefficients */
    /* ---------------------------------------- */
    const int center    = grid.index(i_r, i_theta);
    const double radius = grid.radius(i_r);
    const double theta  = grid.theta(i_theta);

    double coeff_beta, arr, att, art, detDF;
    level_cache.obtainValues(i_r, i_theta, center, radius, theta, coeff_beta, arr, att, art, detDF);

    /* -------------------- */
    /* Node in the interior */
    /* -------------------- */
    if (i_r > 0 && i_r < grid.numberSmootherCircles()) {
        const double h1 = grid.radialSpacing(i_r - 1);
        const double h2 = grid.radialSpacing(i_r);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff1 = 0.5 * (k1 + k2) / h1;
        const double coeff2 = 0.5 * (k1 + k2) / h2;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int center = grid.index(i_r, i_theta);
        const int left   = grid.index(i_r - 1, i_theta);
        const int right  = grid.index(i_r + 1, i_theta);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int top    = grid.index(i_r, i_theta_P1);

        /* Fill result(i,j) */
        result[center] -= (-coeff1 * arr * x[left] /* Left */
                           - coeff2 * arr * x[right]); /* Right */
        /* Fill result(i,j-1) */
        result[bottom] -= (-0.25 * art * x[right] /* Top Right */
                           + 0.25 * art * x[left]); /* Top Left */
        /* Fill result(i,j+1) */
        result[top] -= (+0.25 * art * x[right] /* Bottom Right */
                        - 0.25 * art * x[left]); /* Bottom Left */
    }
    /* -------------------- */
    /* Node on the boundary */
    /* -------------------- */
    else if (i_r == 0) {
        /* ------------------------------------------------ */
        /* Case 1: Dirichlet boundary on the inner boundary */
        /* ------------------------------------------------ */
        if (DirBC_Interior) {
            /* Nothing to be done here */
        }
        else {
            /* ------------------------------------------------------------- */
            /* Case 2: Across origin discretization on the interior boundary */
            /* ------------------------------------------------------------- */
            const double h2 = grid.radialSpacing(i_r);
            const double k1 = grid.angularSpacing(i_theta - 1);
            const double k2 = grid.angularSpacing(i_theta);

            const double coeff2 = 0.5 * (k1 + k2) / h2;

            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

            const int bottom = grid.index(i_r, i_theta_M1);
            const int top    = grid.index(i_r, i_theta_P1);
            const int right  = grid.index(i_r + 1, i_theta);

            /* Fill result(i,j) */
            result[center] -= (-coeff2 * arr * x[right]); /* Right */

            /* Fill result(i,j-1) */
            result[bottom] -= (-0.25 * art * x[right]); /* Top Right */

            /* Fill result(i,j+1) */
            result[top] -= (+0.25 * art * x[right]); /* Bottom Right */
        }
    }
}

template <class LevelCacheType>
static KOKKOS_INLINE_FUNCTION void nodeApplyAscOrthoCircleGiveOutside(int i_r, int i_theta, const PolarGrid& grid,
                                                                      const LevelCacheType& level_cache,
                                                                      bool DirBC_Interior, ConstVector<double>& x,
                                                                      ConstVector<double>& rhs, Vector<double>& result)
{
    assert(0 <= i_r && i_r <= grid.numberSmootherCircles());

    /* ---------------------------------------- */
    /* Compute or retrieve stencil coefficients */
    /* ---------------------------------------- */
    const int center    = grid.index(i_r, i_theta);
    const double radius = grid.radius(i_r);
    const double theta  = grid.theta(i_theta);

    double coeff_beta, arr, att, art, detDF;
    level_cache.obtainValues(i_r, i_theta, center, radius, theta, coeff_beta, arr, att, art, detDF);

    if (0 <= i_r && i_r <= grid.numberSmootherCircles()) {
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int bottom = grid.index(i_r, i_theta_M1);
        const int top    = grid.index(i_r, i_theta_P1);

        /* Fill result(i-1,j) */
        if (i_r > 1 || (i_r == 1 && !DirBC_Interior)) {
            const double h1     = grid.radialSpacing(i_r - 1);
            const double coeff1 = 0.5 * (k1 + k2) / h1;

            const int left = grid.index(i_r - 1, i_theta);

            result[left] -= (-coeff1 * arr * x[center] /* Right */
                             - 0.25 * art * x[top] /* Top Right */
                             + 0.25 * art * x[bottom]); /* Bottom Right */
        }
        /* Fill result(i+1,j) */
        if (i_r < grid.numberSmootherCircles() - 1) {
            const double h2     = grid.radialSpacing(i_r);
            const double coeff2 = 0.5 * (k1 + k2) / h2;

            const int right = grid.index(i_r + 1, i_theta);

            result[right] -= (-coeff2 * arr * x[center] /* Left */
                              + 0.25 * art * x[top] /* Top Left */
                              - 0.25 * art * x[bottom]); /* Bottom Left */
        }
    }
}

template <class LevelCacheType>
static KOKKOS_INLINE_FUNCTION void nodeApplyAscOrthoRadialGiveInside(int i_r, int i_theta, const PolarGrid& grid,
                                                                     const LevelCacheType& level_cache,
                                                                     bool DirBC_Interior, ConstVector<double>& x,
                                                                     ConstVector<double>& rhs, Vector<double>& result)
{
    assert(grid.numberSmootherCircles() - 1 <= i_r && i_r < grid.nr());

    /* ---------------------------------------- */
    /* Compute or retrieve stencil coefficients */
    /* ---------------------------------------- */
    const int center    = grid.index(i_r, i_theta);
    const double radius = grid.radius(i_r);
    const double theta  = grid.theta(i_theta);

    double coeff_beta, arr, att, art, detDF;
    level_cache.obtainValues(i_r, i_theta, center, radius, theta, coeff_beta, arr, att, art, detDF);

    /* Angular neighbours (wrapping always valid) */
    const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
    const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

    const int bottom = grid.index(i_r, i_theta_M1);
    const int top    = grid.index(i_r, i_theta_P1);
    const int left   = grid.index(i_r - 1, i_theta);

    /* h1, k1, k2 are always valid at this point */
    const double h1     = grid.radialSpacing(i_r - 1);
    const double k1     = grid.angularSpacing(i_theta_M1);
    const double k2     = grid.angularSpacing(i_theta);
    const double coeff1 = 0.5 * (k1 + k2) / h1;

    /* ------------------- */
    /* Fill result(center) */
    /* ------------------- */
    if (grid.numberSmootherCircles() <= i_r && i_r < grid.nr() - 1) {
        const double h2     = grid.radialSpacing(i_r); // safe: i_r < nr-1
        const double coeff3 = 0.5 * (h1 + h2) / k1;
        const double coeff4 = 0.5 * (h1 + h2) / k2;
        result[center] -= (-coeff3 * att * x[bottom] /* Bottom */
                           - coeff4 * att * x[top]); /* Top   */
    }

    /* ----------------- */
    /* Fill result(left) */
    /* ----------------- */
    if (grid.numberSmootherCircles() < i_r && i_r < grid.nr()) {
        result[left] -= (-0.25 * art * x[top] /* Top Right    */
                         + 0.25 * art * x[bottom]); /* Bottom Right */
    }

    /* ------------------- */
    /* Fill result(right): */
    /* ------------------- */
    if (grid.numberSmootherCircles() - 1 <= i_r && i_r < grid.nr() - 2) {
        const int right = grid.index(i_r + 1, i_theta);
        result[right] -= (+0.25 * art * x[top] /* Top Left    */
                          - 0.25 * art * x[bottom]); /* Bottom Left */
    }

    /* --------------- */
    /* Symmetry shifts */
    /* --------------- */
    if (i_r == grid.nr() - 2) {
        const double h2     = grid.radialSpacing(i_r);
        const double coeff2 = 0.5 * (k1 + k2) / h2;
        const int right     = grid.index(i_r + 1, i_theta);
        result[center] -= -coeff2 * arr * rhs[right]; /* Right: Symmetry shift! */
    }

    if (i_r == grid.nr() - 1) {
        result[left] -= -coeff1 * arr * rhs[center]; /* Right: Symmetry shift! */
    }
}

template <class LevelCacheType>
static KOKKOS_INLINE_FUNCTION void nodeApplyAscOrthoRadialGiveOutside(int i_r, int i_theta, const PolarGrid& grid,
                                                                      const LevelCacheType& level_cache,
                                                                      bool DirBC_Interior, ConstVector<double>& x,
                                                                      ConstVector<double>& rhs, Vector<double>& result)
{
    assert(grid.numberSmootherCircles() <= i_r && i_r < grid.nr());

    /* ---------------------------------------- */
    /* Compute or retrieve stencil coefficients */
    /* ---------------------------------------- */
    const int center    = grid.index(i_r, i_theta);
    const double radius = grid.radius(i_r);
    const double theta  = grid.theta(i_theta);

    double coeff_beta, arr, att, art, detDF;
    level_cache.obtainValues(i_r, i_theta, center, radius, theta, coeff_beta, arr, att, art, detDF);

    if (grid.numberSmootherCircles() <= i_r && i_r < grid.nr() - 1) {
        double h1 = grid.radialSpacing(i_r - 1);
        double h2 = grid.radialSpacing(i_r);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;

        int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        int left   = grid.index(i_r - 1, i_theta);
        int right  = grid.index(i_r + 1, i_theta);
        int bottom = grid.index(i_r, i_theta_M1);
        int top    = grid.index(i_r, i_theta_P1);

        /* Fill result(i,j-1) */
        result[bottom] -= (-coeff3 * att * x[center] /* Top */
                           - 0.25 * art * x[right] /* Top Right */
                           + 0.25 * art * x[left]); /* Top Left */
        /* Fill result(i,j+1) */
        result[top] -= (-coeff4 * att * x[center] /* Bottom */
                        + 0.25 * art * x[right] /* Bottom Right */
                        - 0.25 * art * x[left]); /* Bottom Left */
    }
}

} // namespace smoother_give

template <class LevelCacheType>
void SmootherGive<LevelCacheType>::applyAscOrthoBlackCircleSection(ConstVector<double> x, ConstVector<double> rhs,
                                                                   Vector<double> temp)
{
    using smoother_give::nodeApplyAscOrthoCircleGiveInside;
    using smoother_give::nodeApplyAscOrthoCircleGiveOutside;

    auto getBatchCount = [](int start, int end, int offset) {
        if (start >= end) {
            return 0;
        }
        return (end - start + offset - 1) / offset;
    };

    const PolarGrid& grid             = Smoother<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = Smoother<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = Smoother<LevelCacheType>::DirBC_Interior_;

    /* ----------------------------------------------- */
    /* 1. Black-Circle update (u_bc):                  */
    /*    A_bc * u_bc = f_bc − A_bc^ortho * u_bc^ortho */
    /* ----------------------------------------------- */
    {
        /* Inside Black Section */
        const int start  = (grid.numberSmootherCircles() % 2 == 0) ? 1 : 0;
        const int end    = grid.numberSmootherCircles();
        const int offset = 2;

        Kokkos::parallel_for(
            "SmootherGive: ApplyAscOrtho (Black Circle - Inside)",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, getBatchCount(start, end, offset)),
            KOKKOS_LAMBDA(const int circle_task) {
                const int i_r = start + circle_task * offset;
                for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                    nodeApplyAscOrthoCircleGiveInside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
                }
            });
        Kokkos::fence();
    }

    {
        /* Outside Black Section (Part 1)*/
        const int start  = (grid.numberSmootherCircles() % 2 == 0) ? 0 : 1;
        const int end    = grid.numberSmootherCircles() + 1;
        const int offset = 4;
        Kokkos::parallel_for(
            "SmootherGive: ApplyAscOrtho (Black Circle - Outside: Part 1)",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, getBatchCount(start, end, offset)),
            KOKKOS_LAMBDA(const int circle_task) {
                const int i_r = start + circle_task * offset;
                for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                    nodeApplyAscOrthoCircleGiveOutside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
                }
            });
        Kokkos::fence();
    }

    {
        /* Outside Black Section (Part 2)*/
        const int start  = (grid.numberSmootherCircles() % 2 == 0) ? 2 : 3;
        const int end    = grid.numberSmootherCircles() + 1;
        const int offset = 4;
        Kokkos::parallel_for(
            "SmootherGive: ApplyAscOrtho (Black Circle - Outside: Part 2)",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, getBatchCount(start, end, offset)),
            KOKKOS_LAMBDA(const int circle_task) {
                const int i_r = start + circle_task * offset;
                for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                    nodeApplyAscOrthoCircleGiveOutside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
                }
            });
        Kokkos::fence();
    }
}

template <class LevelCacheType>
void SmootherGive<LevelCacheType>::applyAscOrthoWhiteCircleSection(ConstVector<double> x, ConstVector<double> rhs,
                                                                   Vector<double> temp)
{
    using smoother_give::nodeApplyAscOrthoCircleGiveInside;
    using smoother_give::nodeApplyAscOrthoCircleGiveOutside;

    auto getBatchCount = [](int start, int end, int offset) {
        if (start >= end) {
            return 0;
        }
        return (end - start + offset - 1) / offset;
    };

    const PolarGrid& grid             = Smoother<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = Smoother<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = Smoother<LevelCacheType>::DirBC_Interior_;

    /* ----------------------------------------------- */
    /* 2. White-Circle update (u_wc):                  */
    /*    A_wc * u_wc = f_wc − A_wc^ortho * u_wc^ortho */
    /* ----------------------------------------------- */

    {
        /* Inside White Section */
        const int start  = (grid.numberSmootherCircles() % 2 == 0) ? 0 : 1;
        const int end    = grid.numberSmootherCircles();
        const int offset = 2;

        Kokkos::parallel_for(
            "SmootherGive: ApplyAscOrtho (White Circle - Inside)",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, getBatchCount(start, end, offset)),
            KOKKOS_LAMBDA(const int circle_task) {
                const int i_r = start + circle_task * offset;
                for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                    nodeApplyAscOrthoCircleGiveInside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
                }
            });
        Kokkos::fence();
    }

    {
        /* Outside White Section (Part 1)*/
        const int start  = (grid.numberSmootherCircles() % 2 == 0) ? 1 : 0;
        const int end    = grid.numberSmootherCircles() + 1;
        const int offset = 4;
        Kokkos::parallel_for(
            "SmootherGive: ApplyAscOrtho (White Circle - Outside: Part 1)",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, getBatchCount(start, end, offset)),
            KOKKOS_LAMBDA(const int circle_task) {
                const int i_r = start + circle_task * offset;
                for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                    nodeApplyAscOrthoCircleGiveOutside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
                }
            });
        Kokkos::fence();
    }

    {
        /* Outside White Section (Part 2)*/
        const int start  = (grid.numberSmootherCircles() % 2 == 0) ? 3 : 2;
        const int end    = grid.numberSmootherCircles() + 1;
        const int offset = 4;
        Kokkos::parallel_for(
            "SmootherGive: ApplyAscOrtho (White Circle - Outside: Part 2)",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, getBatchCount(start, end, offset)),
            KOKKOS_LAMBDA(const int circle_task) {
                const int i_r = start + circle_task * offset;
                for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                    nodeApplyAscOrthoCircleGiveOutside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
                }
            });
        Kokkos::fence();
    }
}

template <class LevelCacheType>
void SmootherGive<LevelCacheType>::applyAscOrthoBlackRadialSection(ConstVector<double> x, ConstVector<double> rhs,
                                                                   Vector<double> temp)
{
    using smoother_give::nodeApplyAscOrthoRadialGiveInside;
    using smoother_give::nodeApplyAscOrthoRadialGiveOutside;

    auto getBatchCount = [](int start, int end, int offset) {
        if (start >= end) {
            return 0;
        }
        return (end - start + offset - 1) / offset;
    };

    const PolarGrid& grid             = Smoother<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = Smoother<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = Smoother<LevelCacheType>::DirBC_Interior_;

    /* ----------------------------------------------- */
    /* 3. Black-Radial update (u_br):                  */
    /*    A_br * u_br = f_br − A_br^ortho * u_br^ortho */
    /* ----------------------------------------------- */

    {
        /* Inside Black Section */
        const int start  = 0;
        const int end    = grid.ntheta();
        const int offset = 2;

        Kokkos::parallel_for(
            "SmootherGive: ApplyAscOrtho (Black Radial - Inside)",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, getBatchCount(start, end, offset)),
            KOKKOS_LAMBDA(const int radial_task) {
                const int i_theta = start + radial_task * offset;
                for (int i_r = grid.numberSmootherCircles() - 1; i_r < grid.nr(); i_r++) {
                    nodeApplyAscOrthoRadialGiveInside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
                }
            });
        Kokkos::fence();
    }

    {
        /* Outside Black Section (Part 1) */
        const int start  = 1;
        const int end    = grid.ntheta();
        const int offset = 4;

        Kokkos::parallel_for(
            "SmootherGive: ApplyAscOrtho (Black Radial - Outside: Part 1)",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, getBatchCount(start, end, offset)),
            KOKKOS_LAMBDA(const int radial_task) {
                const int i_theta = start + radial_task * offset;
                for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                    nodeApplyAscOrthoRadialGiveOutside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
                }
            });
        Kokkos::fence();
    }

    {
        /* Outside Black Section (Part 1) */
        const int start  = 3;
        const int end    = grid.ntheta();
        const int offset = 4;

        Kokkos::parallel_for(
            "SmootherGive: ApplyAscOrtho (Black Radial - Outside: Part 1)",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, getBatchCount(start, end, offset)),
            KOKKOS_LAMBDA(const int radial_task) {
                const int i_theta = start + radial_task * offset;
                for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                    nodeApplyAscOrthoRadialGiveOutside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
                }
            });
        Kokkos::fence();
    }
}

template <class LevelCacheType>
void SmootherGive<LevelCacheType>::applyAscOrthoWhiteRadialSection(ConstVector<double> x, ConstVector<double> rhs,
                                                                   Vector<double> temp)
{
    using smoother_give::nodeApplyAscOrthoRadialGiveInside;
    using smoother_give::nodeApplyAscOrthoRadialGiveOutside;

    auto getBatchCount = [](int start, int end, int offset) {
        if (start >= end) {
            return 0;
        }
        return (end - start + offset - 1) / offset;
    };

    const PolarGrid& grid             = Smoother<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = Smoother<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = Smoother<LevelCacheType>::DirBC_Interior_;

    /* ----------------------------------------------- */
    /* 4. White-Radial update (u_wr):                  */
    /*    A_wr * u_wr = f_wr − A_wr^ortho * u_wr^ortho */
    /* ----------------------------------------------- */

    {
        /* Inside White Section */
        const int start  = 1;
        const int end    = grid.ntheta();
        const int offset = 2;

        Kokkos::parallel_for(
            "SmootherGive: ApplyAscOrtho (White Radial - Inside)",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, getBatchCount(start, end, offset)),
            KOKKOS_LAMBDA(const int radial_task) {
                const int i_theta = start + radial_task * offset;
                for (int i_r = grid.numberSmootherCircles() - 1; i_r < grid.nr(); i_r++) {
                    nodeApplyAscOrthoRadialGiveInside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
                }
            });
        Kokkos::fence();
    }

    {
        /* Outside White Section (Part 1) */
        const int start  = 0;
        const int end    = grid.ntheta();
        const int offset = 4;

        Kokkos::parallel_for(
            "SmootherGive: ApplyAscOrtho (White Radial - Outside: Part 1)",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, getBatchCount(start, end, offset)),
            KOKKOS_LAMBDA(const int radial_task) {
                const int i_theta = start + radial_task * offset;
                for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                    nodeApplyAscOrthoRadialGiveOutside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
                }
            });
        Kokkos::fence();
    }

    {
        /* Outside White Section (Part 1) */
        const int start  = 2;
        const int end    = grid.ntheta();
        const int offset = 4;

        Kokkos::parallel_for(
            "SmootherGive: ApplyAscOrtho (White Radial - Outside: Part 1)",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, getBatchCount(start, end, offset)),
            KOKKOS_LAMBDA(const int radial_task) {
                const int i_theta = start + radial_task * offset;
                for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                    nodeApplyAscOrthoRadialGiveOutside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
                }
            });
        Kokkos::fence();
    }
}
