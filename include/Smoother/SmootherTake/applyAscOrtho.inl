#pragma once

namespace smoother_take
{

static inline void nodeApplyAscOrthoCircleTake(int i_r, int i_theta, const PolarGrid& grid, bool DirBC_Interior,
                                               ConstVector<double>& x, ConstVector<double>& rhs, Vector<double>& result,
                                               ConstVector<double>& arr, ConstVector<double>& att,
                                               ConstVector<double>& art, ConstVector<double>& detDF,
                                               ConstVector<double>& coeff_beta)
{
    assert(i_r >= 0 && i_r < grid.numberSmootherCircles());

    if (i_r > 0 && i_r < grid.numberSmootherCircles()) {
        double h1 = grid.radialSpacing(i_r - 1);
        double h2 = grid.radialSpacing(i_r);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int bottom_left  = grid.index(i_r - 1, i_theta_M1);
        const int left         = grid.index(i_r - 1, i_theta);
        const int top_left     = grid.index(i_r - 1, i_theta_P1);
        const int bottom       = grid.index(i_r, i_theta_M1);
        const int center       = grid.index(i_r, i_theta);
        const int top          = grid.index(i_r, i_theta_P1);
        const int bottom_right = grid.index(i_r + 1, i_theta_M1);
        const int right        = grid.index(i_r + 1, i_theta);
        const int top_right    = grid.index(i_r + 1, i_theta_P1);

        result[center] = rhs[center] - (-coeff1 * (arr[center] + arr[left]) * x[left] /* Left */
                                        - coeff2 * (arr[center] + arr[right]) * x[right] /* Right */

                                        - 0.25 * (art[left] + art[bottom]) * x[bottom_left] /* Bottom Left */
                                        + 0.25 * (art[right] + art[bottom]) * x[bottom_right] /* Bottom Right */
                                        + 0.25 * (art[left] + art[top]) * x[top_left] /* Top Left */
                                        - 0.25 * (art[right] + art[top]) * x[top_right] /* Top Right */
                                       );
    }
    else if (i_r == 0) {
        if (DirBC_Interior) {
            const int center = grid.index(i_r, i_theta);
            result[center]   = rhs[center];
        }
        else {
            /* ------------------------------------------------------------- */
            /* Case 2: Across origin discretization on the interior boundary */
            /* ------------------------------------------------------------- */
            // h1 gets replaced with 2 * R0.
            // (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()/2)).
            // Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil.
            double h1 = 2.0 * grid.radius(0);
            double h2 = grid.radialSpacing(i_r);
            double k1 = grid.angularSpacing(i_theta - 1);
            double k2 = grid.angularSpacing(i_theta);

            double coeff1 = 0.5 * (k1 + k2) / h1;
            double coeff2 = 0.5 * (k1 + k2) / h2;
            double coeff3 = 0.5 * (h1 + h2) / k1;
            double coeff4 = 0.5 * (h1 + h2) / k2;

            const int i_theta_M1     = grid.wrapThetaIndex(i_theta - 1);
            const int i_theta_P1     = grid.wrapThetaIndex(i_theta + 1);
            const int i_theta_Across = grid.wrapThetaIndex(i_theta + grid.ntheta() / 2);

            const int left         = grid.index(i_r, i_theta_Across);
            const int bottom       = grid.index(i_r, i_theta_M1);
            const int center       = grid.index(i_r, i_theta);
            const int top          = grid.index(i_r, i_theta_P1);
            const int bottom_right = grid.index(i_r + 1, i_theta_M1);
            const int right        = grid.index(i_r + 1, i_theta);
            const int top_right    = grid.index(i_r + 1, i_theta_P1);

            result[center] =
                rhs[center] -
                (-coeff2 * (arr[center] + arr[right]) * x[right] /* Right */

                 /* - 0.25 * (art[left] + art[bottom]) * x[bottom_left] // Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
                 + 0.25 * (art[right] + art[bottom]) * x[bottom_right] /* Bottom Right */
                 /* + 0.25 * (art[left] + art[top]) * x[top_left] // Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
                 - 0.25 * (art[right] + art[top]) * x[top_right] /* Top Right */
                );
        }
    }
}

static inline void nodeApplyAscOrthoRadialTake(int i_r, int i_theta, const PolarGrid& grid, bool DirBC_Interior,
                                               ConstVector<double>& x, ConstVector<double>& rhs, Vector<double>& result,
                                               ConstVector<double>& arr, const ConstVector<double>& att,
                                               ConstVector<double>& art, const ConstVector<double>& detDF,
                                               ConstVector<double>& coeff_beta)
{
    assert(i_r >= grid.numberSmootherCircles() && i_r < grid.nr());
    /* -------------------- */
    /* Node in the interior */
    /* -------------------- */
    if (i_r > grid.numberSmootherCircles() && i_r < grid.nr() - 2) {
        double h1 = grid.radialSpacing(i_r - 1);
        double h2 = grid.radialSpacing(i_r);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int bottom_left  = grid.index(i_r - 1, i_theta_M1);
        const int left         = grid.index(i_r - 1, i_theta);
        const int top_left     = grid.index(i_r - 1, i_theta_P1);
        const int bottom       = grid.index(i_r, i_theta_M1);
        const int center       = grid.index(i_r, i_theta);
        const int top          = grid.index(i_r, i_theta_P1);
        const int bottom_right = grid.index(i_r + 1, i_theta_M1);
        const int right        = grid.index(i_r + 1, i_theta);
        const int top_right    = grid.index(i_r + 1, i_theta_P1);

        result[center] = rhs[center] - (-coeff3 * (att[center] + att[bottom]) * x[bottom] /* Bottom */
                                        - coeff4 * (att[center] + att[top]) * x[top] /* Top */

                                        - 0.25 * (art[left] + art[bottom]) * x[bottom_left] /* Bottom Left */
                                        + 0.25 * (art[right] + art[bottom]) * x[bottom_right] /* Bottom Right */
                                        + 0.25 * (art[left] + art[top]) * x[top_left] /* Top Left */
                                        - 0.25 * (art[right] + art[top]) * x[top_right] /* Top Right */
                                       );
    }
    else if (i_r == grid.numberSmootherCircles()) {
        double h1 = grid.radialSpacing(i_r - 1);
        double h2 = grid.radialSpacing(i_r);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int bottom_left  = grid.index(i_r - 1, i_theta_M1);
        const int left         = grid.index(i_r - 1, i_theta);
        const int top_left     = grid.index(i_r - 1, i_theta_P1);
        const int bottom       = grid.index(i_r, i_theta_M1);
        const int center       = grid.index(i_r, i_theta);
        const int top          = grid.index(i_r, i_theta_P1);
        const int bottom_right = grid.index(i_r + 1, i_theta_M1);
        const int right        = grid.index(i_r + 1, i_theta);
        const int top_right    = grid.index(i_r + 1, i_theta_P1);

        result[center] = rhs[center] - (-coeff1 * (arr[center] + arr[left]) * x[left] /* Left */
                                        - coeff3 * (att[center] + att[bottom]) * x[bottom] /* Bottom */
                                        - coeff4 * (att[center] + att[top]) * x[top] /* Top */

                                        - 0.25 * (art[left] + art[bottom]) * x[bottom_left] /* Bottom Left */
                                        + 0.25 * (art[right] + art[bottom]) * x[bottom_right] /* Bottom Right */
                                        + 0.25 * (art[left] + art[top]) * x[top_left] /* Top Left */
                                        - 0.25 * (art[right] + art[top]) * x[top_right] /* Top Right */
                                       );
    }
    else if (i_r == grid.nr() - 2) {
        double h1 = grid.radialSpacing(i_r - 1);
        double h2 = grid.radialSpacing(i_r);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int bottom_left  = grid.index(i_r - 1, i_theta_M1);
        const int left         = grid.index(i_r - 1, i_theta);
        const int top_left     = grid.index(i_r - 1, i_theta_P1);
        const int bottom       = grid.index(i_r, i_theta_M1);
        const int center       = grid.index(i_r, i_theta);
        const int top          = grid.index(i_r, i_theta_P1);
        const int bottom_right = grid.index(i_r + 1, i_theta_M1);
        const int right        = grid.index(i_r + 1, i_theta);
        const int top_right    = grid.index(i_r + 1, i_theta_P1);

        /* "Right" is part of the radial Asc smoother matrices, */
        /* but is shifted over to the rhs to make the radial Asc smoother matrices symmetric. */
        /* Note that the circle Asc smoother matrices are symmetric by default. */
        /* Note that rhs[right] contains the correct boundary value of u_D. */
        result[center] = rhs[center] - (-coeff2 * (arr[center] + arr[right]) * rhs[right] /* Right */
                                        - coeff3 * (att[center] + att[bottom]) * x[bottom] /* Bottom */
                                        - coeff4 * (att[center] + att[top]) * x[top] /* Top */

                                        - 0.25 * (art[left] + art[bottom]) * x[bottom_left] /* Bottom Left */
                                        + 0.25 * (art[right] + art[bottom]) * x[bottom_right] /* Bottom Right */
                                        + 0.25 * (art[left] + art[top]) * x[top_left] /* Top Left */
                                        - 0.25 * (art[right] + art[top]) * x[top_right] /* Top Right */
                                       );
    }
    else if (i_r == grid.nr() - 1) {
        const int center = grid.index(i_r, i_theta);
        result[center]   = rhs[center];
    }
}

} // namespace smoother_take

template <class LevelCacheType>
void SmootherTake<LevelCacheType>::applyAscOrthoBlackCircleSection(ConstVector<double> x, ConstVector<double> rhs,
                                                                   Vector<double> temp)
{
    using smoother_take::nodeApplyAscOrthoCircleTake;

    const PolarGrid& grid             = Smoother<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = Smoother<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = Smoother<LevelCacheType>::DirBC_Interior_;
    const int num_omp_threads         = Smoother<LevelCacheType>::num_omp_threads_;

    assert(level_cache.cacheDensityProfileCoefficients());
    assert(level_cache.cacheDomainGeometry());

    ConstVector<double> arr        = level_cache.arr();
    ConstVector<double> att        = level_cache.att();
    ConstVector<double> art        = level_cache.art();
    ConstVector<double> detDF      = level_cache.detDF();
    ConstVector<double> coeff_beta = level_cache.coeff_beta();

    /* The outer most circle next to the radial section is defined to be black. */
    const int start_black_circles = (grid.numberSmootherCircles() % 2 == 0) ? 1 : 0;
    const int num_black_circles   = (grid.numberSmootherCircles() - start_black_circles + 1) / 2;

    Kokkos::parallel_for(
        "Smoother Take: ApplyAscOrtho (Black Circular)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>( // Rank of the index space
            {0, 0}, // Starting point of the index space
            {num_black_circles, grid.ntheta()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_CLASS_LAMBDA(const int circle_task, const int i_theta) {
            int i_r = start_black_circles + circle_task * 2;
            nodeApplyAscOrthoCircleTake(i_r, i_theta, grid, DirBC_Interior, x, rhs, temp, arr, att, art, detDF,
                                        coeff_beta);
        });

    Kokkos::fence();
}

template <class LevelCacheType>
void SmootherTake<LevelCacheType>::applyAscOrthoWhiteCircleSection(ConstVector<double> x, ConstVector<double> rhs,
                                                                   Vector<double> temp)
{
    using smoother_take::nodeApplyAscOrthoCircleTake;

    const PolarGrid& grid             = Smoother<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = Smoother<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = Smoother<LevelCacheType>::DirBC_Interior_;
    const int num_omp_threads         = Smoother<LevelCacheType>::num_omp_threads_;

    assert(level_cache.cacheDensityProfileCoefficients());
    assert(level_cache.cacheDomainGeometry());

    ConstVector<double> arr        = level_cache.arr();
    ConstVector<double> att        = level_cache.att();
    ConstVector<double> art        = level_cache.art();
    ConstVector<double> detDF      = level_cache.detDF();
    ConstVector<double> coeff_beta = level_cache.coeff_beta();

    /* The outer most circle next to the radial section is defined to be black. */
    const int start_white_circles = (grid.numberSmootherCircles() % 2 == 0) ? 0 : 1;
    const int num_white_circles   = (grid.numberSmootherCircles() - start_white_circles + 1) / 2;

    Kokkos::parallel_for(
        "Smoother Take: ApplyAscOrtho (White Circular)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>( // Rank of the index space
            {0, 0}, // Starting point of the index space
            {num_white_circles, grid.ntheta()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_CLASS_LAMBDA(const int circle_task, const int i_theta) {
            const int i_r = start_white_circles + circle_task * 2;
            nodeApplyAscOrthoCircleTake(i_r, i_theta, grid, DirBC_Interior, x, rhs, temp, arr, att, art, detDF,
                                        coeff_beta);
        });

    Kokkos::fence();
}

template <class LevelCacheType>
void SmootherTake<LevelCacheType>::applyAscOrthoBlackRadialSection(ConstVector<double> x, ConstVector<double> rhs,
                                                                   Vector<double> temp)
{
    using smoother_take::nodeApplyAscOrthoRadialTake;

    const PolarGrid& grid             = Smoother<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = Smoother<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = Smoother<LevelCacheType>::DirBC_Interior_;
    const int num_omp_threads         = Smoother<LevelCacheType>::num_omp_threads_;

    assert(level_cache.cacheDensityProfileCoefficients());
    assert(level_cache.cacheDomainGeometry());

    ConstVector<double> arr        = level_cache.arr();
    ConstVector<double> att        = level_cache.att();
    ConstVector<double> art        = level_cache.art();
    ConstVector<double> detDF      = level_cache.detDF();
    ConstVector<double> coeff_beta = level_cache.coeff_beta();

    assert(grid.ntheta() % 2 == 0);
    const int start_black_radials    = 0;
    const int num_black_radial_lines = grid.ntheta() / 2;

    Kokkos::parallel_for(
        "Smoother Take: ApplyAscOrtho (Black Radial)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>( // Rank of the index space
            {0, grid.numberSmootherCircles()}, // Starting point of the index space
            {num_black_radial_lines, grid.nr()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_CLASS_LAMBDA(const int radial_task, const int i_r) {
            const int i_theta = start_black_radials + radial_task * 2;
            nodeApplyAscOrthoRadialTake(i_r, i_theta, grid, DirBC_Interior, x, rhs, temp, arr, att, art, detDF,
                                        coeff_beta);
        });

    Kokkos::fence();
}

template <class LevelCacheType>
void SmootherTake<LevelCacheType>::applyAscOrthoWhiteRadialSection(ConstVector<double> x, ConstVector<double> rhs,
                                                                   Vector<double> temp)
{
    using smoother_take::nodeApplyAscOrthoRadialTake;

    const PolarGrid& grid             = Smoother<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = Smoother<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = Smoother<LevelCacheType>::DirBC_Interior_;
    const int num_omp_threads         = Smoother<LevelCacheType>::num_omp_threads_;

    assert(level_cache.cacheDensityProfileCoefficients());
    assert(level_cache.cacheDomainGeometry());

    ConstVector<double> arr        = level_cache.arr();
    ConstVector<double> att        = level_cache.att();
    ConstVector<double> art        = level_cache.art();
    ConstVector<double> detDF      = level_cache.detDF();
    ConstVector<double> coeff_beta = level_cache.coeff_beta();

    assert(grid.ntheta() % 2 == 0);
    const int start_white_radials    = 1;
    const int num_white_radial_lines = grid.ntheta() / 2;

    Kokkos::parallel_for(
        "Smoother Take: ApplyAscOrtho (White Radial)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>( // Rank of the index space
            {0, grid.numberSmootherCircles()}, // Starting point of the index space
            {num_white_radial_lines, grid.nr()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_CLASS_LAMBDA(const int radial_task, const int i_r) {
            const int i_theta = start_white_radials + radial_task * 2;
            nodeApplyAscOrthoRadialTake(i_r, i_theta, grid, DirBC_Interior, x, rhs, temp, arr, att, art, detDF,
                                        coeff_beta);
        });

    Kokkos::fence();
}
