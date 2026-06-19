#pragma once

namespace extrapolated_smoother_give
{

template <typename LevelCacheType>
static KOKKOS_INLINE_FUNCTION void
nodeApplyAscOrthoCircleGiveInside(const int i_r, const int i_theta, const PolarGrid& grid,
                                  const LevelCacheType& level_cache, const bool DirBC_Interior, ConstVector<double>& x,
                                  ConstVector<double>& rhs, Vector<double>& result)
{
    KOKKOS_ASSERT(i_r >= 0 && i_r < grid.numberSmootherCircles());

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
        const double coeff3 = 0.5 * (h1 + h2) / k1;
        const double coeff4 = 0.5 * (h1 + h2) / k2;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int left   = grid.index(i_r - 1, i_theta);
        const int right  = grid.index(i_r + 1, i_theta);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int top    = grid.index(i_r, i_theta_P1);

        if (i_r & 1) {
            /* i_r % 2 == 1 and i_theta % 2 == 1 */
            /* | x | o | x | */
            /* |   |   |   | */
            /* | o | O | o | */
            /* |   |   |   | */
            /* | x | o | x | */
            /* or */
            /* i_r % 2 == 1 and i_theta % 2 == 0 */
            /* | o | o | o | */
            /* |   |   |   | */
            /* | x | O | x | */
            /* |   |   |   | */
            /* | o | o | o | */

            /* Fill result(i,j) */
            result[center] -= (-coeff1 * arr * x[left] /* Left */
                               - coeff2 * arr * x[right] /* Right */
            ); /* Fill result(i,j-1) */
            result[bottom] -= (-0.25 * art * x[right] /* Top Right */
                               + 0.25 * art * x[left]); /* Top Left */
            /* Fill result(i,j+1) */
            result[top] -= (+0.25 * art * x[right] /* Bottom Right */
                            - 0.25 * art * x[left]); /* Bottom Left */
        }
        else {
            if (i_theta & 1) {
                /* i_r % 2 == 0 and i_theta % 2 == 1 */
                /* | o | x | o | */
                /* |   |   |   | */
                /* | o | O | o | */
                /* |   |   |   | */
                /* | o | x | o | */

                /* Fill result(i,j) */
                result[center] -= (-coeff1 * arr * x[left] /* Left */
                                   - coeff2 * arr * x[right] /* Right */
                                   - coeff3 * att * x[bottom] /* Bottom */
                                   - coeff4 * att * x[top] /* Top */
                );
            }
            else {
                /* i_r % 2 == 0 and i_theta % 2 == 0 */
                /* | o | o | o | */
                /* |   |   |   | */
                /* | o | X | o | */
                /* |   |   |   | */
                /* | o | o | o | */

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
    }
    /* -------------------- */
    /* Node on the boundary */
    /* -------------------- */
    else if (i_r == 0) {
        /* ------------------------------------------------ */
        /* Case 1: Dirichlet boundary on the inner boundary */
        /* ------------------------------------------------ */
        if (DirBC_Interior) {

            if (i_theta & 1) {
                /* i_theta % 2 == 1 */
                /* || x | o | x | */
                /* ||   |   |   | */
                /* || O | o | o | */
                /* ||   |   |   | */
                /* || x | o | x | */

                /* Nothing to do! */
            }
            else {
                /* i_theta % 2 == 0 */
                /* || o | o | o | */
                /* ||   |   |   | */
                /* || X | o | x | */
                /* ||   |   |   | */
                /* || o | o | o | */

                /* Nothing to do! */
            }
        }
        else {
            /* ------------------------------------------------------------- */
            /* Case 2: Across origin discretization on the interior boundary */
            /* ------------------------------------------------------------- */
            // h1 gets replaced with 2 * R0.
            // (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()/2)).
            // Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil.
            const double h1 = 2.0 * grid.radius(0);
            const double h2 = grid.radialSpacing(i_r);
            const double k1 = grid.angularSpacing(i_theta - 1);
            const double k2 = grid.angularSpacing(i_theta);

            const double coeff2 = 0.5 * (k1 + k2) / h2;
            const double coeff3 = 0.5 * (h1 + h2) / k1;
            const double coeff4 = 0.5 * (h1 + h2) / k2;

            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);
            //const int i_theta_Across = grid.wrapThetaIndex(i_theta + grid.ntheta() / 2);

            //const int left   = grid.index(i_r, i_theta_Across);
            const int bottom = grid.index(i_r, i_theta_M1);
            const int top    = grid.index(i_r, i_theta_P1);
            const int right  = grid.index(i_r + 1, i_theta);

            if (i_theta & 1) {
                /* i_theta % 2 == 1 */
                /* -| x | o | x | */
                /* -|   |   |   | */
                /* -| O | o | o | */
                /* -|   |   |   | */
                /* -| x | o | x | */

                /* Fill result(i,j) */
                result[center] -= (
                    /* - coeff1 * arr * x[left] // Left: Not in Asc_ortho */
                    -coeff2 * arr * x[right] /* Right */
                    - coeff3 * att * x[bottom] /* Bottom */
                    - coeff4 * att * x[top] /* Top */
                );
            }
            else {
                /* i_theta % 2 == 0 */
                /* -| o | o | o | */
                /* -|   |   |   | */
                /* -| X | o | x | */
                /* -|   |   |   | */
                /* -| o | o | o | */

                /* Fill result(i,j-1) */
                result[bottom] -= (-coeff3 * att * x[center] /* Top */
                                   - 0.25 * art * x[right]); /* Top Right */
                /*  + 0.25 * art * x[grid.index(i_r-1,i_theta)]; // Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
                /* Fill result(i,j+1) */
                result[top] -= (-coeff4 * att * x[center] /* Bottom */
                                + 0.25 * art * x[right]); /* Bottom Right */
                /*  - 0.25 * art * x[grid.index(i_r-1,i_theta)]; // Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
            }
        }
    }
}

template <typename LevelCacheType>
static KOKKOS_INLINE_FUNCTION void
nodeApplyAscOrthoCircleGiveOutside(const int i_r, const int i_theta, const PolarGrid& grid,
                                   const LevelCacheType& level_cache, const bool DirBC_Interior, ConstVector<double>& x,
                                   ConstVector<double>& rhs, Vector<double>& result)
{
    KOKKOS_ASSERT(i_r >= 0 && i_r <= grid.numberSmootherCircles());

    bool give_left  = false;
    bool give_right = false;

    /* -------------------- */
    /* Node in the interior */
    /* -------------------- */
    if (i_r > 0 && i_r < grid.numberSmootherCircles()) {
        if (i_r & 1) {
            if (i_theta & 1) {
                /* i_r % 2 == 1 and i_theta % 2 == 1 */
                /* | x | o | x | */
                /* |   |   |   | */
                /* | o | O | o | */
                /* |   |   |   | */
                /* | x | o | x | */

                if (i_r > 1 || !DirBC_Interior) {
                    give_left = true;
                }

                if (i_r < grid.numberSmootherCircles() - 1) {
                    give_right = true;
                }
            }
            else {
                /* i_r % 2 == 1 and i_theta % 2 == 0 */
                /* | o | o | o | */
                /* |   |   |   | */
                /* | x | O | x | */
                /* |   |   |   | */
                /* | o | o | o | */

                /* Nothing to do! */
            }
        }
        else {
            if (i_theta & 1) {
                /* i_r % 2 == 0 and i_theta % 2 == 1 */
                /* | o | x | o | */
                /* |   |   |   | */
                /* | o | O | o | */
                /* |   |   |   | */
                /* | o | x | o | */

                if (i_r > 1 || !DirBC_Interior) {
                    give_left = true;
                }

                if (i_r < grid.numberSmootherCircles() - 1) {
                    give_right = true;
                }
            }
            else {
                /* i_r % 2 == 0 and i_theta % 2 == 0 */
                /* | o | o | o | */
                /* |   |   |   | */
                /* | o | X | o | */
                /* |   |   |   | */
                /* | o | o | o | */

                if (i_r > 1 || !DirBC_Interior) {
                    give_left = true;
                }

                if (i_r < grid.numberSmootherCircles() - 1) {
                    give_right = true;
                }
            }
        }
    }
    /* -------------------- */
    /* Node on the boundary */
    /* -------------------- */
    else if (i_r == 0) {
        give_right = true;
    }
    /* ----------------------------- */
    /* Node next to circular section */
    /* ----------------------------- */
    else if (i_r == grid.numberSmootherCircles()) {
        /* i_theta % 2 == 1 and i_r % 2 == 1 */
        /* | x | o | x || o   x   o   x  */
        /* |   |   |   || -------------- */
        /* | o | o | o || O   o   o   o  */
        /* |   |   |   || -------------- */
        /* | x | o | x || o   x   o   x  */
        /* -> Give Left */

        /* i_theta % 2 == 1 and i_r % 2 == 0 */
        /* | o | x | o || x   o   x   o  */
        /* |   |   |   || -------------- */
        /* | o | o | o || O   o   o   o  */
        /* |   |   |   || -------------- */
        /* | o | x | o || x   o   x   o  */
        /* -> Give Left */

        /* i_theta % 2 == 0 and i_r % 2 == 1 */
        /* | o | o | o || o   o   o   o  */
        /* |   |   |   || -------------- */
        /* | x | o | x || O   x   o   x  */
        /* |   |   |   || -------------- */
        /* | o | o | o || o   o   o   o  */
        /* -> Don't give to the Left! */

        /* i_theta % 2 == 0 and i_r % 2 == 0 */
        /* | o | o | o || o   o   o   o  */
        /* |   |   |   || -------------- */
        /* | o | x | o || X   o   x   o  */
        /* |   |   |   || -------------- */
        /* | o | o | o || o   o   o   o  */
        /* -> Give Left */

        if (i_theta & 1 || !(i_r & 1)) {
            give_left = true;
        }
    }

    /* ---------------------------------------- */
    /* Compute or retrieve stencil coefficients */
    /* ---------------------------------------- */
    const int center    = grid.index(i_r, i_theta);
    const double radius = grid.radius(i_r);
    const double theta  = grid.theta(i_theta);

    double coeff_beta, arr, att, art, detDF;
    level_cache.obtainValues(i_r, i_theta, center, radius, theta, coeff_beta, arr, att, art, detDF);

    const double k1 = grid.angularSpacing(i_theta - 1);
    const double k2 = grid.angularSpacing(i_theta);

    const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
    const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

    const int bottom = grid.index(i_r, i_theta_M1);
    const int top    = grid.index(i_r, i_theta_P1);

    if (give_left) {
        const double h1     = grid.radialSpacing(i_r - 1);
        const double coeff1 = 0.5 * (k1 + k2) / h1;
        const int left      = grid.index(i_r - 1, i_theta);

        /* Fill result(i-1,j) */
        result[left] -= (-coeff1 * arr * x[center] /* Right */
                         - 0.25 * art * x[top] /* Top Right */
                         + 0.25 * art * x[bottom]); /* Bottom Right */
    }

    if (give_right) {
        const double h2     = grid.radialSpacing(i_r);
        const double coeff2 = 0.5 * (k1 + k2) / h2;
        const int right     = grid.index(i_r + 1, i_theta);

        /* Fill result(i+1,j) */
        result[right] -= (-coeff2 * arr * x[center] /* Left */
                          + 0.25 * art * x[top] /* Top Left */
                          - 0.25 * art * x[bottom]); /* Bottom Left */
    }
}

template <typename LevelCacheType>
static KOKKOS_INLINE_FUNCTION void
nodeApplyAscOrthoRadialGiveInside(const int i_r, const int i_theta, const PolarGrid& grid,
                                  const LevelCacheType& level_cache, const bool DirBC_Interior, ConstVector<double>& x,
                                  ConstVector<double>& rhs, Vector<double>& result)
{
    KOKKOS_ASSERT(i_r >= grid.numberSmootherCircles() - 1 && i_r < grid.nr());

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
    if (i_r > grid.numberSmootherCircles() && i_r < grid.nr() - 2) {
        const double h1 = grid.radialSpacing(i_r - 1);
        const double h2 = grid.radialSpacing(i_r);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff1 = 0.5 * (k1 + k2) / h1;
        const double coeff2 = 0.5 * (k1 + k2) / h2;
        const double coeff3 = 0.5 * (h1 + h2) / k1;
        const double coeff4 = 0.5 * (h1 + h2) / k2;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int left   = grid.index(i_r - 1, i_theta);
        const int right  = grid.index(i_r + 1, i_theta);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int top    = grid.index(i_r, i_theta_P1);

        if (i_theta & 1) {

            /* i_theta % 2 == 1 and i_r % 2 == 1 */
            /* ---------- */
            /* x   o   x  */
            /* ---------- */
            /* o   O   o  */
            /* ---------- */
            /* x   o   x  */
            /* ---------- */
            /* or */
            /* i_theta % 2 == 1 and i_r % 2 == 0 */
            /* ---------- */
            /* o   x   o  */
            /* ---------- */
            /* o   O   o  */
            /* ---------- */
            /* o   x   o  */
            /* ---------- */

            /* Fill result(i,j) */
            result[center] -= (-coeff3 * att * x[bottom] /* Bottom */
                               - coeff4 * att * x[top]); /* Top */
            /* Fill result(i-1,j) */
            result[left] -= (-0.25 * art * x[top] /* Top Right */
                             + 0.25 * art * x[bottom]); /* Bottom Right */
            /* Fill result(i+1,j) */
            result[right] -= (+0.25 * art * x[top] /* Top Left */
                              - 0.25 * art * x[bottom]); /* Bottom Left */
        }
        else {
            if (i_r & 1) {
                /* i_theta % 2 == 0 and i_r % 2 == 1 */
                /* ---------- */
                /* o   o   o  */
                /* ---------- */
                /* x   O   x  */
                /* ---------- */
                /* o   o   o  */
                /* ---------- */

                /* Fill result(i,j) */
                result[center] -= (-coeff1 * arr * x[left] /* Left */
                                   - coeff2 * arr * x[right] /* Right */
                                   - coeff3 * att * x[bottom] /* Bottom */
                                   - coeff4 * att * x[top]); /* Top */
            }
            else {
                /* i_theta % 2 == 0 and i_r % 2 == 0 */
                /* ---------- */
                /* o   o   o  */
                /* ---------- */
                /* o   X   o  */
                /* ---------- */
                /* o   o   o  */
                /* ---------- */

                /* Fill result(i-1,j) */
                result[left] -= (-coeff1 * arr * x[center] /* Right */
                                 - 0.25 * art * x[top] /* Top Right */
                                 + 0.25 * art * x[bottom]); /* Bottom Right */
                /* Fill result(i+1,j) */
                result[right] -= (-coeff2 * arr * x[center] /* Left */
                                  + 0.25 * art * x[top] /* Top Left */
                                  - 0.25 * art * x[bottom]); /* Bottom Left */
            }
        }
    }
    else if (i_r == grid.numberSmootherCircles() - 1) {
        const double h2 = grid.radialSpacing(i_r);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff2 = 0.5 * (k1 + k2) / h2;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int right  = grid.index(i_r + 1, i_theta);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int top    = grid.index(i_r, i_theta_P1);

        /* Dont give to the right when this case occurs! */
        /* i_theta % 2 = 0 and i_r % 2 == 1 */
        /* | o | o | o || o   o   o   o  */
        /* |   |   |   || -------------- */
        /* | o | x | O || x   o   x   o  */
        /* |   |   |   || -------------- */
        /* | o | o | o || o   o   o   o  */
        if ((!(i_r & 1) || (i_theta & 1))) {
            /* Fill result(i+1,j) */
            result[right] -= (-coeff2 * arr * x[center] /* Left */
                              + 0.25 * art * x[top] /* Top Left */
                              - 0.25 * art * x[bottom]); /* Bottom Left */
        }
    }
    else if (i_r == grid.numberSmootherCircles()) {
        const double h1 = grid.radialSpacing(i_r - 1);
        const double h2 = grid.radialSpacing(i_r);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff1 = 0.5 * (k1 + k2) / h1;
        const double coeff2 = 0.5 * (k1 + k2) / h2;
        const double coeff3 = 0.5 * (h1 + h2) / k1;
        const double coeff4 = 0.5 * (h1 + h2) / k2;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int left   = grid.index(i_r - 1, i_theta);
        const int right  = grid.index(i_r + 1, i_theta);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int top    = grid.index(i_r, i_theta_P1);

        if (i_theta & 1) {
            if (i_r & 1) {
                /* i_theta % 2 == 1 and i_r % 2 == 1 */
                /* | x | o | x || o   x   o   x  */
                /* |   |   |   || -------------- */
                /* | o | o | o || O   o   o   o  */
                /* |   |   |   || -------------- */
                /* | x | o | x || o   x   o   x  */

                /* Fill result(i,j) */
                result[center] -= (-coeff1 * arr * x[left] /* Left */
                                   - coeff3 * att * x[bottom] /* Bottom */
                                   - coeff4 * att * x[top] /* Top */
                );
                /* Fill result(i+1,j) */
                result[right] -= (+0.25 * art * x[top] /* Top Left */
                                  - 0.25 * art * x[bottom]); /* Bottom Left */
            }
            else {
                /* i_theta % 2 == 1 and i_r % 2 == 0 */
                /* | o | x | o || x   o   x   o  */
                /* |   |   |   || -------------- */
                /* | o | o | o || O   o   o   o  */
                /* |   |   |   || -------------- */
                /* | o | x | o || x   o   x   o  */

                /* Fill result(i,j) */
                result[center] -= (-coeff1 * arr * x[left] /* Left */
                                   - coeff3 * att * x[bottom] /* Bottom */
                                   - coeff4 * att * x[top] /* Top */
                );
                /* Fill result(i+1,j) */
                result[right] -= (+0.25 * art * x[top] /* Top Left */
                                  - 0.25 * art * x[bottom]); /* Bottom Left */
            }
        }
        else {
            if (i_r & 1) {
                /* i_theta % 2 == 0 and i_r % 2 == 1 */
                /* | o | o | o || o   o   o   o  */
                /* |   |   |   || -------------- */
                /* | x | o | x || O   x   o   x  */
                /* |   |   |   || -------------- */
                /* | o | o | o || o   o   o   o  */

                /* Fill result(i,j) */
                result[center] -= (-coeff1 * arr * x[left] /* Left */
                                   - coeff2 * arr * x[right] /* Right */
                                   - coeff3 * att * x[bottom] /* Bottom */
                                   - coeff4 * att * x[top] /* Top */
                );
            }
            else {
                /* i_theta % 2 == 0 and i_r % 2 == 0 */
                /* | o | o | o || o   o   o   o  */
                /* |   |   |   || -------------- */
                /* | o | x | o || X   o   x   o  */
                /* |   |   |   || -------------- */
                /* | o | o | o || o   o   o   o  */

                /* Fill result(i+1,j) */
                result[right] -= (-coeff2 * arr * x[center] /* Left */
                                  + 0.25 * art * x[top] /* Top Left */
                                  - 0.25 * art * x[bottom]); /* Bottom Left */
            }
        }
    }
    else if (i_r == grid.nr() - 2) {
        KOKKOS_ASSERT(i_r & 1);

        const double h1 = grid.radialSpacing(i_r - 1);
        const double h2 = grid.radialSpacing(i_r);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff1 = 0.5 * (k1 + k2) / h1;
        const double coeff2 = 0.5 * (k1 + k2) / h2;
        const double coeff3 = 0.5 * (h1 + h2) / k1;
        const double coeff4 = 0.5 * (h1 + h2) / k2;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int left   = grid.index(i_r - 1, i_theta);
        const int right  = grid.index(i_r + 1, i_theta);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int top    = grid.index(i_r, i_theta_P1);

        if (i_theta & 1) {
            /* i_theta % 2 == 1 */
            /* ---------------|| */
            /* o   x   o   x  || */
            /* ---------------|| */
            /* o   o   O   o  || */
            /* ---------------|| */
            /* o   x   o   x  || */
            /* ---------------|| */

            /* Fill result(i,j) */
            result[center] -= (-coeff3 * att * x[bottom] /* Bottom */
                               - coeff4 * att * x[top] /* Top */
            );
            /* Fill result(i-1,j) */
            result[left] -= (-0.25 * art * x[top] /* Top Right */
                             + 0.25 * art * x[bottom]); /* Bottom Right */

            /* "Right" is part of the radial Asc smoother matrices, */
            /* but is shifted over to the rhs to make the radial Asc smoother matrices symmetric. */
            /* Note that the circle Asc smoother matrices are symmetric by default. */
            /* Note that rhs[right] contains the correct boundary value of u_D. */
            result[center] -= -coeff2 * arr * rhs[right]; /* Right: Symmetry shift! */
        }
        else {
            /* ---------------|| */
            /* o   o   o   o  || */
            /* ---------------|| */
            /* o   x   O   x  || */
            /* ---------------|| */
            /* o   o   o   o  || */
            /* ---------------|| */

            /* Fill result(i,j) */
            result[center] -= (-coeff1 * arr * x[left] /* Left */
                               - coeff2 * arr * x[right] /* Right */
                               - coeff3 * att * x[bottom] /* Bottom */
                               - coeff4 * att * x[top] /* Top */
            );
        }
    }
    else if (i_r == grid.nr() - 1) {
        KOKKOS_ASSERT(!(i_r & 1));

        const double h1 = grid.radialSpacing(i_r - 1);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff1 = 0.5 * (k1 + k2) / h1;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int left   = grid.index(i_r - 1, i_theta);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int top    = grid.index(i_r, i_theta_P1);

        if (i_theta & 1) {
            /* i_theta % 2 == 1 */
            /* -----------|| */
            /* x   o   x  || */
            /* -----------|| */
            /* o   o   O  || */
            /* -----------|| */
            /* x   o   x  || */
            /* -----------|| */

            /* Fill result(i-1,j) */
            result[left] -= (-0.25 * art * x[top] /* Top Right */
                             + 0.25 * art * x[bottom] /* Bottom Right */
            );
            /* "Right" is part of the radial Asc smoother matrices, */
            /* but is shifted over to the rhs to make the radial Asc smoother matrices symmetric. */
            /* Note that the circle Asc smoother matrices are symmetric by default. */
            /* Note that rhs[center] contains the correct boundary value of u_D. */
            result[left] -= (-coeff1 * arr * rhs[center] /* Right */
            );
        }
        else {
            /* -----------|| */
            /* o   o   o  || */
            /* -----------|| */
            /* x   o   X  || */
            /* -----------|| */
            /* o   o   o  || */
            /* -----------|| */

            /* Fill result(i-1,j) */
            result[left] -= (-coeff1 * arr * x[center] /* Right */
                             - 0.25 * art * x[top] /* Top Right */
                             + 0.25 * art * x[bottom]); /* Bottom Right */
        }
    }
}

template <typename LevelCacheType>
static KOKKOS_INLINE_FUNCTION void
nodeApplyAscOrthoRadialGiveOutside(const int i_r, const int i_theta, const PolarGrid& grid,
                                   const LevelCacheType& level_cache, const bool DirBC_Interior, ConstVector<double>& x,
                                   ConstVector<double>& rhs, Vector<double>& result)
{
    KOKKOS_ASSERT(i_r >= grid.numberSmootherCircles() && i_r < grid.nr());

    if (i_r == grid.nr() - 1) {
        return;
    }

    bool give_bottom = false;
    bool give_top    = false;

    /* -------------------- */
    /* Node in the interior */
    /* -------------------- */
    if (i_r > grid.numberSmootherCircles() && i_r < grid.nr() - 2) {
        if (i_theta & 1) {
            if (i_r & 1) {
                /* i_theta % 2 == 1 and i_r % 2 == 1 */
                /* ---------- */
                /* x   o   x  */
                /* ---------- */
                /* o   O   o  */
                /* ---------- */
                /* x   o   x  */
                /* ---------- */

                give_bottom = true;
                give_top    = true;
            }
            else {
                /* i_theta % 2 == 1 and i_r % 2 == 0 */
                /* ---------- */
                /* o   x   o  */
                /* ---------- */
                /* o   O   o  */
                /* ---------- */
                /* o   x   o  */
                /* ---------- */

                /* Nothing to do! */
            }
        }
        else {
            if (i_r & 1) {
                /* i_theta % 2 == 0 and i_r % 2 == 1 */
                /* ---------- */
                /* o   o   o  */
                /* ---------- */
                /* x   O   x  */
                /* ---------- */
                /* o   o   o  */
                /* ---------- */

                give_bottom = true;
                give_top    = true;
            }
            else {
                /* i_theta % 2 == 0 and i_r % 2 == 0 */
                /* ---------- */
                /* o   o   o  */
                /* ---------- */
                /* o   X   o  */
                /* ---------- */
                /* o   o   o  */
                /* ---------- */

                give_bottom = true;
                give_top    = true;
            }
        }
    }
    else if (i_r == grid.numberSmootherCircles()) {
        /* Dont give to bottom and up when this case occurs! */
        /* i_theta % 2 == 1 and i_r % 2 == 0 */
        /* | o | x | o || x   o   x   o  */
        /* |   |   |   || -------------- */
        /* | o | o | o || O   o   o   o  */
        /* |   |   |   || -------------- */
        /* | o | x | o || x   o   x   o  */
        if (i_r & 1 || !(i_theta & 1)) {
            give_bottom = true;
            give_top    = true;
        }
    }
    else if (i_r == grid.nr() - 2) {
        KOKKOS_ASSERT(i_r & 1);

        give_bottom = true;
        give_top    = true;
    }

    /* ---------------------------------------- */
    /* Compute or retrieve stencil coefficients */
    /* ---------------------------------------- */
    const int center    = grid.index(i_r, i_theta);
    const double radius = grid.radius(i_r);
    const double theta  = grid.theta(i_theta);

    double coeff_beta, arr, att, art, detDF;
    level_cache.obtainValues(i_r, i_theta, center, radius, theta, coeff_beta, arr, att, art, detDF);

    const double h1 = grid.radialSpacing(i_r - 1);
    const double h2 = grid.radialSpacing(i_r);

    const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
    const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

    const int left  = grid.index(i_r - 1, i_theta);
    const int right = grid.index(i_r + 1, i_theta);

    if (give_bottom) {
        const double k1     = grid.angularSpacing(i_theta - 1);
        const double coeff3 = 0.5 * (h1 + h2) / k1;
        const int bottom    = grid.index(i_r, i_theta_M1);

        /* Fill result(i,j-1) */
        result[bottom] -= (-coeff3 * att * x[center] /* Top */
                           - 0.25 * art * x[right] /* Top Right */
                           + 0.25 * art * x[left]); /* Top Left */
    }

    if (give_top) {
        const double k2     = grid.angularSpacing(i_theta);
        const double coeff4 = 0.5 * (h1 + h2) / k2;
        const int top       = grid.index(i_r, i_theta_P1);
        /* Fill result(i,j+1) */
        result[top] -= (-coeff4 * att * x[center] /* Bottom */
                        + 0.25 * art * x[right] /* Bottom Right */
                        - 0.25 * art * x[left]); /* Bottom Left */
    }
}

} // namespace extrapolated_smoother_give

template <class LevelCacheType>
void ExtrapolatedSmootherGive<LevelCacheType>::applyAscOrthoBlackCircleSection(ConstVector<double> x,
                                                                               ConstVector<double> rhs,
                                                                               Vector<double> temp)
{
    using extrapolated_smoother_give::nodeApplyAscOrthoCircleGiveInside;
    using extrapolated_smoother_give::nodeApplyAscOrthoCircleGiveOutside;

    auto getBatchCount = [](int start, int end, int offset) {
        if (start >= end) {
            return 0;
        }
        return (end - start + offset - 1) / offset;
    };

    const PolarGrid& grid             = ExtrapolatedSmoother<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = ExtrapolatedSmoother<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = ExtrapolatedSmoother<LevelCacheType>::DirBC_Interior_;

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
            "ExtrapolatedSmootherGive: ApplyAscOrtho (Black Circle - Inside)",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, getBatchCount(start, end, offset)),
            KOKKOS_LAMBDA(const int circle_task) {
                const int i_r = start + circle_task * offset;
                // Serial loop to avoid race conditions
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
            "ExtrapolatedSmootherGive: ApplyAscOrtho (Black Circle - Outside: Part 1)",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
                {0, 0}, // Starting point of the index space
                {getBatchCount(start, end, offset), grid.ntheta()} // Ending point of the index space
                ),
            // Kokkos lambda function to execute for each point in the index space
            KOKKOS_LAMBDA(const int circle_task, const int i_theta) {
                const int i_r = start + circle_task * offset;
                nodeApplyAscOrthoCircleGiveOutside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
            });

        Kokkos::fence();
    }
    {
        /* Outside Black Section (Part 2)*/
        const int start  = (grid.numberSmootherCircles() % 2 == 0) ? 2 : 3;
        const int end    = grid.numberSmootherCircles() + 1;
        const int offset = 4;

        Kokkos::parallel_for(
            "ExtrapolatedSmootherGive: ApplyAscOrtho (Black Circle - Outside: Part 2)",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
                {0, 0}, // Starting point of the index space
                {getBatchCount(start, end, offset), grid.ntheta()} // Ending point of the index space
                ),
            // Kokkos lambda function to execute for each point in the index space
            KOKKOS_LAMBDA(const int circle_task, const int i_theta) {
                const int i_r = start + circle_task * offset;
                nodeApplyAscOrthoCircleGiveOutside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
            });

        Kokkos::fence();
    }
}

template <class LevelCacheType>
void ExtrapolatedSmootherGive<LevelCacheType>::applyAscOrthoWhiteCircleSection(ConstVector<double> x,
                                                                               ConstVector<double> rhs,
                                                                               Vector<double> temp)
{
    using extrapolated_smoother_give::nodeApplyAscOrthoCircleGiveInside;
    using extrapolated_smoother_give::nodeApplyAscOrthoCircleGiveOutside;

    auto getBatchCount = [](int start, int end, int offset) {
        if (start >= end) {
            return 0;
        }
        return (end - start + offset - 1) / offset;
    };

    const PolarGrid& grid             = ExtrapolatedSmoother<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = ExtrapolatedSmoother<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = ExtrapolatedSmoother<LevelCacheType>::DirBC_Interior_;

    /* ----------------------------------------------- */
    /* 2. White-Circle update (u_wc):                  */
    /*    A_wc * u_wc = f_wc − A_wc^ortho * u_wc^ortho */
    /* ----------------------------------------------- */

    {
        /* Inside White Section */
        const int start  = (grid.numberSmootherCircles() % 2 == 0) ? 0 : 1;
        const int end    = grid.numberSmootherCircles() - 1;
        const int offset = 2;

        Kokkos::parallel_for(
            "ExtrapolatedSmootherGive: ApplyAscOrtho (White Circle - Inside)",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, getBatchCount(start, end, offset)),
            KOKKOS_LAMBDA(const int circle_task) {
                const int i_r = start + circle_task * offset;
                // Serial loop to avoid race conditions
                for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                    nodeApplyAscOrthoCircleGiveInside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
                }
            });

        Kokkos::fence();
    }
    {
        /* Outside White Section (Part 1)*/
        const int start  = (grid.numberSmootherCircles() % 2 == 0) ? 1 : 0;
        const int end    = grid.numberSmootherCircles();
        const int offset = 4;

        Kokkos::parallel_for(
            "ExtrapolatedSmootherGive: ApplyAscOrtho (White Circle - Outside: Part 1)",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
                {0, 0}, // Starting point of the index space
                {getBatchCount(start, end, offset), grid.ntheta()} // Ending point of the index space
                ),
            // Kokkos lambda function to execute for each point in the index space
            KOKKOS_LAMBDA(const int circle_task, const int i_theta) {
                const int i_r = start + circle_task * offset;
                nodeApplyAscOrthoCircleGiveOutside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
            });

        Kokkos::fence();
    }
    {
        /* Outside White Section (Part 2)*/
        const int start  = (grid.numberSmootherCircles() % 2 == 0) ? 3 : 2;
        const int end    = grid.numberSmootherCircles();
        const int offset = 4;

        Kokkos::parallel_for(
            "ExtrapolatedSmootherGive: ApplyAscOrtho (White Circle - Outside: Part 2)",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
                {0, 0}, // Starting point of the index space
                {getBatchCount(start, end, offset), grid.ntheta()} // Ending point of the index space
                ),
            // Kokkos lambda function to execute for each point in the index space
            KOKKOS_LAMBDA(const int circle_task, const int i_theta) {
                const int i_r = start + circle_task * offset;
                nodeApplyAscOrthoCircleGiveOutside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
            });

        Kokkos::fence();
    }
}

template <class LevelCacheType>
void ExtrapolatedSmootherGive<LevelCacheType>::applyAscOrthoBlackRadialSection(ConstVector<double> x,
                                                                               ConstVector<double> rhs,
                                                                               Vector<double> temp)
{
    using extrapolated_smoother_give::nodeApplyAscOrthoRadialGiveInside;
    using extrapolated_smoother_give::nodeApplyAscOrthoRadialGiveOutside;

    auto getBatchCount = [](int start, int end, int offset) {
        if (start >= end) {
            return 0;
        }
        return (end - start + offset - 1) / offset;
    };

    const PolarGrid& grid             = ExtrapolatedSmoother<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = ExtrapolatedSmoother<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = ExtrapolatedSmoother<LevelCacheType>::DirBC_Interior_;

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
            "ExtrapolatedSmootherGive: ApplyAscOrtho (Black Radial - Inside)",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, getBatchCount(start, end, offset)),
            KOKKOS_LAMBDA(const int radial_task) {
                const int i_theta = start + radial_task * offset;
                // Serial loop to avoid race conditions
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
            "ExtrapolatedSmootherGive: ApplyAscOrtho (Black Radial - Outside: Part 1)",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
                {0, grid.numberSmootherCircles()}, // Starting point of the index space
                {getBatchCount(start, end, offset), grid.nr()} // Ending point of the index space
                ),
            // Kokkos lambda function to execute for each point in the index space
            KOKKOS_LAMBDA(const int radial_task, const int i_r) {
                const int i_theta = start + radial_task * offset;
                nodeApplyAscOrthoRadialGiveOutside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
            });

        Kokkos::fence();
    }
    {
        /* Outside Black Section (Part 1) */
        const int start  = 3;
        const int end    = grid.ntheta();
        const int offset = 4;

        Kokkos::parallel_for(
            "ExtrapolatedSmootherGive: ApplyAscOrtho (Black Radial - Outside: Part 2)",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
                {0, grid.numberSmootherCircles()}, // Starting point of the index space
                {getBatchCount(start, end, offset), grid.nr()} // Ending point of the index space
                ),
            // Kokkos lambda function to execute for each point in the index space
            KOKKOS_LAMBDA(const int radial_task, const int i_r) {
                const int i_theta = start + radial_task * offset;
                nodeApplyAscOrthoRadialGiveOutside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
            });

        Kokkos::fence();
    }
}

template <class LevelCacheType>
void ExtrapolatedSmootherGive<LevelCacheType>::applyAscOrthoWhiteRadialSection(ConstVector<double> x,
                                                                               ConstVector<double> rhs,
                                                                               Vector<double> temp)
{
    using extrapolated_smoother_give::nodeApplyAscOrthoRadialGiveInside;
    using extrapolated_smoother_give::nodeApplyAscOrthoRadialGiveOutside;

    auto getBatchCount = [](int start, int end, int offset) {
        if (start >= end) {
            return 0;
        }
        return (end - start + offset - 1) / offset;
    };

    const PolarGrid& grid             = ExtrapolatedSmoother<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = ExtrapolatedSmoother<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = ExtrapolatedSmoother<LevelCacheType>::DirBC_Interior_;

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
            "ExtrapolatedSmootherGive: ApplyAscOrtho (White Radial - Inside)",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, getBatchCount(start, end, offset)),
            KOKKOS_LAMBDA(const int radial_task) {
                const int i_theta = start + radial_task * offset;
                // Serial loop to avoid race conditions
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
            "ExtrapolatedSmootherGive: ApplyAscOrtho (White Radial - Outside: Part 1)",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
                {0, grid.numberSmootherCircles()}, // Starting point of the index space
                {getBatchCount(start, end, offset), grid.nr()} // Ending point of the index space
                ),
            // Kokkos lambda function to execute for each point in the index space
            KOKKOS_LAMBDA(const int radial_task, const int i_r) {
                const int i_theta = start + radial_task * offset;
                nodeApplyAscOrthoRadialGiveOutside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
            });

        Kokkos::fence();
    }
    {
        /* Outside White Section (Part 1) */
        const int start  = 2;
        const int end    = grid.ntheta();
        const int offset = 4;

        Kokkos::parallel_for(
            "ExtrapolatedSmootherGive: ApplyAscOrtho (White Radial - Outside: Part 2)",
            Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
                {0, grid.numberSmootherCircles()}, // Starting point of the index space
                {getBatchCount(start, end, offset), grid.nr()} // Ending point of the index space
                ),
            // Kokkos lambda function to execute for each point in the index space
            KOKKOS_LAMBDA(const int radial_task, const int i_r) {
                const int i_theta = start + radial_task * offset;
                nodeApplyAscOrthoRadialGiveOutside(i_r, i_theta, grid, level_cache, DirBC_Interior, x, rhs, temp);
            });

        Kokkos::fence();
    }
}
