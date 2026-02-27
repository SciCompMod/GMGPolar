#include "../../../include/Smoother/SmootherTake/smootherTake.h"

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

void SmootherTake::applyAscOrthoCircleSection(int i_r, ConstVector<double> x, ConstVector<double> rhs,
                                              Vector<double> temp)
{
    assert(i_r >= 0 && i_r < grid_.numberSmootherCircles());

    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    ConstVector<double> arr        = level_cache_.arr();
    ConstVector<double> att        = level_cache_.att();
    ConstVector<double> art        = level_cache_.art();
    ConstVector<double> detDF      = level_cache_.detDF();
    ConstVector<double> coeff_beta = level_cache_.coeff_beta();

    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        nodeApplyAscOrthoCircleTake(i_r, i_theta, grid_, DirBC_Interior_, x, rhs, temp, arr, att, art, detDF,
                                    coeff_beta);
    }
}

void SmootherTake::applyAscOrthoRadialSection(int i_theta, ConstVector<double> x, ConstVector<double> rhs,
                                              Vector<double> temp)
{
    assert(i_theta >= 0 && i_theta < grid_.ntheta());

    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    ConstVector<double> arr        = level_cache_.arr();
    ConstVector<double> att        = level_cache_.att();
    ConstVector<double> art        = level_cache_.art();
    ConstVector<double> detDF      = level_cache_.detDF();
    ConstVector<double> coeff_beta = level_cache_.coeff_beta();

    for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
        nodeApplyAscOrthoRadialTake(i_r, i_theta, grid_, DirBC_Interior_, x, rhs, temp, arr, att, art, detDF,
                                    coeff_beta);
    }
}