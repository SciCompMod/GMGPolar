#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherGive/extrapolatedSmootherGive.h"

#include "../../../include/common/geometry_helper.h"

static inline void nodeApplyAscOrthoCircleGive(int i_r, int i_theta, const PolarGrid& grid, bool DirBC_Interior,
                                               SmootherColor smoother_color, ConstVector<double>& x,
                                               ConstVector<double>& rhs, Vector<double>& result, double arr, double att,
                                               double art, double detDF, double coeff_beta)
{
    assert(i_r >= 0 && i_r <= grid.numberSmootherCircles());

    bool isOddNumberSmootherCircles = (grid.numberSmootherCircles() & 1);
    bool isOddRadialIndex           = (i_r & 1);

    SmootherColor node_color =
        (isOddNumberSmootherCircles == isOddRadialIndex) ? SmootherColor::White : SmootherColor::Black;

    /* -------------------- */
    /* Node in the interior */
    /* -------------------- */
    if (i_r > 0 && i_r < grid.numberSmootherCircles()) {
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

        int center = grid.index(i_r, i_theta);
        int left   = grid.index(i_r - 1, i_theta);
        int right  = grid.index(i_r + 1, i_theta);
        int bottom = grid.index(i_r, i_theta_M1);
        int top    = grid.index(i_r, i_theta_P1);

        /* -------------------- */
        /* Inside Section Parts */
        if (node_color == smoother_color) {
            if (i_r & 1) {
                if (i_theta & 1) {
                    /* i_r % 2 == 1 and i_theta % 2 == 1 */
                    /* | x | o | x | */
                    /* |   |   |   | */
                    /* | o | O | o | */
                    /* |   |   |   | */
                    /* | x | o | x | */

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
                    /* i_r % 2 == 1 and i_theta % 2 == 0 */
                    /* | o | o | o | */
                    /* |   |   |   | */
                    /* | x | O | x | */
                    /* |   |   |   | */
                    /* | o | o | o | */

                    /* Fill result(i,j) */
                    result[center] -= (-coeff1 * arr * x[left] /* Left */
                                       - coeff2 * arr * x[right] /* Right */
                    );
                    /* Fill result(i,j-1) */
                    result[bottom] -= (-0.25 * art * x[right] /* Top Right */
                                       + 0.25 * art * x[left]); /* Top Left */
                    /* Fill result(i,j+1) */
                    result[top] -= (+0.25 * art * x[right] /* Bottom Right */
                                    - 0.25 * art * x[left]); /* Bottom Left */
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
        /* --------------------- */
        /* Outside Section Parts */
        else if (node_color != smoother_color) {
            if (i_r & 1) {
                if (i_theta & 1) {
                    /* i_r % 2 == 1 and i_theta % 2 == 1 */
                    /* | x | o | x | */
                    /* |   |   |   | */
                    /* | o | O | o | */
                    /* |   |   |   | */
                    /* | x | o | x | */

                    /* Fill result(i-1,j) */
                    if (i_r > 1 || !DirBC_Interior) {
                        result[left] -= (-coeff1 * arr * x[center] /* Right */
                                         - 0.25 * art * x[top] /* Top Right */
                                         + 0.25 * art * x[bottom]); /* Bottom Right */
                    } /* Fill result(i+1,j) */
                    if (i_r < grid.numberSmootherCircles() - 1) {
                        result[right] -= (-coeff2 * arr * x[center] /* Left */
                                          + 0.25 * art * x[top] /* Top Left */
                                          - 0.25 * art * x[bottom]); /* Bottom Left */
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

                    /* Fill result(i-1,j) */
                    if (i_r > 1 || !DirBC_Interior) {
                        result[left] -= (-coeff1 * arr * x[center] /* Right */
                                         - 0.25 * art * x[top] /* Top Right */
                                         + 0.25 * art * x[bottom]); /* Bottom Right */
                    }
                    /* Fill result(i+1,j) */
                    if (i_r < grid.numberSmootherCircles() - 1) {
                        result[right] -= (-coeff2 * arr * x[center] /* Left */
                                          + 0.25 * art * x[top] /* Top Left */
                                          - 0.25 * art * x[bottom]); /* Bottom Left */
                    }
                }
                else {
                    /* i_r % 2 == 0 and i_theta % 2 == 0 */
                    /* | o | o | o | */
                    /* |   |   |   | */
                    /* | o | X | o | */
                    /* |   |   |   | */
                    /* | o | o | o | */

                    /* Fill result(i-1,j) */
                    if (i_r > 1 || !DirBC_Interior) {
                        result[left] -= (-coeff1 * arr * x[center] /* Right */
                                         - 0.25 * art * x[top] /* Top Right */
                                         + 0.25 * art * x[bottom]); /* Bottom Right */
                    } /* Fill result(i+1,j) */
                    if (i_r < grid.numberSmootherCircles() - 1) {
                        result[right] -= (-coeff2 * arr * x[center] /* Left */
                                          + 0.25 * art * x[top] /* Top Left */
                                          - 0.25 * art * x[bottom]); /* Bottom Left */
                    }
                }
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
            double h2 = grid.radialSpacing(i_r);
            double k1 = grid.angularSpacing(i_theta - 1);
            double k2 = grid.angularSpacing(i_theta);

            double coeff2 = 0.5 * (k1 + k2) / h2;

            int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
            int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

            int center = grid.index(i_r, i_theta);
            int right  = grid.index(i_r + 1, i_theta);
            int bottom = grid.index(i_r, i_theta_M1);
            int top    = grid.index(i_r, i_theta_P1);

            /* -------------------- */
            /* Inside Section Parts */
            if (node_color == smoother_color) {
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
            /* --------------------- */
            /* Outside Section Parts */
            else if (node_color != smoother_color) {
                /* Fill result(i+1,j) */
                result[right] -= (-coeff2 * arr * x[center] /* Left */
                                  + 0.25 * art * x[top] /* Top Left */
                                  - 0.25 * art * x[bottom]); /* Bottom Left */
            }
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

            const int left   = grid.index(i_r, i_theta_Across);
            const int bottom = grid.index(i_r, i_theta_M1);
            const int center = grid.index(i_r, i_theta);
            const int top    = grid.index(i_r, i_theta_P1);
            const int right  = grid.index(i_r + 1, i_theta);

            /* -------------------- */
            /* Inside Section Parts */
            if (node_color == smoother_color) {
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
            /* --------------------- */
            /* Outside Section Parts */
            else if (node_color != smoother_color) {
                /* Fill result(i+1,j) */
                result[right] -= (-coeff2 * arr * x[center] /* Left */
                                  + 0.25 * art * x[top] /* Top Left */
                                  - 0.25 * art * x[bottom]); /* Bottom Left */
            }
        }
    }
    /* ----------------------------- */
    /* Node next to circular section */
    /* ----------------------------- */
    else if (i_r == grid.numberSmootherCircles()) {
        assert(node_color == SmootherColor::White);

        if (smoother_color == SmootherColor::Black) {
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

            int center = grid.index(i_r, i_theta);
            int left   = grid.index(i_r - 1, i_theta);
            int bottom = grid.index(i_r, i_theta_M1);
            int top    = grid.index(i_r, i_theta_P1);

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
                /* --------------------- */
                /* Outside Section Parts */

                /* Fill result(i-1,j) */
                result[left] -= (-coeff1 * arr * x[center] /* Right */
                                 - 0.25 * art * x[top] /* Top Right */
                                 + 0.25 * art * x[bottom]); /* Bottom Right */
            }
        }
    }
}

static inline void nodeApplyAscOrthoRadialGive(int i_r, int i_theta, const PolarGrid& grid, bool DirBC_Interior,
                                               SmootherColor smoother_color, ConstVector<double>& x,
                                               ConstVector<double>& rhs, Vector<double>& result, double arr, double att,
                                               double art, double detDF, double coeff_beta)
{
    assert(i_r >= grid.numberSmootherCircles() - 1 && i_r < grid.nr());

    SmootherColor node_color = (i_theta & 1) ? SmootherColor::White : SmootherColor::Black;

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

        int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        int center = grid.index(i_r, i_theta);
        int left   = grid.index(i_r - 1, i_theta);
        int right  = grid.index(i_r + 1, i_theta);
        int bottom = grid.index(i_r, i_theta_M1);
        int top    = grid.index(i_r, i_theta_P1);

        /* -------------------- */
        /* Inside Section Parts */
        if (node_color == smoother_color) {
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

                    /* Fill result(i,j) */
                    result[center] -= (-coeff3 * att * x[bottom] /* Bottom */
                                       - coeff4 * att * x[top] /* Top */
                    );
                    /* Fill result(i-1,j) */
                    result[left] -= (-0.25 * art * x[top] /* Top Right */
                                     + 0.25 * art * x[bottom]); /* Bottom Right */
                    /* Fill result(i+1,j) */
                    result[right] -= (+0.25 * art * x[top] /* Top Left */
                                      - 0.25 * art * x[bottom]); /* Bottom Left */
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

                    /* Fill result(i,j) */
                    result[center] -= (-coeff3 * att * x[bottom] /* Bottom */
                                       - coeff4 * att * x[top] /* Top */
                    );
                    /* Fill result(i-1,j) */
                    result[left] -= (-0.25 * art * x[top] /* Top Right */
                                     + 0.25 * art * x[bottom]); /* Bottom Right */
                    /* Fill result(i+1,j) */
                    result[right] -= (+0.25 * art * x[top] /* Top Left */
                                      - 0.25 * art * x[bottom]); /* Bottom Left */
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

                    /* Fill result(i,j) */
                    result[center] -= (-coeff1 * arr * x[left] /* Left */
                                       - coeff2 * arr * x[right] /* Right */
                                       - coeff3 * att * x[bottom] /* Bottom */
                                       - coeff4 * att * x[top] /* Top */
                    );
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
        /* --------------------- */
        /* Outside Section Parts */
        else if (node_color != smoother_color) {
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

                    /* Fill result(i,j-1) */
                    result[bottom] -= (-coeff3 * att * x[center] /* Top */
                                       - 0.25 * art * x[right] /* Top Right */
                                       + 0.25 * art * x[left]); /* Top Left */
                    /* Fill result(i,j+1) */
                    result[top] -= (-coeff4 * att * x[center] /* Bottom */
                                    + 0.25 * art * x[right] /* Bottom Right */
                                    - 0.25 * art * x[left]); /* Bottom Left */
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

                    /* Fill result(i,j-1) */
                    result[bottom] -= (-coeff3 * att * x[center] /* Top */
                                       - 0.25 * art * x[right] /* Top Right */
                                       + 0.25 * art * x[left]); /* Top Left */
                    /* Fill result(i,j+1) */
                    result[top] -= (-coeff4 * att * x[center] /* Bottom */
                                    + 0.25 * art * x[right] /* Bottom Right */
                                    - 0.25 * art * x[left]); /* Bottom Left */
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
    }
    else if (i_r == grid.numberSmootherCircles() - 1) {
        double h1 = grid.radialSpacing(i_r - 1);
        double h2 = grid.radialSpacing(i_r);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff2 = 0.5 * (k1 + k2) / h2;

        int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        int center = grid.index(i_r, i_theta);
        int right  = grid.index(i_r + 1, i_theta);
        int bottom = grid.index(i_r, i_theta_M1);
        int top    = grid.index(i_r, i_theta_P1);

        /* -------------------- */
        /* Inside Section Parts */
        if (node_color == smoother_color) {
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
        /* --------------------- */
        /* Outside Section Parts */
        else if (node_color != smoother_color) {
            /* Nothing to be done here */
        }
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

        int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        int center = grid.index(i_r, i_theta);
        int left   = grid.index(i_r - 1, i_theta);
        int right  = grid.index(i_r + 1, i_theta);
        int bottom = grid.index(i_r, i_theta_M1);
        int top    = grid.index(i_r, i_theta_P1);

        /* -------------------- */
        /* Inside Section Parts */
        if (node_color == smoother_color) {
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
        /* --------------------- */
        /* Outside Section Parts */
        else if (node_color != smoother_color) {
            /* Dont give to bottom and up when this case occurs! */
            /* i_theta % 2 == 1 and i_r % 2 == 0 */
            /* | o | x | o || x   o   x   o  */
            /* |   |   |   || -------------- */
            /* | o | o | o || O   o   o   o  */
            /* |   |   |   || -------------- */
            /* | o | x | o || x   o   x   o  */
            if (i_r & 1 || !(i_theta & 1)) {
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
    else if (i_r == grid.nr() - 2) {
        assert(i_r & 1);

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

        int center = grid.index(i_r, i_theta);
        int left   = grid.index(i_r - 1, i_theta);
        int right  = grid.index(i_r + 1, i_theta);
        int bottom = grid.index(i_r, i_theta_M1);
        int top    = grid.index(i_r, i_theta_P1);

        /* -------------------- */
        /* Inside Section Parts */
        if (node_color == smoother_color) {
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
                result[center] -= /* Right: Symmetry shift! */
                    -coeff2 * arr * rhs[right];
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
        /* --------------------- */
        /* Outside Section Parts */
        else if (node_color != smoother_color) {
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
    else if (i_r == grid.nr() - 1) {
        assert(!(i_r & 1));

        double h1 = grid.radialSpacing(i_r - 1);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;

        int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        int center = grid.index(i_r, i_theta);
        int left   = grid.index(i_r - 1, i_theta);
        int bottom = grid.index(i_r, i_theta_M1);
        int top    = grid.index(i_r, i_theta_P1);

        /* -------------------- */
        /* Inside Section Parts */
        if (node_color == smoother_color) {
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
        /* --------------------- */
        /* Outside Section Parts */
        else if (node_color != smoother_color) {
            /* Nothing to be done here */
        }
    }
}

void ExtrapolatedSmootherGive::applyAscOrthoCircleSection(const int i_r, const SmootherColor smoother_color,
                                                          ConstVector<double> x, ConstVector<double> rhs,
                                                          Vector<double> temp)
{
    assert(i_r >= 0 && i_r < grid_.numberSmootherCircles() + 1);

    const double r = grid_.radius(i_r);

    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        const double theta = grid_.theta(i_theta);
        const int index    = grid_.index(i_r, i_theta);

        double coeff_beta, arr, att, art, detDF;
        level_cache_.obtainValues(i_r, i_theta, index, r, theta, coeff_beta, arr, att, art, detDF);

        // Apply Asc Ortho at the current node
        nodeApplyAscOrthoCircleGive(i_r, i_theta, grid_, DirBC_Interior_, smoother_color, x, rhs, temp, arr, att, art,
                                    detDF, coeff_beta);
    }
}

void ExtrapolatedSmootherGive::applyAscOrthoRadialSection(const int i_theta, const SmootherColor smoother_color,
                                                          ConstVector<double> x, ConstVector<double> rhs,
                                                          Vector<double> temp)
{
    const double theta = grid_.theta(i_theta);

    /* !!! i_r = grid_.numberSmootherCircles()-1 !!! */
    for (int i_r = grid_.numberSmootherCircles() - 1; i_r < grid_.nr(); i_r++) {
        const double r  = grid_.radius(i_r);
        const int index = grid_.index(i_r, i_theta);

        double coeff_beta, arr, att, art, detDF;
        level_cache_.obtainValues(i_r, i_theta, index, r, theta, coeff_beta, arr, att, art, detDF);

        // Apply Asc Ortho at the current node
        nodeApplyAscOrthoRadialGive(i_r, i_theta, grid_, DirBC_Interior_, smoother_color, x, rhs, temp, arr, att, art,
                                    detDF, coeff_beta);
    }
}

void ExtrapolatedSmootherGive::solveCircleSection(const int i_r, Vector<double> x, Vector<double> temp,
                                                  Vector<double> solver_storage_1, Vector<double> solver_storage_2)
{
    const int start = grid_.index(i_r, 0);
    const int end   = start + grid_.ntheta();
    if (i_r == 0) {
#ifdef GMGPOLAR_USE_MUMPS
        inner_boundary_mumps_solver_.job    = JOB_COMPUTE_SOLUTION;
        inner_boundary_mumps_solver_.nrhs   = 1; // single rhs vector
        inner_boundary_mumps_solver_.nz_rhs = grid_.ntheta(); // non-zeros in rhs
        inner_boundary_mumps_solver_.rhs    = temp.data() + start;
        inner_boundary_mumps_solver_.lrhs   = grid_.ntheta(); // leading dimension of rhs
        dmumps_c(&inner_boundary_mumps_solver_);
        if (inner_boundary_mumps_solver_.info[0] != 0) {
            std::cerr << "Error solving the system: " << inner_boundary_mumps_solver_.info[0] << std::endl;
        }
#else
        inner_boundary_lu_solver_.solveInPlace(temp.data());
#endif
    }
    else {
        if (i_r & 1) {
            circle_tridiagonal_solver_[i_r / 2].solveInPlace(temp.data() + start, solver_storage_1.data(),
                                                             solver_storage_2.data());
        }
        else {
            circle_diagonal_solver_[i_r / 2].solveInPlace(temp.data() + start);
        }
    }
    // Move updated values to x
    Kokkos::deep_copy(Kokkos::subview(x, Kokkos::make_pair(start, end)),
                      Kokkos::subview(temp, Kokkos::make_pair(start, end)));
}

void ExtrapolatedSmootherGive::solveRadialSection(const int i_theta, Vector<double> x, Vector<double> temp,
                                                  Vector<double> solver_storage)
{
    const int start = grid_.index(grid_.numberSmootherCircles(), i_theta);
    const int end   = start + grid_.lengthSmootherRadial();
    if (i_theta & 1) {
        radial_tridiagonal_solver_[i_theta / 2].solveInPlace(temp.data() + start, solver_storage.data());
    }
    else {
        radial_diagonal_solver_[i_theta / 2].solveInPlace(temp.data() + start);
    }
    // Move updated values to x
    Kokkos::deep_copy(Kokkos::subview(x, Kokkos::make_pair(start, end)),
                      Kokkos::subview(temp, Kokkos::make_pair(start, end)));
}

// Quick overview:
// Step 1: temp = rhs - Asc_ortho(x)
// Step 2: Solve in place
// Step 3: Update x

/* ------------------ */
/* Sequential Version */
/* ------------------ */

void ExtrapolatedSmootherGive::extrapolatedSmoothingSequential(Vector<double> x, ConstVector<double> rhs,
                                                               Vector<double> temp)
{
    assert(x.size() == rhs.size());
    assert(temp.size() == rhs.size());

    for (int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
        for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
            const int index = grid_.index(i_r, i_theta);
            temp[index]     = (i_r & 1 || i_theta & 1) ? rhs[index] : x[index];
        }
    }

    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
            const int index = grid_.index(i_r, i_theta);
            temp[index]     = (i_r & 1 || i_theta & 1) ? rhs[index] : x[index];
        }
    }

    /* Single-threaded execution */
    Vector<double> circle_solver_storage_1("circle_solver_storage1", grid_.ntheta());
    Vector<double> circle_solver_storage_2("circle_solver_storage1", grid_.ntheta());
    Vector<double> radial_solver_storage("radial_solver_storage", grid_.lengthSmootherRadial());

    /* ---------------------------- */
    /* ------ CIRCLE SECTION ------ */

    /* The outer most circle next to the radial section is defined to be black. */
    /* Priority: Black -> White. */
    for (int i_r = 0; i_r < grid_.numberSmootherCircles() + 1; i_r++) {
        applyAscOrthoCircleSection(i_r, SmootherColor::Black, x, rhs, temp);
    }
    const int start_black_circles = (grid_.numberSmootherCircles() % 2 == 0) ? 1 : 0;
    for (int i_r = start_black_circles; i_r < grid_.numberSmootherCircles(); i_r += 2) {
        solveCircleSection(i_r, x, temp, circle_solver_storage_1, circle_solver_storage_2);
    }
    for (int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
        applyAscOrthoCircleSection(i_r, SmootherColor::White, x, rhs, temp);
    }
    const int start_white_circles = (grid_.numberSmootherCircles() % 2 == 0) ? 0 : 1;
    for (int i_r = start_white_circles; i_r < grid_.numberSmootherCircles(); i_r += 2) {
        solveCircleSection(i_r, x, temp, circle_solver_storage_1, circle_solver_storage_2);
    }
    /* ---------------------------- */
    /* ------ RADIAL SECTION ------ */
    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        applyAscOrthoRadialSection(i_theta, SmootherColor::Black, x, rhs, temp);
    }
    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta += 2) {
        solveRadialSection(i_theta, x, temp, radial_solver_storage);
    }
    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        applyAscOrthoRadialSection(i_theta, SmootherColor::White, x, rhs, temp);
    }
    for (int i_theta = 1; i_theta < grid_.ntheta(); i_theta += 2) {
        solveRadialSection(i_theta, x, temp, radial_solver_storage);
    }
}

/* ------------------------------------ */
/* Parallelization Version 1: For Loops */
/* ------------------------------------ */

// clang-format off
void ExtrapolatedSmootherGive::extrapolatedSmoothingForLoop(Vector<double> x, ConstVector<double> rhs,
                                                            Vector<double> temp)
{
    assert(x.size() == rhs.size());
    assert(temp.size() == rhs.size());

    if (num_omp_threads_ == 1) {
        extrapolatedSmoothingSequential(x, rhs, temp);
    }
    else {
        #pragma omp parallel num_threads(num_omp_threads_)
        {
            #pragma omp for nowait
            for (int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
                for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
                    const int index = grid_.index(i_r, i_theta);
                    temp[index]     = (i_r & 1 || i_theta & 1) ? rhs[index] : x[index];
                }
            }
            #pragma omp for
            for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
                for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
                    const int index = grid_.index(i_r, i_theta);
                    temp[index]     = (i_r & 1 || i_theta & 1) ? rhs[index] : x[index];
                }
            }
        }

        /* Multi-threaded execution */
        const int num_circle_tasks = grid_.numberSmootherCircles();
        const int num_radial_tasks = grid_.ntheta();

        #pragma omp parallel
        {
            Vector<double> circle_solver_storage_1("circle_solver_storage_1",grid_.ntheta());
            Vector<double> circle_solver_storage_2("circle_solver_storage_2",grid_.ntheta());
            Vector<double> radial_solver_storage("radial_solver_storage",grid_.lengthSmootherRadial());

            /* ---------------------------- */
            /* ------ CIRCLE SECTION ------ */
            /* ---------------------------- */

            /* ---------------------------- */
            /* Asc ortho Black Circle Tasks */

            /* Inside Black Section */
            #pragma omp for
            for (int circle_task = 0; circle_task < num_circle_tasks; circle_task += 2) {
                int i_r = num_circle_tasks - circle_task - 1;
                applyAscOrthoCircleSection(i_r, SmootherColor::Black, x, rhs, temp);
            }
            /* Outside Black Section (Part 1)*/
            #pragma omp for
            for (int circle_task = -1; circle_task < num_circle_tasks; circle_task += 4) {
                int i_r = num_circle_tasks - circle_task - 1;
                applyAscOrthoCircleSection(i_r, SmootherColor::Black, x, rhs, temp);
            }
            /* Outside Black Section (Part 2)*/
            #pragma omp for
            for (int circle_task = 1; circle_task < num_circle_tasks; circle_task += 4) {
                int i_r = num_circle_tasks - circle_task - 1;
                applyAscOrthoCircleSection(i_r, SmootherColor::Black, x, rhs, temp);
            }

            /* Black Circle Smoother */
            #pragma omp for
            for (int circle_task = 0; circle_task < num_circle_tasks; circle_task += 2) {
                int i_r = num_circle_tasks - circle_task - 1;
                solveCircleSection(i_r, x, temp, circle_solver_storage_1, circle_solver_storage_2);
            }

            /* ---------------------------- */
            /* Asc ortho White Circle Tasks */
            /* Inside White Section */
            #pragma omp for nowait
            for (int circle_task = 1; circle_task < num_circle_tasks; circle_task += 2) {
                int i_r = num_circle_tasks - circle_task - 1;
                applyAscOrthoCircleSection(i_r, SmootherColor::White, x, rhs, temp);
            }
            /* ---------------------------- */
            /* Asc ortho Black Radial Tasks */
            /* Inside Black Section */
            #pragma omp for
            for (int radial_task = 0; radial_task < num_radial_tasks; radial_task += 2) {
                int i_theta = radial_task;
                applyAscOrthoRadialSection(i_theta, SmootherColor::Black, x, rhs, temp);
            }

            /* ---------------------------- */
            /* Asc ortho White Circle Tasks */
            /* Outside White Section (Part 1)*/
            #pragma omp for nowait
            for (int circle_task = 0; circle_task < num_circle_tasks; circle_task += 4) {
                int i_r = num_circle_tasks - circle_task - 1;
                applyAscOrthoCircleSection(i_r, SmootherColor::White, x, rhs, temp);
            }
            /* ---------------------------- */
            /* Asc ortho Black Radial Tasks */
            /* Outside Black Section (Part 1) */
            #pragma omp for
            for (int radial_task = 1; radial_task < num_radial_tasks; radial_task += 4) {
                int i_theta = radial_task;
                applyAscOrthoRadialSection(i_theta, SmootherColor::Black, x, rhs, temp);
            }

            /* ---------------------------- */
            /* Asc ortho White Circle Tasks */
            /* Outside White Section (Part 2)*/
            #pragma omp for nowait
            for (int circle_task = 2; circle_task < num_circle_tasks; circle_task += 4) {
                int i_r = num_circle_tasks - circle_task - 1;
                applyAscOrthoCircleSection(i_r, SmootherColor::White, x, rhs, temp);
            }
            /* ---------------------------- */
            /* Asc ortho Black Radial Tasks */
            /* Outside Black Section (Part 2) */
            #pragma omp for
            for (int radial_task = 3; radial_task < num_radial_tasks; radial_task += 4) {
                int i_theta = radial_task;
                applyAscOrthoRadialSection(i_theta, SmootherColor::Black, x, rhs, temp);
            }

            /* White Circle Smoother */
            #pragma omp for nowait
            for (int circle_task = 1; circle_task < num_circle_tasks; circle_task += 2) {
                int i_r = num_circle_tasks - circle_task - 1;
                solveCircleSection(i_r, x, temp, circle_solver_storage_1, circle_solver_storage_2);
            }
            /* Black Radial Smoother */
            #pragma omp for
            for (int radial_task = 0; radial_task < num_radial_tasks; radial_task += 2) {
                int i_theta = radial_task;
                solveRadialSection(i_theta, x, temp, radial_solver_storage);
            }

            /* ---------------------------- */
            /* Asc ortho White Circle Tasks */

            /* Inside White Section */
            #pragma omp for
            for (int radial_task = 1; radial_task < num_radial_tasks; radial_task += 2) {
                int i_theta = radial_task;
                applyAscOrthoRadialSection(i_theta, SmootherColor::White, x, rhs, temp);
            }
            /* Outside White Section (Part 1) */
            #pragma omp for
            for (int radial_task = 0; radial_task < num_radial_tasks; radial_task += 4) {
                int i_theta = radial_task;
                applyAscOrthoRadialSection(i_theta, SmootherColor::White, x, rhs, temp);
            }
            /* Outside White Section (Part 2) */
            #pragma omp for
            for (int radial_task = 2; radial_task < num_radial_tasks; radial_task += 4) {
                int i_theta = radial_task;
                applyAscOrthoRadialSection(i_theta, SmootherColor::White, x, rhs, temp);
            }

            /* White Radial Smoother */
            #pragma omp for
            for (int radial_task = 1; radial_task < num_radial_tasks; radial_task += 2) {
                int i_theta = radial_task;
                solveRadialSection(i_theta, x, temp, radial_solver_storage);
            }
        }
    }
}
