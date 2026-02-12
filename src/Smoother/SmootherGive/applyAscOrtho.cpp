#include "../../../include/Smoother/SmootherGive/smootherGive.h"

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
        /* --------------------- */
        /* Outside Section Parts */
        else if (node_color != smoother_color) {
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
                /* Nothing to be done here */
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
                /* Fill result(i,j) */
                result[center] -= (
                    /* - coeff1 * arr * x[left] // Left: Not in Asc_ortho */
                    -coeff2 * arr * x[right]); /* Right */
                /* Fill result(i,j-1) */
                result[bottom] -= (-0.25 * art * x[right]); /* Top Right */
                /*  + 0.25 * art * x[left]; // Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
                /* Fill result(i,j+1) */
                result[top] -= (+0.25 * art * x[right]); /* Bottom Right */
                /*  - 0.25 * art * x[left]; // Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
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
            double k1 = grid.angularSpacing(i_theta - 1);
            double k2 = grid.angularSpacing(i_theta);

            double coeff1 = 0.5 * (k1 + k2) / h1;

            int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
            int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

            int center = grid.index(i_r, i_theta);
            int left   = grid.index(i_r - 1, i_theta);
            int bottom = grid.index(i_r, i_theta_M1);
            int top    = grid.index(i_r, i_theta_P1);

            /* --------------------- */
            /* Outside Section Parts */
            /* Fill result(i-1,j) */
            result[left] -= (-coeff1 * arr * x[center] /* Right */
                             - 0.25 * art * x[top] /* Top Right */
                             + 0.25 * art * x[bottom]); /* Bottom Right */
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
            /* Fill result(i,j) */
            result[center] -= (-coeff3 * att * x[bottom] /* Bottom */
                               - coeff4 * att * x[top] /* Top */
            ); /* Fill result(i-1,j) */
            result[left] -= (-0.25 * art * x[top] /* Top Right */
                             + 0.25 * art * x[bottom]); /* Bottom Right */
            /* Fill result(i+1,j) */
            result[right] -= (+0.25 * art * x[top] /* Top Left */
                              - 0.25 * art * x[bottom]); /* Bottom Left */
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

    /* -------------------------------- */
    /* Node left next to radial section */
    /* -------------------------------- */
    else if (i_r == grid.numberSmootherCircles() - 1) {
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
        int right  = grid.index(i_r + 1, i_theta);
        int bottom = grid.index(i_r, i_theta_M1);
        int top    = grid.index(i_r, i_theta_P1);

        /* -------------------- */
        /* Inside Section Parts */
        if (node_color == smoother_color) {
            /* Fill result(i+1,j) */
            result[right] -= (-coeff2 * arr * x[center] /* Left */
                              + 0.25 * art * x[top] /* Top Left */
                              - 0.25 * art * x[bottom]); /* Bottom Left */
        }
        /* --------------------- */
        /* Outside Section Parts */
        else if (node_color != smoother_color) {
            /* Nothing to be done here */
        }
    }
    /* --------------------------- */
    /* Node next to circle section */
    /* --------------------------- */
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
            /* Fill result(i,j) */
            result[center] -= (-coeff1 * arr * x[left] /* Left */
                               - coeff3 * att * x[bottom] /* Bottom */
                               - coeff4 * att * x[top] /* Top */
            );
            /* Fill result(i+1,j) */
            result[right] -= (+0.25 * art * x[top] /* Top Left */
                              - 0.25 * art * x[bottom]); /* Bottom Left */
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
    /* ------------------------------- */
    /* Node next to dirichlet boundary */
    /* ------------------------------- */
    else if (i_r == grid.nr() - 2) {
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
            /* Fill result(i,j) */
            result[center] -= (-coeff3 * att * x[bottom] /* Bottom */
                               - coeff4 * att * x[top] /* Top */
            );
            /* Fill result(i-1,j) */
            result[left] -= (-0.25 * art * x[top] /* Top Right */
                             + 0.25 * art * x[bottom]); /* Bottom Right */

            // "Right" is part of the radial Asc smoother matrices,
            // but is shifted over to the rhs to make the radial Asc smoother matrices symmetric.
            // Note that the circle Asc smoother matrices are symmetric by default.
            // Note that rhs[right] contains the correct boundary value of u_D.
            result[center] -= -coeff2 * arr * rhs[right]; /* Right: Symmetry shift! */
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
    /* ------------------------------ */
    /* Node on the dirichlet boundary */
    /* ------------------------------ */
    else if (i_r == grid.nr() - 1) {
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
            /* Fill result(i-1,j) */
            result[left] -= (-0.25 * art * x[top] /* Top Right */
                             + 0.25 * art * x[bottom] /* Bottom Right */
            );
            // "Right" is part of the radial Asc smoother matrices,
            // but is shifted over to the rhs to make the radial Asc smoother matrices symmetric.
            // Note that the circle Asc smoother matrices are symmetric by default.
            // Note that rhs[center] contains the correct boundary value of u_D.
            result[left] -= (-coeff1 * arr * rhs[center] /* Right */
            );
        }
        /* --------------------- */
        /* Outside Section Parts */
        else if (node_color != smoother_color) {
            /* Nothing to be done here */
        }
    }
}

void SmootherGive::applyAscOrthoCircleSection(const int i_r, const SmootherColor smoother_color, ConstVector<double> x,
                                              ConstVector<double> rhs, Vector<double> temp)
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

void SmootherGive::applyAscOrthoRadialSection(const int i_theta, const SmootherColor smoother_color,
                                              ConstVector<double> x, ConstVector<double> rhs, Vector<double> temp)
{
    const double theta = grid_.theta(i_theta);

    /* We need to obtain left contributions from the circular section for AscOrtho. */
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