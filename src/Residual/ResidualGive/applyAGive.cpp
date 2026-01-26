#include "../../../include/Residual/ResidualGive/residualGive.h"

inline void node_apply_a_give(int i_r, int i_theta, double r, double theta, const PolarGrid& grid,
                              const LevelCache& level_cache, bool DirBC_Interior, Vector<double>& result,
                              ConstVector<double>& x)
{
    /* ---------------------------------------- */
    /* Compute or retrieve stencil coefficients */
    /* ---------------------------------------- */
    const int center = grid.index(i_r, i_theta);
    double coeff_beta, arr, att, art, detDF;
    level_cache.obtainValues(i_r, i_theta, center, r, theta, coeff_beta, arr, att, art, detDF);

    /* -------------------- */
    /* Node in the interior */
    /* -------------------- */
    if (i_r > 1 && i_r < grid.nr() - 2) {
        double h1 = grid.radialSpacing(i_r - 1);
        double h2 = grid.radialSpacing(i_r);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;

        const int left   = grid.index(i_r - 1, i_theta);
        const int right  = grid.index(i_r + 1, i_theta);
        const int bottom = grid.index(i_r, i_theta - 1);
        const int top    = grid.index(i_r, i_theta + 1);

        /* Fill result(i,j) */
        result[center] -= (0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF) * x[center] /* beta_{i,j} */
                           - coeff1 * arr * x[left] /* Left */
                           - coeff2 * arr * x[right] /* Right */
                           - coeff3 * att * x[bottom] /* Bottom */
                           - coeff4 * att * x[top] /* Top */
                           + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) *
                                 x[center]) /* Center: (Left, Right, Bottom, Top) */;
        /* Fill result(i-1,j) */
        result[left] -= (-coeff1 * arr * x[center] /* Right */
                         + coeff1 * arr * x[left] /* Center: (Right) */
                         - 0.25 * art * x[top] /* Top Right */
                         + 0.25 * art * x[bottom]); /* Bottom Right */
        /* Fill result(i+1,j) */
        result[right] -= (-coeff2 * arr * x[center] /* Left */
                          + coeff2 * arr * x[right] /* Center: (Left) */
                          + 0.25 * art * x[top] /* Top Left */
                          - 0.25 * art * x[bottom]); /* Bottom Left */
        /* Fill result(i,j-1) */
        result[bottom] -= (-coeff3 * att * x[center] /* Top */
                           + coeff3 * att * x[bottom] /* Center: (Top) */
                           - 0.25 * art * x[right] /* Top Right */
                           + 0.25 * art * x[left]); /* Top Left */
        /* Fill result(i,j+1) */
        result[top] -= (-coeff4 * att * x[center] /* Bottom */
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
            result[center] -= x[center];

            /* Give value to the interior nodes! */
            double h2 = grid.radialSpacing(i_r);
            double k1 = grid.angularSpacing(i_theta - 1);
            double k2 = grid.angularSpacing(i_theta);

            double coeff2 = 0.5 * (k1 + k2) / h2;

            const int right  = grid.index(i_r + 1, i_theta);
            const int bottom = grid.index(i_r, i_theta - 1);
            const int top    = grid.index(i_r, i_theta + 1);

            /* Fill result(i+1,j) */
            result[right] -= (-coeff2 * arr * x[center] /* Left */
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
            double h1 = 2.0 * grid.radius(0);
            double h2 = grid.radialSpacing(i_r);
            double k1 = grid.angularSpacing(i_theta - 1);
            double k2 = grid.angularSpacing(i_theta);

            double coeff1 = 0.5 * (k1 + k2) / h1;
            double coeff2 = 0.5 * (k1 + k2) / h2;
            double coeff3 = 0.5 * (h1 + h2) / k1;
            double coeff4 = 0.5 * (h1 + h2) / k2;

            const int left   = grid.index(i_r, i_theta + (grid.ntheta() / 2));
            const int right  = grid.index(i_r + 1, i_theta);
            const int bottom = grid.index(i_r, i_theta - 1);
            const int top    = grid.index(i_r, i_theta + 1);

            /* Fill result(i,j) */
            result[center] -= (0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF) * x[center] /* beta_{i,j} */
                               - coeff1 * arr * x[left] /* Left */
                               - coeff2 * arr * x[right] /* Right */
                               - coeff3 * att * x[bottom] /* Bottom */
                               - coeff4 * att * x[top] /* Top */
                               + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) *
                                     x[center]); /* Center: (Left, Right, Bottom, Top) */
            /* Fill result(i-1,j) */
            /* From view the view of the across origin node, the directions are roatated by 180 degrees in the stencil! */
            result[left] -= (-coeff1 * arr * x[center] /* Right -> Left */
                             + coeff1 * arr * x[left]); /* Center: (Right) -> Center: (Left)*/
            /* + 0.25 * art * x[top]; // Top Right -> Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
            /* - 0.25 * art * x[bottom]; // Bottom Right -> Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
            /* Fill result(i+1,j) */
            result[right] -= (-coeff2 * arr * x[center] /* Left */
                              + coeff2 * arr * x[right] /* Center: (Left) */
                              + 0.25 * art * x[top] /* Top Left */
                              - 0.25 * art * x[bottom]); /* Bottom Left */
            /* Fill result(i,j-1) */
            result[bottom] -= (-coeff3 * att * x[center] /* Top */
                               + coeff3 * att * x[bottom] /* Center: (Top) */
                               - 0.25 * art * x[right]); /* Top Right */
            /* + 0.25 * art * x[left]; // Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
            /* Fill result(i,j+1) */
            result[top] -= (-coeff4 * att * x[center] /* Bottom */
                            + coeff4 * att * x[top] /* Center: (Bottom) */
                            + 0.25 * art * x[right]); /* Bottom Right */
            /* - 0.25 * art * x[left]; // Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
        }
    }
    /* ------------------------------- */
    /* Node next to the inner boundary */
    /* ------------------------------- */
    else if (i_r == 1) {
        double h1 = grid.radialSpacing(i_r - 1);
        double h2 = grid.radialSpacing(i_r);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;

        const int left   = grid.index(i_r - 1, i_theta);
        const int right  = grid.index(i_r + 1, i_theta);
        const int bottom = grid.index(i_r, i_theta - 1);
        const int top    = grid.index(i_r, i_theta + 1);

        /* Fill result(i,j) */
        result[center] -= (0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF) * x[center] /* beta_{i,j} */
                           - coeff1 * arr * x[left] /* Left */
                           - coeff2 * arr * x[right] /* Right */
                           - coeff3 * att * x[bottom] /* Bottom */
                           - coeff4 * att * x[top] /* Top */
                           + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) *
                                 x[center]); /* Center: (Left, Right, Bottom, Top) */
        /* Fill result(i-1,j) */
        if (!DirBC_Interior) { /* Don't give to the inner dirichlet boundary! */
            result[left] -= (-coeff1 * arr * x[center] /* Right */
                             + coeff1 * arr * x[left] /* Center: (Right) */
                             - 0.25 * art * x[top] /* Top Right */
                             + 0.25 * art * x[bottom]); /* Bottom Right */
        }
        /* Fill result(i+1,j) */
        result[right] -= (-coeff2 * arr * x[center] /* Left */
                          + coeff2 * arr * x[right] /* Center: (Left) */
                          + 0.25 * art * x[top] /* Top Left */
                          - 0.25 * art * x[bottom]); /* Bottom Left */
        /* Fill result(i,j-1) */
        result[bottom] -= (-coeff3 * att * x[center] /* Top */
                           + coeff3 * att * x[bottom] /* Center: (Top) */
                           - 0.25 * art * x[right] /* Top Right */
                           + 0.25 * art * x[left]); /* Top Left */
        /* Fill result(i,j+1) */
        result[top] -= (-coeff4 * att * x[center] /* Bottom */
                        + coeff4 * att * x[top] /* Center: (Bottom) */
                        + 0.25 * art * x[right] /* Bottom Right */
                        - 0.25 * art * x[left]); /* Bottom Left */
    }
    /* ------------------------------- */
    /* Node next to the outer boundary */
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

        const int left   = grid.index(i_r - 1, i_theta);
        const int right  = grid.index(i_r + 1, i_theta);
        const int bottom = grid.index(i_r, i_theta - 1);
        const int top    = grid.index(i_r, i_theta + 1);

        /* Fill result(i,j) */
        result[center] -= (0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF) * x[center] /* beta_{i,j} */
                           - coeff1 * arr * x[left] /* Left */
                           - coeff2 * arr * x[right] /* Right */
                           - coeff3 * att * x[bottom] /* Bottom */
                           - coeff4 * att * x[top] /* Top */
                           + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) *
                                 x[center]); /* Center: (Left, Right, Bottom, Top) */
        /* Fill result(i-1,j) */
        result[left] -= (-coeff1 * arr * x[center] /* Right */
                         + coeff1 * arr * x[left] /* Center: (Right) */
                         - 0.25 * art * x[top] /* Top Right */
                         + 0.25 * art * x[bottom]); /* Bottom Right */
        /* Don't give to the outer dirichlet boundary! */
        /* Fill result(i+1,j) */
        /* result[grid.index(right] -= ( */
        /*     - coeff2 * arr * x[center] // Left */
        /*     + coeff2 * arr * x[right] // Center: (Left) */
        /*     + 0.25 * art * x[top] // Top Left */
        /*     - 0.25 * art * x[bottom]); // Bottom Left */
        /* Fill result(i,j-1) */
        result[bottom] -= (-coeff3 * att * x[center] /* Top */
                           + coeff3 * att * x[bottom] /* Center: (Top) */
                           - 0.25 * art * x[right] /* Top Right */
                           + 0.25 * art * x[left]); /* Top Left */
        /* Fill result(i,j+1) */
        result[top] -= (-coeff4 * att * x[center] /* Bottom */
                        + coeff4 * att * x[top] /* Center: (Bottom) */
                        + 0.25 * art * x[right] /* Bottom Right */
                        - 0.25 * art * x[left]); /* Bottom Left */
    }
    /* ----------------------------- */
    /* Node on to the outer boundary */
    /* ----------------------------- */
    else if (i_r == grid.nr() - 1) {
        /* Fill result of (i,j) */
        result[center] -= x[center];

        /* Give value to the interior nodes! */
        double h1 = grid.radialSpacing(i_r - 1);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;

        const int left   = grid.index(i_r - 1, i_theta);
        const int bottom = grid.index(i_r, i_theta - 1);
        const int top    = grid.index(i_r, i_theta + 1);

        /* Fill result(i-1,j) */
        result[left] -= (-coeff1 * arr * x[center] /* Right */
                         + coeff1 * arr * x[left] /* Center: (Right) */
                         - 0.25 * art * x[top] /* Top Right */
                         + 0.25 * art * x[bottom]); /* Bottom Right */
    }
}

void ResidualGive::applyCircleSection(const int i_r, Vector<double> result, ConstVector<double> x) const
{
    const double r = grid_.radius(i_r);
    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        const double theta = grid_.theta(i_theta);
        node_apply_a_give(i_r, i_theta, r, theta, grid_, level_cache_, DirBC_Interior_, result, x);
    }
}

void ResidualGive::applyRadialSection(const int i_theta, Vector<double> result, ConstVector<double> x) const
{
    const double theta = grid_.theta(i_theta);
    for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
        const double r = grid_.radius(i_r);
        node_apply_a_give(i_r, i_theta, r, theta, grid_, level_cache_, DirBC_Interior_, result, x);
    }
}
