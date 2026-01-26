#include "../../../include/Residual/ResidualGive/residualGive.h"

#include "../../../include/common/geometry_helper.h"

inline void node_apply_a_give(int i_r, int i_theta, double r, double theta, const PolarGrid& grid,
                              const LevelCache& level_cache, bool DirBC_Interior, Vector<double>& result,
                              ConstVector<double>& x)
{
    const int global_index = grid.index(i_r, i_theta);
    double coeff_beta, arr, att, art, detDF;
    level_cache.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art, detDF);
    /* -------------------- */
    /* Node in the interior */
    /* -------------------- */
    if (i_r > 1 && i_r < grid.nr() - 2) {
        double h1     = grid.radialSpacing(i_r - 1);
        double h2     = grid.radialSpacing(i_r);
        double k1     = grid.angularSpacing(i_theta - 1);
        double k2     = grid.angularSpacing(i_theta);
        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;
        /* Fill result(i,j) */
        result[grid.index(i_r, i_theta)] -=
            (0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF) * x[grid.index(i_r, i_theta)] /* beta_{i,j} */
             - coeff1 * arr * x[grid.index(i_r - 1, i_theta)] /* Left */
             - coeff2 * arr * x[grid.index(i_r + 1, i_theta)] /* Right */
             - coeff3 * att * x[grid.index(i_r, i_theta - 1)] /* Bottom */
             - coeff4 * att * x[grid.index(i_r, i_theta + 1)] /* Top */ /* Center: (Left, Right, Bottom, Top) */
             + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) * x[grid.index(i_r, i_theta)]);
        /* Fill result(i-1,j) */
        result[grid.index(i_r - 1, i_theta)] -= (-coeff1 * arr * x[grid.index(i_r, i_theta)] /* Right */
                                                 + coeff1 * arr * x[grid.index(i_r - 1, i_theta)] /* Center: (Right) */
                                                 - 0.25 * art * x[grid.index(i_r, i_theta + 1)] /* Top Right */
                                                 + 0.25 * art * x[grid.index(i_r, i_theta - 1)]); /* Bottom Right */
        /* Fill result(i+1,j) */
        result[grid.index(i_r + 1, i_theta)] -= (-coeff2 * arr * x[grid.index(i_r, i_theta)] /* Left */
                                                 + coeff2 * arr * x[grid.index(i_r + 1, i_theta)] /* Center: (Left) */
                                                 + 0.25 * art * x[grid.index(i_r, i_theta + 1)] /* Top Left */
                                                 - 0.25 * art * x[grid.index(i_r, i_theta - 1)]); /* Bottom Left */
        /* Fill result(i,j-1) */
        result[grid.index(i_r, i_theta - 1)] -= (-coeff3 * att * x[grid.index(i_r, i_theta)] /* Top */
                                                 + coeff3 * att * x[grid.index(i_r, i_theta - 1)] /* Center: (Top) */
                                                 - 0.25 * art * x[grid.index(i_r + 1, i_theta)] /* Top Right */
                                                 + 0.25 * art * x[grid.index(i_r - 1, i_theta)]); /* Top Left */
        /* Fill result(i,j+1) */
        result[grid.index(i_r, i_theta + 1)] -= (-coeff4 * att * x[grid.index(i_r, i_theta)] /* Bottom */
                                                 + coeff4 * att * x[grid.index(i_r, i_theta + 1)] /* Center: (Bottom) */
                                                 + 0.25 * art * x[grid.index(i_r + 1, i_theta)] /* Bottom Right */
                                                 - 0.25 * art * x[grid.index(i_r - 1, i_theta)]); /* Bottom Left */
        /* -------------------------- */
        /* Node on the inner boundary */
        /* -------------------------- */
    }
    else if (i_r == 0) {
        /* ------------------------------------------------ */
        /* Case 1: Dirichlet boundary on the inner boundary */
        /* ------------------------------------------------ */
        if (DirBC_Interior) {
            /* Fill result(i,j) */
            result[grid.index(i_r, i_theta)] -= x[grid.index(i_r, i_theta)];
            /* Give value to the interior nodes! */
            double h2     = grid.radialSpacing(i_r);
            double k1     = grid.angularSpacing(i_theta - 1);
            double k2     = grid.angularSpacing(i_theta);
            double coeff2 = 0.5 * (k1 + k2) / h2;
            /* Fill result(i+1,j) */
            result[grid.index(i_r + 1, i_theta)] -=
                (-coeff2 * arr * x[grid.index(i_r, i_theta)] /* Left */
                 + coeff2 * arr * x[grid.index(i_r + 1, i_theta)] /* Center: (Left) */
                 + 0.25 * art * x[grid.index(i_r, i_theta + 1)] /* Top Left */
                 - 0.25 * art * x[grid.index(i_r, i_theta - 1)]); /* Bottom Left */
        }
        else {
            /* ------------------------------------------------------------- */
            /* Case 2: Across origin discretization on the interior boundary */
            /* ------------------------------------------------------------- */
            /* h1 gets replaced with 2 * R0. */
            /* (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()/2)). */
            /* Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil. */
            double h1     = 2.0 * grid.radius(0);
            double h2     = grid.radialSpacing(i_r);
            double k1     = grid.angularSpacing(i_theta - 1);
            double k2     = grid.angularSpacing(i_theta);
            double coeff1 = 0.5 * (k1 + k2) / h1;
            double coeff2 = 0.5 * (k1 + k2) / h2;
            double coeff3 = 0.5 * (h1 + h2) / k1;
            double coeff4 = 0.5 * (h1 + h2) / k2;
            /* Fill result(i,j) */
            result[grid.index(i_r, i_theta)] -=
                (0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF) * x[grid.index(i_r, i_theta)] /* beta_{i,j} */
                 - coeff1 * arr * x[grid.index(i_r, i_theta + (grid.ntheta() / 2))] /* Left */
                 - coeff2 * arr * x[grid.index(i_r + 1, i_theta)] /* Right */
                 - coeff3 * att * x[grid.index(i_r, i_theta - 1)] /* Bottom */
                 - coeff4 * att * x[grid.index(i_r, i_theta + 1)] /* Top */
                 + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) *
                       x[grid.index(i_r, i_theta)]); /* Center: (Left, Right, Bottom, Top) */
            /* Fill result(i-1,j) */
            /* From view the view of the across origin node, the directions are roatated by 180 degrees in the stencil! */
            result[grid.index(i_r, i_theta + (grid.ntheta() / 2))] -=
                (-coeff1 * arr * x[grid.index(i_r, i_theta)] /* Right -> Left */
                 + coeff1 * arr *
                       x[grid.index(i_r, i_theta + (grid.ntheta() / 2))]); /* Center: (Right) -> Center: (Left)*/
            /*  + 0.25 * art * x[grid.index(i_r,i_theta+1)]; // Top Right -> Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
            /*  - 0.25 * art * x[grid.index(i_r,i_theta-1)]; // Bottom Right -> Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
            /* Fill result(i+1,j) */
            result[grid.index(i_r + 1, i_theta)] -=
                (-coeff2 * arr * x[grid.index(i_r, i_theta)] /* Left */
                 + coeff2 * arr * x[grid.index(i_r + 1, i_theta)] /* Center: (Left) */
                 + 0.25 * art * x[grid.index(i_r, i_theta + 1)] /* Top Left */
                 - 0.25 * art * x[grid.index(i_r, i_theta - 1)]); /* Bottom Left */
            /* Fill result(i,j-1) */
            result[grid.index(i_r, i_theta - 1)] -=
                (-coeff3 * att * x[grid.index(i_r, i_theta)] /* Top */
                 + coeff3 * att * x[grid.index(i_r, i_theta - 1)] /* Center: (Top) */
                 - 0.25 * art * x[grid.index(i_r + 1, i_theta)]); /* Top Right */
            /*  + 0.25 * art * x[grid.index(i_r-1,i_theta)]; // Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
            /* Fill result(i,j+1) */
            result[grid.index(i_r, i_theta + 1)] -=
                (-coeff4 * att * x[grid.index(i_r, i_theta)] /* Bottom */
                 + coeff4 * att * x[grid.index(i_r, i_theta + 1)] /* Center: (Bottom) */
                 + 0.25 * art * x[grid.index(i_r + 1, i_theta)]); /* Bottom Right */
            /*  - 0.25 * art * x[grid.index(i_r-1,i_theta)]; // Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
        }
        /* ------------------------------- */
        /* Node next to the inner boundary */
        /* ------------------------------- */
    }
    else if (i_r == 1) {
        double h1     = grid.radialSpacing(i_r - 1);
        double h2     = grid.radialSpacing(i_r);
        double k1     = grid.angularSpacing(i_theta - 1);
        double k2     = grid.angularSpacing(i_theta);
        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;
        /* Fill result(i,j) */
        result[grid.index(i_r, i_theta)] -=
            (0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF) * x[grid.index(i_r, i_theta)] /* beta_{i,j} */
             - coeff1 * arr * x[grid.index(i_r - 1, i_theta)] /* Left */
             - coeff2 * arr * x[grid.index(i_r + 1, i_theta)] /* Right */
             - coeff3 * att * x[grid.index(i_r, i_theta - 1)] /* Bottom */
             - coeff4 * att * x[grid.index(i_r, i_theta + 1)] /* Top */
             + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) *
                   x[grid.index(i_r, i_theta)]); /* Center: (Left, Right, Bottom, Top) */
        /* Fill result(i-1,j) */
        if (!DirBC_Interior) { /* Don't give to the inner dirichlet boundary! */
            result[grid.index(i_r - 1, i_theta)] -=
                (-coeff1 * arr * x[grid.index(i_r, i_theta)] /* Right */
                 + coeff1 * arr * x[grid.index(i_r - 1, i_theta)] /* Center: (Right) */
                 - 0.25 * art * x[grid.index(i_r, i_theta + 1)] /* Top Right */
                 + 0.25 * art * x[grid.index(i_r, i_theta - 1)]); /* Bottom Right */
        }
        /* Fill result(i+1,j) */
        result[grid.index(i_r + 1, i_theta)] -= (-coeff2 * arr * x[grid.index(i_r, i_theta)] /* Left */
                                                 + coeff2 * arr * x[grid.index(i_r + 1, i_theta)] /* Center: (Left) */
                                                 + 0.25 * art * x[grid.index(i_r, i_theta + 1)] /* Top Left */
                                                 - 0.25 * art * x[grid.index(i_r, i_theta - 1)]); /* Bottom Left */
        /* Fill result(i,j-1) */
        result[grid.index(i_r, i_theta - 1)] -= (-coeff3 * att * x[grid.index(i_r, i_theta)] /* Top */
                                                 + coeff3 * att * x[grid.index(i_r, i_theta - 1)] /* Center: (Top) */
                                                 - 0.25 * art * x[grid.index(i_r + 1, i_theta)] /* Top Right */
                                                 + 0.25 * art * x[grid.index(i_r - 1, i_theta)]); /* Top Left */
        /* Fill result(i,j+1) */
        result[grid.index(i_r, i_theta + 1)] -= (-coeff4 * att * x[grid.index(i_r, i_theta)] /* Bottom */
                                                 + coeff4 * att * x[grid.index(i_r, i_theta + 1)] /* Center: (Bottom) */
                                                 + 0.25 * art * x[grid.index(i_r + 1, i_theta)] /* Bottom Right */
                                                 - 0.25 * art * x[grid.index(i_r - 1, i_theta)]); /* Bottom Left */
        /* ------------------------------- */
        /* Node next to the outer boundary */
        /* ------------------------------- */
    }
    else if (i_r == grid.nr() - 2) {
        double h1     = grid.radialSpacing(i_r - 1);
        double h2     = grid.radialSpacing(i_r);
        double k1     = grid.angularSpacing(i_theta - 1);
        double k2     = grid.angularSpacing(i_theta);
        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;
        /* Fill result(i,j) */
        result[grid.index(i_r, i_theta)] -=
            (0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF) * x[grid.index(i_r, i_theta)] /* beta_{i,j} */
             - coeff1 * arr * x[grid.index(i_r - 1, i_theta)] /* Left */
             - coeff2 * arr * x[grid.index(i_r + 1, i_theta)] /* Right */
             - coeff3 * att * x[grid.index(i_r, i_theta - 1)] /* Bottom */
             - coeff4 * att * x[grid.index(i_r, i_theta + 1)] /* Top */
             + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) *
                   x[grid.index(i_r, i_theta)]); /* Center: (Left, Right, Bottom, Top) */
        /* Fill result(i-1,j) */
        result[grid.index(i_r - 1, i_theta)] -= (-coeff1 * arr * x[grid.index(i_r, i_theta)] /* Right */
                                                 + coeff1 * arr * x[grid.index(i_r - 1, i_theta)] /* Center: (Right) */
                                                 - 0.25 * art * x[grid.index(i_r, i_theta + 1)] /* Top Right */
                                                 + 0.25 * art * x[grid.index(i_r, i_theta - 1)]); /* Bottom Right */
        /* Don't give to the outer dirichlet boundary! */
        /* Fill result(i+1,j) */
        /* result[grid.index(i_r+1,i_theta)] -= ( */
        /*     - coeff2 * arr * x[grid.index(i_r,i_theta)] // Left */
        /*     + coeff2 * arr * x[grid.index(i_r+1,i_theta)] // Center: (Left) */
        /*     + 0.25 * art * x[grid.index(i_r,i_theta+1)] // Top Left */
        /*     - 0.25 * art * x[grid.index(i_r,i_theta-1)] ); // Bottom Left */
        /* Fill result(i,j-1) */
        result[grid.index(i_r, i_theta - 1)] -= (-coeff3 * att * x[grid.index(i_r, i_theta)] /* Top */
                                                 + coeff3 * att * x[grid.index(i_r, i_theta - 1)] /* Center: (Top) */
                                                 - 0.25 * art * x[grid.index(i_r + 1, i_theta)] /* Top Right */
                                                 + 0.25 * art * x[grid.index(i_r - 1, i_theta)]); /* Top Left */
        /* Fill result(i,j+1) */
        result[grid.index(i_r, i_theta + 1)] -= (-coeff4 * att * x[grid.index(i_r, i_theta)] /* Bottom */
                                                 + coeff4 * att * x[grid.index(i_r, i_theta + 1)] /* Center: (Bottom) */
                                                 + 0.25 * art * x[grid.index(i_r + 1, i_theta)] /* Bottom Right */
                                                 - 0.25 * art * x[grid.index(i_r - 1, i_theta)]); /* Bottom Left */
        /* ----------------------------- */
        /* Node on to the outer boundary */
        /* ----------------------------- */
    }
    else if (i_r == grid.nr() - 1) {
        /* Fill result of (i,j) */
        result[grid.index(i_r, i_theta)] -= x[grid.index(i_r, i_theta)];
        /* Give value to the interior nodes! */
        double h1     = grid.radialSpacing(i_r - 1);
        double k1     = grid.angularSpacing(i_theta - 1);
        double k2     = grid.angularSpacing(i_theta);
        double coeff1 = 0.5 * (k1 + k2) / h1;
        /* Fill result(i-1,j) */
        result[grid.index(i_r - 1, i_theta)] -= (-coeff1 * arr * x[grid.index(i_r, i_theta)] /* Right */
                                                 + coeff1 * arr * x[grid.index(i_r - 1, i_theta)] /* Center: (Right) */
                                                 - 0.25 * art * x[grid.index(i_r, i_theta + 1)] /* Top Right */
                                                 + 0.25 * art * x[grid.index(i_r, i_theta - 1)]); /* Bottom Right */
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
