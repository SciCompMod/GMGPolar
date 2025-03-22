#include "../../../include/Smoother/SmootherGive/smootherGive.h"

#define COMPUTE_JACOBIAN_ELEMENTS(domain_geometry, r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDF)  \
    do {                                                                                                               \
        /* Calculate the elements of the Jacobian matrix for the transformation mapping */                             \
        /* The Jacobian matrix is: */                                                                                  \
        /* [Jrr, Jrt] */                                                                                               \
        /* [Jtr, Jtt] */                                                                                               \
        const double Jrr = domain_geometry.dFx_dr(r, theta, sin_theta, cos_theta);                                     \
        const double Jtr = domain_geometry.dFy_dr(r, theta, sin_theta, cos_theta);                                     \
        const double Jrt = domain_geometry.dFx_dt(r, theta, sin_theta, cos_theta);                                     \
        const double Jtt = domain_geometry.dFy_dt(r, theta, sin_theta, cos_theta);                                     \
        /* Compute the determinant of the Jacobian matrix */                                                           \
        detDF = Jrr * Jtt - Jrt * Jtr;                                                                                 \
        /* Compute the elements of the symmetric matrix: */                                                            \
        /* 0.5 * alpha * DF^{-1} * DF^{-T} * |det(DF)| */                                                              \
        /* which is represented by: */                                                                                 \
        /* [arr, 0.5*art] */                                                                                           \
        /* [0.5*atr, att] */                                                                                           \
        arr = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / fabs(detDF);                                               \
        att = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / fabs(detDF);                                               \
        art = (-Jtt * Jtr - Jrt * Jrr) * coeff_alpha / fabs(detDF);                                                    \
        /* Note that the inverse Jacobian matrix DF^{-1} is: */                                                        \
        /* 1.0 / det(DF) *   */                                                                                        \
        /* [Jtt, -Jrt] */                                                                                              \
        /* [-Jtr, Jrr] */                                                                                              \
    } while (0)

#define NODE_APPLY_ASC_ORTHO_CIRCLE_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, grid, DirBC_Interior,                     \
                                         smoother_color, x, rhs, temp, arr, att, art, detDF, coeff_beta)                         \
    do {                                                                                                                         \
        assert(i_r >= 0 && i_r <= grid_.numberSmootherCircles());                                                                \
        bool isOddNumberSmootherCircles = (grid.numberSmootherCircles() & 1);                                                    \
        bool isOddRadialIndex           = (i_r & 1);                                                                             \
        SmootherColor node_color =                                                                                               \
            (isOddNumberSmootherCircles == isOddRadialIndex) ? SmootherColor::White : SmootherColor::Black;                      \
        /* -------------------- */                                                                                               \
        /* Node in the interior */                                                                                               \
        /* -------------------- */                                                                                               \
        if (i_r > 0 && i_r < grid.numberSmootherCircles()) {                                                                     \
            double h1     = grid.radialSpacing(i_r - 1);                                                                         \
            double h2     = grid.radialSpacing(i_r);                                                                             \
            double k1     = grid.angularSpacing(i_theta - 1);                                                                    \
            double k2     = grid.angularSpacing(i_theta);                                                                        \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                                \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                                \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                                \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                                \
            /* -------------------- */                                                                                           \
            /* Inside Section Parts */                                                                                           \
            if (node_color == smoother_color) {                                                                                  \
                /* Fill temp(i,j) */                                                                                             \
                temp[grid.index(i_r, i_theta)] -= (-coeff1 * arr * x[grid.index(i_r - 1, i_theta)] /* Left */                    \
                                                   - coeff2 * arr * x[grid.index(i_r + 1, i_theta)] /* Right */                  \
                );                                                                                                               \
                /* Fill temp(i,j-1) */                                                                                           \
                temp[grid.index(i_r, i_theta - 1)] -= (-0.25 * art * x[grid.index(i_r + 1, i_theta)] /* Top Right */             \
                                                       + 0.25 * art * x[grid.index(i_r - 1, i_theta)]); /* Top Left */           \
                /* Fill temp(i,j+1) */                                                                                           \
                temp[grid.index(i_r, i_theta + 1)] -=                                                                            \
                    (+0.25 * art * x[grid.index(i_r + 1, i_theta)] /* Bottom Right */                                            \
                     - 0.25 * art * x[grid.index(i_r - 1, i_theta)]); /* Bottom Left */                                          \
            }                                                                                                                    \
            /* --------------------- */                                                                                          \
            /* Outside Section Parts */                                                                                          \
            else if (node_color != smoother_color) {                                                                             \
                /* Fill temp(i-1,j) */                                                                                           \
                if (i_r > 1 || !DirBC_Interior) {                                                                                \
                    temp[grid.index(i_r - 1, i_theta)] -=                                                                        \
                        (-coeff1 * arr * x[grid.index(i_r, i_theta)] /* Right */                                                 \
                         - 0.25 * art * x[grid.index(i_r, i_theta + 1)] /* Top Right */                                          \
                         + 0.25 * art * x[grid.index(i_r, i_theta - 1)]); /* Bottom Right */                                     \
                }                                                                                                                \
                /* Fill temp(i+1,j) */                                                                                           \
                if (i_r < grid.numberSmootherCircles() - 1) {                                                                    \
                    temp[grid.index(i_r + 1, i_theta)] -=                                                                        \
                        (-coeff2 * arr * x[grid.index(i_r, i_theta)] /* Left */                                                  \
                         + 0.25 * art * x[grid.index(i_r, i_theta + 1)] /* Top Left */                                           \
                         - 0.25 * art * x[grid.index(i_r, i_theta - 1)]); /* Bottom Left */                                      \
                }                                                                                                                \
            }                                                                                                                    \
        }                                                                                                                        \
        /* -------------------- */                                                                                               \
        /* Node on the boundary */                                                                                               \
        /* -------------------- */                                                                                               \
        else if (i_r == 0) {                                                                                                     \
            /* ------------------------------------------------ */                                                               \
            /* Case 1: Dirichlet boundary on the inner boundary */                                                               \
            /* ------------------------------------------------ */                                                               \
            if (DirBC_Interior) {                                                                                                \
                double h2     = grid.radialSpacing(i_r);                                                                         \
                double k1     = grid.angularSpacing(i_theta - 1);                                                                \
                double k2     = grid.angularSpacing(i_theta);                                                                    \
                double coeff2 = 0.5 * (k1 + k2) / h2;                                                                            \
                /* -------------------- */                                                                                       \
                /* Inside Section Parts */                                                                                       \
                if (node_color == smoother_color) {                                                                              \
                }                                                                                                                \
                /* --------------------- */                                                                                      \
                /* Outside Section Parts */                                                                                      \
                else if (node_color != smoother_color) {                                                                         \
                    /* Fill temp(i+1,j) */                                                                                       \
                    temp[grid.index(i_r + 1, i_theta)] -=                                                                        \
                        (-coeff2 * arr * x[grid.index(i_r, i_theta)] /* Left */                                                  \
                         + 0.25 * art * x[grid.index(i_r, i_theta + 1)] /* Top Left */                                           \
                         - 0.25 * art * x[grid.index(i_r, i_theta - 1)]); /* Bottom Left */                                      \
                }                                                                                                                \
            }                                                                                                                    \
            else {                                                                                                               \
                /* ------------------------------------------------------------- */                                              \
                /* Case 2: Across origin discretization on the interior boundary */                                              \
                /* ------------------------------------------------------------- */                                              \
                /* h1 gets replaced with 2 * R0. */                                                                              \
                /* (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()/2)). */                                     \
                /* Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil. */              \
                double h1     = 2.0 * grid.radius(0);                                                                            \
                double h2     = grid.radialSpacing(i_r);                                                                         \
                double k1     = grid.angularSpacing(i_theta - 1);                                                                \
                double k2     = grid.angularSpacing(i_theta);                                                                    \
                double coeff1 = 0.5 * (k1 + k2) / h1;                                                                            \
                double coeff2 = 0.5 * (k1 + k2) / h2;                                                                            \
                double coeff3 = 0.5 * (h1 + h2) / k1;                                                                            \
                double coeff4 = 0.5 * (h1 + h2) / k2;                                                                            \
                /* -------------------- */                                                                                       \
                /* Inside Section Parts */                                                                                       \
                if (node_color == smoother_color) {                                                                              \
                    /* Fill temp(i,j) */                                                                                         \
                    temp[grid.index(i_r, i_theta)] -=                                                                            \
                        (/* - coeff1 * arr * x[grid.index(i_r, i_theta + (grid.ntheta()/2))] // Left: Not in Asc_ortho */        \
                         -coeff2 * arr * x[grid.index(i_r + 1, i_theta)] /* Right */                                             \
                        );                                                                                                       \
                    /* Fill temp(i,j-1) */                                                                                       \
                    temp[grid.index(i_r, i_theta - 1)] -=                                                                        \
                        (-0.25 * art * x[grid.index(i_r + 1, i_theta)]); /* Top Right */                                         \
                    /*  + 0.25 * art * x[grid.index(i_r-1,i_theta)]; // Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */    \
                    /* Fill temp(i,j+1) */                                                                                       \
                    temp[grid.index(i_r, i_theta + 1)] -=                                                                        \
                        (+0.25 * art * x[grid.index(i_r + 1, i_theta)]); /* Bottom Right */                                      \
                    /*  - 0.25 * art * x[grid.index(i_r-1,i_theta)]; // Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
                }                                                                                                                \
                /* --------------------- */                                                                                      \
                /* Outside Section Parts */                                                                                      \
                else if (node_color != smoother_color) {                                                                         \
                    /* Fill temp(i+1,j) */                                                                                       \
                    temp[grid.index(i_r + 1, i_theta)] -=                                                                        \
                        (-coeff2 * arr * x[grid.index(i_r, i_theta)] /* Left */                                                  \
                         + 0.25 * art * x[grid.index(i_r, i_theta + 1)] /* Top Left */                                           \
                         - 0.25 * art * x[grid.index(i_r, i_theta - 1)]); /* Bottom Left */                                      \
                }                                                                                                                \
            }                                                                                                                    \
        }                                                                                                                        \
        /* ----------------------------- */                                                                                      \
        /* Node next to circular section */                                                                                      \
        /* ----------------------------- */                                                                                      \
        else if (i_r == grid.numberSmootherCircles()) {                                                                          \
            assert(node_color == SmootherColor::White);                                                                          \
            if (smoother_color == SmootherColor::Black) {                                                                        \
                double h1     = grid.radialSpacing(i_r - 1);                                                                     \
                double k1     = grid.angularSpacing(i_theta - 1);                                                                \
                double k2     = grid.angularSpacing(i_theta);                                                                    \
                double coeff1 = 0.5 * (k1 + k2) / h1;                                                                            \
                /* --------------------- */                                                                                      \
                /* Outside Section Parts */                                                                                      \
                /* Fill temp(i-1,j) */                                                                                           \
                temp[grid.index(i_r - 1, i_theta)] -=                                                                            \
                    (-coeff1 * arr * x[grid.index(i_r, i_theta)] /* Right */                                                     \
                     - 0.25 * art * x[grid.index(i_r, i_theta + 1)] /* Top Right */                                              \
                     + 0.25 * art * x[grid.index(i_r, i_theta - 1)]); /* Bottom Right */                                         \
            }                                                                                                                    \
        }                                                                                                                        \
    } while (0)

#define NODE_APPLY_ASC_ORTHO_RADIAL_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, grid, DirBC_Interior,           \
                                         smoother_color, x, rhs, temp, arr, att, art, detDF, coeff_beta)               \
    do {                                                                                                               \
        assert(i_r >= grid.numberSmootherCircles() - 1 && i_r < grid.nr());                                            \
        SmootherColor node_color = (i_theta & 1) ? SmootherColor::White : SmootherColor::Black;                        \
        /* -------------------- */                                                                                     \
        /* Node in the interior */                                                                                     \
        /* -------------------- */                                                                                     \
        if (i_r > grid.numberSmootherCircles() && i_r < grid.nr() - 2) {                                               \
            double h1     = grid.radialSpacing(i_r - 1);                                                               \
            double h2     = grid.radialSpacing(i_r);                                                                   \
            double k1     = grid.angularSpacing(i_theta - 1);                                                          \
            double k2     = grid.angularSpacing(i_theta);                                                              \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                      \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                      \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                      \
            /* -------------------- */                                                                                 \
            /* Inside Section Parts */                                                                                 \
            if (node_color == smoother_color) {                                                                        \
                /* Fill temp(i,j) */                                                                                   \
                temp[grid.index(i_r, i_theta)] -= (-coeff3 * att * x[grid.index(i_r, i_theta - 1)] /* Bottom */        \
                                                   - coeff4 * att * x[grid.index(i_r, i_theta + 1)] /* Top */          \
                );                                                                                                     \
                /* Fill temp(i-1,j) */                                                                                 \
                temp[grid.index(i_r - 1, i_theta)] -=                                                                  \
                    (-0.25 * art * x[grid.index(i_r, i_theta + 1)] /* Top Right */                                     \
                     + 0.25 * art * x[grid.index(i_r, i_theta - 1)]); /* Bottom Right */                               \
                /* Fill temp(i+1,j) */                                                                                 \
                temp[grid.index(i_r + 1, i_theta)] -=                                                                  \
                    (+0.25 * art * x[grid.index(i_r, i_theta + 1)] /* Top Left */                                      \
                     - 0.25 * art * x[grid.index(i_r, i_theta - 1)]); /* Bottom Left */                                \
            }                                                                                                          \
            /* --------------------- */                                                                                \
            /* Outside Section Parts */                                                                                \
            else if (node_color != smoother_color) {                                                                   \
                /* Fill temp(i,j-1) */                                                                                 \
                temp[grid.index(i_r, i_theta - 1)] -= (-coeff3 * att * x[grid.index(i_r, i_theta)] /* Top */           \
                                                       - 0.25 * art * x[grid.index(i_r + 1, i_theta)] /* Top Right */  \
                                                       + 0.25 * art * x[grid.index(i_r - 1, i_theta)]); /* Top Left */ \
                /* Fill temp(i,j+1) */                                                                                 \
                temp[grid.index(i_r, i_theta + 1)] -=                                                                  \
                    (-coeff4 * att * x[grid.index(i_r, i_theta)] /* Bottom */                                          \
                     + 0.25 * art * x[grid.index(i_r + 1, i_theta)] /* Bottom Right */                                 \
                     - 0.25 * art * x[grid.index(i_r - 1, i_theta)]); /* Bottom Left */                                \
            }                                                                                                          \
        }                                                                                                              \
        else if (i_r == grid.numberSmootherCircles() - 1) {                                                            \
            double h1     = grid.radialSpacing(i_r - 1);                                                               \
            double h2     = grid.radialSpacing(i_r);                                                                   \
            double k1     = grid.angularSpacing(i_theta - 1);                                                          \
            double k2     = grid.angularSpacing(i_theta);                                                              \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                      \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                      \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                      \
            /* -------------------- */                                                                                 \
            /* Inside Section Parts */                                                                                 \
            if (node_color == smoother_color) {                                                                        \
                /* Fill temp(i+1,j) */                                                                                 \
                temp[grid.index(i_r + 1, i_theta)] -=                                                                  \
                    (-coeff2 * arr * x[grid.index(i_r, i_theta)] /* Left */                                            \
                     + 0.25 * art * x[grid.index(i_r, i_theta + 1)] /* Top Left */                                     \
                     - 0.25 * art * x[grid.index(i_r, i_theta - 1)]); /* Bottom Left */                                \
            }                                                                                                          \
            /* --------------------- */                                                                                \
            /* Outside Section Parts */                                                                                \
            else if (node_color != smoother_color) {                                                                   \
                /* Nothing to be done here */                                                                          \
            }                                                                                                          \
        }                                                                                                              \
        else if (i_r == grid.numberSmootherCircles()) {                                                                \
            double h1     = grid.radialSpacing(i_r - 1);                                                               \
            double h2     = grid.radialSpacing(i_r);                                                                   \
            double k1     = grid.angularSpacing(i_theta - 1);                                                          \
            double k2     = grid.angularSpacing(i_theta);                                                              \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                      \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                      \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                      \
            /* -------------------- */                                                                                 \
            /* Inside Section Parts */                                                                                 \
            if (node_color == smoother_color) {                                                                        \
                /* Fill temp(i,j) */                                                                                   \
                temp[grid.index(i_r, i_theta)] -= (-coeff1 * arr * x[grid.index(i_r - 1, i_theta)] /* Left */          \
                                                   - coeff3 * att * x[grid.index(i_r, i_theta - 1)] /* Bottom */       \
                                                   - coeff4 * att * x[grid.index(i_r, i_theta + 1)] /* Top */          \
                );                                                                                                     \
                /* Fill temp(i+1,j) */                                                                                 \
                temp[grid.index(i_r + 1, i_theta)] -=                                                                  \
                    (+0.25 * art * x[grid.index(i_r, i_theta + 1)] /* Top Left */                                      \
                     - 0.25 * art * x[grid.index(i_r, i_theta - 1)]); /* Bottom Left */                                \
            }                                                                                                          \
            /* --------------------- */                                                                                \
            /* Outside Section Parts */                                                                                \
            else if (node_color != smoother_color) {                                                                   \
                /* Fill temp(i,j-1) */                                                                                 \
                temp[grid.index(i_r, i_theta - 1)] -= (-coeff3 * att * x[grid.index(i_r, i_theta)] /* Top */           \
                                                       - 0.25 * art * x[grid.index(i_r + 1, i_theta)] /* Top Right */  \
                                                       + 0.25 * art * x[grid.index(i_r - 1, i_theta)]); /* Top Left */ \
                /* Fill temp(i,j+1) */                                                                                 \
                temp[grid.index(i_r, i_theta + 1)] -=                                                                  \
                    (-coeff4 * att * x[grid.index(i_r, i_theta)] /* Bottom */                                          \
                     + 0.25 * art * x[grid.index(i_r + 1, i_theta)] /* Bottom Right */                                 \
                     - 0.25 * art * x[grid.index(i_r - 1, i_theta)]); /* Bottom Left */                                \
            }                                                                                                          \
        }                                                                                                              \
        else if (i_r == grid.nr() - 2) {                                                                               \
            double h1     = grid.radialSpacing(i_r - 1);                                                               \
            double h2     = grid.radialSpacing(i_r);                                                                   \
            double k1     = grid.angularSpacing(i_theta - 1);                                                          \
            double k2     = grid.angularSpacing(i_theta);                                                              \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                      \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                      \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                      \
            /* -------------------- */                                                                                 \
            /* Inside Section Parts */                                                                                 \
            if (node_color == smoother_color) {                                                                        \
                /* Fill temp(i,j) */                                                                                   \
                temp[grid.index(i_r, i_theta)] -= (-coeff3 * att * x[grid.index(i_r, i_theta - 1)] /* Bottom */        \
                                                   - coeff4 * att * x[grid.index(i_r, i_theta + 1)] /* Top */          \
                );                                                                                                     \
                /* Fill temp(i-1,j) */                                                                                 \
                temp[grid.index(i_r - 1, i_theta)] -=                                                                  \
                    (-0.25 * art * x[grid.index(i_r, i_theta + 1)] /* Top Right */                                     \
                     + 0.25 * art * x[grid.index(i_r, i_theta - 1)]); /* Bottom Right */                               \
                                                                                                                       \
                /* "Right" is part of the radial Asc smoother matrices, */                                             \
                /* but is shifted over to the rhs to make the radial Asc smoother matrices symmetric. */               \
                /* Note that the circle Asc smoother matrices are symmetric by default. */                             \
                temp[grid.index(i_r, i_theta)] -= /* Right: Symmetry shift! */                                         \
                    -coeff2 * arr * rhs[grid.index(i_r + 1, i_theta)];                                                 \
            }                                                                                                          \
            /* --------------------- */                                                                                \
            /* Outside Section Parts */                                                                                \
            else if (node_color != smoother_color) {                                                                   \
                /* Fill temp(i,j-1) */                                                                                 \
                temp[grid.index(i_r, i_theta - 1)] -= (-coeff3 * att * x[grid.index(i_r, i_theta)] /* Top */           \
                                                       - 0.25 * art * x[grid.index(i_r + 1, i_theta)] /* Top Right */  \
                                                       + 0.25 * art * x[grid.index(i_r - 1, i_theta)]); /* Top Left */ \
                /* Fill temp(i,j+1) */                                                                                 \
                temp[grid.index(i_r, i_theta + 1)] -=                                                                  \
                    (-coeff4 * att * x[grid.index(i_r, i_theta)] /* Bottom */                                          \
                     + 0.25 * art * x[grid.index(i_r + 1, i_theta)] /* Bottom Right */                                 \
                     - 0.25 * art * x[grid.index(i_r - 1, i_theta)]); /* Bottom Left */                                \
            }                                                                                                          \
        }                                                                                                              \
        else if (i_r == grid.nr() - 1) {                                                                               \
            double h1     = grid.radialSpacing(i_r - 1);                                                               \
            double k1     = grid.angularSpacing(i_theta - 1);                                                          \
            double k2     = grid.angularSpacing(i_theta);                                                              \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            /* -------------------- */                                                                                 \
            /* Inside Section Parts */                                                                                 \
            if (node_color == smoother_color) {                                                                        \
                /* Fill temp(i-1,j) */                                                                                 \
                temp[grid.index(i_r - 1, i_theta)] -=                                                                  \
                    (-0.25 * art * x[grid.index(i_r, i_theta + 1)] /* Top Right */                                     \
                     + 0.25 * art * x[grid.index(i_r, i_theta - 1)] /* Bottom Right */                                 \
                    );                                                                                                 \
                /* "Right" is part of the radial Asc smoother matrices, */                                             \
                /* but is shifted over to the rhs to make the radial Asc smoother matrices symmetric. */               \
                /* Note that the circle Asc smoother matrices are symmetric by default. */                             \
                temp[grid.index(i_r - 1, i_theta)] -= (-coeff1 * arr * rhs[grid.index(i_r, i_theta)] /* Right */       \
                );                                                                                                     \
            }                                                                                                          \
            /* --------------------- */                                                                                \
            /* Outside Section Parts */                                                                                \
            else if (node_color != smoother_color) {                                                                   \
                /* Nothing to be done here */                                                                          \
            }                                                                                                          \
        }                                                                                                              \
    } while (0)

void SmootherGive::applyAscOrthoCircleSection(const int i_r, const SmootherColor smoother_color,
                                              const Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp)
{
    assert(i_r >= 0 && i_r < grid_.numberSmootherCircles() + 1);

    const auto& sin_theta_cache = level_cache_.sin_theta();
    const auto& cos_theta_cache = level_cache_.cos_theta();

    const double r = grid_.radius(i_r);

    double coeff_beta;
    if (level_cache_.cacheDensityProfileCoefficients()) {
        coeff_beta = level_cache_.coeff_beta()[i_r];
    }
    else {
        coeff_beta = density_profile_coefficients_.beta(r);
    }

    double coeff_alpha;
    if (!level_cache_.cacheDomainGeometry()) {
        if (level_cache_.cacheDensityProfileCoefficients()) {
            coeff_alpha = level_cache_.coeff_alpha()[i_r];
        }
        else {
            coeff_alpha = density_profile_coefficients_.alpha(r);
        }
    }

    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        const double theta     = grid_.theta(i_theta);
        const double sin_theta = sin_theta_cache[i_theta];
        const double cos_theta = cos_theta_cache[i_theta];

        /* Compute arr, att, art, detDF value at the current node */
        double arr, att, art, detDF;
        if (level_cache_.cacheDomainGeometry()) {
            const int index = grid_.index(i_r, i_theta);
            arr             = level_cache_.arr()[index];
            att             = level_cache_.att()[index];
            art             = level_cache_.art()[index];
            detDF           = level_cache_.detDF()[index];
        }
        else {
            COMPUTE_JACOBIAN_ELEMENTS(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art,
                                      detDF);
        }

        // Apply Asc Ortho at the current node
        NODE_APPLY_ASC_ORTHO_CIRCLE_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, grid_, DirBC_Interior_,
                                         smoother_color, x, rhs, temp, arr, att, art, detDF, coeff_beta);
    }
}

void SmootherGive::applyAscOrthoRadialSection(const int i_theta, const SmootherColor smoother_color,
                                              const Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp)
{
    const auto& sin_theta_cache = level_cache_.sin_theta();
    const auto& cos_theta_cache = level_cache_.cos_theta();

    const double theta     = grid_.theta(i_theta);
    const double sin_theta = sin_theta_cache[i_theta];
    const double cos_theta = cos_theta_cache[i_theta];

    /* !!! i_r = grid_.numberSmootherCircles()-1 !!! */
    for (int i_r = grid_.numberSmootherCircles() - 1; i_r < grid_.nr(); i_r++) {
        const double r = grid_.radius(i_r);

        double coeff_beta;
        if (level_cache_.cacheDensityProfileCoefficients()) {
            coeff_beta = level_cache_.coeff_beta()[i_r];
        }
        else {
            coeff_beta = density_profile_coefficients_.beta(r);
        }

        double coeff_alpha;
        if (!level_cache_.cacheDomainGeometry()) {
            if (level_cache_.cacheDensityProfileCoefficients()) {
                coeff_alpha = level_cache_.coeff_alpha()[i_r];
            }
            else {
                coeff_alpha = density_profile_coefficients_.alpha(r);
            }
        }

        /* Compute arr, att, art, detDF value at the current node */
        double arr, att, art, detDF;
        if (level_cache_.cacheDomainGeometry()) {
            const int index = grid_.index(i_r, i_theta);
            arr             = level_cache_.arr()[index];
            att             = level_cache_.att()[index];
            art             = level_cache_.art()[index];
            detDF           = level_cache_.detDF()[index];
        }
        else {
            COMPUTE_JACOBIAN_ELEMENTS(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art,
                                      detDF);
        }

        // Apply Asc Ortho at the current node
        NODE_APPLY_ASC_ORTHO_RADIAL_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, grid_, DirBC_Interior_,
                                         smoother_color, x, rhs, temp, arr, att, art, detDF, coeff_beta);
    }
}

void SmootherGive::solveCircleSection(const int i_r, Vector<double>& x, Vector<double>& temp,
                                      Vector<double>& solver_storage_1, Vector<double>& solver_storage_2)
{
    const int start = grid_.index(i_r, 0);
    const int end   = start + grid_.ntheta();
    if (i_r == 0) {
        inner_boundary_mumps_solver_.job    = JOB_COMPUTE_SOLUTION;
        inner_boundary_mumps_solver_.nrhs   = 1; // single rhs vector
        inner_boundary_mumps_solver_.nz_rhs = grid_.ntheta(); // non-zeros in rhs
        inner_boundary_mumps_solver_.rhs    = temp.begin() + start;
        inner_boundary_mumps_solver_.lrhs   = grid_.ntheta(); // leading dimension of rhs
        dmumps_c(&inner_boundary_mumps_solver_);
        if (inner_boundary_mumps_solver_.info[0] != 0) {
            std::cerr << "Error solving the system: " << inner_boundary_mumps_solver_.info[0] << std::endl;
        }
    }
    else {
        circle_tridiagonal_solver_[i_r].solveInPlace(temp.begin() + start, solver_storage_1.begin(),
                                                     solver_storage_2.begin());
    }
    // Move updated values to x
    std::move(temp.begin() + start, temp.begin() + end, x.begin() + start);
}

void SmootherGive::solveRadialSection(const int i_theta, Vector<double>& x, Vector<double>& temp,
                                      Vector<double>& solver_storage)
{
    const int start = grid_.index(grid_.numberSmootherCircles(), i_theta);
    const int end   = start + grid_.lengthSmootherRadial();

    radial_tridiagonal_solver_[i_theta].solveInPlace(temp.begin() + start, solver_storage.begin());
    // Move updated values to x
    std::move(temp.begin() + start, temp.begin() + end, x.begin() + start);
}

/* ------------------ */
/* Sequential Version */
/* ------------------ */

void SmootherGive::smoothingSequential(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp)
{
    assert(x.size() == rhs.size());
    assert(temp.size() == rhs.size());

    temp = rhs;

    /* Single-threaded execution */
    Vector<double> circle_solver_storage_1(grid_.ntheta());
    Vector<double> circle_solver_storage_2(grid_.ntheta());
    Vector<double> radial_solver_storage(grid_.lengthSmootherRadial());

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
void SmootherGive::smoothingForLoop(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp)
{
    assert(x.size() == rhs.size());
    assert(temp.size() == rhs.size());

    omp_set_num_threads(num_omp_threads_);

    if (omp_get_max_threads() == 1) {
        smoothingSequential(x, rhs, temp);
    }
    else {
        temp = rhs;

        /* Multi-threaded execution */
        const int num_circle_tasks = grid_.numberSmootherCircles();
        const int num_radial_tasks = grid_.ntheta();

        #pragma omp parallel
        {
            Vector<double> circle_solver_storage_1(grid_.ntheta());
            Vector<double> circle_solver_storage_2(grid_.ntheta());
            Vector<double> radial_solver_storage(grid_.lengthSmootherRadial());

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
// clang-format on
