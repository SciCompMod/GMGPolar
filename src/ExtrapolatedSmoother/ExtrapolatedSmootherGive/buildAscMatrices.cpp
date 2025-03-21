#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherGive/extrapolatedSmootherGive.h"

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
        /*  [arr, 0.5*art]  */                                                                                         \
        /*  [0.5*atr, att]  */                                                                                         \
        arr = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / fabs(detDF);                                               \
        att = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / fabs(detDF);                                               \
        art = (-Jtt * Jtr - Jrt * Jrr) * coeff_alpha / fabs(detDF);                                                    \
        /* Note that the inverse Jacobian matrix DF^{-1} is: */                                                        \
        /*  1.0 / det(DF) *  */                                                                                        \
        /*  [Jtt, -Jrt]      */                                                                                        \
        /*  [-Jtr, Jrr]      */                                                                                        \
    } while (0)

#define NODE_BUILD_SMOOTHER_GIVE(i_r, i_theta, grid, DirBC_Interior, inner_boundary_circle_matrix,                     \
                                 circle_diagonal_solver, circle_tridiagonal_solver, radial_diagonal_solver,            \
                                 radial_tridiagonal_solver)                                                            \
    do {                                                                                                               \
        assert(i_r >= 0 && i_r < grid.nr());                                                                           \
        assert(i_theta >= 0 && i_theta < grid.ntheta());                                                               \
                                                                                                                       \
        const int numberSmootherCircles = grid.numberSmootherCircles();                                                \
        const int lengthSmootherRadial  = grid.lengthSmootherRadial();                                                 \
                                                                                                                       \
        assert(numberSmootherCircles >= 3);                                                                            \
        assert(lengthSmootherRadial >= 3);                                                                             \
                                                                                                                       \
        int row, column;                                                                                               \
        double value;                                                                                                  \
        /* ------------------------------------------ */                                                               \
        /* Node in the interior of the Circle Section */                                                               \
        /* ------------------------------------------ */                                                               \
        if (i_r > 0 && i_r < numberSmootherCircles - 1) {                                                              \
            double h1     = grid.radialSpacing(i_r - 1);                                                               \
            double h2     = grid.radialSpacing(i_r);                                                                   \
            double k1     = grid.angularSpacing(i_theta - 1);                                                          \
            double k2     = grid.angularSpacing(i_theta);                                                              \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                      \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                      \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                      \
                                                                                                                       \
            int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                         \
            int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                         \
                                                                                                                       \
            int center_index = i_theta;                                                                                \
            int left_index   = i_theta;                                                                                \
            int right_index  = i_theta;                                                                                \
            int bottom_index = i_theta_M1;                                                                             \
            int top_index    = i_theta_P1;                                                                             \
            /* -------------------------- */                                                                           \
            /* Cyclic Tridiagonal Section */                                                                           \
            /* i_r % 2 == 1               */                                                                           \
            if (i_r & 1) {                                                                                             \
                /* i_theta % 2 == 1 */                                                                                 \
                /* | X | O | X | */                                                                                    \
                /* |   |   |   | */                                                                                    \
                /* | 0 | Õ | O | */                                                                                    \
                /* |   |   |   | */                                                                                    \
                /* | X | O | X | */                                                                                    \
                /* or */                                                                                               \
                /* i_theta % 2 == 0 */                                                                                 \
                /* | O | O | O | */                                                                                    \
                /* |   |   |   | */                                                                                    \
                /* | X | Õ | X | */                                                                                    \
                /* |   |   |   | */                                                                                    \
                /* | O | O | O | */                                                                                    \
                                                                                                                       \
                auto& center_matrix = circle_tridiagonal_solver[i_r / 2];                                              \
                auto& left_matrix   = circle_diagonal_solver[(i_r - 1) / 2];                                           \
                auto& right_matrix  = circle_diagonal_solver[(i_r + 1) / 2];                                           \
                                                                                                                       \
                /* Fill matrix row of (i,j) */                                                                         \
                row    = center_index;                                                                                 \
                column = center_index;                                                                                 \
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */             \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                row    = center_index;                                                                                 \
                column = bottom_index;                                                                                 \
                value  = -coeff3 * att; /* Bottom */                                                                   \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                row    = center_index;                                                                                 \
                column = top_index;                                                                                    \
                value  = -coeff4 * att; /* Top */                                                                      \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                row    = center_index;                                                                                 \
                column = center_index;                                                                                 \
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */   \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                /* Fill matrix row of (i,j-1) */                                                                       \
                row    = bottom_index;                                                                                 \
                column = center_index;                                                                                 \
                value  = -coeff3 * att; /* Top */                                                                      \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                row    = bottom_index;                                                                                 \
                column = bottom_index;                                                                                 \
                value  = coeff3 * att; /* Center: (Top) */                                                             \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                /* Fill matrix row of (i,j+1) */                                                                       \
                row    = top_index;                                                                                    \
                column = center_index;                                                                                 \
                value  = -coeff4 * att; /* Bottom */                                                                   \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                row    = top_index;                                                                                    \
                column = top_index;                                                                                    \
                value  = coeff4 * att; /* Center: (Bottom) */                                                          \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                if (i_theta & 1) {                                                                                     \
                    /* i_theta % 2 == 1 */                                                                             \
                    /* | X | O | X | */                                                                                \
                    /* |   |   |   | */                                                                                \
                    /* | 0 | Õ | O | */                                                                                \
                    /* |   |   |   | */                                                                                \
                    /* | X | O | X | */                                                                                \
                                                                                                                       \
                    /* Fill matrix row of (i-1,j) */                                                                   \
                    if (i_r == 1) {                                                                                    \
                        /* Only in the case of AcrossOrigin */                                                         \
                        if (!DirBC_Interior) {                                                                         \
                            const Stencil& LeftStencil = getStencil(i_r - 1, i_theta);                                 \
                            int left_nz_index          = getCircleAscIndex(i_r - 1, i_theta);                          \
                            int nz_index               = left_nz_index + LeftStencil[StencilPosition::Center];         \
                            inner_boundary_circle_matrix.row_index(nz_index) = left_index;                         \
                            inner_boundary_circle_matrix.col_index(nz_index) = left_index;                         \
                            inner_boundary_circle_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) */        \
                        }                                                                                              \
                    }                                                                                                  \
                    else {                                                                                             \
                        row    = left_index;                                                                           \
                        column = left_index;                                                                           \
                        value  = coeff1 * arr; /* Center: (Right) */                                                   \
                        left_matrix.diagonal(row) += value;                                                            \
                    }                                                                                                  \
                                                                                                                       \
                    /* Fill matrix row of (i+1,j) */                                                                   \
                    row    = right_index;                                                                              \
                    column = right_index;                                                                              \
                    value  = coeff2 * arr; /* Center: (Left) */                                                        \
                    right_matrix.diagonal(row) += value;                                                               \
                }                                                                                                      \
            }                                                                                                          \
            /* ---------------- */                                                                                     \
            /* Diagonal Section */                                                                                     \
            /* i_r % 2 == 0     */                                                                                     \
            else {                                                                                                     \
                /* i_theta % 2 == 1 */                                                                                 \
                /* | O | X | O | */                                                                                    \
                /* |   |   |   | */                                                                                    \
                /* | O | Õ | O | */                                                                                    \
                /* |   |   |   | */                                                                                    \
                /* | O | X | O | */                                                                                    \
                /* or */                                                                                               \
                /* i_theta % 2 == 0 */                                                                                 \
                /* | O | O | O | */                                                                                    \
                /* |   |   |   | */                                                                                    \
                /* | O | X̃ | O | */                                                                                    \
                /* |   |   |   | */                                                                                    \
                /* | O | O | O | */                                                                                    \
                                                                                                                       \
                auto& center_matrix = circle_diagonal_solver[i_r / 2];                                                 \
                auto& left_matrix   = circle_tridiagonal_solver[(i_r - 1) / 2];                                        \
                auto& right_matrix  = circle_tridiagonal_solver[(i_r + 1) / 2];                                        \
                                                                                                                       \
                if (i_theta & 1) { /* i_theta % 2 == 1 */                                                              \
                    /* Fill matrix row of (i,j) */                                                                     \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */         \
                    center_matrix.diagonal(row) += value;                                                              \
                                                                                                                       \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value =                                                                                            \
                        (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */    \
                    center_matrix.diagonal(row) += value;                                                              \
                }                                                                                                      \
                else { /* i_theta % 2 == 0 */                                                                          \
                    /* Fill matrix row of (i,j) */                                                                     \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    center_matrix.diagonal(row) += 1.0;                                                                \
                                                                                                                       \
                    /* Fill matrix row of (i,j-1) */                                                                   \
                    row    = bottom_index;                                                                             \
                    column = bottom_index;                                                                             \
                    value  = coeff3 * att; /* Center: (Top) */                                                         \
                    center_matrix.diagonal(row) += value;                                                              \
                                                                                                                       \
                    /* Fill matrix row of (i,j+1) */                                                                   \
                    row    = top_index;                                                                                \
                    column = top_index;                                                                                \
                    value  = coeff4 * att; /* Center: (Bottom) */                                                      \
                    center_matrix.diagonal(row) += value;                                                              \
                }                                                                                                      \
                /* Fill matrix row of (i-1,j) */                                                                       \
                row    = left_index;                                                                                   \
                column = left_index;                                                                                   \
                value  = coeff1 * arr; /* Center: (Right) */                                                           \
                if (row == column)                                                                                     \
                    left_matrix.main_diagonal(row) += value;                                                           \
                else if (row == column - 1)                                                                            \
                    left_matrix.sub_diagonal(row) += value;                                                            \
                else if (row == 0 && column == left_matrix.columns() - 1)                                              \
                    left_matrix.cyclic_corner_element() += value;                                                      \
                                                                                                                       \
                /* Fill matrix row of (i+1,j) */                                                                       \
                row    = right_index;                                                                                  \
                column = right_index;                                                                                  \
                value  = coeff2 * arr; /* Center: (Left) */                                                            \
                if (row == column)                                                                                     \
                    right_matrix.main_diagonal(row) += value;                                                          \
                else if (row == column - 1)                                                                            \
                    right_matrix.sub_diagonal(row) += value;                                                           \
                else if (row == 0 && column == right_matrix.columns() - 1)                                             \
                    right_matrix.cyclic_corner_element() += value;                                                     \
            }                                                                                                          \
        }                                                                                                              \
        /* ------------------------------------------ */                                                               \
        /* Node in the interior of the Radial Section */                                                               \
        /* ------------------------------------------ */                                                               \
        else if (i_r > numberSmootherCircles && i_r < grid.nr() - 2) {                                                 \
            double h1     = grid.radialSpacing(i_r - 1);                                                               \
            double h2     = grid.radialSpacing(i_r);                                                                   \
            double k1     = grid.angularSpacing(i_theta - 1);                                                          \
            double k2     = grid.angularSpacing(i_theta);                                                              \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                      \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                      \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                      \
                                                                                                                       \
            int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                         \
            int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                         \
                                                                                                                       \
            int center_index = i_r - numberSmootherCircles;                                                            \
            int left_index   = i_r - numberSmootherCircles - 1;                                                        \
            int right_index  = i_r - numberSmootherCircles + 1;                                                        \
            int bottom_index = i_r - numberSmootherCircles;                                                            \
            int top_index    = i_r - numberSmootherCircles;                                                            \
            /* ------------------- */                                                                                  \
            /* Tridiagonal Section */                                                                                  \
            /* i_theta % 2 == 1    */                                                                                  \
            if (i_theta & 1) {                                                                                         \
                /* i_r % 2 == 1 */                                                                                     \
                /* ---------- */                                                                                       \
                /* X   O   X  */                                                                                       \
                /* ---------- */                                                                                       \
                /* O   Õ   O  */                                                                                       \
                /* ---------- */                                                                                       \
                /* X   O   X  */                                                                                       \
                /* ---------- */                                                                                       \
                /* or */                                                                                               \
                /* i_r % 2 == 0 */                                                                                     \
                /* ---------- */                                                                                       \
                /* O   X   O  */                                                                                       \
                /* ---------- */                                                                                       \
                /* O   Õ   O  */                                                                                       \
                /* ---------- */                                                                                       \
                /* O   X   O  */                                                                                       \
                /* ---------- */                                                                                       \
                                                                                                                       \
                auto& center_matrix = radial_tridiagonal_solver[i_theta / 2];                                          \
                auto& bottom_matrix = radial_diagonal_solver[i_theta_M1 / 2];                                          \
                auto& top_matrix    = radial_diagonal_solver[i_theta_P1 / 2];                                          \
                                                                                                                       \
                /* Fill matrix row of (i,j) */                                                                         \
                row    = center_index;                                                                                 \
                column = center_index;                                                                                 \
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */             \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                row    = center_index;                                                                                 \
                column = left_index;                                                                                   \
                value  = -coeff1 * arr; /* Left */                                                                     \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                row    = center_index;                                                                                 \
                column = right_index;                                                                                  \
                value  = -coeff2 * arr; /* Right */                                                                    \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                row    = center_index;                                                                                 \
                column = center_index;                                                                                 \
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */   \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                /* Fill matrix row of (i-1,j) */                                                                       \
                row    = left_index;                                                                                   \
                column = center_index;                                                                                 \
                value  = -coeff1 * arr; /* Right */                                                                    \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                row    = left_index;                                                                                   \
                column = left_index;                                                                                   \
                value  = coeff1 * arr; /* Center: (Right) */                                                           \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                /* Fill matrix row of (i+1,j) */                                                                       \
                row    = right_index;                                                                                  \
                column = center_index;                                                                                 \
                value  = -coeff2 * arr; /* Left */                                                                     \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                row    = right_index;                                                                                  \
                column = right_index;                                                                                  \
                value  = coeff2 * arr; /* Center: (Left) */                                                            \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                if (i_r & 1) { /* i_r % 2 == 1 */                                                                      \
                    /* ---------- */                                                                                   \
                    /* X   O   X  */                                                                                   \
                    /* ---------- */                                                                                   \
                    /* O   Õ   O  */                                                                                   \
                    /* ---------- */                                                                                   \
                    /* X   O   X  */                                                                                   \
                    /* ---------- */                                                                                   \
                    /* Fill matrix row of (i,j-1) */                                                                   \
                    row    = bottom_index;                                                                             \
                    column = bottom_index;                                                                             \
                    value  = coeff3 * att; /* Center: (Top) */                                                         \
                    bottom_matrix.diagonal(row) += value;                                                              \
                                                                                                                       \
                    /* Fill matrix row of (i,j+1) */                                                                   \
                    row    = top_index;                                                                                \
                    column = top_index;                                                                                \
                    value  = coeff4 * att; /* Center: (Bottom) */                                                      \
                    top_matrix.diagonal(row) += value;                                                                 \
                }                                                                                                      \
            }                                                                                                          \
            /* ---------------- */                                                                                     \
            /* Diagonal Section */                                                                                     \
            /* i_theta % 2 == 0 */                                                                                     \
            else {                                                                                                     \
                /* i_r % 2 == 1 */                                                                                     \
                /* ---------- */                                                                                       \
                /* O   O   O  */                                                                                       \
                /* ---------- */                                                                                       \
                /* X   Õ   X  */                                                                                       \
                /* ---------- */                                                                                       \
                /* O   O   O  */                                                                                       \
                /* ---------- */                                                                                       \
                /* or */                                                                                               \
                /* i_r % 2 == 0 */                                                                                     \
                /* ---------- */                                                                                       \
                /* O   O   O  */                                                                                       \
                /* ---------- */                                                                                       \
                /* O   X̃   O  */                                                                                       \
                /* ---------- */                                                                                       \
                /* O   O   O  */                                                                                       \
                /* ---------- */                                                                                       \
                                                                                                                       \
                auto& center_matrix = radial_diagonal_solver[i_theta / 2];                                             \
                auto& bottom_matrix = radial_tridiagonal_solver[i_theta_M1 / 2];                                       \
                auto& top_matrix    = radial_tridiagonal_solver[i_theta_P1 / 2];                                       \
                if (i_r & 1) { /* i_r % 2 == 1 */                                                                      \
                    /* Fill matrix row of (i,j) */                                                                     \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */         \
                    center_matrix.diagonal(row) += value;                                                              \
                                                                                                                       \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value =                                                                                            \
                        (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */    \
                    center_matrix.diagonal(row) += value;                                                              \
                }                                                                                                      \
                else { /* i_r % 2 == 0 */                                                                              \
                    /* Fill matrix row of (i,j) */                                                                     \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    center_matrix.diagonal(row) += 1.0;                                                                \
                                                                                                                       \
                    /* Fill matrix row of (i-1,j) */                                                                   \
                    row    = left_index;                                                                               \
                    column = left_index;                                                                               \
                    value  = coeff1 * arr; /* Center: (Right) */                                                       \
                    center_matrix.diagonal(row) += value;                                                              \
                                                                                                                       \
                    /* Fill matrix row of (i+1,j) */                                                                   \
                    row    = right_index;                                                                              \
                    column = right_index;                                                                              \
                    value  = coeff2 * arr; /* Center: (Left) */                                                        \
                    center_matrix.diagonal(row) += value;                                                              \
                }                                                                                                      \
                /* Fill matrix row of (i,j-1) */                                                                       \
                row    = bottom_index;                                                                                 \
                column = bottom_index;                                                                                 \
                value  = coeff3 * att; /* Center: (Top) */                                                             \
                if (row == column)                                                                                     \
                    bottom_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    bottom_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == bottom_matrix.columns() - 1)                                            \
                    bottom_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                /* Fill matrix row of (i,j+1) */                                                                       \
                row    = top_index;                                                                                    \
                column = top_index;                                                                                    \
                value  = coeff4 * att; /* Center: (Bottom) */                                                          \
                if (row == column)                                                                                     \
                    top_matrix.main_diagonal(row) += value;                                                            \
                else if (row == column - 1)                                                                            \
                    top_matrix.sub_diagonal(row) += value;                                                             \
                else if (row == 0 && column == top_matrix.columns() - 1)                                               \
                    top_matrix.cyclic_corner_element() += value;                                                       \
            }                                                                                                          \
        }                                                                                                              \
        /* ------------------------------------------ */                                                               \
        /* Circle Section: Node in the inner boundary */                                                               \
        /* ------------------------------------------ */                                                               \
        else if (i_r == 0) {                                                                                           \
            /* ------------------------------------------------ */                                                     \
            /* Case 1: Dirichlet boundary on the inner boundary */                                                     \
            /* ------------------------------------------------ */                                                     \
            if (DirBC_Interior) {                                                                                      \
                /* Fill result(i,j) */                                                                                 \
                double h2     = grid.radialSpacing(i_r);                                                               \
                double k1     = grid.angularSpacing(i_theta - 1);                                                      \
                double k2     = grid.angularSpacing(i_theta);                                                          \
                double coeff2 = 0.5 * (k1 + k2) / h2;                                                                  \
                                                                                                                       \
                int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                     \
                int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                     \
                                                                                                                       \
                auto& center_matrix = inner_boundary_circle_matrix;                                                    \
                auto& right_matrix  = circle_tridiagonal_solver[(i_r + 1) / 2];                                        \
                                                                                                                       \
                int center_index = i_theta;                                                                            \
                int right_index  = i_theta;                                                                            \
                int bottom_index = i_theta_M1;                                                                         \
                int top_index    = i_theta_P1;                                                                         \
                                                                                                                       \
                /* Fill matrix row of (i,j) */                                                                         \
                const Stencil& CenterStencil      = getStencil(i_r, i_theta);                                          \
                int center_nz_index               = getCircleAscIndex(i_r, i_theta);                                   \
                int nz_index                      = center_nz_index + CenterStencil[StencilPosition::Center];          \
                center_matrix.row_index(nz_index) = center_index;                                                  \
                center_matrix.col_index(nz_index) = center_index;                                                  \
                center_matrix.value(nz_index) += 1.0;                                                                  \
                                                                                                                       \
                /* Fill matrix row of (i+1,j) */                                                                       \
                row    = right_index;                                                                                  \
                column = right_index;                                                                                  \
                value  = coeff2 * arr; /* Center: (Left) */                                                            \
                if (row == column)                                                                                     \
                    right_matrix.main_diagonal(row) += value;                                                          \
                else if (row == column - 1)                                                                            \
                    right_matrix.sub_diagonal(row) += value;                                                           \
                else if (row == 0 && column == right_matrix.columns() - 1)                                             \
                    right_matrix.cyclic_corner_element() += value;                                                     \
            }                                                                                                          \
            else {                                                                                                     \
                /* ------------------------------------------------------------- */                                    \
                /* Case 2: Across origin discretization on the interior boundary */                                    \
                /* ------------------------------------------------------------- */                                    \
                /* h1 gets replaced with 2 * R0. */                                                                    \
                /* (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()/2)). */                           \
                /* Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil. */    \
                double h1     = 2.0 * grid.radius(0);                                                                  \
                double h2     = grid.radialSpacing(i_r);                                                               \
                double k1     = grid.angularSpacing(i_theta - 1);                                                      \
                double k2     = grid.angularSpacing(i_theta);                                                          \
                double coeff1 = 0.5 * (k1 + k2) / h1;                                                                  \
                double coeff2 = 0.5 * (k1 + k2) / h2;                                                                  \
                double coeff3 = 0.5 * (h1 + h2) / k1;                                                                  \
                double coeff4 = 0.5 * (h1 + h2) / k2;                                                                  \
                                                                                                                       \
                const int i_theta_M1           = grid.wrapThetaIndex(i_theta - 1);                                     \
                const int i_theta_P1           = grid.wrapThetaIndex(i_theta + 1);                                     \
                const int i_theta_AcrossOrigin = grid.wrapThetaIndex(i_theta + (grid.ntheta() / 2));                   \
                                                                                                                       \
                const int center_index = i_theta;                                                                      \
                const int left_index   = i_theta_AcrossOrigin;                                                         \
                const int right_index  = i_theta;                                                                      \
                const int bottom_index = i_theta_M1;                                                                   \
                const int top_index    = i_theta_P1;                                                                   \
                                                                                                                       \
                const int center_nz_index = getCircleAscIndex(i_r, i_theta);                                           \
                const int bottom_nz_index = getCircleAscIndex(i_r, i_theta_M1);                                        \
                const int top_nz_index    = getCircleAscIndex(i_r, i_theta_P1);                                        \
                const int left_nz_index   = getCircleAscIndex(i_r, i_theta_AcrossOrigin);                              \
                                                                                                                       \
                int nz_index;                                                                                          \
                const Stencil& CenterStencil = getStencil(i_r, i_theta);                                               \
                                                                                                                       \
                if (i_theta & 1) {                                                                                     \
                    /* i_theta % 2 == 1 */                                                                             \
                    /* -| X | O | X | */                                                                               \
                    /* -|   |   |   | */                                                                               \
                    /* -| Õ | O | O | */                                                                               \
                    /* -|   |   |   | */                                                                               \
                    /* -| X | O | X | */                                                                               \
                                                                                                                       \
                    auto& center_matrix = inner_boundary_circle_matrix;                                                \
                    auto& right_matrix  = circle_tridiagonal_solver[(i_r + 1) / 2];                                    \
                    auto& left_matrix   = inner_boundary_circle_matrix;                                                \
                    /* Fill matrix row of (i,j) */                                                                     \
                    nz_index                          = center_nz_index + CenterStencil[StencilPosition::Center];      \
                    center_matrix.row_index(nz_index) = center_index;                                              \
                    center_matrix.col_index(nz_index) = center_index;                                              \
                    center_matrix.value(nz_index) +=                                                                   \
                        0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */                      \
                                                                                                                       \
                    nz_index                          = center_nz_index + CenterStencil[StencilPosition::Left];        \
                    center_matrix.row_index(nz_index) = center_index;                                              \
                    center_matrix.col_index(nz_index) = left_index;                                                \
                    center_matrix.value(nz_index) += -coeff1 * arr; /* Left */                                         \
                                                                                                                       \
                    nz_index                          = center_nz_index + CenterStencil[StencilPosition::Center];      \
                    center_matrix.row_index(nz_index) = center_index;                                              \
                    center_matrix.col_index(nz_index) = center_index;                                              \
                    /* Center: (Left, Right, Bottom, Top) */                                                           \
                    center_matrix.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att;                \
                                                                                                                       \
                    /* Fill matrix row of (i-1,j) */                                                                   \
                    /* From view the view of the across origin node, */                                                \
                    /* the directions are roatated by 180 degrees in the stencil! */                                   \
                    const Stencil& LeftStencil = CenterStencil;                                                        \
                                                                                                                       \
                    nz_index                        = left_nz_index + LeftStencil[StencilPosition::Left];              \
                    left_matrix.row_index(nz_index) = left_index;                                                  \
                    left_matrix.col_index(nz_index) = center_index;                                                \
                    left_matrix.value(nz_index) += -coeff1 * arr; /* Right -> Left*/                                   \
                                                                                                                       \
                    nz_index                        = left_nz_index + LeftStencil[StencilPosition::Center];            \
                    left_matrix.row_index(nz_index) = left_index;                                                  \
                    left_matrix.col_index(nz_index) = left_index;                                                  \
                    left_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) -> Center: (Left) */               \
                                                                                                                       \
                    /* Fill matrix row of (i+1,j) */                                                                   \
                    row    = right_index;                                                                              \
                    column = right_index;                                                                              \
                    value  = coeff2 * arr; /* Center: (Left) */                                                        \
                    if (row == column)                                                                                 \
                        right_matrix.main_diagonal(row) += value;                                                      \
                    else if (row == column - 1)                                                                        \
                        right_matrix.sub_diagonal(row) += value;                                                       \
                    else if (row == 0 && column == right_matrix.columns() - 1)                                         \
                        right_matrix.cyclic_corner_element() += value;                                                 \
                }                                                                                                      \
                else {                                                                                                 \
                    /* i_theta % 2 == 0 */                                                                             \
                    /* -| O | O | O | */                                                                               \
                    /* -|   |   |   | */                                                                               \
                    /* -| X̃ | O | X | */                                                                               \
                    /* -|   |   |   | */                                                                               \
                    /* -| O | O | O | */                                                                               \
                                                                                                                       \
                    auto& center_matrix = inner_boundary_circle_matrix;                                                \
                    auto& right_matrix  = circle_tridiagonal_solver[(i_r + 1) / 2];                                    \
                    auto& left_matrix   = inner_boundary_circle_matrix;                                                \
                    /* Fill matrix row of (i,j) */                                                                     \
                    nz_index                          = center_nz_index + CenterStencil[StencilPosition::Center];      \
                    center_matrix.row_index(nz_index) = center_index;                                              \
                    center_matrix.col_index(nz_index) = center_index;                                              \
                    center_matrix.value(nz_index) += 1.0;                                                              \
                    /* Fill matrix row of (i,j-1) */                                                                   \
                    const Stencil& BottomStencil = CenterStencil;                                                      \
                                                                                                                       \
                    nz_index                          = bottom_nz_index + BottomStencil[StencilPosition::Center];      \
                    center_matrix.row_index(nz_index) = bottom_index;                                              \
                    center_matrix.col_index(nz_index) = bottom_index;                                              \
                    center_matrix.value(nz_index) += coeff3 * att; /* Center: (Top) */                                 \
                                                                                                                       \
                    /* Fill matrix row of (i,j+1) */                                                                   \
                    const Stencil& TopStencil = CenterStencil;                                                         \
                                                                                                                       \
                    nz_index                          = top_nz_index + TopStencil[StencilPosition::Center];            \
                    center_matrix.row_index(nz_index) = top_index;                                                 \
                    center_matrix.col_index(nz_index) = top_index;                                                 \
                    center_matrix.value(nz_index) += coeff4 * att; /* Center: (Bottom) */                              \
                                                                                                                       \
                    /* Fill matrix row of (i+1,j) */                                                                   \
                    row    = right_index;                                                                              \
                    column = right_index;                                                                              \
                    value  = coeff2 * arr; /* Center: (Left) */                                                        \
                    if (row == column)                                                                                 \
                        right_matrix.main_diagonal(row) += value;                                                      \
                    else if (row == column - 1)                                                                        \
                        right_matrix.sub_diagonal(row) += value;                                                       \
                    else if (row == 0 && column == right_matrix.columns() - 1)                                         \
                        right_matrix.cyclic_corner_element() += value;                                                 \
                }                                                                                                      \
            }                                                                                                          \
        }                                                                                                              \
        /* ------------------------------------------- */                                                              \
        /* Circle Section: Node next to radial section */                                                              \
        /* ------------------------------------------- */                                                              \
        else if (i_r == numberSmootherCircles - 1) {                                                                   \
            assert(i_r > 1);                                                                                           \
                                                                                                                       \
            double h1     = grid.radialSpacing(i_r - 1);                                                               \
            double h2     = grid.radialSpacing(i_r);                                                                   \
            double k1     = grid.angularSpacing(i_theta - 1);                                                          \
            double k2     = grid.angularSpacing(i_theta);                                                              \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                      \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                      \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                      \
                                                                                                                       \
            int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                         \
            int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                         \
                                                                                                                       \
            int center_index = i_theta;                                                                                \
            int left_index   = i_theta;                                                                                \
            int right_index  = 0;                                                                                      \
            int bottom_index = i_theta_M1;                                                                             \
            int top_index    = i_theta_P1;                                                                             \
                                                                                                                       \
            if (i_r & 1) {                                                                                             \
                if (i_theta & 1) {                                                                                     \
                    /* i_r % 2 == 1 and i_theta % 2 == 1 */                                                            \
                    /* | O | X | O || X   O   X   O  */                                                                \
                    /* |   |   |   || -------------- */                                                                \
                    /* | 0 | O | Õ || O   O   O   O  */                                                                \
                    /* |   |   |   || -------------- */                                                                \
                    /* | O | X | O || X   O   X   O  */                                                                \
                                                                                                                       \
                    auto& center_matrix = circle_tridiagonal_solver[i_r / 2];                                          \
                    auto& left_matrix   = circle_diagonal_solver[(i_r - 1) / 2];                                       \
                    auto& right_matrix  = radial_tridiagonal_solver[i_theta / 2];                                      \
                                                                                                                       \
                    /* Fill matrix row of (i,j) */                                                                     \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */         \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    row    = center_index;                                                                             \
                    column = bottom_index;                                                                             \
                    value  = -coeff3 * att; /* Bottom */                                                               \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    row    = center_index;                                                                             \
                    column = top_index;                                                                                \
                    value  = -coeff4 * att; /* Top */                                                                  \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value =                                                                                            \
                        (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */    \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    /* Fill matrix row of (i,j-1) */                                                                   \
                    row    = bottom_index;                                                                             \
                    column = center_index;                                                                             \
                    value  = -coeff3 * att; /* Top */                                                                  \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    row    = bottom_index;                                                                             \
                    column = bottom_index;                                                                             \
                    value  = coeff3 * att; /* Center: (Top) */                                                         \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    /* Fill matrix row of (i,j+1) */                                                                   \
                    row    = top_index;                                                                                \
                    column = center_index;                                                                             \
                    value  = -coeff4 * att; /* Bottom */                                                               \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    row    = top_index;                                                                                \
                    column = top_index;                                                                                \
                    value  = coeff4 * att; /* Center: (Bottom) */                                                      \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    /* Fill matrix row of (i-1,j) */                                                                   \
                    row    = left_index;                                                                               \
                    column = left_index;                                                                               \
                    value  = coeff1 * arr; /* Center: (Right) */                                                       \
                    left_matrix.diagonal(row) += value;                                                                \
                                                                                                                       \
                    /* Fill matrix row of (i+1,j) */                                                                   \
                    row    = right_index;                                                                              \
                    column = right_index;                                                                              \
                    value  = coeff2 * arr; /* Center: (Left) */                                                        \
                    if (row == column)                                                                                 \
                        right_matrix.main_diagonal(row) += value;                                                      \
                    else if (row == column - 1)                                                                        \
                        right_matrix.sub_diagonal(row) += value;                                                       \
                    else if (row == 0 && column == right_matrix.columns() - 1)                                         \
                        right_matrix.cyclic_corner_element() += value;                                                 \
                }                                                                                                      \
                else {                                                                                                 \
                    /* i_r % 2 == 1 and i_theta % 2 == 0 */                                                            \
                    /* | O | O | O || O   O   O   O  */                                                                \
                    /* |   |   |   || -------------- */                                                                \
                    /* | 0 | X | Õ || X   O   X   O  */                                                                \
                    /* |   |   |   || -------------- */                                                                \
                    /* | O | O | O || O   O   O   O  */                                                                \
                                                                                                                       \
                    auto& center_matrix = circle_tridiagonal_solver[i_r / 2];                                          \
                    auto& left_matrix   = circle_diagonal_solver[(i_r - 1) / 2];                                       \
                    auto& right_matrix  = radial_diagonal_solver[i_theta / 2];                                         \
                                                                                                                       \
                    /* Fill matrix row of (i,j) */                                                                     \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */         \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    row    = center_index;                                                                             \
                    column = bottom_index;                                                                             \
                    value  = -coeff3 * att; /* Bottom */                                                               \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    row    = center_index;                                                                             \
                    column = top_index;                                                                                \
                    value  = -coeff4 * att; /* Top */                                                                  \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value =                                                                                            \
                        (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */    \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    /* Fill matrix row of (i,j-1) */                                                                   \
                    row    = bottom_index;                                                                             \
                    column = center_index;                                                                             \
                    value  = -coeff3 * att; /* Top */                                                                  \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    row    = bottom_index;                                                                             \
                    column = bottom_index;                                                                             \
                    value  = coeff3 * att; /* Center: (Top) */                                                         \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    /* Fill matrix row of (i,j+1) */                                                                   \
                    row    = top_index;                                                                                \
                    column = center_index;                                                                             \
                    value  = -coeff4 * att; /* Bottom */                                                               \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    row    = top_index;                                                                                \
                    column = top_index;                                                                                \
                    value  = coeff4 * att; /* Center: (Bottom) */                                                      \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                }                                                                                                      \
            }                                                                                                          \
            else {                                                                                                     \
                if (i_theta & 1) {                                                                                     \
                    /* i_r % 2 == 0 and i_theta % 2 == 1 */                                                            \
                    /* | X | O | X || O   X   O   X  */                                                                \
                    /* |   |   |   || -------------- */                                                                \
                    /* | 0 | O | Õ || O   O   O   O  */                                                                \
                    /* |   |   |   || -------------- */                                                                \
                    /* | X | O | X || O   X   O   X  */                                                                \
                                                                                                                       \
                    auto& center_matrix = circle_diagonal_solver[i_r / 2];                                             \
                    auto& left_matrix   = circle_tridiagonal_solver[(i_r - 1) / 2];                                    \
                    auto& right_matrix  = radial_tridiagonal_solver[i_theta / 2];                                      \
                                                                                                                       \
                    /* Fill matrix row of (i,j) */                                                                     \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */         \
                    center_matrix.diagonal(row) += value;                                                              \
                                                                                                                       \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value =                                                                                            \
                        (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */    \
                    center_matrix.diagonal(row) += value;                                                              \
                                                                                                                       \
                    /* Fill matrix row of (i-1,j) */                                                                   \
                    row    = left_index;                                                                               \
                    column = left_index;                                                                               \
                    value  = coeff1 * arr; /* Center: (Right) */                                                       \
                    if (row == column)                                                                                 \
                        left_matrix.main_diagonal(row) += value;                                                       \
                    else if (row == column - 1)                                                                        \
                        left_matrix.sub_diagonal(row) += value;                                                        \
                    else if (row == 0 && column == left_matrix.columns() - 1)                                          \
                        left_matrix.cyclic_corner_element() += value;                                                  \
                                                                                                                       \
                    /* Fill matrix row of (i+1,j) */                                                                   \
                    row    = right_index;                                                                              \
                    column = right_index;                                                                              \
                    value  = coeff2 * arr; /* Center: (Left) */                                                        \
                    if (row == column)                                                                                 \
                        right_matrix.main_diagonal(row) += value;                                                      \
                    else if (row == column - 1)                                                                        \
                        right_matrix.sub_diagonal(row) += value;                                                       \
                    else if (row == 0 && column == right_matrix.columns() - 1)                                         \
                        right_matrix.cyclic_corner_element() += value;                                                 \
                }                                                                                                      \
                else {                                                                                                 \
                    /* i_r % 2 == 0 and i_theta % 2 == 0 */                                                            \
                    /* | O | O | O || O   O   O   O  */                                                                \
                    /* |   |   |   || -------------- */                                                                \
                    /* | X | O | X̃ || O   X   O   X  */                                                                \
                    /* |   |   |   || -------------- */                                                                \
                    /* | O | O | O || O   O   O   O  */                                                                \
                                                                                                                       \
                    auto& center_matrix = circle_diagonal_solver[i_r / 2];                                             \
                    auto& left_matrix   = circle_tridiagonal_solver[(i_r - 1) / 2];                                    \
                    auto& right_matrix  = radial_diagonal_solver[i_theta / 2];                                         \
                                                                                                                       \
                    /* Fill matrix row of (i,j) */                                                                     \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    center_matrix.diagonal(row) += 1.0;                                                                \
                                                                                                                       \
                    /* Fill matrix row of (i,j-1) */                                                                   \
                    row    = bottom_index;                                                                             \
                    column = bottom_index;                                                                             \
                    value  = coeff3 * att; /* Center: (Top) */                                                         \
                    center_matrix.diagonal(row) += value;                                                              \
                                                                                                                       \
                    /* Fill matrix row of (i,j+1) */                                                                   \
                    row    = top_index;                                                                                \
                    column = top_index;                                                                                \
                    value  = coeff4 * att; /* Center: (Bottom) */                                                      \
                    center_matrix.diagonal(row) += value;                                                              \
                                                                                                                       \
                    /* Fill matrix row of (i-1,j) */                                                                   \
                    row    = left_index;                                                                               \
                    column = left_index;                                                                               \
                    value  = coeff1 * arr; /* Center: (Right) */                                                       \
                    if (row == column)                                                                                 \
                        left_matrix.main_diagonal(row) += value;                                                       \
                    else if (row == column - 1)                                                                        \
                        left_matrix.sub_diagonal(row) += value;                                                        \
                    else if (row == 0 && column == left_matrix.columns() - 1)                                          \
                        left_matrix.cyclic_corner_element() += value;                                                  \
                                                                                                                       \
                    /* Fill matrix row of (i+1,j) */                                                                   \
                    row    = right_index;                                                                              \
                    column = right_index;                                                                              \
                    value  = coeff2 * arr; /* Center: (Left) */                                                        \
                    right_matrix.diagonal(row) += value;                                                               \
                }                                                                                                      \
            }                                                                                                          \
        }                                                                                                              \
        /* --------------------------------------------- */                                                            \
        /* Radial Section: Node next to circular section */                                                            \
        /* --------------------------------------------- */                                                            \
        else if (i_r == numberSmootherCircles) {                                                                       \
            double h1     = grid.radialSpacing(i_r - 1);                                                               \
            double h2     = grid.radialSpacing(i_r);                                                                   \
            double k1     = grid.angularSpacing(i_theta - 1);                                                          \
            double k2     = grid.angularSpacing(i_theta);                                                              \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                      \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                      \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                      \
                                                                                                                       \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                   \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                   \
                                                                                                                       \
            const int center_index = i_r - numberSmootherCircles;                                                      \
            const int left_index   = i_theta;                                                                          \
            const int right_index  = i_r - numberSmootherCircles + 1;                                                  \
            const int bottom_index = i_r - numberSmootherCircles;                                                      \
            const int top_index    = i_r - numberSmootherCircles;                                                      \
                                                                                                                       \
            if (i_theta & 1) {                                                                                         \
                if (i_r & 1) {                                                                                         \
                    /* i_theta % 2 == 1 and i_r % 2 == 1 */                                                            \
                    /* | X | O | X || O   X   O   X  */                                                                \
                    /* |   |   |   || -------------- */                                                                \
                    /* | 0 | O | O || Õ   O   O   O  */                                                                \
                    /* |   |   |   || -------------- */                                                                \
                    /* | X | O | X || O   X   O   X  */                                                                \
                                                                                                                       \
                    auto& center_matrix = radial_tridiagonal_solver[i_theta / 2];                                      \
                    auto& bottom_matrix = radial_diagonal_solver[i_theta_M1 / 2];                                      \
                    auto& top_matrix    = radial_diagonal_solver[i_theta_P1 / 2];                                      \
                    auto& left_matrix   = circle_diagonal_solver[(i_r - 1) / 2];                                       \
                                                                                                                       \
                    /* Fill matrix row of (i,j) */                                                                     \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */         \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    row    = center_index;                                                                             \
                    column = right_index;                                                                              \
                    value  = -coeff2 * arr; /* Right */                                                                \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value =                                                                                            \
                        (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */    \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    /* Fill matrix row of (i-1,j) */                                                                   \
                    row    = left_index;                                                                               \
                    column = left_index;                                                                               \
                    value  = coeff1 * arr; /* Center: (Right) */                                                       \
                    left_matrix.diagonal(row) += value;                                                                \
                                                                                                                       \
                    /* Fill matrix row of (i+1,j) */                                                                   \
                    row    = right_index;                                                                              \
                    column = center_index;                                                                             \
                    value  = -coeff2 * arr; /* Left */                                                                 \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    row    = right_index;                                                                              \
                    column = right_index;                                                                              \
                    value  = coeff2 * arr; /* Center: (Left) */                                                        \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    /* Fill matrix row of (i,j-1) */                                                                   \
                    row    = bottom_index;                                                                             \
                    column = bottom_index;                                                                             \
                    value  = coeff3 * att; /* Center: (Top) */                                                         \
                    bottom_matrix.diagonal(row) += value;                                                              \
                                                                                                                       \
                    /* Fill matrix row of (i,j+1) */                                                                   \
                    row    = top_index;                                                                                \
                    column = top_index;                                                                                \
                    value  = coeff4 * att; /* Center: (Bottom) */                                                      \
                    top_matrix.diagonal(row) += value;                                                                 \
                }                                                                                                      \
                else {                                                                                                 \
                    /* i_theta % 2 == 1 and i_r % 2 == 0 */                                                            \
                    /* | O | X | O || X   O   X   O  */                                                                \
                    /* |   |   |   || -------------- */                                                                \
                    /* | 0 | O | O || Õ   O   O   O  */                                                                \
                    /* |   |   |   || -------------- */                                                                \
                    /* | O | X | O || X   O   X   O  */                                                                \
                                                                                                                       \
                    auto& center_matrix = radial_tridiagonal_solver[i_theta / 2];                                      \
                    auto& bottom_matrix = radial_diagonal_solver[i_theta_M1 / 2];                                      \
                    auto& top_matrix    = radial_diagonal_solver[i_theta_P1 / 2];                                      \
                    auto& left_matrix   = circle_tridiagonal_solver[(i_r - 1) / 2];                                    \
                                                                                                                       \
                    /* Fill matrix row of (i,j) */                                                                     \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */         \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    row    = center_index;                                                                             \
                    column = right_index;                                                                              \
                    value  = -coeff2 * arr; /* Right */                                                                \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value =                                                                                            \
                        (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */    \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    /* Fill matrix row of (i-1,j) */                                                                   \
                    row    = left_index;                                                                               \
                    column = left_index;                                                                               \
                    value  = coeff1 * arr; /* Center: (Right) */                                                       \
                    if (row == column)                                                                                 \
                        left_matrix.main_diagonal(row) += value;                                                       \
                    else if (row == column - 1)                                                                        \
                        left_matrix.sub_diagonal(row) += value;                                                        \
                    else if (row == 0 && column == left_matrix.columns() - 1)                                          \
                        left_matrix.cyclic_corner_element() += value;                                                  \
                                                                                                                       \
                    /* Fill matrix row of (i+1,j) */                                                                   \
                    row    = right_index;                                                                              \
                    column = center_index;                                                                             \
                    value  = -coeff2 * arr; /* Left */                                                                 \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    row    = right_index;                                                                              \
                    column = right_index;                                                                              \
                    value  = coeff2 * arr; /* Center: (Left) */                                                        \
                    if (row == column)                                                                                 \
                        center_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        center_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == center_matrix.columns() - 1)                                        \
                        center_matrix.cyclic_corner_element() += value;                                                \
                }                                                                                                      \
            }                                                                                                          \
            else {                                                                                                     \
                if (i_r & 1) {                                                                                         \
                    /* i_theta % 2 == 0 and i_r % 2 == 1 */                                                            \
                    /* | O | O | O || O   O   O   O  */                                                                \
                    /* |   |   |   || -------------- */                                                                \
                    /* | X | O | X || Õ   X   O   X  */                                                                \
                    /* |   |   |   || -------------- */                                                                \
                    /* | O | O | O || O   O   O   O  */                                                                \
                                                                                                                       \
                    auto& center_matrix = radial_diagonal_solver[i_theta / 2];                                         \
                    auto& bottom_matrix = radial_tridiagonal_solver[i_theta_M1 / 2];                                   \
                    auto& top_matrix    = radial_tridiagonal_solver[i_theta_P1 / 2];                                   \
                    auto& left_matrix   = circle_diagonal_solver[(i_r - 1) / 2];                                       \
                                                                                                                       \
                    /* Fill matrix row of (i,j) */                                                                     \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */         \
                    center_matrix.diagonal(row) += value;                                                              \
                                                                                                                       \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value =                                                                                            \
                        (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */    \
                    center_matrix.diagonal(row) += value;                                                              \
                                                                                                                       \
                    /* Fill matrix row of (i,j-1) */                                                                   \
                    row    = bottom_index;                                                                             \
                    column = bottom_index;                                                                             \
                    value  = coeff3 * att; /* Center: (Top) */                                                         \
                    if (row == column)                                                                                 \
                        bottom_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        bottom_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == bottom_matrix.columns() - 1)                                        \
                        bottom_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    /* Fill matrix row of (i,j+1) */                                                                   \
                    row    = top_index;                                                                                \
                    column = top_index;                                                                                \
                    value  = coeff4 * att; /* Center: (Bottom) */                                                      \
                    if (row == column)                                                                                 \
                        top_matrix.main_diagonal(row) += value;                                                        \
                    else if (row == column - 1)                                                                        \
                        top_matrix.sub_diagonal(row) += value;                                                         \
                    else if (row == 0 && column == top_matrix.columns() - 1)                                           \
                        top_matrix.cyclic_corner_element() += value;                                                   \
                }                                                                                                      \
                else {                                                                                                 \
                    /* i_theta % 2 == 0 and i_r % 2 == 0 */                                                            \
                    /* | O | O | O || O   O   O   O  */                                                                \
                    /* |   |   |   || -------------- */                                                                \
                    /* | O | X | O || X̃   O   X   O  */                                                                \
                    /* |   |   |   || -------------- */                                                                \
                    /* | O | O | O || O   O   O   O  */                                                                \
                                                                                                                       \
                    auto& center_matrix = radial_diagonal_solver[i_theta / 2];                                         \
                    auto& bottom_matrix = radial_tridiagonal_solver[i_theta_M1 / 2];                                   \
                    auto& top_matrix    = radial_tridiagonal_solver[i_theta_P1 / 2];                                   \
                    auto& left_matrix   = circle_tridiagonal_solver[(i_r - 1) / 2];                                    \
                                                                                                                       \
                    /* Fill matrix row of (i,j) */                                                                     \
                                                                                                                       \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    center_matrix.diagonal(row) += 1.0;                                                                \
                                                                                                                       \
                    /* Fill matrix row of (i-1,j) */                                                                   \
                    row    = left_index;                                                                               \
                    column = left_index;                                                                               \
                    value  = coeff1 * arr; /* Center: (Right) */                                                       \
                    if (row == column)                                                                                 \
                        left_matrix.main_diagonal(row) += value;                                                       \
                    else if (row == column - 1)                                                                        \
                        left_matrix.sub_diagonal(row) += value;                                                        \
                    else if (row == 0 && column == left_matrix.columns() - 1)                                          \
                        left_matrix.cyclic_corner_element() += value;                                                  \
                                                                                                                       \
                    /* Fill matrix row of (i+1,j) */                                                                   \
                    row    = right_index;                                                                              \
                    column = right_index;                                                                              \
                    value  = coeff2 * arr; /* Center: (Left) */                                                        \
                    center_matrix.diagonal(row) += value;                                                              \
                                                                                                                       \
                    /* Fill matrix row of (i,j-1) */                                                                   \
                    row    = bottom_index;                                                                             \
                    column = bottom_index;                                                                             \
                    value  = coeff3 * att; /* Center: (Top) */                                                         \
                    if (row == column)                                                                                 \
                        bottom_matrix.main_diagonal(row) += value;                                                     \
                    else if (row == column - 1)                                                                        \
                        bottom_matrix.sub_diagonal(row) += value;                                                      \
                    else if (row == 0 && column == bottom_matrix.columns() - 1)                                        \
                        bottom_matrix.cyclic_corner_element() += value;                                                \
                                                                                                                       \
                    /* Fill matrix row of (i,j+1) */                                                                   \
                    row    = top_index;                                                                                \
                    column = top_index;                                                                                \
                    value  = coeff4 * att; /* Center: (Bottom) */                                                      \
                    if (row == column)                                                                                 \
                        top_matrix.main_diagonal(row) += value;                                                        \
                    else if (row == column - 1)                                                                        \
                        top_matrix.sub_diagonal(row) += value;                                                         \
                    else if (row == 0 && column == top_matrix.columns() - 1)                                           \
                        top_matrix.cyclic_corner_element() += value;                                                   \
                }                                                                                                      \
            }                                                                                                          \
        }                                                                                                              \
        /* ------------------------------------------- */                                                              \
        /* Radial Section: Node next to outer boundary */                                                              \
        /* ------------------------------------------- */                                                              \
        else if (i_r == grid.nr() - 2) {                                                                               \
            assert(i_r % 2 == 1);                                                                                      \
                                                                                                                       \
            double h1     = grid.radialSpacing(i_r - 1);                                                               \
            double h2     = grid.radialSpacing(i_r);                                                                   \
            double k1     = grid.angularSpacing(i_theta - 1);                                                          \
            double k2     = grid.angularSpacing(i_theta);                                                              \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                      \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                      \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                      \
                                                                                                                       \
            int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                         \
            int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                         \
                                                                                                                       \
            int center_index = i_r - numberSmootherCircles;                                                            \
            int left_index   = i_r - numberSmootherCircles - 1;                                                        \
            int right_index  = i_r - numberSmootherCircles + 1;                                                        \
            int bottom_index = i_r - numberSmootherCircles;                                                            \
            int top_index    = i_r - numberSmootherCircles;                                                            \
                                                                                                                       \
            if (i_theta & 1) {                                                                                         \
                /* i_theta % 2 == 1 */                                                                                 \
                /* ---------------|| */                                                                                \
                /* O   X   O   X  || */                                                                                \
                /* ---------------|| */                                                                                \
                /* O   O   Õ   O  || */                                                                                \
                /* ---------------|| */                                                                                \
                /* O   X   O   X  || */                                                                                \
                /* ---------------|| */                                                                                \
                auto& center_matrix = radial_tridiagonal_solver[i_theta / 2];                                          \
                auto& bottom_matrix = radial_diagonal_solver[i_theta_M1 / 2];                                          \
                auto& top_matrix    = radial_diagonal_solver[i_theta_P1 / 2];                                          \
                /* Fill matrix row of (i,j) */                                                                         \
                row    = center_index;                                                                                 \
                column = center_index;                                                                                 \
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */             \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                row    = center_index;                                                                                 \
                column = left_index;                                                                                   \
                value  = -coeff1 * arr; /* Left */                                                                     \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                /* Remark: Right is not included here due to the symmetry shift */                                     \
                                                                                                                       \
                row    = center_index;                                                                                 \
                column = center_index;                                                                                 \
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */   \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                /* Fill matrix row of (i-1,j) */                                                                       \
                row    = left_index;                                                                                   \
                column = center_index;                                                                                 \
                value  = -coeff1 * arr; /* Right */                                                                    \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                row    = left_index;                                                                                   \
                column = left_index;                                                                                   \
                value  = coeff1 * arr; /* Center: (Right) */                                                           \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                /* Fill matrix row of (i+1,j) */                                                                       \
                /* Nothing to be done */                                                                               \
                                                                                                                       \
                /* Fill matrix row of (i,j-1) */                                                                       \
                row    = bottom_index;                                                                                 \
                column = bottom_index;                                                                                 \
                value  = coeff3 * att; /* Center: (Top) */                                                             \
                bottom_matrix.diagonal(row) += value;                                                                  \
                                                                                                                       \
                /* Fill matrix row of (i,j+1) */                                                                       \
                row    = top_index;                                                                                    \
                column = top_index;                                                                                    \
                value  = coeff4 * att; /* Center: (Bottom) */                                                          \
                top_matrix.diagonal(row) += value;                                                                     \
            }                                                                                                          \
            else {                                                                                                     \
                /* i_theta % 2 == 0 */                                                                                 \
                /* ---------------|| */                                                                                \
                /* O   O   O   O  || */                                                                                \
                /* ---------------|| */                                                                                \
                /* O   X   Õ   X  || */                                                                                \
                /* ---------------|| */                                                                                \
                /* O   O   O   O  || */                                                                                \
                /* ---------------|| */                                                                                \
                auto& center_matrix = radial_diagonal_solver[i_theta / 2];                                             \
                auto& bottom_matrix = radial_tridiagonal_solver[i_theta_M1 / 2];                                       \
                auto& top_matrix    = radial_tridiagonal_solver[i_theta_P1 / 2];                                       \
                /* Fill matrix row of (i,j) */                                                                         \
                row    = center_index;                                                                                 \
                column = center_index;                                                                                 \
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */             \
                center_matrix.diagonal(row) += value;                                                                  \
                                                                                                                       \
                row    = center_index;                                                                                 \
                column = center_index;                                                                                 \
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */   \
                center_matrix.diagonal(row) += value;                                                                  \
                                                                                                                       \
                /* Fill matrix row of (i,j-1) */                                                                       \
                row    = bottom_index;                                                                                 \
                column = bottom_index;                                                                                 \
                value  = coeff3 * att; /* Center: (Top) */                                                             \
                if (row == column)                                                                                     \
                    bottom_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    bottom_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == bottom_matrix.columns() - 1)                                            \
                    bottom_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                /* Fill matrix row of (i,j+1) */                                                                       \
                row    = top_index;                                                                                    \
                column = top_index;                                                                                    \
                value  = coeff4 * att; /* Center: (Bottom) */                                                          \
                if (row == column)                                                                                     \
                    top_matrix.main_diagonal(row) += value;                                                            \
                else if (row == column - 1)                                                                            \
                    top_matrix.sub_diagonal(row) += value;                                                             \
                else if (row == 0 && column == top_matrix.columns() - 1)                                               \
                    top_matrix.cyclic_corner_element() += value;                                                       \
            }                                                                                                          \
        }                                                                                                              \
        /* ------------------------------------------ */                                                               \
        /* Radial Section: Node on the outer boundary */                                                               \
        /* ------------------------------------------ */                                                               \
        else if (i_r == grid.nr() - 1) {                                                                               \
            assert(!i_r % 2 == 0);                                                                                     \
                                                                                                                       \
            double h1     = grid.radialSpacing(i_r - 1);                                                               \
            double k1     = grid.angularSpacing(i_theta - 1);                                                          \
            double k2     = grid.angularSpacing(i_theta);                                                              \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
                                                                                                                       \
            int center_index = i_r - numberSmootherCircles;                                                            \
            int left_index   = i_r - numberSmootherCircles - 1;                                                        \
                                                                                                                       \
            if (i_theta & 1) {                                                                                         \
                /* i_theta % 2 == 1 */                                                                                 \
                /* -----------|| */                                                                                    \
                /* X   O   X  || */                                                                                    \
                /* -----------|| */                                                                                    \
                /* O   O   Õ  || */                                                                                    \
                /* -----------|| */                                                                                    \
                /* X   O   X  || */                                                                                    \
                /* -----------|| */                                                                                    \
                auto& center_matrix = radial_tridiagonal_solver[i_theta / 2];                                          \
                                                                                                                       \
                /* Fill matrix row of (i,j) */                                                                         \
                row    = center_index;                                                                                 \
                column = center_index;                                                                                 \
                value  = 1.0;                                                                                          \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
                                                                                                                       \
                /* Fill matrix row of (i-1,j) */                                                                       \
                row    = left_index;                                                                                   \
                column = left_index;                                                                                   \
                value  = coeff1 * arr; /* Center: (Right) */                                                           \
                if (row == column)                                                                                     \
                    center_matrix.main_diagonal(row) += value;                                                         \
                else if (row == column - 1)                                                                            \
                    center_matrix.sub_diagonal(row) += value;                                                          \
                else if (row == 0 && column == center_matrix.columns() - 1)                                            \
                    center_matrix.cyclic_corner_element() += value;                                                    \
            }                                                                                                          \
            else {                                                                                                     \
                /* i_theta % 2 == 0 */                                                                                 \
                /* -----------|| */                                                                                    \
                /* O   O   O  || */                                                                                    \
                /* -----------|| */                                                                                    \
                /* X   O   X̃  || */                                                                                    \
                /* -----------|| */                                                                                    \
                /* O   O   O  || */                                                                                    \
                /* -----------|| */                                                                                    \
                auto& center_matrix = radial_diagonal_solver[i_theta / 2];                                             \
                                                                                                                       \
                /* Fill matrix row of (i,j) */                                                                         \
                row    = center_index;                                                                                 \
                column = center_index;                                                                                 \
                center_matrix.diagonal(row) += 1.0;                                                                    \
                                                                                                                       \
                /* Fill matrix row of (i-1,j) */                                                                       \
                row    = left_index;                                                                                   \
                column = left_index;                                                                                   \
                value  = coeff1 * arr; /* Center: (Right) */                                                           \
                center_matrix.diagonal(row) += value;                                                                  \
            }                                                                                                          \
        }                                                                                                              \
    } while (0)

void ExtrapolatedSmootherGive::buildAscCircleSection(const int i_r)
{
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

        // Build Asc at the current node
        NODE_BUILD_SMOOTHER_GIVE(i_r, i_theta, grid_, DirBC_Interior_, inner_boundary_circle_matrix_,
                                 circle_diagonal_solver_, circle_tridiagonal_solver_, radial_diagonal_solver_,
                                 radial_tridiagonal_solver_);
    }
}

void ExtrapolatedSmootherGive::buildAscRadialSection(const int i_theta)
{
    const auto& sin_theta_cache = level_cache_.sin_theta();
    const auto& cos_theta_cache = level_cache_.cos_theta();

    const double theta     = grid_.theta(i_theta);
    const double sin_theta = sin_theta_cache[i_theta];
    const double cos_theta = cos_theta_cache[i_theta];

    for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
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

        // Build Asc at the current node
        NODE_BUILD_SMOOTHER_GIVE(i_r, i_theta, grid_, DirBC_Interior_, inner_boundary_circle_matrix_,
                                 circle_diagonal_solver_, circle_tridiagonal_solver_, radial_diagonal_solver_,
                                 radial_tridiagonal_solver_);
    }
}

void ExtrapolatedSmootherGive::buildAscMatrices()
{
    omp_set_num_threads(num_omp_threads_);

    /* -------------------------------------- */
    /* Part 1: Allocate Asc Smoother matrices */
    /* -------------------------------------- */

    const int number_smoother_circles = grid_.numberSmootherCircles();
    const int length_smoother_radial  = grid_.lengthSmootherRadial();

    const int num_circle_nodes = grid_.ntheta();
    circle_tridiagonal_solver_.resize(number_smoother_circles / 2);
    circle_diagonal_solver_.resize(number_smoother_circles - number_smoother_circles / 2);

    assert((grid_.ntheta() / 2) % 2 == 0);
    const int num_radial_nodes = length_smoother_radial;
    radial_tridiagonal_solver_.resize(grid_.ntheta() / 2);
    radial_diagonal_solver_.resize(grid_.ntheta() / 2);

// Remark: circle_diagonal_solver_[0] is undefnied.
// Use inner_boundary_circle_matrix_ instead.
#pragma omp parallel if (grid_.numberOfNodes() > 10'000)
    {
// ---------------- //
// Circular Section //
#pragma omp for nowait
        for (int circle_Asc_index = 0; circle_Asc_index < number_smoother_circles; circle_Asc_index++) {

            /* Inner boundary circle */
            if (circle_Asc_index == 0) {
                // Although the matrix is symmetric, we need to store all its entries, so we disable the symmetry.
                const int nnz                 = getNonZeroCountCircleAsc(circle_Asc_index);
                inner_boundary_circle_matrix_ = SparseMatrixCOO<double>(num_circle_nodes, num_circle_nodes, nnz);
                inner_boundary_circle_matrix_.is_symmetric(false);
                for (int i = 0; i < nnz; i++) {
                    inner_boundary_circle_matrix_.value(i) = 0.0;
                }
            }

            /* Interior Circle Section */
            else {
                if (circle_Asc_index & 1) {
                    const int circle_tridiagonal_solver_index = circle_Asc_index / 2;
                    auto& solver_matrix = circle_tridiagonal_solver_[circle_tridiagonal_solver_index];

                    solver_matrix = SymmetricTridiagonalSolver<double>(num_circle_nodes);
                    solver_matrix.is_cyclic(true);
                    solver_matrix.cyclic_corner_element() = 0.0;

                    for (int i = 0; i < num_circle_nodes; i++) {
                        solver_matrix.main_diagonal(i) = 0.0;
                        if (i < num_circle_nodes - 1) {
                            solver_matrix.sub_diagonal(i) = 0.0;
                        }
                    }
                }
                else {
                    const int circle_diagonal_solver_index = circle_Asc_index / 2;
                    auto& solver_matrix                    = circle_diagonal_solver_[circle_diagonal_solver_index];

                    solver_matrix = DiagonalSolver<double>(num_circle_nodes);

                    for (int i = 0; i < num_circle_nodes; i++) {
                        solver_matrix.diagonal(i) = 0.0;
                    }
                }
            }
        }

// -------------- //
// Radial Section //
#pragma omp for nowait
        for (int radial_Asc_index = 0; radial_Asc_index < grid_.ntheta(); radial_Asc_index++) {

            if (radial_Asc_index & 1) {
                const int radial_tridiagonal_solver_index = radial_Asc_index / 2;
                auto& solver_matrix                       = radial_tridiagonal_solver_[radial_tridiagonal_solver_index];

                solver_matrix = SymmetricTridiagonalSolver<double>(num_radial_nodes);
                solver_matrix.is_cyclic(false);

                for (int i = 0; i < num_radial_nodes; i++) {
                    solver_matrix.main_diagonal(i) = 0.0;
                    if (i < num_radial_nodes - 1) {
                        solver_matrix.sub_diagonal(i) = 0.0;
                    }
                }
            }
            else {
                const int radial_diagonal_solver_index = radial_Asc_index / 2;
                auto& solver_matrix                    = radial_diagonal_solver_[radial_diagonal_solver_index];

                solver_matrix = DiagonalSolver<double>(num_radial_nodes);

                for (int i = 0; i < num_radial_nodes; i++) {
                    solver_matrix.diagonal(i) = 0.0;
                }
            }
        }
    }

    /* ---------------------------------- */
    /* Part 2: Fill Asc Smoother matrices */
    /* ---------------------------------- */

    bool use_simple_parallelism = true; // Fastest: true

    if (omp_get_max_threads() == 1) {
        /* Single-threaded execution */
        for (int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
            buildAscCircleSection(i_r);
        }
        for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
            buildAscRadialSection(i_theta);
        }
    }
    else {
        if (use_simple_parallelism) {
            /*  Multi-threaded execution: For Loops */
            const int num_circle_tasks        = grid_.numberSmootherCircles();
            const int additional_radial_tasks = grid_.ntheta() % 3;
            const int num_radial_tasks        = grid_.ntheta() - additional_radial_tasks;

#pragma omp parallel
            {
#pragma omp for
                for (int circle_task = 0; circle_task < num_circle_tasks; circle_task += 3) {
                    int i_r = grid_.numberSmootherCircles() - circle_task - 1;
                    buildAscCircleSection(i_r);
                }
#pragma omp for
                for (int circle_task = 1; circle_task < num_circle_tasks; circle_task += 3) {
                    int i_r = grid_.numberSmootherCircles() - circle_task - 1;
                    buildAscCircleSection(i_r);
                }
#pragma omp for nowait
                for (int circle_task = 2; circle_task < num_circle_tasks; circle_task += 3) {
                    int i_r = grid_.numberSmootherCircles() - circle_task - 1;
                    buildAscCircleSection(i_r);
                }

#pragma omp for
                for (int radial_task = 0; radial_task < num_radial_tasks; radial_task += 3) {
                    if (radial_task > 0) {
                        int i_theta = radial_task + additional_radial_tasks;
                        buildAscRadialSection(i_theta);
                    }
                    else {
                        if (additional_radial_tasks == 0) {
                            buildAscRadialSection(0);
                        }
                        else if (additional_radial_tasks >= 1) {
                            buildAscRadialSection(0);
                            buildAscRadialSection(1);
                        }
                    }
                }
#pragma omp for
                for (int radial_task = 1; radial_task < num_radial_tasks; radial_task += 3) {
                    if (radial_task > 1) {
                        int i_theta = radial_task + additional_radial_tasks;
                        buildAscRadialSection(i_theta);
                    }
                    else {
                        if (additional_radial_tasks == 0) {
                            buildAscRadialSection(1);
                        }
                        else if (additional_radial_tasks == 1) {
                            buildAscRadialSection(2);
                        }
                        else if (additional_radial_tasks == 2) {
                            buildAscRadialSection(2);
                            buildAscRadialSection(3);
                        }
                    }
                }
#pragma omp for
                for (int radial_task = 2; radial_task < num_radial_tasks; radial_task += 3) {
                    int i_theta = radial_task + additional_radial_tasks;
                    buildAscRadialSection(i_theta);
                }
            }
        }
        else {
            /*  Multi-threaded execution: Task Dependencies */
            const int num_circle_tasks        = grid_.numberSmootherCircles();
            const int additional_radial_tasks = grid_.ntheta() % 3;
            const int num_radial_tasks        = grid_.ntheta() - additional_radial_tasks;

            assert(num_circle_tasks >= 2);
            assert(num_radial_tasks >= 3 && num_radial_tasks % 3 == 0);

            /* Make sure to deallocate at the end */
            const int boundary_margin = 2; // Additional space to ensure safe access
            int* circle_dep           = new int[num_circle_tasks + boundary_margin];
            int* radial_dep           = new int[num_radial_tasks];

#pragma omp parallel
            {
#pragma omp single
                {
                    /* ------------ */
                    /* Circle Tasks */
                    /* ------------ */

                    /* Mod 0 Circles */
                    for (int circle_task = 0; circle_task < num_circle_tasks; circle_task += 3) {
#pragma omp task depend(out : circle_dep[circle_task])
                        {
                            const int i_r = grid_.numberSmootherCircles() - circle_task - 1;
                            buildAscCircleSection(i_r);
                        }
                    }
                    /* Mod 2 Circles */
                    for (int circle_task = 1; circle_task < num_circle_tasks; circle_task += 3) {
#pragma omp task depend(out : circle_dep[circle_task])                                                                 \
    depend(in : circle_dep[circle_task - 1], circle_dep[circle_task + 2])
                        {
                            const int i_r = grid_.numberSmootherCircles() - circle_task - 1;
                            buildAscCircleSection(i_r);
                        }
                    }
                    /* Mod 2 Circles */
                    for (int circle_task = 2; circle_task < num_circle_tasks; circle_task += 3) {
#pragma omp task depend(out : circle_dep[circle_task])                                                                 \
    depend(in : circle_dep[circle_task - 1], circle_dep[circle_task + 2])
                        {
                            const int i_r = grid_.numberSmootherCircles() - circle_task - 1;
                            buildAscCircleSection(i_r);
                        }
                    }

                    /* ------------ */
                    /* Radial Tasks */
                    /* ------------ */

                    /* Mod 0 Radials */
                    for (int radial_task = 0; radial_task < num_radial_tasks; radial_task += 3) {
#pragma omp task depend(out : radial_dep[radial_task]) depend(in : circle_dep[1])
                        {
                            if (radial_task > 2) {
                                const int i_theta = radial_task + additional_radial_tasks;
                                buildAscRadialSection(i_theta);
                            }
                            else {
                                if (additional_radial_tasks == 0) {
                                    buildAscRadialSection(0);
                                }
                                else if (additional_radial_tasks >= 1) {
                                    buildAscRadialSection(0);
                                    buildAscRadialSection(1);
                                }
                            }
                        }
                    }
                    /* Mod 1 Radials */
                    for (int radial_task = 1; radial_task < num_radial_tasks; radial_task += 3) {
#pragma omp task depend(out : radial_dep[radial_task])                                                                 \
    depend(in : radial_dep[radial_task - 1], radial_dep[(radial_task + 2) % num_radial_tasks])
                        {
                            if (radial_task > 2) {
                                int i_theta = radial_task + additional_radial_tasks;
                                buildAscRadialSection(i_theta);
                            }
                            else {
                                if (additional_radial_tasks == 0) {
                                    buildAscRadialSection(1);
                                }
                                else if (additional_radial_tasks == 1) {
                                    buildAscRadialSection(2);
                                }
                                else if (additional_radial_tasks == 2) {
                                    buildAscRadialSection(2);
                                    buildAscRadialSection(3);
                                }
                            }
                        }
                    }
                    /* Mod 2 Radials */
                    for (int radial_task = 2; radial_task < num_radial_tasks; radial_task += 3) {
#pragma omp task depend(out : radial_dep[radial_task])                                                                 \
    depend(in : radial_dep[radial_task - 1], radial_dep[(radial_task + 2) % num_radial_tasks])
                        {
                            int i_theta = radial_task + additional_radial_tasks;
                            buildAscRadialSection(i_theta);
                        }
                    }
                }
            }
            delete[] circle_dep;
            delete[] radial_dep;
        }
    }

    /* ------------------------------------------------------------------- */
    /* Part 3: Convert inner_boundary_circle_matrix_ to a symmetric matrix */
    /* ------------------------------------------------------------------- */

    SparseMatrixCOO<double> full_matrix = std::move(inner_boundary_circle_matrix_);

    const int nnz           = full_matrix.non_zero_size();
    const int numRows       = full_matrix.rows();
    const int numColumns    = full_matrix.columns();
    const int symmetric_nnz = nnz - (nnz - numRows) / 2;

    inner_boundary_circle_matrix_ = SparseMatrixCOO<double>(numRows, numColumns, symmetric_nnz);
    inner_boundary_circle_matrix_.is_symmetric(true);

    int current_nz = 0;
    for (int nz_index = 0; nz_index < full_matrix.non_zero_size(); nz_index++) {
        int current_row = full_matrix.row_index(nz_index);
        int current_col = full_matrix.col_index(nz_index);
        if (current_row <= current_col) {
            inner_boundary_circle_matrix_.row_index(current_nz) = current_row;
            inner_boundary_circle_matrix_.col_index(current_nz) = current_col;
            inner_boundary_circle_matrix_.value(current_nz)     = std::move(full_matrix.value(nz_index));
            current_nz++;
        }
    }
}