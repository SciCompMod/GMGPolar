#include "../../../include/DirectSolver/DirectSolverTake/directSolverTake.h"

#ifdef GMGPOLAR_USE_MUMPS

    #define UPDATE_MATRIX_ELEMENT(matrix, ptr, offset, row, col, val)                                                  \
        do {                                                                                                           \
            matrix.row_index(ptr + offset) = row;                                                                      \
            matrix.col_index(ptr + offset) = col;                                                                      \
            matrix.value(ptr + offset)     = val;                                                                      \
        } while (0)

    #define NODE_BUILD_SOLVER_MATRIX_TAKE(i_r, i_theta, grid, DirBC_Interior, solver_matrix, arr, att, art, detDF,      \
                                          coeff_beta)                                                                   \
        do {                                                                                                            \
            int ptr, offset;                                                                                            \
            int row, col;                                                                                               \
            double val;                                                                                                 \
            /* -------------------- */                                                                                  \
            /* Node in the interior */                                                                                  \
            /* -------------------- */                                                                                  \
            if (i_r > 1 && i_r < grid.nr() - 2) {                                                                       \
                const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                \
                const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                \
                                                                                                                        \
                const double h1     = grid.radialSpacing(i_r - 1);                                                      \
                const double h2     = grid.radialSpacing(i_r);                                                          \
                const double k1     = grid.angularSpacing(i_theta_M1);                                                  \
                const double k2     = grid.angularSpacing(i_theta);                                                     \
                const double coeff1 = 0.5 * (k1 + k2) / h1;                                                             \
                const double coeff2 = 0.5 * (k1 + k2) / h2;                                                             \
                const double coeff3 = 0.5 * (h1 + h2) / k1;                                                             \
                const double coeff4 = 0.5 * (h1 + h2) / k2;                                                             \
                                                                                                                        \
                const int center_nz_index = getSolverMatrixIndex(i_r, i_theta);                                         \
                                                                                                                        \
                const int center_index       = grid.index(i_r, i_theta);                                                \
                const int left_index         = grid.index(i_r - 1, i_theta);                                            \
                const int right_index        = grid.index(i_r + 1, i_theta);                                            \
                const int bottom_index       = grid.index(i_r, i_theta_M1);                                             \
                const int top_index          = grid.index(i_r, i_theta_P1);                                             \
                const int bottom_left_index  = grid.index(i_r - 1, i_theta_M1);                                         \
                const int bottom_right_index = grid.index(i_r + 1, i_theta_M1);                                         \
                const int top_left_index     = grid.index(i_r - 1, i_theta_P1);                                         \
                const int top_right_index    = grid.index(i_r + 1, i_theta_P1);                                         \
                                                                                                                        \
                const double left_value   = -coeff1 * (arr[center_index] + arr[left_index]); /* Left */                 \
                const double right_value  = -coeff2 * (arr[center_index] + arr[right_index]); /* Right */               \
                const double bottom_value = -coeff3 * (att[center_index] + att[bottom_index]); /* Bottom */             \
                const double top_value    = -coeff4 * (att[center_index] + att[top_index]); /* Top */                   \
                                                                                                                        \
                const double center_value = (+0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center_index] *                 \
                                                 fabs(detDF[center_index]) /* beta_{i,j} */                             \
                                             - left_value /* Center: (Left) */                                          \
                                             - right_value /* Center: (Right) */                                        \
                                             - bottom_value /* Center: (Bottom) */                                      \
                                             - top_value /* Center: (Top) */                                            \
                );                                                                                                      \
                                                                                                                        \
                const double bottom_left_value  = -0.25 * (art[left_index] + art[bottom_index]); /* Bottom Left */      \
                const double bottom_right_value = +0.25 * (art[right_index] + art[bottom_index]); /* Bottom Right */    \
                const double top_left_value     = +0.25 * (art[left_index] + art[top_index]); /* Top Left */            \
                const double top_right_value    = -0.25 * (art[right_index] + art[top_index]); /* Top Right */          \
                                                                                                                        \
                /* Fill matrix row of (i,j) */                                                                          \
                row = center_index;                                                                                     \
                ptr = center_nz_index;                                                                                  \
                                                                                                                        \
                const Stencil& CenterStencil = getStencil(i_r);                                                         \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::Center];                                                        \
                col    = center_index;                                                                                  \
                val    = center_value;                                                                                  \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::Left];                                                          \
                col    = left_index;                                                                                    \
                val    = left_value;                                                                                    \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::Right];                                                         \
                col    = right_index;                                                                                   \
                val    = right_value;                                                                                   \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::Bottom];                                                        \
                col    = bottom_index;                                                                                  \
                val    = bottom_value;                                                                                  \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::Top];                                                           \
                col    = top_index;                                                                                     \
                val    = top_value;                                                                                     \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::BottomLeft];                                                    \
                col    = bottom_left_index;                                                                             \
                val    = bottom_left_value;                                                                             \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::BottomRight];                                                   \
                col    = bottom_right_index;                                                                            \
                val    = bottom_right_value;                                                                            \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::TopLeft];                                                       \
                col    = top_left_index;                                                                                \
                val    = top_left_value;                                                                                \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::TopRight];                                                      \
                col    = top_right_index;                                                                               \
                val    = top_right_value;                                                                               \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
            }                                                                                                           \
            /* -------------------------- */                                                                            \
            /* Node on the inner boundary */                                                                            \
            /* -------------------------- */                                                                            \
            else if (i_r == 0) {                                                                                        \
                /* ------------------------------------------------ */                                                  \
                /* Case 1: Dirichlet boundary on the inner boundary */                                                  \
                /* ------------------------------------------------ */                                                  \
                if (DirBC_Interior) {                                                                                   \
                    const int center_nz_index = getSolverMatrixIndex(i_r, i_theta);                                     \
                                                                                                                        \
                    const int center_index = grid.index(i_r, i_theta);                                                  \
                                                                                                                        \
                    /* Fill matrix row of (i,j) */                                                                      \
                    row = center_index;                                                                                 \
                    ptr = center_nz_index;                                                                              \
                                                                                                                        \
                    const Stencil& CenterStencil = getStencil(i_r);                                                     \
                                                                                                                        \
                    offset = CenterStencil[StencilPosition::Center];                                                    \
                    col    = center_index;                                                                              \
                    val    = 1.0;                                                                                       \
                    UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                   \
                }                                                                                                       \
                else {                                                                                                  \
                    /* ------------------------------------------------------------- */                                 \
                    /* Case 2: Across origin discretization on the interior boundary */                                 \
                    /* ------------------------------------------------------------- */                                 \
                    /* h1 gets replaced with 2 * R0. */                                                                 \
                    /* (i_r-1,i_theta) gets replaced with (i_r, i_theta + grid.ntheta()/2). */                          \
                    /* Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil. */ \
                    const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                            \
                    const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                            \
                                                                                                                        \
                    assert(grid_.ntheta() % 2 == 0);                                                                    \
                    const int i_theta_AcrossOrigin = grid.wrapThetaIndex(i_theta + grid.ntheta() / 2);                  \
                                                                                                                        \
                    double h1 = 2.0 * grid.radius(0);                                                                   \
                    double h2 = grid.radialSpacing(i_r);                                                                \
                    double k1 = grid.angularSpacing(i_theta_M1);                                                        \
                    double k2 = grid.angularSpacing(i_theta);                                                           \
                                                                                                                        \
                    const double coeff1 = 0.5 * (k1 + k2) / h1;                                                         \
                    const double coeff2 = 0.5 * (k1 + k2) / h2;                                                         \
                    const double coeff3 = 0.5 * (h1 + h2) / k1;                                                         \
                    const double coeff4 = 0.5 * (h1 + h2) / k2;                                                         \
                                                                                                                        \
                    const int center_nz_index = getSolverMatrixIndex(i_r, i_theta);                                     \
                                                                                                                        \
                    const int center_index       = grid.index(i_r, i_theta);                                            \
                    const int left_index         = grid.index(i_r, i_theta_AcrossOrigin);                               \
                    const int right_index        = grid.index(i_r + 1, i_theta);                                        \
                    const int bottom_index       = grid.index(i_r, i_theta_M1);                                         \
                    const int top_index          = grid.index(i_r, i_theta_P1);                                         \
                    const int bottom_right_index = grid.index(i_r + 1, i_theta_M1);                                     \
                    const int top_right_index    = grid.index(i_r + 1, i_theta_P1);                                     \
                                                                                                                        \
                    const double left_value   = -coeff1 * (arr[center_index] + arr[left_index]); /* Left */             \
                    const double right_value  = -coeff2 * (arr[center_index] + arr[right_index]); /* Right */           \
                    const double bottom_value = -coeff3 * (att[center_index] + att[bottom_index]); /* Bottom */         \
                    const double top_value    = -coeff4 * (att[center_index] + att[top_index]); /* Top */               \
                                                                                                                        \
                    const double center_value = (+0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center_index] *             \
                                                     fabs(detDF[center_index]) /* beta_{i,j} */                         \
                                                 - left_value /* Center: (Left) */                                      \
                                                 - right_value /* Center: (Right) */                                    \
                                                 - bottom_value /* Center: (Bottom) */                                  \
                                                 - top_value /* Center: (Top) */                                        \
                    );                                                                                                  \
                                                                                                                        \
                    const double bottom_right_value =                                                                   \
                        +0.25 * (art[right_index] + art[bottom_index]); /* Bottom Right */                              \
                    const double top_right_value = -0.25 * (art[right_index] + art[top_index]); /* Top Right */         \
                                                                                                                        \
                    /* Fill matrix row of (i,j) */                                                                      \
                    row = center_index;                                                                                 \
                    ptr = center_nz_index;                                                                              \
                                                                                                                        \
                    const Stencil& CenterStencil = getStencil(i_r);                                                     \
                                                                                                                        \
                    offset = CenterStencil[StencilPosition::Center];                                                    \
                    col    = center_index;                                                                              \
                    val    = center_value;                                                                              \
                    UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                   \
                                                                                                                        \
                    offset = CenterStencil[StencilPosition::Left];                                                      \
                    col    = left_index;                                                                                \
                    val    = left_value;                                                                                \
                    UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                   \
                                                                                                                        \
                    offset = CenterStencil[StencilPosition::Right];                                                     \
                    col    = right_index;                                                                               \
                    val    = right_value;                                                                               \
                    UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                   \
                                                                                                                        \
                    offset = CenterStencil[StencilPosition::Bottom];                                                    \
                    col    = bottom_index;                                                                              \
                    val    = bottom_value;                                                                              \
                    UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                   \
                                                                                                                        \
                    offset = CenterStencil[StencilPosition::Top];                                                       \
                    col    = top_index;                                                                                 \
                    val    = top_value;                                                                                 \
                    UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                   \
                                                                                                                        \
                    /* BottomLeft: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */                                          \
                                                                                                                        \
                    offset = CenterStencil[StencilPosition::BottomRight];                                               \
                    col    = bottom_right_index;                                                                        \
                    val    = bottom_right_value;                                                                        \
                    UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                   \
                                                                                                                        \
                    /* TopLeft: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */                                             \
                                                                                                                        \
                    offset = CenterStencil[StencilPosition::TopRight];                                                  \
                    col    = top_right_index;                                                                           \
                    val    = top_right_value;                                                                           \
                    UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                   \
                }                                                                                                       \
            }                                                                                                           \
            /* ------------------------------- */                                                                       \
            /* Node next to the inner boundary */                                                                       \
            /* ------------------------------- */                                                                       \
            else if (i_r == 1) {                                                                                        \
                const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                \
                const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                \
                                                                                                                        \
                const double h1     = grid.radialSpacing(i_r - 1);                                                      \
                const double h2     = grid.radialSpacing(i_r);                                                          \
                const double k1     = grid.angularSpacing(i_theta_M1);                                                  \
                const double k2     = grid.angularSpacing(i_theta);                                                     \
                const double coeff1 = 0.5 * (k1 + k2) / h1;                                                             \
                const double coeff2 = 0.5 * (k1 + k2) / h2;                                                             \
                const double coeff3 = 0.5 * (h1 + h2) / k1;                                                             \
                const double coeff4 = 0.5 * (h1 + h2) / k2;                                                             \
                                                                                                                        \
                const int center_nz_index = getSolverMatrixIndex(i_r, i_theta);                                         \
                                                                                                                        \
                const int center_index       = grid.index(i_r, i_theta);                                                \
                const int left_index         = grid.index(i_r - 1, i_theta);                                            \
                const int right_index        = grid.index(i_r + 1, i_theta);                                            \
                const int bottom_index       = grid.index(i_r, i_theta_M1);                                             \
                const int top_index          = grid.index(i_r, i_theta_P1);                                             \
                const int bottom_left_index  = grid.index(i_r - 1, i_theta_M1);                                         \
                const int bottom_right_index = grid.index(i_r + 1, i_theta_M1);                                         \
                const int top_left_index     = grid.index(i_r - 1, i_theta_P1);                                         \
                const int top_right_index    = grid.index(i_r + 1, i_theta_P1);                                         \
                                                                                                                        \
                const double left_value   = -coeff1 * (arr[center_index] + arr[left_index]); /* Left */                 \
                const double right_value  = -coeff2 * (arr[center_index] + arr[right_index]); /* Right */               \
                const double bottom_value = -coeff3 * (att[center_index] + att[bottom_index]); /* Bottom */             \
                const double top_value    = -coeff4 * (att[center_index] + att[top_index]); /* Top */                   \
                                                                                                                        \
                const double center_value = (+0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center_index] *                 \
                                                 fabs(detDF[center_index]) /* beta_{i,j} */                             \
                                             - left_value /* Center: (Left) */                                          \
                                             - right_value /* Center: (Right) */                                        \
                                             - bottom_value /* Center: (Bottom) */                                      \
                                             - top_value /* Center: (Top) */                                            \
                );                                                                                                      \
                                                                                                                        \
                const double bottom_left_value  = -0.25 * (art[left_index] + art[bottom_index]); /* Bottom Left */      \
                const double bottom_right_value = +0.25 * (art[right_index] + art[bottom_index]); /* Bottom Right */    \
                const double top_left_value     = +0.25 * (art[left_index] + art[top_index]); /* Top Left */            \
                const double top_right_value    = -0.25 * (art[right_index] + art[top_index]); /* Top Right */          \
                                                                                                                        \
                /* Fill matrix row of (i,j) */                                                                          \
                row = center_index;                                                                                     \
                ptr = center_nz_index;                                                                                  \
                                                                                                                        \
                const Stencil& CenterStencil = getStencil(i_r);                                                         \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::Center];                                                        \
                col    = center_index;                                                                                  \
                val    = center_value;                                                                                  \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
                                                                                                                        \
                /* REMOVED: Moved to the right hand side to make the matrix symmetric */                                \
                if (!DirBC_Interior) {                                                                                  \
                    offset = CenterStencil[StencilPosition::Left];                                                      \
                    col    = left_index;                                                                                \
                    val    = left_value;                                                                                \
                    UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                   \
                }                                                                                                       \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::Right];                                                         \
                col    = right_index;                                                                                   \
                val    = right_value;                                                                                   \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::Bottom];                                                        \
                col    = bottom_index;                                                                                  \
                val    = bottom_value;                                                                                  \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::Top];                                                           \
                col    = top_index;                                                                                     \
                val    = top_value;                                                                                     \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
                                                                                                                        \
                /* REMOVED: Moved to the right hand side to make the matrix symmetric */                                \
                if (!DirBC_Interior) {                                                                                  \
                    offset = CenterStencil[StencilPosition::BottomLeft];                                                \
                    col    = bottom_left_index;                                                                         \
                    val    = bottom_left_value;                                                                         \
                    UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                   \
                }                                                                                                       \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::BottomRight];                                                   \
                col    = bottom_right_index;                                                                            \
                val    = bottom_right_value;                                                                            \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
                                                                                                                        \
                /* REMOVED: Moved to the right hand side to make the matrix symmetric */                                \
                if (!DirBC_Interior) {                                                                                  \
                    offset = CenterStencil[StencilPosition::TopLeft];                                                   \
                    col    = top_left_index;                                                                            \
                    val    = top_left_value;                                                                            \
                    UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                   \
                }                                                                                                       \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::TopRight];                                                      \
                col    = top_right_index;                                                                               \
                val    = top_right_value;                                                                               \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
            }                                                                                                           \
            /* ------------------------------- */                                                                       \
            /* Node next to the outer boundary */                                                                       \
            /* ------------------------------- */                                                                       \
            else if (i_r == grid.nr() - 2) {                                                                            \
                const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                \
                const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                \
                                                                                                                        \
                const double h1 = grid.radialSpacing(i_r - 1);                                                          \
                const double h2 = grid.radialSpacing(i_r);                                                              \
                const double k1 = grid.angularSpacing(i_theta_M1);                                                      \
                const double k2 = grid.angularSpacing(i_theta);                                                         \
                                                                                                                        \
                const double coeff1 = 0.5 * (k1 + k2) / h1;                                                             \
                const double coeff2 = 0.5 * (k1 + k2) / h2;                                                             \
                const double coeff3 = 0.5 * (h1 + h2) / k1;                                                             \
                const double coeff4 = 0.5 * (h1 + h2) / k2;                                                             \
                                                                                                                        \
                const int center_nz_index = getSolverMatrixIndex(i_r, i_theta);                                         \
                                                                                                                        \
                const int center_index       = grid.index(i_r, i_theta);                                                \
                const int left_index         = grid.index(i_r - 1, i_theta);                                            \
                const int right_index        = grid.index(i_r + 1, i_theta);                                            \
                const int bottom_index       = grid.index(i_r, i_theta_M1);                                             \
                const int top_index          = grid.index(i_r, i_theta_P1);                                             \
                const int bottom_left_index  = grid.index(i_r - 1, i_theta_M1);                                         \
                const int bottom_right_index = grid.index(i_r + 1, i_theta_M1);                                         \
                const int top_left_index     = grid.index(i_r - 1, i_theta_P1);                                         \
                const int top_right_index    = grid.index(i_r + 1, i_theta_P1);                                         \
                                                                                                                        \
                const double left_value   = -coeff1 * (arr[center_index] + arr[left_index]); /* Left */                 \
                const double right_value  = -coeff2 * (arr[center_index] + arr[right_index]); /* Right */               \
                const double bottom_value = -coeff3 * (att[center_index] + att[bottom_index]); /* Bottom */             \
                const double top_value    = -coeff4 * (att[center_index] + att[top_index]); /* Top */                   \
                                                                                                                        \
                const double center_value = (+0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center_index] *                 \
                                                 fabs(detDF[center_index]) /* beta_{i,j} */                             \
                                             - left_value /* Center: (Left) */                                          \
                                             - right_value /* Center: (Right) */                                        \
                                             - bottom_value /* Center: (Bottom) */                                      \
                                             - top_value /* Center: (Top) */                                            \
                );                                                                                                      \
                                                                                                                        \
                const double bottom_left_value  = -0.25 * (art[left_index] + art[bottom_index]); /* Bottom Left */      \
                const double bottom_right_value = +0.25 * (art[right_index] + art[bottom_index]); /* Bottom Right */    \
                const double top_left_value     = +0.25 * (art[left_index] + art[top_index]); /* Top Left */            \
                const double top_right_value    = -0.25 * (art[right_index] + art[top_index]); /* Top Right */          \
                                                                                                                        \
                /* Fill matrix row of (i,j) */                                                                          \
                row = center_index;                                                                                     \
                ptr = center_nz_index;                                                                                  \
                                                                                                                        \
                const Stencil& CenterStencil = getStencil(i_r);                                                         \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::Center];                                                        \
                col    = center_index;                                                                                  \
                val    = center_value;                                                                                  \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::Left];                                                          \
                col    = left_index;                                                                                    \
                val    = left_value;                                                                                    \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
                                                                                                                        \
                /* Right REMOVED: Moved to the right hand side to make the matrix symmetric */                          \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::Bottom];                                                        \
                col    = bottom_index;                                                                                  \
                val    = bottom_value;                                                                                  \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::Top];                                                           \
                col    = top_index;                                                                                     \
                val    = top_value;                                                                                     \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::BottomLeft];                                                    \
                col    = bottom_left_index;                                                                             \
                val    = bottom_left_value;                                                                             \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
                                                                                                                        \
                /* BottomRight REMOVED: Moved to the right hand side to make the matrix symmetric */                    \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::TopLeft];                                                       \
                col    = top_left_index;                                                                                \
                val    = top_left_value;                                                                                \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
                                                                                                                        \
                /* TopRight REMOVED: Moved to the right hand side to make the matrix symmetric */                       \
            }                                                                                                           \
            /* ------------------------------------ */                                                                  \
            /* Node on the outer dirichlet boundary */                                                                  \
            /* ------------------------------------ */                                                                  \
            else if (i_r == grid.nr() - 1) {                                                                            \
                const int center_nz_index = getSolverMatrixIndex(i_r, i_theta);                                         \
                                                                                                                        \
                const int center_index = grid.index(i_r, i_theta);                                                      \
                                                                                                                        \
                /* Fill matrix row of (i,j) */                                                                          \
                row = center_index;                                                                                     \
                ptr = center_nz_index;                                                                                  \
                                                                                                                        \
                const Stencil& CenterStencil = getStencil(i_r);                                                         \
                                                                                                                        \
                offset = CenterStencil[StencilPosition::Center];                                                        \
                col    = center_index;                                                                                  \
                val    = 1.0;                                                                                           \
                UPDATE_MATRIX_ELEMENT(solver_matrix, ptr, offset, row, col, val);                                       \
            }                                                                                                           \
        } while (0)

void DirectSolverTake::buildSolverMatrixCircleSection(const int i_r, SparseMatrixCOO<double>& solver_matrix)
{
    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    const auto& arr        = level_cache_.arr();
    const auto& att        = level_cache_.att();
    const auto& art        = level_cache_.art();
    const auto& detDF      = level_cache_.detDF();
    const auto& coeff_beta = level_cache_.coeff_beta();

    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        // Build solver matrix at the current node
        NODE_BUILD_SOLVER_MATRIX_TAKE(i_r, i_theta, grid_, DirBC_Interior_, solver_matrix, arr, att, art, detDF,
                                      coeff_beta);
    }
}

void DirectSolverTake::buildSolverMatrixRadialSection(const int i_theta, SparseMatrixCOO<double>& solver_matrix)
{
    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    const auto& arr        = level_cache_.arr();
    const auto& att        = level_cache_.att();
    const auto& art        = level_cache_.art();
    const auto& detDF      = level_cache_.detDF();
    const auto& coeff_beta = level_cache_.coeff_beta();

    for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
        // Build solver matrix at the current node
        NODE_BUILD_SOLVER_MATRIX_TAKE(i_r, i_theta, grid_, DirBC_Interior_, solver_matrix, arr, att, art, detDF,
                                      coeff_beta);
    }
}

// clang-format off

/* ------------------------------------------------------------------------ */
/* If the indexing is not smoother-based, please adjust the access patterns */
SparseMatrixCOO<double> DirectSolverTake::buildSolverMatrix()
{
    const int n   = grid_.numberOfNodes();
    const int nnz = getNonZeroCountSolverMatrix();

    // Although the matrix is symmetric, we need to store all its entries, so we disable the symmetry.
    SparseMatrixCOO<double> solver_matrix(n, n, nnz);
    solver_matrix.is_symmetric(false);

    if (num_omp_threads_ == 1) {
        /* Single-threaded execution */
        for (int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
            buildSolverMatrixCircleSection(i_r, solver_matrix);
        }
        for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
            buildSolverMatrixRadialSection(i_theta, solver_matrix);
        }
    }
    else {
        /* Multi-threaded execution */
        #pragma omp parallel num_threads(num_omp_threads_)
        {
            /* Circle Section */
            #pragma omp for nowait
            for (int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
                buildSolverMatrixCircleSection(i_r, solver_matrix);
            }
            /* Radial Section */
            #pragma omp for nowait
            for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
                buildSolverMatrixRadialSection(i_theta, solver_matrix);
            }
        }
    }

    /* Mumps: In the case of symmetric matrices, only half of the matrix should be provided. */
    const bool construct_symmetric = true;

    if (!construct_symmetric) {
        return solver_matrix;
    }

    /* Only store the upper tridiagonal entries of the symmetric solver_matrix */
    const int symmetric_nnz = nnz - (nnz - n) / 2;
    SparseMatrixCOO<double> symmetric_solver_matrix(n, n, symmetric_nnz);
    symmetric_solver_matrix.is_symmetric(true);

    int current_nz = 0;
    for (int nz_index = 0; nz_index < nnz; nz_index++) {
        const int row = solver_matrix.row_index(nz_index);
        const int col = solver_matrix.col_index(nz_index);
        if (row <= col) {
            symmetric_solver_matrix.row_index(current_nz) = row;
            symmetric_solver_matrix.col_index(current_nz) = col;
            symmetric_solver_matrix.value(current_nz)     = std::move(solver_matrix.value(nz_index));
            current_nz++;
        }
    }

    return symmetric_solver_matrix;
}
// clang-format on

#endif
