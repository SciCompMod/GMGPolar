#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTake/extrapolatedSmootherTake.h"

/* Tridiagonal matrices */
#define UPDATE_TRIDIAGONAL_ELEMENT(matrix, row, column, value)                                                         \
    do {                                                                                                               \
        if (row == column)                                                                                             \
            matrix.main_diagonal(row) = value;                                                                         \
        else if (row == column - 1)                                                                                    \
            matrix.sub_diagonal(row) = value;                                                                          \
        else if (row == 0 && column == matrix.columns() - 1)                                                           \
            matrix.cyclic_corner_element() = value;                                                                    \
    } while (0)

/* Diagonal matrices */
#define UPDATE_DIAGONAL_ELEMENT(matrix, row, column, value)                                                            \
    do {                                                                                                               \
        matrix.diagonal(row) = value;                                                                                  \
    } while (0)

/* Inner Boundary COO/CSR matrix */
#ifdef GMGPOLAR_USE_MUMPS
    #define COO_CSR_UPDATE(matrix, ptr, offset, row, col, val)                                                         \
        do {                                                                                                           \
            matrix.row_index(ptr + offset) = row;                                                                      \
            matrix.col_index(ptr + offset) = col;                                                                      \
            matrix.value(ptr + offset)     = val;                                                                      \
        } while (0)
#else
    #define COO_CSR_UPDATE(matrix, ptr, offset, row, col, val)                                                         \
        do {                                                                                                           \
            matrix.row_nz_index(row, offset) = col;                                                                    \
            matrix.row_nz_entry(row, offset) = val;                                                                    \
        } while (0)
#endif

#define NODE_BUILD_SMOOTHER_TAKE(i_r, i_theta, grid, DirBC_Interior, inner_boundary_circle_matrix,                     \
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
        int ptr, offset;                                                                                               \
        int row, column, col;                                                                                          \
        double value, val;                                                                                             \
        /* ------------------------------------------ */                                                               \
        /* Node in the interior of the Circle Section */                                                               \
        /* ------------------------------------------ */                                                               \
        if (i_r > 0 && i_r < numberSmootherCircles) { /* i_r = numberSmootherCircles-1 is included here! */            \
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
            const int left   = grid.index(i_r - 1, i_theta);                                                           \
            const int bottom = grid.index(i_r, i_theta_M1);                                                            \
            const int center = grid.index(i_r, i_theta);                                                               \
            const int top    = grid.index(i_r, i_theta_P1);                                                            \
            const int right  = grid.index(i_r + 1, i_theta);                                                           \
                                                                                                                       \
            int center_index = i_theta;                                                                                \
            int bottom_index = i_theta_M1;                                                                             \
            int top_index    = i_theta_P1;                                                                             \
            /* -------------------------- */                                                                           \
            /* Cyclic Tridiagonal Section */                                                                           \
            /* i_r % 2 == 1               */                                                                           \
            if (i_r & 1) {                                                                                             \
                /* i_theta % 2 == 1 */                                                                                 \
                /* | X | O | X | */                                                                                    \
                /* |   |   |   | */                                                                                    \
                /* | 0 | Õ | O | */                                                                                   \
                /* |   |   |   | */                                                                                    \
                /* | X | O | X | */                                                                                    \
                /* or */                                                                                               \
                /* i_theta % 2 == 0 */                                                                                 \
                /* | O | O | O | */                                                                                    \
                /* |   |   |   | */                                                                                    \
                /* | X | Õ | X | */                                                                                   \
                /* |   |   |   | */                                                                                    \
                /* | O | O | O | */                                                                                    \
                                                                                                                       \
                auto& matrix = circle_tridiagonal_solver[i_r / 2];                                                     \
                                                                                                                       \
                /* Center: (Left, Right, Bottom, Top) */                                                               \
                row    = center_index;                                                                                 \
                column = center_index;                                                                                 \
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center]) +                        \
                        coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +                     \
                        coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);                      \
                UPDATE_TRIDIAGONAL_ELEMENT(matrix, row, column, value);                                                \
                                                                                                                       \
                /* Bottom */                                                                                           \
                row    = center_index;                                                                                 \
                column = bottom_index;                                                                                 \
                value  = -coeff3 * (att[center] + att[bottom]);                                                        \
                UPDATE_TRIDIAGONAL_ELEMENT(matrix, row, column, value);                                                \
                                                                                                                       \
                /* Top */                                                                                              \
                row    = center_index;                                                                                 \
                column = top_index;                                                                                    \
                value  = -coeff4 * (att[center] + att[top]);                                                           \
                UPDATE_TRIDIAGONAL_ELEMENT(matrix, row, column, value);                                                \
            }                                                                                                          \
            /* ---------------- */                                                                                     \
            /* Diagonal Section */                                                                                     \
            /* i_r % 2 == 0     */                                                                                     \
            else {                                                                                                     \
                /* i_theta % 2 == 1 */                                                                                 \
                /* | O | X | O | */                                                                                    \
                /* |   |   |   | */                                                                                    \
                /* | O | Õ | O | */                                                                                   \
                /* |   |   |   | */                                                                                    \
                /* | O | X | O | */                                                                                    \
                /* or */                                                                                               \
                /* i_theta % 2 == 0 */                                                                                 \
                /* | O | O | O | */                                                                                    \
                /* |   |   |   | */                                                                                    \
                /* | O | X̃ | O | */                                                                                  \
                /* |   |   |   | */                                                                                    \
                /* | O | O | O | */                                                                                    \
                                                                                                                       \
                auto& matrix = circle_diagonal_solver[i_r / 2];                                                        \
                                                                                                                       \
                if (i_theta & 1) { /* i_theta % 2 == 1 */                                                              \
                    /* Center: (Left, Right, Bottom, Top) */                                                           \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center]) +                    \
                            coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +                 \
                            coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);                  \
                    UPDATE_DIAGONAL_ELEMENT(matrix, row, column, value);                                               \
                }                                                                                                      \
                else { /* i_theta % 2 == 0 */                                                                          \
                    /* Center: Coarse */                                                                               \
                    auto& matrix     = circle_diagonal_solver[i_r / 2];                                                \
                    int center_index = i_theta;                                                                        \
                                                                                                                       \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value  = 1.0;                                                                                      \
                    UPDATE_DIAGONAL_ELEMENT(matrix, row, column, value);                                               \
                }                                                                                                      \
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
                auto& matrix              = inner_boundary_circle_matrix;                                              \
                const int center_index    = i_theta;                                                                   \
                const int center_nz_index = getCircleAscIndex(i_r, i_theta);                                           \
                                                                                                                       \
                /* Fill matrix row of (i,j) */                                                                         \
                row = center_index;                                                                                    \
                ptr = center_nz_index;                                                                                 \
                                                                                                                       \
                const Stencil& CenterStencil = getStencil(i_r, i_theta);                                               \
                                                                                                                       \
                offset = CenterStencil[StencilPosition::Center];                                                       \
                col    = center_index;                                                                                 \
                val    = 1.0;                                                                                          \
                COO_CSR_UPDATE(matrix, ptr, offset, row, col, val);                                                    \
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
                auto& matrix = inner_boundary_circle_matrix;                                                           \
                                                                                                                       \
                int nz_index;                                                                                          \
                const Stencil& CenterStencil = getStencil(i_r, i_theta);                                               \
                                                                                                                       \
                if (i_theta & 1) {                                                                                     \
                    /* i_theta % 2 == 1 */                                                                             \
                    /* -| X | O | X | */                                                                               \
                    /* -|   |   |   | */                                                                               \
                    /* -| Õ | O | O | */                                                                              \
                    /* -|   |   |   | */                                                                               \
                    /* -| X | O | X | */                                                                               \
                                                                                                                       \
                    const int left   = grid.index(i_r, i_theta_AcrossOrigin);                                          \
                    const int bottom = grid.index(i_r, i_theta_M1);                                                    \
                    const int center = grid.index(i_r, i_theta);                                                       \
                    const int top    = grid.index(i_r, i_theta_P1);                                                    \
                    const int right  = grid.index(i_r + 1, i_theta);                                                   \
                                                                                                                       \
                    const double center_value =                                                                        \
                        0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center]) +                         \
                        coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +                     \
                        coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);                      \
                    const double left_value = -coeff1 * (arr[center] + arr[left]);                                     \
                                                                                                                       \
                    /* Fill matrix row of (i,j) */                                                                     \
                    row = center_index;                                                                                \
                    ptr = center_nz_index;                                                                             \
                                                                                                                       \
                    const Stencil& CenterStencil = getStencil(i_r, i_theta);                                           \
                                                                                                                       \
                    offset = CenterStencil[StencilPosition::Center];                                                   \
                    col    = center_index;                                                                             \
                    val    = center_value;                                                                             \
                    COO_CSR_UPDATE(matrix, ptr, offset, row, col, val);                                                \
                                                                                                                       \
                    offset = CenterStencil[StencilPosition::Left];                                                     \
                    col    = left_index;                                                                               \
                    val    = left_value;                                                                               \
                    COO_CSR_UPDATE(matrix, ptr, offset, row, col, val);                                                \
                }                                                                                                      \
                else {                                                                                                 \
                    /* i_theta % 2 == 0 */                                                                             \
                    /* -| O | O | O | */                                                                               \
                    /* -|   |   |   | */                                                                               \
                    /* -| X̃ | O | X | */                                                                             \
                    /* -|   |   |   | */                                                                               \
                    /* -| O | O | O | */                                                                               \
                                                                                                                       \
                    /* Fill matrix row of (i,j) */                                                                     \
                    row = center_index;                                                                                \
                    ptr = center_nz_index;                                                                             \
                                                                                                                       \
                    const Stencil& CenterStencil = getStencil(i_r, i_theta);                                           \
                                                                                                                       \
                    offset = CenterStencil[StencilPosition::Center];                                                   \
                    col    = center_index;                                                                             \
                    val    = 1.0;                                                                                      \
                    COO_CSR_UPDATE(matrix, ptr, offset, row, col, val);                                                \
                }                                                                                                      \
            }                                                                                                          \
        }                                                                                                              \
        /* ------------------------------------------ */                                                               \
        /* Node in the interior of the Radial Section */                                                               \
        /* ------------------------------------------ */                                                               \
        else if (i_r > numberSmootherCircles && i_r < grid.nr() - 2) {                                                 \
            double h1 = grid.radialSpacing(i_r - 1);                                                                   \
            double h2 = grid.radialSpacing(i_r);                                                                       \
            double k1 = grid.angularSpacing(i_theta - 1);                                                              \
            double k2 = grid.angularSpacing(i_theta);                                                                  \
                                                                                                                       \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                      \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                      \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                      \
                                                                                                                       \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                   \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                   \
                                                                                                                       \
            const int left   = grid.index(i_r - 1, i_theta);                                                           \
            const int bottom = grid.index(i_r, i_theta_M1);                                                            \
            const int center = grid.index(i_r, i_theta);                                                               \
            const int top    = grid.index(i_r, i_theta_P1);                                                            \
            const int right  = grid.index(i_r + 1, i_theta);                                                           \
                                                                                                                       \
            const int center_index = i_r - numberSmootherCircles;                                                      \
            const int left_index   = i_r - numberSmootherCircles - 1;                                                  \
            const int right_index  = i_r - numberSmootherCircles + 1;                                                  \
            /* ------------------- */                                                                                  \
            /* Tridiagonal Section */                                                                                  \
            /* i_theta % 2 == 1    */                                                                                  \
            if (i_theta & 1) {                                                                                         \
                /* i_r % 2 == 1 */                                                                                     \
                /* ---------- */                                                                                       \
                /* X   O   X  */                                                                                       \
                /* ---------- */                                                                                       \
                /* O   Õ   O  */                                                                                      \
                /* ---------- */                                                                                       \
                /* X   O   X  */                                                                                       \
                /* ---------- */                                                                                       \
                /* or */                                                                                               \
                /* i_r % 2 == 0 */                                                                                     \
                /* ---------- */                                                                                       \
                /* O   X   O  */                                                                                       \
                /* ---------- */                                                                                       \
                /* O   Õ   O  */                                                                                      \
                /* ---------- */                                                                                       \
                /* O   X   O  */                                                                                       \
                /* ---------- */                                                                                       \
                                                                                                                       \
                auto& matrix = radial_tridiagonal_solver[i_theta / 2];                                                 \
                                                                                                                       \
                /* Center: (Left, Right, Bottom, Top) */                                                               \
                row    = center_index;                                                                                 \
                column = center_index;                                                                                 \
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center]) +                        \
                        coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +                     \
                        coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);                      \
                UPDATE_TRIDIAGONAL_ELEMENT(matrix, row, column, value);                                                \
                                                                                                                       \
                /* Left */                                                                                             \
                row    = center_index;                                                                                 \
                column = left_index;                                                                                   \
                value  = -coeff1 * (arr[center] + arr[left]);                                                          \
                UPDATE_TRIDIAGONAL_ELEMENT(matrix, row, column, value);                                                \
                                                                                                                       \
                /* Right */                                                                                            \
                row    = center_index;                                                                                 \
                column = right_index;                                                                                  \
                value  = -coeff2 * (arr[center] + arr[right]);                                                         \
                UPDATE_TRIDIAGONAL_ELEMENT(matrix, row, column, value);                                                \
            }                                                                                                          \
            /* ---------------- */                                                                                     \
            /* Diagonal Section */                                                                                     \
            /* i_theta % 2 == 0 */                                                                                     \
            else {                                                                                                     \
                /* i_r % 2 == 1 */                                                                                     \
                /* ---------- */                                                                                       \
                /* O   O   O  */                                                                                       \
                /* ---------- */                                                                                       \
                /* X   Õ   X  */                                                                                      \
                /* ---------- */                                                                                       \
                /* O   O   O  */                                                                                       \
                /* ---------- */                                                                                       \
                /* or */                                                                                               \
                /* i_r % 2 == 0 */                                                                                     \
                /* ---------- */                                                                                       \
                /* O   O   O  */                                                                                       \
                /* ---------- */                                                                                       \
                /* O   X̃   O  */                                                                                     \
                /* ---------- */                                                                                       \
                /* O   O   O  */                                                                                       \
                /* ---------- */                                                                                       \
                                                                                                                       \
                auto& matrix = radial_diagonal_solver[i_theta / 2];                                                    \
                                                                                                                       \
                if (i_r & 1) { /* i_r % 2 == 1 */                                                                      \
                    /* Center: (Left, Right, Bottom, Top) */                                                           \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center]) +                    \
                            coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +                 \
                            coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);                  \
                    UPDATE_DIAGONAL_ELEMENT(matrix, row, column, value);                                               \
                }                                                                                                      \
                else { /* i_r % 2 == 0 */                                                                              \
                    /* Center: Coarse */                                                                               \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value  = 1.0;                                                                                      \
                    UPDATE_DIAGONAL_ELEMENT(matrix, row, column, value);                                               \
                }                                                                                                      \
            }                                                                                                          \
        }                                                                                                              \
        /* --------------------------------------------- */                                                            \
        /* Radial Section: Node next to circular section */                                                            \
        /* --------------------------------------------- */                                                            \
        else if (i_r == numberSmootherCircles) {                                                                       \
                                                                                                                       \
            double h1 = grid.radialSpacing(i_r - 1);                                                                   \
            double h2 = grid.radialSpacing(i_r);                                                                       \
            double k1 = grid.angularSpacing(i_theta - 1);                                                              \
            double k2 = grid.angularSpacing(i_theta);                                                                  \
                                                                                                                       \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                      \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                      \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                      \
                                                                                                                       \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                   \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                   \
                                                                                                                       \
            const int left   = grid.index(i_r - 1, i_theta);                                                           \
            const int bottom = grid.index(i_r, i_theta_M1);                                                            \
            const int center = grid.index(i_r, i_theta);                                                               \
            const int top    = grid.index(i_r, i_theta_P1);                                                            \
            const int right  = grid.index(i_r + 1, i_theta);                                                           \
                                                                                                                       \
            const int center_index = i_r - numberSmootherCircles;                                                      \
            const int right_index  = i_r - numberSmootherCircles + 1;                                                  \
                                                                                                                       \
            if (i_theta & 1) {                                                                                         \
                /* i_theta % 2 == 1 and i_r % 2 == 1 */                                                                \
                /* | X | O | X || O   X   O   X  */                                                                    \
                /* |   |   |   || -------------- */                                                                    \
                /* | 0 | O | O || Õ   O   O   O  */                                                                   \
                /* |   |   |   || -------------- */                                                                    \
                /* | X | O | X || O   X   O   X  */                                                                    \
                /* or */                                                                                               \
                /* i_theta % 2 == 1 and i_r % 2 == 0 */                                                                \
                /* | O | X | O || X   O   X   O  */                                                                    \
                /* |   |   |   || -------------- */                                                                    \
                /* | 0 | O | O || Õ   O   O   O  */                                                                   \
                /* |   |   |   || -------------- */                                                                    \
                /* | O | X | O || X   O   X   O  */                                                                    \
                                                                                                                       \
                auto& matrix = radial_tridiagonal_solver[i_theta / 2];                                                 \
                                                                                                                       \
                /* Center: (Left, Right, Bottom, Top) */                                                               \
                row    = center_index;                                                                                 \
                column = center_index;                                                                                 \
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center]) +                        \
                        coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +                     \
                        coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);                      \
                UPDATE_TRIDIAGONAL_ELEMENT(matrix, row, column, value);                                                \
                                                                                                                       \
                /* Right */                                                                                            \
                row    = center_index;                                                                                 \
                column = right_index;                                                                                  \
                value  = -coeff2 * (arr[center] + arr[right]);                                                         \
                UPDATE_TRIDIAGONAL_ELEMENT(matrix, row, column, value);                                                \
            }                                                                                                          \
            else {                                                                                                     \
                                                                                                                       \
                auto& matrix = radial_diagonal_solver[i_theta / 2];                                                    \
                                                                                                                       \
                if (i_r & 1) {                                                                                         \
                    /* i_theta % 2 == 0 and i_r % 2 == 1 */                                                            \
                    /* | O | O | O || O   O   O   O  */                                                                \
                    /* |   |   |   || -------------- */                                                                \
                    /* | X | O | X || Õ   X   O   X  */                                                               \
                    /* |   |   |   || -------------- */                                                                \
                    /* | O | O | O || O   O   O   O  */                                                                \
                                                                                                                       \
                    /* Center: (Left, Right, Bottom, Top) */                                                           \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center]) +                    \
                            coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +                 \
                            coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);                  \
                    UPDATE_DIAGONAL_ELEMENT(matrix, row, column, value);                                               \
                }                                                                                                      \
                else {                                                                                                 \
                    /* i_theta % 2 == 0 and i_r % 2 == 0 */                                                            \
                    /* | O | O | O || O   O   O   O  */                                                                \
                    /* |   |   |   || -------------- */                                                                \
                    /* | O | X | O || X̃   O   X   O  */                                                              \
                    /* |   |   |   || -------------- */                                                                \
                    /* | O | O | O || O   O   O   O  */                                                                \
                    /* Center: Coarse */                                                                               \
                    row    = center_index;                                                                             \
                    column = center_index;                                                                             \
                    value  = 1.0;                                                                                      \
                    UPDATE_DIAGONAL_ELEMENT(matrix, row, column, value);                                               \
                }                                                                                                      \
            }                                                                                                          \
        }                                                                                                              \
        /* ------------------------------------------- */                                                              \
        /* Radial Section: Node next to outer boundary */                                                              \
        /* ------------------------------------------- */                                                              \
        else if (i_r == grid.nr() - 2) {                                                                               \
            assert(i_r % 2 == 1);                                                                                      \
                                                                                                                       \
            double h1 = grid.radialSpacing(i_r - 1);                                                                   \
            double h2 = grid.radialSpacing(i_r);                                                                       \
            double k1 = grid.angularSpacing(i_theta - 1);                                                              \
            double k2 = grid.angularSpacing(i_theta);                                                                  \
                                                                                                                       \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                      \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                      \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                      \
                                                                                                                       \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                   \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                   \
                                                                                                                       \
            const int left   = grid.index(i_r - 1, i_theta);                                                           \
            const int bottom = grid.index(i_r, i_theta_M1);                                                            \
            const int center = grid.index(i_r, i_theta);                                                               \
            const int top    = grid.index(i_r, i_theta_P1);                                                            \
            const int right  = grid.index(i_r + 1, i_theta);                                                           \
                                                                                                                       \
            const int center_index = i_r - numberSmootherCircles;                                                      \
            const int left_index   = i_r - numberSmootherCircles - 1;                                                  \
            const int right_index  = i_r - numberSmootherCircles + 1;                                                  \
                                                                                                                       \
            if (i_theta & 1) {                                                                                         \
                /* i_theta % 2 == 1 */                                                                                 \
                /* ---------------|| */                                                                                \
                /* O   X   O   X  || */                                                                                \
                /* ---------------|| */                                                                                \
                /* O   O   Õ   O  || */                                                                               \
                /* ---------------|| */                                                                                \
                /* O   X   O   X  || */                                                                                \
                /* ---------------|| */                                                                                \
                                                                                                                       \
                auto& matrix = radial_tridiagonal_solver[i_theta / 2];                                                 \
                                                                                                                       \
                /* Center: (Left, Right, Bottom, Top) */                                                               \
                row    = center_index;                                                                                 \
                column = center_index;                                                                                 \
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center]) +                        \
                        coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +                     \
                        coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);                      \
                UPDATE_TRIDIAGONAL_ELEMENT(matrix, row, column, value);                                                \
                                                                                                                       \
                /* Left */                                                                                             \
                row    = center_index;                                                                                 \
                column = left_index;                                                                                   \
                value  = -coeff1 * (arr[center] + arr[left]);                                                          \
                UPDATE_TRIDIAGONAL_ELEMENT(matrix, row, column, value);                                                \
                /* Right */                                                                                            \
                row    = center_index;                                                                                 \
                column = right_index;                                                                                  \
                value  = 0.0; /* Make tridiagonal matrix symmetric */                                                  \
                UPDATE_TRIDIAGONAL_ELEMENT(matrix, row, column, value);                                                \
            }                                                                                                          \
            else {                                                                                                     \
                /* i_theta % 2 == 0 */                                                                                 \
                /* ---------------|| */                                                                                \
                /* O   O   O   O  || */                                                                                \
                /* ---------------|| */                                                                                \
                /* O   X   Õ   X  || */                                                                               \
                /* ---------------|| */                                                                                \
                /* O   O   O   O  || */                                                                                \
                /* ---------------|| */                                                                                \
                                                                                                                       \
                auto& matrix = radial_diagonal_solver[i_theta / 2];                                                    \
                                                                                                                       \
                /* Center: (Left, Right, Bottom, Top) */                                                               \
                row    = center_index;                                                                                 \
                column = center_index;                                                                                 \
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center]) +                        \
                        coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +                     \
                        coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);                      \
                UPDATE_DIAGONAL_ELEMENT(matrix, row, column, value);                                                   \
            }                                                                                                          \
        }                                                                                                              \
        /* ------------------------------------------ */                                                               \
        /* Radial Section: Node on the outer boundary */                                                               \
        /* ------------------------------------------ */                                                               \
        else if (i_r == grid.nr() - 1) {                                                                               \
            assert(!i_r % 2 == 0);                                                                                     \
                                                                                                                       \
            int center_index = i_r - numberSmootherCircles;                                                            \
            int left_index   = i_r - numberSmootherCircles - 1;                                                        \
                                                                                                                       \
            if (i_theta & 1) {                                                                                         \
                /* i_theta % 2 == 1 */                                                                                 \
                /* -----------|| */                                                                                    \
                /* X   O   X  || */                                                                                    \
                /* -----------|| */                                                                                    \
                /* O   O   Õ  || */                                                                                   \
                /* -----------|| */                                                                                    \
                /* X   O   X  || */                                                                                    \
                /* -----------|| */                                                                                    \
                                                                                                                       \
                auto& matrix = radial_tridiagonal_solver[i_theta / 2];                                                 \
                                                                                                                       \
                /* Fill matrix row of (i,j) */                                                                         \
                row    = center_index;                                                                                 \
                column = center_index;                                                                                 \
                value  = 1.0;                                                                                          \
                UPDATE_TRIDIAGONAL_ELEMENT(matrix, row, column, value);                                                \
                                                                                                                       \
                row    = center_index;                                                                                 \
                column = left_index;                                                                                   \
                value  = 0.0; /* Make tridiagonal matrix symmetric */                                                  \
                UPDATE_TRIDIAGONAL_ELEMENT(matrix, row, column, value);                                                \
            }                                                                                                          \
            else {                                                                                                     \
                /* i_theta % 2 == 0 */                                                                                 \
                /* -----------|| */                                                                                    \
                /* O   O   O  || */                                                                                    \
                /* -----------|| */                                                                                    \
                /* X   O   X̃  || */                                                                                  \
                /* -----------|| */                                                                                    \
                /* O   O   O  || */                                                                                    \
                /* -----------|| */                                                                                    \
                                                                                                                       \
                auto& matrix = radial_diagonal_solver[i_theta / 2];                                                    \
                                                                                                                       \
                /* Fill matrix row of (i,j) */                                                                         \
                row    = center_index;                                                                                 \
                column = center_index;                                                                                 \
                value  = 1.0;                                                                                          \
                UPDATE_DIAGONAL_ELEMENT(matrix, row, column, value);                                                   \
            }                                                                                                          \
        }                                                                                                              \
    } while (0)

void ExtrapolatedSmootherTake::buildAscCircleSection(const int i_r)
{
    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    const auto& arr        = level_cache_.arr();
    const auto& att        = level_cache_.att();
    const auto& art        = level_cache_.art();
    const auto& detDF      = level_cache_.detDF();
    const auto& coeff_beta = level_cache_.coeff_beta();

    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        // Build Asc at the current node
        NODE_BUILD_SMOOTHER_TAKE(i_r, i_theta, grid_, DirBC_Interior_, inner_boundary_circle_matrix_,
                                 circle_diagonal_solver_, circle_tridiagonal_solver_, radial_diagonal_solver_,
                                 radial_tridiagonal_solver_);
    }
}

void ExtrapolatedSmootherTake::buildAscRadialSection(const int i_theta)
{
    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    const auto& arr        = level_cache_.arr();
    const auto& att        = level_cache_.att();
    const auto& art        = level_cache_.art();
    const auto& detDF      = level_cache_.detDF();
    const auto& coeff_beta = level_cache_.coeff_beta();

    for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
        // Build Asc at the current node
        NODE_BUILD_SMOOTHER_TAKE(i_r, i_theta, grid_, DirBC_Interior_, inner_boundary_circle_matrix_,
                                 circle_diagonal_solver_, circle_tridiagonal_solver_, radial_diagonal_solver_,
                                 radial_tridiagonal_solver_);
    }
}

/* clang-format off */
void ExtrapolatedSmootherTake::buildAscMatrices()
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
                #ifdef GMGPOLAR_USE_MUMPS
                // Although the matrix is symmetric, we need to store all its entries, so we disable the symmetry.
                const int nnz                 = getNonZeroCountCircleAsc(circle_Asc_index);
                inner_boundary_circle_matrix_ = SparseMatrixCOO<double>(num_circle_nodes, num_circle_nodes, nnz);
                inner_boundary_circle_matrix_.is_symmetric(false);
                #else
                std::function<int(int)> nnz_per_row = [&](int i_theta) {
                    if(DirBC_Interior_) return 1;
                    else return i_theta % 2 == 0 ? 1 : 2;
                };
                inner_boundary_circle_matrix_ = SparseMatrixCSR<double>(num_circle_nodes, num_circle_nodes, nnz_per_row);
                #endif
            }

            /* Interior Circle Section */
            else {
                if (circle_Asc_index & 1) {
                    const int circle_tridiagonal_solver_index = circle_Asc_index / 2;
                    auto& solver_matrix = circle_tridiagonal_solver_[circle_tridiagonal_solver_index];
                    solver_matrix       = SymmetricTridiagonalSolver<double>(num_circle_nodes);
                    solver_matrix.is_cyclic(true);
                }
                else {
                    const int circle_diagonal_solver_index = circle_Asc_index / 2;
                    auto& solver_matrix                    = circle_diagonal_solver_[circle_diagonal_solver_index];
                    solver_matrix                          = DiagonalSolver<double>(num_circle_nodes);
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
                solver_matrix                             = SymmetricTridiagonalSolver<double>(num_radial_nodes);
                solver_matrix.is_cyclic(false);
            }
            else {
                const int radial_diagonal_solver_index = radial_Asc_index / 2;
                auto& solver_matrix                    = radial_diagonal_solver_[radial_diagonal_solver_index];
                solver_matrix                          = DiagonalSolver<double>(num_radial_nodes);
            }
        }
    }

    /* ---------------------------------- */
    /* Part 2: Fill Asc Smoother matrices */
    /* ---------------------------------- */

    #pragma omp parallel
    {
        #pragma omp for nowait
        for (int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
            buildAscCircleSection(i_r);
        }

        #pragma omp for nowait
        for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
            buildAscRadialSection(i_theta);
        }
    }

    #ifdef GMGPOLAR_USE_MUMPS
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
    #endif
}
/* clang-format on */