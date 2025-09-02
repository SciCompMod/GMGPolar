#include "../../../include/Smoother/SmootherGive/smootherGive.h"

#include "../../../include/common/geometry_helper.h"

/* Tridiagonal matrices */
#define UPDATE_MATRIX_ELEMENT(matrix, row, column, value)                                                              \
    do {                                                                                                               \
        if (row == column)                                                                                             \
            matrix.main_diagonal(row) += value;                                                                        \
        else if (row == column - 1)                                                                                    \
            matrix.sub_diagonal(row) += value;                                                                         \
        else if (row == 0 && column == matrix.columns() - 1)                                                           \
            matrix.cyclic_corner_element() += value;                                                                   \
    } while (0)

/* Inner Boundary COO/CSR matrix */
#ifdef GMGPOLAR_USE_MUMPS
    #define COO_CSR_UPDATE(matrix, ptr, offset, row, col, val)                                                         \
        do {                                                                                                           \
            matrix.row_index(ptr + offset) = row;                                                                      \
            matrix.col_index(ptr + offset) = col;                                                                      \
            matrix.value(ptr + offset) += val;                                                                         \
        } while (0)
#else
    #define COO_CSR_UPDATE(matrix, ptr, offset, row, col, val)                                                         \
        do {                                                                                                           \
            matrix.row_nz_index(row, offset) = col;                                                                    \
            matrix.row_nz_entry(row, offset) += val;                                                                   \
        } while (0)
#endif

#define NODE_BUILD_SMOOTHER_GIVE(i_r, i_theta, grid, DirBC_Interior, inner_boundary_circle_matrix,                     \
                                 circle_tridiagonal_solver, radial_tridiagonal_solver)                                 \
    do {                                                                                                               \
        assert(i_r >= 0 && i_r < grid.nr());                                                                           \
        assert(i_theta >= 0 && i_theta < grid.ntheta());                                                               \
                                                                                                                       \
        const int numberSmootherCircles = grid.numberSmootherCircles();                                                \
        const int lengthSmootherRadial  = grid.lengthSmootherRadial();                                                 \
                                                                                                                       \
        assert(numberSmootherCircles >= 2);                                                                            \
        assert(lengthSmootherRadial >= 3);                                                                             \
                                                                                                                       \
        int ptr, offset;                                                                                               \
        int row, column, col;                                                                                          \
        double value, val;                                                                                             \
                                                                                                                       \
        /* ------------------------------------------ */                                                               \
        /* Node in the interior of the Circle Section */                                                               \
        /* ------------------------------------------ */                                                               \
        if (i_r > 0 && i_r < numberSmootherCircles) {                                                                  \
            const double h1     = grid.radialSpacing(i_r - 1);                                                         \
            const double h2     = grid.radialSpacing(i_r);                                                             \
            const double k1     = grid.angularSpacing(i_theta - 1);                                                    \
            const double k2     = grid.angularSpacing(i_theta);                                                        \
            const double coeff1 = 0.5 * (k1 + k2) / h1;                                                                \
            const double coeff2 = 0.5 * (k1 + k2) / h2;                                                                \
            const double coeff3 = 0.5 * (h1 + h2) / k1;                                                                \
            const double coeff4 = 0.5 * (h1 + h2) / k2;                                                                \
                                                                                                                       \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                   \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                   \
                                                                                                                       \
            const int center_index = i_theta;                                                                          \
            const int left_index   = i_theta;                                                                          \
            const int right_index  = (i_r + 1 == numberSmootherCircles) ? 0 : i_theta;                                 \
            const int bottom_index = i_theta_M1;                                                                       \
            const int top_index    = i_theta_P1;                                                                       \
                                                                                                                       \
            /* Visualization of the sourrounding tridiagonal matrices. */                                              \
            /* left_matrix, center_matrix, right_matrix */                                                             \
            /* | O | O | O | */                                                                                        \
            /* |   |   |   | */                                                                                        \
            /* | O | Õ | O | */                                                                                        \
            /* |   |   |   | */                                                                                        \
            /* | O | O | O | */                                                                                        \
            /* or */                                                                                                   \
            /* left_matrix, right_matrix */                                                                            \
            /* | O | O | O || O   O   O   O  */                                                                        \
            /* |   |   |   || -------------- */                                                                        \
            /* | O | O | Õ || O   O   O   O  <- right_matrix */                                                        \
            /* |   |   |   || -------------- */                                                                        \
            /* | O | O | O || O   O   O   O  */                                                                        \
            auto& left_matrix   = circle_tridiagonal_solver[i_r - 1];                                                  \
            auto& center_matrix = circle_tridiagonal_solver[i_r];                                                      \
            auto& right_matrix  = (i_r + 1 == numberSmootherCircles) ? radial_tridiagonal_solver[i_theta]              \
                                                                     : circle_tridiagonal_solver[i_r + 1];             \
                                                                                                                       \
            /* Fill matrix row of (i,j) */                                                                             \
            row    = center_index;                                                                                     \
            column = center_index;                                                                                     \
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */                 \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            row    = center_index;                                                                                     \
            column = bottom_index;                                                                                     \
            value  = -coeff3 * att; /* Bottom */                                                                       \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            row    = center_index;                                                                                     \
            column = top_index;                                                                                        \
            value  = -coeff4 * att; /* Top */                                                                          \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            row    = center_index;                                                                                     \
            column = center_index;                                                                                     \
            value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */       \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            /* Fill matrix row of (i,j-1) */                                                                           \
            row    = bottom_index;                                                                                     \
            column = center_index;                                                                                     \
            value  = -coeff3 * att; /* Top */                                                                          \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            row    = bottom_index;                                                                                     \
            column = bottom_index;                                                                                     \
            value  = coeff3 * att; /* Center: (Top) */                                                                 \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            /* Fill matrix row of (i,j+1) */                                                                           \
            row    = top_index;                                                                                        \
            column = center_index;                                                                                     \
            value  = -coeff4 * att; /* Bottom */                                                                       \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            row    = top_index;                                                                                        \
            column = top_index;                                                                                        \
            value  = coeff4 * att; /* Center: (Bottom) */                                                              \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            /* Fill matrix row of (i-1,j) */                                                                           \
            if (!DirBC_Interior && i_r == 1) {                                                                         \
                row = left_index;                                                                                      \
                ptr = getCircleAscIndex(i_r - 1, i_theta);                                                             \
                                                                                                                       \
                const Stencil& LeftStencil = getStencil(i_r - 1);                                                      \
                                                                                                                       \
                offset = LeftStencil[StencilPosition::Center];                                                         \
                col    = left_index;                                                                                   \
                val    = +coeff1 * arr; /* Center: (Right) */                                                          \
                COO_CSR_UPDATE(inner_boundary_circle_matrix, ptr, offset, row, col, val);                              \
            }                                                                                                          \
            if (i_r > 1) {                                                                                             \
                row    = left_index;                                                                                   \
                column = left_index;                                                                                   \
                value  = coeff1 * arr; /* Center: (Right) */                                                           \
                UPDATE_MATRIX_ELEMENT(left_matrix, row, column, value);                                                \
            }                                                                                                          \
                                                                                                                       \
            /* Fill matrix row of (i+1,j) */                                                                           \
            row    = right_index;                                                                                      \
            column = right_index;                                                                                      \
            value  = coeff2 * arr; /* Center: (Left) */                                                                \
            UPDATE_MATRIX_ELEMENT(right_matrix, row, column, value);                                                   \
        }                                                                                                              \
        /* ------------------------------------------ */                                                               \
        /* Node in the interior of the Radial Section */                                                               \
        /* ------------------------------------------ */                                                               \
        else if (i_r > numberSmootherCircles && i_r < grid.nr() - 2) {                                                 \
            const double h1     = grid.radialSpacing(i_r - 1);                                                         \
            const double h2     = grid.radialSpacing(i_r);                                                             \
            const double k1     = grid.angularSpacing(i_theta - 1);                                                    \
            const double k2     = grid.angularSpacing(i_theta);                                                        \
            const double coeff1 = 0.5 * (k1 + k2) / h1;                                                                \
            const double coeff2 = 0.5 * (k1 + k2) / h2;                                                                \
            const double coeff3 = 0.5 * (h1 + h2) / k1;                                                                \
            const double coeff4 = 0.5 * (h1 + h2) / k2;                                                                \
                                                                                                                       \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                   \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                   \
                                                                                                                       \
            /* ---------- */                                                                                           \
            /* O   O   O  <- top_matrix */                                                                             \
            /* ---------- */                                                                                           \
            /* O   Õ   O  <- center_matrix */                                                                          \
            /* ---------- */                                                                                           \
            /* O   O   O  <- bottom_matrix */                                                                          \
            /* ---------- */                                                                                           \
            auto& bottom_matrix = radial_tridiagonal_solver[i_theta_M1];                                               \
            auto& center_matrix = radial_tridiagonal_solver[i_theta];                                                  \
            auto& top_matrix    = radial_tridiagonal_solver[i_theta_P1];                                               \
                                                                                                                       \
            const int center_index = i_r - numberSmootherCircles;                                                      \
            const int left_index   = i_r - numberSmootherCircles - 1;                                                  \
            const int right_index  = i_r - numberSmootherCircles + 1;                                                  \
            const int bottom_index = i_r - numberSmootherCircles;                                                      \
            const int top_index    = i_r - numberSmootherCircles;                                                      \
                                                                                                                       \
            /* Fill matrix row of (i,j) */                                                                             \
            row    = center_index;                                                                                     \
            column = center_index;                                                                                     \
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */                 \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            row    = center_index;                                                                                     \
            column = left_index;                                                                                       \
            value  = -coeff1 * arr; /* Left */                                                                         \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            row    = center_index;                                                                                     \
            column = right_index;                                                                                      \
            value  = -coeff2 * arr; /* Right */                                                                        \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            row    = center_index;                                                                                     \
            column = center_index;                                                                                     \
            value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */       \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            /* Fill matrix row of (i-1,j) */                                                                           \
            row    = left_index;                                                                                       \
            column = center_index;                                                                                     \
            value  = -coeff1 * arr; /* Right */                                                                        \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            row    = left_index;                                                                                       \
            column = left_index;                                                                                       \
            value  = coeff1 * arr; /* Center: (Right) */                                                               \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            /* Fill matrix row of (i+1,j) */                                                                           \
            row    = right_index;                                                                                      \
            column = center_index;                                                                                     \
            value  = -coeff2 * arr; /* Left */                                                                         \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            row    = right_index;                                                                                      \
            column = right_index;                                                                                      \
            value  = coeff2 * arr; /* Center: (Left) */                                                                \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            /* Fill matrix row of (i,j-1) */                                                                           \
            row    = bottom_index;                                                                                     \
            column = bottom_index;                                                                                     \
            value  = coeff3 * att; /* Center: (Top) */                                                                 \
            UPDATE_MATRIX_ELEMENT(bottom_matrix, row, column, value);                                                  \
                                                                                                                       \
            /* Fill matrix row of (i,j+1) */                                                                           \
            row    = top_index;                                                                                        \
            column = top_index;                                                                                        \
            value  = coeff4 * att; /* Center: (Bottom) */                                                              \
            UPDATE_MATRIX_ELEMENT(top_matrix, row, column, value);                                                     \
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
                const double h2     = grid.radialSpacing(i_r);                                                         \
                const double k1     = grid.angularSpacing(i_theta - 1);                                                \
                const double k2     = grid.angularSpacing(i_theta);                                                    \
                const double coeff2 = 0.5 * (k1 + k2) / h2;                                                            \
                                                                                                                       \
                const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                               \
                const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                               \
                                                                                                                       \
                auto& center_matrix = inner_boundary_circle_matrix;                                                    \
                auto& right_matrix  = circle_tridiagonal_solver[i_r + 1];                                              \
                                                                                                                       \
                const int center_index = i_theta;                                                                      \
                const int right_index  = i_theta;                                                                      \
                const int bottom_index = i_theta_M1;                                                                   \
                const int top_index    = i_theta_P1;                                                                   \
                                                                                                                       \
                /* Fill matrix row of (i,j) */                                                                         \
                row = center_index;                                                                                    \
                ptr = getCircleAscIndex(i_r, i_theta);                                                                 \
                                                                                                                       \
                const Stencil& CenterStencil = getStencil(i_r);                                                        \
                                                                                                                       \
                offset = CenterStencil[StencilPosition::Center];                                                       \
                col    = center_index;                                                                                 \
                val    = 1.0;                                                                                          \
                COO_CSR_UPDATE(inner_boundary_circle_matrix, ptr, offset, row, col, val);                              \
                                                                                                                       \
                /* Fill matrix row of (i+1,j) */                                                                       \
                row    = right_index;                                                                                  \
                column = right_index;                                                                                  \
                value  = coeff2 * arr; /* Center: (Left) */                                                            \
                UPDATE_MATRIX_ELEMENT(right_matrix, row, column, value);                                               \
            }                                                                                                          \
            else {                                                                                                     \
                /* ------------------------------------------------------------- */                                    \
                /* Case 2: Across origin discretization on the interior boundary */                                    \
                /* ------------------------------------------------------------- */                                    \
                /* h1 gets replaced with 2 * R0. */                                                                    \
                /* (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()>>1)). */                          \
                /* Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil. */    \
                const double h1     = 2 * grid.radius(0);                                                              \
                const double h2     = grid.radialSpacing(i_r);                                                         \
                const double k1     = grid.angularSpacing(i_theta - 1);                                                \
                const double k2     = grid.angularSpacing(i_theta);                                                    \
                const double coeff1 = 0.5 * (k1 + k2) / h1;                                                            \
                const double coeff2 = 0.5 * (k1 + k2) / h2;                                                            \
                const double coeff3 = 0.5 * (h1 + h2) / k1;                                                            \
                const double coeff4 = 0.5 * (h1 + h2) / k2;                                                            \
                                                                                                                       \
                /* left_matrix (across-the origin), center_matrix, right_matrix */                                     \
                /* -| X | O | X | */                                                                                   \
                /* -|   |   |   | */                                                                                   \
                /* -| Õ | O | O | */                                                                                   \
                /* -|   |   |   | */                                                                                   \
                /* -| X | O | X | */                                                                                   \
                auto& center_matrix = inner_boundary_circle_matrix;                                                    \
                auto& right_matrix  = circle_tridiagonal_solver[i_r + 1];                                              \
                auto& left_matrix   = inner_boundary_circle_matrix;                                                    \
                                                                                                                       \
                const int i_theta_M1           = grid.wrapThetaIndex(i_theta - 1);                                     \
                const int i_theta_P1           = grid.wrapThetaIndex(i_theta + 1);                                     \
                const int i_theta_AcrossOrigin = grid.wrapThetaIndex(i_theta + grid.ntheta() / 2);                     \
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
                /* Fill matrix row of (i,j) */                                                                         \
                row = center_index;                                                                                    \
                ptr = center_nz_index;                                                                                 \
                                                                                                                       \
                const Stencil& CenterStencil = getStencil(i_r);                                                        \
                                                                                                                       \
                offset = CenterStencil[StencilPosition::Center];                                                       \
                col    = center_index;                                                                                 \
                val    = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */                     \
                COO_CSR_UPDATE(inner_boundary_circle_matrix, ptr, offset, row, col, val);                              \
                                                                                                                       \
                offset = CenterStencil[StencilPosition::Left];                                                         \
                col    = left_index;                                                                                   \
                val    = -coeff1 * arr; /* Left */                                                                     \
                COO_CSR_UPDATE(inner_boundary_circle_matrix, ptr, offset, row, col, val);                              \
                                                                                                                       \
                offset = CenterStencil[StencilPosition::Bottom];                                                       \
                col    = bottom_index;                                                                                 \
                val    = -coeff3 * att; /* Bottom */                                                                   \
                COO_CSR_UPDATE(inner_boundary_circle_matrix, ptr, offset, row, col, val);                              \
                                                                                                                       \
                offset = CenterStencil[StencilPosition::Top];                                                          \
                col    = top_index;                                                                                    \
                val    = -coeff4 * att; /* Top */                                                                      \
                COO_CSR_UPDATE(inner_boundary_circle_matrix, ptr, offset, row, col, val);                              \
                                                                                                                       \
                offset = CenterStencil[StencilPosition::Center];                                                       \
                col    = center_index;                                                                                 \
                val    = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */   \
                COO_CSR_UPDATE(inner_boundary_circle_matrix, ptr, offset, row, col, val);                              \
                                                                                                                       \
                /* Fill matrix row of (i-1,j) */                                                                       \
                /* From view the view of the across origin node, */                                                    \
                /* the directions are roatated by 180 degrees in the stencil! */                                       \
                row = left_index;                                                                                      \
                ptr = left_nz_index;                                                                                   \
                                                                                                                       \
                const Stencil& LeftStencil = CenterStencil;                                                            \
                                                                                                                       \
                offset = LeftStencil[StencilPosition::Left];                                                           \
                col    = center_index;                                                                                 \
                val    = -coeff1 * arr; /* Right -> Left*/                                                             \
                COO_CSR_UPDATE(inner_boundary_circle_matrix, ptr, offset, row, col, val);                              \
                                                                                                                       \
                offset = LeftStencil[StencilPosition::Center];                                                         \
                col    = left_index;                                                                                   \
                val    = +coeff1 * arr; /* Center: (Right) -> Center: (Left) */                                        \
                COO_CSR_UPDATE(inner_boundary_circle_matrix, ptr, offset, row, col, val);                              \
                                                                                                                       \
                /* Top Right -> Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */                               \
                                                                                                                       \
                /* Bottom Right -> Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */                               \
                                                                                                                       \
                /* Fill matrix row of (i+1,j) */                                                                       \
                row    = right_index;                                                                                  \
                column = right_index;                                                                                  \
                value  = coeff2 * arr; /* Center: (Left) */                                                            \
                UPDATE_MATRIX_ELEMENT(right_matrix, row, column, value);                                               \
                                                                                                                       \
                /* Fill matrix row of (i,j-1) */                                                                       \
                row = bottom_index;                                                                                    \
                ptr = bottom_nz_index;                                                                                 \
                                                                                                                       \
                const Stencil& BottomStencil = CenterStencil;                                                          \
                                                                                                                       \
                offset = BottomStencil[StencilPosition::Top];                                                          \
                col    = center_index;                                                                                 \
                val    = -coeff3 * att; /* Top */                                                                      \
                COO_CSR_UPDATE(inner_boundary_circle_matrix, ptr, offset, row, col, val);                              \
                                                                                                                       \
                offset = BottomStencil[StencilPosition::Center];                                                       \
                col    = bottom_index;                                                                                 \
                val    = +coeff3 * att; /* Center: (Top) */                                                            \
                COO_CSR_UPDATE(inner_boundary_circle_matrix, ptr, offset, row, col, val);                              \
                                                                                                                       \
                /* TopLeft: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */                                                \
                                                                                                                       \
                /* Fill matrix row of (i,j+1) */                                                                       \
                row = top_index;                                                                                       \
                ptr = top_nz_index;                                                                                    \
                                                                                                                       \
                const Stencil& TopStencil = CenterStencil;                                                             \
                                                                                                                       \
                offset = TopStencil[StencilPosition::Bottom];                                                          \
                col    = center_index;                                                                                 \
                val    = -coeff4 * att; /* Bottom */                                                                   \
                COO_CSR_UPDATE(inner_boundary_circle_matrix, ptr, offset, row, col, val);                              \
                                                                                                                       \
                offset = TopStencil[StencilPosition::Center];                                                          \
                col    = top_index;                                                                                    \
                val    = +coeff4 * att; /* Center: (Bottom) */                                                         \
                COO_CSR_UPDATE(inner_boundary_circle_matrix, ptr, offset, row, col, val);                              \
                                                                                                                       \
                /* BottomLeft: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */                                             \
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
            /* | O | O | O || O   O   O   O  <- top_matrix */                                                          \
            /* |   |   |   || -------------- */                                                                        \
            /* | O | O | O || Õ   O   O   O  <- center_matrix */                                                       \
            /* |   |   |   || -------------- */                                                                        \
            /* | O | O | O || O   O   O   O  <- bottom_matrix */                                                       \
            auto& bottom_matrix = radial_tridiagonal_solver[i_theta_M1];                                               \
            auto& center_matrix = radial_tridiagonal_solver[i_theta];                                                  \
            auto& top_matrix    = radial_tridiagonal_solver[i_theta_P1];                                               \
            auto& left_matrix   = circle_tridiagonal_solver[i_r - 1];                                                  \
                                                                                                                       \
            const int center_index = i_r - numberSmootherCircles;                                                      \
            const int left_index   = i_theta;                                                                          \
            const int right_index  = i_r - numberSmootherCircles + 1;                                                  \
            const int bottom_index = i_r - numberSmootherCircles;                                                      \
            const int top_index    = i_r - numberSmootherCircles;                                                      \
                                                                                                                       \
            /* Fill matrix row of (i,j) */                                                                             \
            row    = center_index;                                                                                     \
            column = center_index;                                                                                     \
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */                 \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            row    = center_index;                                                                                     \
            column = right_index;                                                                                      \
            value  = -coeff2 * arr; /* Right */                                                                        \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            row    = center_index;                                                                                     \
            column = center_index;                                                                                     \
            value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */       \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            /* Fill matrix row of (i-1,j) */                                                                           \
            row    = left_index;                                                                                       \
            column = left_index;                                                                                       \
            value  = coeff1 * arr; /* Center: (Right) */                                                               \
            UPDATE_MATRIX_ELEMENT(left_matrix, row, column, value);                                                    \
                                                                                                                       \
            /* Fill matrix row of (i+1,j) */                                                                           \
            row    = right_index;                                                                                      \
            column = center_index;                                                                                     \
            value  = -coeff2 * arr; /* Left */                                                                         \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            row    = right_index;                                                                                      \
            column = right_index;                                                                                      \
            value  = coeff2 * arr; /* Center: (Left) */                                                                \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            /* Fill matrix row of (i,j-1) */                                                                           \
            row    = bottom_index;                                                                                     \
            column = bottom_index;                                                                                     \
            value  = coeff3 * att; /* Center: (Top) */                                                                 \
            UPDATE_MATRIX_ELEMENT(bottom_matrix, row, column, value);                                                  \
                                                                                                                       \
            /* Fill matrix row of (i,j+1) */                                                                           \
            row    = top_index;                                                                                        \
            column = top_index;                                                                                        \
            value  = coeff4 * att; /* Center: (Bottom) */                                                              \
            UPDATE_MATRIX_ELEMENT(top_matrix, row, column, value);                                                     \
        }                                                                                                              \
        /* ------------------------------------------- */                                                              \
        /* Radial Section: Node next to outer boundary */                                                              \
        /* ------------------------------------------- */                                                              \
        else if (i_r == grid.nr() - 2) {                                                                               \
            const double h1     = grid.radialSpacing(i_r - 1);                                                         \
            const double h2     = grid.radialSpacing(i_r);                                                             \
            const double k1     = grid.angularSpacing(i_theta - 1);                                                    \
            const double k2     = grid.angularSpacing(i_theta);                                                        \
            const double coeff1 = 0.5 * (k1 + k2) / h1;                                                                \
            const double coeff2 = 0.5 * (k1 + k2) / h2;                                                                \
            const double coeff3 = 0.5 * (h1 + h2) / k1;                                                                \
            const double coeff4 = 0.5 * (h1 + h2) / k2;                                                                \
                                                                                                                       \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                   \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                   \
                                                                                                                       \
            /* ---------------|| */                                                                                    \
            /* O   O   O   O  || <- top_matrix */                                                                      \
            /* ---------------|| */                                                                                    \
            /* O   O   Õ   O  || <- center_matrix */                                                                   \
            /* ---------------|| */                                                                                    \
            /* O   O   O   O  || <- bottom_matrix */                                                                   \
            /* ---------------|| */                                                                                    \
            auto& bottom_matrix = radial_tridiagonal_solver[i_theta_M1];                                               \
            auto& center_matrix = radial_tridiagonal_solver[i_theta];                                                  \
            auto& top_matrix    = radial_tridiagonal_solver[i_theta_P1];                                               \
                                                                                                                       \
            const int center_index = i_r - numberSmootherCircles;                                                      \
            const int left_index   = i_r - numberSmootherCircles - 1;                                                  \
            const int right_index  = i_r - numberSmootherCircles + 1;                                                  \
            const int bottom_index = i_r - numberSmootherCircles;                                                      \
            const int top_index    = i_r - numberSmootherCircles;                                                      \
                                                                                                                       \
            /* ---------------------------- */                                                                         \
            /* Give values to center matrix */                                                                         \
            /* ---------------------------- */                                                                         \
            /* Fill matrix row of (i,j) */                                                                             \
            row    = center_index;                                                                                     \
            column = center_index;                                                                                     \
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */                 \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            row    = center_index;                                                                                     \
            column = left_index;                                                                                       \
            value  = -coeff1 * arr; /* Left */                                                                         \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            row    = center_index;                                                                                     \
            column = center_index;                                                                                     \
            value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */       \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            /* Fill matrix row of (i-1,j) */                                                                           \
            row    = left_index;                                                                                       \
            column = center_index;                                                                                     \
            value  = -coeff1 * arr; /* Right */                                                                        \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            row    = left_index;                                                                                       \
            column = left_index;                                                                                       \
            value  = coeff1 * arr; /* Center: (Right) */                                                               \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            /* Fill matrix row of (i,j-1) */                                                                           \
            row    = bottom_index;                                                                                     \
            column = bottom_index;                                                                                     \
            value  = coeff3 * att; /* Center: (Top) */                                                                 \
            UPDATE_MATRIX_ELEMENT(bottom_matrix, row, column, value);                                                  \
                                                                                                                       \
            /* Fill matrix row of (i,j+1) */                                                                           \
            row    = top_index;                                                                                        \
            column = top_index;                                                                                        \
            value  = coeff4 * att; /* Center: (Bottom) */                                                              \
            UPDATE_MATRIX_ELEMENT(top_matrix, row, column, value);                                                     \
        }                                                                                                              \
        /* ------------------------------------------ */                                                               \
        /* Radial Section: Node on the outer boundary */                                                               \
        /* ------------------------------------------ */                                                               \
        else if (i_r == grid.nr() - 1) {                                                                               \
            double h1     = grid.radialSpacing(i_r - 1);                                                               \
            double k1     = grid.angularSpacing(i_theta - 1);                                                          \
            double k2     = grid.angularSpacing(i_theta);                                                              \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
                                                                                                                       \
            /* -----------|| */                                                                                        \
            /* O   O   O  || */                                                                                        \
            /* -----------|| */                                                                                        \
            /* O   O   Õ  || <- center_matrix*/                                                                        \
            /* -----------|| */                                                                                        \
            /* O   O   O  || */                                                                                        \
            /* -----------|| */                                                                                        \
            auto& center_matrix = radial_tridiagonal_solver[i_theta];                                                  \
                                                                                                                       \
            const int center_index = i_r - numberSmootherCircles;                                                      \
            const int left_index   = i_r - numberSmootherCircles - 1;                                                  \
                                                                                                                       \
            /* Fill matrix row of (i,j) */                                                                             \
            row    = center_index;                                                                                     \
            column = center_index;                                                                                     \
            value  = 1.0;                                                                                              \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
                                                                                                                       \
            /* Fill matrix row of (i-1,j) */                                                                           \
            row    = left_index;                                                                                       \
            column = left_index;                                                                                       \
            value  = coeff1 * arr; /* Center: (Right) */                                                               \
            UPDATE_MATRIX_ELEMENT(center_matrix, row, column, value);                                                  \
        }                                                                                                              \
    } while (0)

void SmootherGive::buildAscCircleSection(const int i_r)
{
    const double r = grid_.radius(i_r);
    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        const int global_index = grid_.index(i_r, i_theta);
        const double theta     = grid_.theta(i_theta);

        double sin_theta, cos_theta;
        double coeff_beta, arr, att, art, detDF;
        level_cache_.obtainValues(i_r, i_theta, global_index, r, theta, sin_theta, cos_theta, coeff_beta, arr, att, art,
                                  detDF);

        // Build Asc at the current node
        NODE_BUILD_SMOOTHER_GIVE(i_r, i_theta, grid_, DirBC_Interior_, inner_boundary_circle_matrix_,
                                 circle_tridiagonal_solver_, radial_tridiagonal_solver_);
    }
}

void SmootherGive::buildAscRadialSection(const int i_theta)
{
    const double theta = grid_.theta(i_theta);
    for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
        const int global_index = grid_.index(i_r, i_theta);
        const double r         = grid_.radius(i_r);

        double sin_theta, cos_theta;
        double coeff_beta, arr, att, art, detDF;
        level_cache_.obtainValues(i_r, i_theta, global_index, r, theta, sin_theta, cos_theta, coeff_beta, arr, att, art,
                                  detDF);

        // Build Asc at the current node
        NODE_BUILD_SMOOTHER_GIVE(i_r, i_theta, grid_, DirBC_Interior_, inner_boundary_circle_matrix_,
                                 circle_tridiagonal_solver_, radial_tridiagonal_solver_);
    }
}

// clang-format off
void SmootherGive::buildAscMatrices()
{
    /* -------------------------------------- */
    /* Part 1: Allocate Asc Smoother matrices */
    /* -------------------------------------- */

    const int number_smoother_circles = grid_.numberSmootherCircles();
    const int length_smoother_radial  = grid_.lengthSmootherRadial();

    const int num_circle_nodes = grid_.ntheta();
    circle_tridiagonal_solver_.resize(number_smoother_circles);

    const int num_radial_nodes = length_smoother_radial;
    radial_tridiagonal_solver_.resize(grid_.ntheta());

    // Remark: circle_tridiagonal_solver_[0] is unitialized.
    // Please use inner_boundary_circle_matrix_ instead!
    #pragma omp parallel num_threads(num_omp_threads_) if (grid_.numberOfNodes() > 10'000)
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
                for (int i = 0; i < nnz; i++) {
                    inner_boundary_circle_matrix_.value(i) = 0.0;
                }
                #else
                std::function<int(int)> nnz_per_row = [&](int i_theta) {
                    return DirBC_Interior_? 1 : 4;
                };
                inner_boundary_circle_matrix_ = SparseMatrixCSR<double>(num_circle_nodes, num_circle_nodes, nnz_per_row);
                for (int i = 0; i < inner_boundary_circle_matrix_.non_zero_size(); i++) {
                    inner_boundary_circle_matrix_.values_data()[i] = 0.0;
                }
                #endif
            }

            /* Interior Circle Section */
            else {
                auto& solverMatrix = circle_tridiagonal_solver_[circle_Asc_index];

                solverMatrix = SymmetricTridiagonalSolver<double>(num_circle_nodes);
                solverMatrix.is_cyclic(true);
                solverMatrix.cyclic_corner_element() = 0.0;

                for (int i = 0; i < num_circle_nodes; i++) {
                    solverMatrix.main_diagonal(i) = 0.0;
                    if (i < num_circle_nodes - 1) {
                        solverMatrix.sub_diagonal(i) = 0.0;
                    }
                }
            }
        }

        // -------------- //
        // Radial Section //
        #pragma omp for nowait
        for (int radial_Asc_index = 0; radial_Asc_index < grid_.ntheta(); radial_Asc_index++) {
            auto& solverMatrix = radial_tridiagonal_solver_[radial_Asc_index];

            solverMatrix = SymmetricTridiagonalSolver<double>(num_radial_nodes);
            solverMatrix.is_cyclic(false);

            for (int i = 0; i < num_radial_nodes; i++) {
                solverMatrix.main_diagonal(i) = 0.0;
                if (i < num_radial_nodes - 1) {
                    solverMatrix.sub_diagonal(i) = 0.0;
                }
            }
        }
    }

    /* ---------------------------------- */
    /* Part 2: Fill Asc Smoother matrices */
    /* ---------------------------------- */

    if (num_omp_threads_ == 1) {
        /* Single-threaded execution */
        for (int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
            buildAscCircleSection(i_r);
        }
        for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
            buildAscRadialSection(i_theta);
        }
    }
    else {
        /*  Multi-threaded execution: For Loops */
        const int num_circle_tasks        = grid_.numberSmootherCircles();
        const int additional_radial_tasks = grid_.ntheta() % 3;
        const int num_radial_tasks        = grid_.ntheta() - additional_radial_tasks;

        #pragma omp parallel num_threads(num_omp_threads_)
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

    #ifdef GMGPOLAR_USE_MUMPS
    /* ------------------------------------------------------------------ */
    /* Part 3: Convert inner_boundary_circle_matrix to a symmetric matrix */
    /* ------------------------------------------------------------------ */

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
// clang-format on
