#include "../../../include/DirectSolver/DirectSolverGive/directSolverGive.h"

#include "../../../include/common/geometry_helper.h"

#define NODE_BUILD_SOLVER_MATRIX_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, grid, DirBC_Interior,                      \
                                      solver_matrix, arr, att, art, detDF, coeff_beta)                                         \
    do {                                                                                                                       \
        /* -------------------- */                                                                                             \
        /* Node in the interior */                                                                                             \
        /* -------------------- */                                                                                             \
        if (i_r > 1 && i_r < grid.nr() - 2) {                                                                                  \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                           \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                           \
                                                                                                                               \
            const double h1     = grid.radialSpacing(i_r - 1);                                                                 \
            const double h2     = grid.radialSpacing(i_r);                                                                     \
            const double k1     = grid.angularSpacing(i_theta_M1);                                                             \
            const double k2     = grid.angularSpacing(i_theta);                                                                \
            const double coeff1 = 0.5 * (k1 + k2) / h1;                                                                        \
            const double coeff2 = 0.5 * (k1 + k2) / h2;                                                                        \
            const double coeff3 = 0.5 * (h1 + h2) / k1;                                                                        \
            const double coeff4 = 0.5 * (h1 + h2) / k2;                                                                        \
                                                                                                                               \
            const int center_nz_index = getSolverMatrixIndex(i_r, i_theta);                                                    \
            const int left_nz_index   = getSolverMatrixIndex(i_r - 1, i_theta);                                                \
            const int right_nz_index  = getSolverMatrixIndex(i_r + 1, i_theta);                                                \
            const int bottom_nz_index = getSolverMatrixIndex(i_r, i_theta_M1);                                                 \
            const int top_nz_index    = getSolverMatrixIndex(i_r, i_theta_P1);                                                 \
                                                                                                                               \
            int nz_index; /* Current non_zero index in solver_matrix */                                                        \
                                                                                                                               \
            const int center_index = grid.index(i_r, i_theta);                                                                 \
            const int left_index   = grid.index(i_r - 1, i_theta);                                                             \
            const int right_index  = grid.index(i_r + 1, i_theta);                                                             \
            const int bottom_index = grid.index(i_r, i_theta_M1);                                                              \
            const int top_index    = grid.index(i_r, i_theta_P1);                                                              \
                                                                                                                               \
            /* Fill matrix row of (i,j) */                                                                                     \
            const Stencil& CenterStencil = getStencil(i_r);                                                                    \
                                                                                                                               \
            nz_index                          = center_nz_index + CenterStencil[StencilPosition::Center];                          \
            solver_matrix.row_index(nz_index) = center_index;                                                              \
            solver_matrix.col_index(nz_index) = center_index;                                                              \
            solver_matrix.value(nz_index) += 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */         \
                                                                                                                               \
            nz_index                          = center_nz_index + CenterStencil[StencilPosition::Left];                            \
            solver_matrix.row_index(nz_index) = center_index;                                                              \
            solver_matrix.col_index(nz_index) = left_index;                                                                \
            solver_matrix.value(nz_index) += -coeff1 * arr; /* Left */                                                         \
                                                                                                                               \
            nz_index                          = center_nz_index + CenterStencil[StencilPosition::Right];                           \
            solver_matrix.row_index(nz_index) = center_index;                                                              \
            solver_matrix.col_index(nz_index) = right_index;                                                               \
            solver_matrix.value(nz_index) += -coeff2 * arr; /* Right */                                                        \
                                                                                                                               \
            nz_index                          = center_nz_index + CenterStencil[StencilPosition::Bottom];                          \
            solver_matrix.row_index(nz_index) = center_index;                                                              \
            solver_matrix.col_index(nz_index) = bottom_index;                                                              \
            solver_matrix.value(nz_index) += -coeff3 * att; /* Bottom */                                                       \
                                                                                                                               \
            nz_index                          = center_nz_index + CenterStencil[StencilPosition::Top];                             \
            solver_matrix.row_index(nz_index) = center_index;                                                              \
            solver_matrix.col_index(nz_index) = top_index;                                                                 \
            solver_matrix.value(nz_index) += -coeff4 * att; /* Top */                                                          \
                                                                                                                               \
            nz_index                          = center_nz_index + CenterStencil[StencilPosition::Center];                          \
            solver_matrix.row_index(nz_index) = center_index;                                                              \
            solver_matrix.col_index(nz_index) = center_index;                                                              \
            /* Center: (Left, Right, Bottom, Top) */                                                                           \
            solver_matrix.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att;                                \
                                                                                                                               \
            /* Fill matrix row of (i-1,j) */                                                                                   \
            const Stencil& LeftStencil = getStencil(i_r - 1);                                                                  \
                                                                                                                               \
            nz_index                          = left_nz_index + LeftStencil[StencilPosition::Right];                               \
            solver_matrix.row_index(nz_index) = left_index;                                                                \
            solver_matrix.col_index(nz_index) = center_index;                                                              \
            solver_matrix.value(nz_index) += -coeff1 * arr; /* Right */                                                        \
                                                                                                                               \
            nz_index                          = left_nz_index + LeftStencil[StencilPosition::Center];                              \
            solver_matrix.row_index(nz_index) = left_index;                                                                \
            solver_matrix.col_index(nz_index) = left_index;                                                                \
            solver_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) */                                               \
                                                                                                                               \
            nz_index                          = left_nz_index + LeftStencil[StencilPosition::TopRight];                            \
            solver_matrix.row_index(nz_index) = left_index;                                                                \
            solver_matrix.col_index(nz_index) = top_index;                                                                 \
            solver_matrix.value(nz_index) += -0.25 * art; /* Top Right */                                                      \
                                                                                                                               \
            nz_index                          = left_nz_index + LeftStencil[StencilPosition::BottomRight];                         \
            solver_matrix.row_index(nz_index) = left_index;                                                                \
            solver_matrix.col_index(nz_index) = bottom_index;                                                              \
            solver_matrix.value(nz_index) += 0.25 * art; /* Bottom Right */                                                    \
                                                                                                                               \
            /* Fill matrix row of (i+1,j) */                                                                                   \
            const Stencil& RightStencil = getStencil(i_r + 1);                                                                 \
                                                                                                                               \
            nz_index                          = right_nz_index + RightStencil[StencilPosition::Left];                              \
            solver_matrix.row_index(nz_index) = right_index;                                                               \
            solver_matrix.col_index(nz_index) = center_index;                                                              \
            solver_matrix.value(nz_index) += -coeff2 * arr; /* Left */                                                         \
                                                                                                                               \
            nz_index                          = right_nz_index + RightStencil[StencilPosition::Center];                            \
            solver_matrix.row_index(nz_index) = right_index;                                                               \
            solver_matrix.col_index(nz_index) = right_index;                                                               \
            solver_matrix.value(nz_index) += coeff2 * arr; /* Center: (Left) */                                                \
                                                                                                                               \
            nz_index                          = right_nz_index + RightStencil[StencilPosition::TopLeft];                           \
            solver_matrix.row_index(nz_index) = right_index;                                                               \
            solver_matrix.col_index(nz_index) = top_index;                                                                 \
            solver_matrix.value(nz_index) += 0.25 * art; /* Top Left */                                                        \
                                                                                                                               \
            nz_index                          = right_nz_index + RightStencil[StencilPosition::BottomLeft];                        \
            solver_matrix.row_index(nz_index) = right_index;                                                               \
            solver_matrix.col_index(nz_index) = bottom_index;                                                              \
            solver_matrix.value(nz_index) += -0.25 * art; /* Bottom Left */                                                    \
                                                                                                                               \
            /* Fill matrix row of (i,j-1) */                                                                                   \
            const Stencil& BottomStencil = CenterStencil;                                                                      \
                                                                                                                               \
            nz_index                          = bottom_nz_index + BottomStencil[StencilPosition::Top];                             \
            solver_matrix.row_index(nz_index) = bottom_index;                                                              \
            solver_matrix.col_index(nz_index) = center_index;                                                              \
            solver_matrix.value(nz_index) += -coeff3 * att; /* Top */                                                          \
                                                                                                                               \
            nz_index                          = bottom_nz_index + BottomStencil[StencilPosition::Center];                          \
            solver_matrix.row_index(nz_index) = bottom_index;                                                              \
            solver_matrix.col_index(nz_index) = bottom_index;                                                              \
            solver_matrix.value(nz_index) += coeff3 * att; /* Center: (Top) */                                                 \
                                                                                                                               \
            nz_index                          = bottom_nz_index + BottomStencil[StencilPosition::TopRight];                        \
            solver_matrix.row_index(nz_index) = bottom_index;                                                              \
            solver_matrix.col_index(nz_index) = right_index;                                                               \
            solver_matrix.value(nz_index) += -0.25 * art; /* Top Right */                                                      \
                                                                                                                               \
            nz_index                          = bottom_nz_index + BottomStencil[StencilPosition::TopLeft];                         \
            solver_matrix.row_index(nz_index) = bottom_index;                                                              \
            solver_matrix.col_index(nz_index) = left_index;                                                                \
            solver_matrix.value(nz_index) += 0.25 * art; /* Top Left */                                                        \
                                                                                                                               \
            /* Fill matrix row of (i,j+1) */                                                                                   \
            const Stencil& TopStencil = CenterStencil;                                                                         \
                                                                                                                               \
            nz_index                          = top_nz_index + TopStencil[StencilPosition::Bottom];                                \
            solver_matrix.row_index(nz_index) = top_index;                                                                 \
            solver_matrix.col_index(nz_index) = center_index;                                                              \
            solver_matrix.value(nz_index) += -coeff4 * att; /* Bottom */                                                       \
                                                                                                                               \
            nz_index                          = top_nz_index + TopStencil[StencilPosition::Center];                                \
            solver_matrix.row_index(nz_index) = top_index;                                                                 \
            solver_matrix.col_index(nz_index) = top_index;                                                                 \
            solver_matrix.value(nz_index) += coeff4 * att; /* Center: (Bottom) */                                              \
                                                                                                                               \
            nz_index                          = top_nz_index + TopStencil[StencilPosition::BottomRight];                           \
            solver_matrix.row_index(nz_index) = top_index;                                                                 \
            solver_matrix.col_index(nz_index) = right_index;                                                               \
            solver_matrix.value(nz_index) += 0.25 * art; /* Bottom Right */                                                    \
                                                                                                                               \
            nz_index                          = top_nz_index + TopStencil[StencilPosition::BottomLeft];                            \
            solver_matrix.row_index(nz_index) = top_index;                                                                 \
            solver_matrix.col_index(nz_index) = left_index;                                                                \
            solver_matrix.value(nz_index) += -0.25 * art; /* Bottom Left */                                                    \
                                                                                                                               \
            /* -------------------------- */                                                                                   \
            /* Node on the inner boundary */                                                                                   \
            /* -------------------------- */                                                                                   \
        }                                                                                                                      \
        else if (i_r == 0) {                                                                                                   \
            /* ------------------------------------------------ */                                                             \
            /* Case 1: Dirichlet boundary on the inner boundary */                                                             \
            /* ------------------------------------------------ */                                                             \
            if (DirBC_Interior) {                                                                                              \
                const double h2     = grid.radialSpacing(i_r);                                                                 \
                const double k1     = grid.angularSpacing(i_theta - 1);                                                        \
                const double k2     = grid.angularSpacing(i_theta);                                                            \
                const double coeff2 = 0.5 * (k1 + k2) / h2;                                                                    \
                                                                                                                               \
                const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                       \
                const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                       \
                                                                                                                               \
                const int center_nz_index = getSolverMatrixIndex(i_r, i_theta);                                                \
                const int right_nz_index  = getSolverMatrixIndex(i_r + 1, i_theta);                                            \
                                                                                                                               \
                int nz_index; /* Current non_zero index in solver_matrix */                                                    \
                                                                                                                               \
                const int center_index = grid.index(i_r, i_theta);                                                             \
                const int right_index  = grid.index(i_r + 1, i_theta);                                                         \
                const int bottom_index = grid.index(i_r, i_theta_M1);                                                          \
                const int top_index    = grid.index(i_r, i_theta_P1);                                                          \
                                                                                                                               \
                /* Fill matrix row of (i,j) */                                                                                 \
                const Stencil& CenterStencil = getStencil(i_r);                                                                \
                                                                                                                               \
                nz_index                          = center_nz_index + CenterStencil[StencilPosition::Center];                      \
                solver_matrix.row_index(nz_index) = center_index;                                                          \
                solver_matrix.col_index(nz_index) = center_index;                                                          \
                solver_matrix.value(nz_index) += 1.0;                                                                          \
                                                                                                                               \
                /* Fill matrix row of (i+1,j) */                                                                               \
                const Stencil& RightStencil = getStencil(i_r + 1);                                                             \
                                                                                                                               \
                /* REMOVED: Moved to the right hand side to make the matrix symmetric */                                       \
                /* nz_index = right_nz_index + RightStencil[StencilPosition::Left]; */                                             \
                /* solver_matrix.row_index(nz_index) = right_index; */                                                     \
                /* solver_matrix.col_index(nz_index) = center_index; */                                                    \
                /* solver_matrix.value(nz_index) += - coeff2 * arr; // Left */                                                 \
                                                                                                                               \
                nz_index                          = right_nz_index + RightStencil[StencilPosition::Center];                        \
                solver_matrix.row_index(nz_index) = right_index;                                                           \
                solver_matrix.col_index(nz_index) = right_index;                                                           \
                solver_matrix.value(nz_index) += coeff2 * arr; /* Center: (Left) */                                            \
                                                                                                                               \
                /* REMOVED: Moved to the right hand side to make the matrix symmetric */                                       \
                /* nz_index = right_nz_index + RightStencil[StencilPosition::TopLeft]; */                                          \
                /* solver_matrix.row_index(nz_index) = right_index; */                                                     \
                /* solver_matrix.col_index(nz_index) = top_index; */                                                       \
                /* solver_matrix.value(nz_index) += 0.25 * art; // Top Left */                                                 \
                                                                                                                               \
                /* REMOVED: Moved to the right hand side to make the matrix symmetric */                                       \
                /* nz_index = right_nz_index + RightStencil[StencilPosition::BottomLeft]; */                                       \
                /* solver_matrix.row_index(nz_index) = right_index; */                                                     \
                /* solver_matrix.col_index(nz_index) = bottom_index; */                                                    \
                /* solver_matrix.value(nz_index) += - 0.25 * art; // Bottom Left */                                            \
            }                                                                                                                  \
            else {                                                                                                             \
                /* ------------------------------------------------------------- */                                            \
                /* Case 2: Across origin discretization on the interior boundary */                                            \
                /* ------------------------------------------------------------- */                                            \
                /* h1 gets replaced with 2 * R0. */                                                                            \
                /* (i_r-1,i_theta) gets replaced with (i_r, i_theta + grid.ntheta()/2). */                                     \
                /* Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil. */            \
                const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                       \
                const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                       \
                                                                                                                               \
                assert(grid_.ntheta() % 2 == 0);                                                                               \
                const int i_theta_AcrossOrigin = grid.wrapThetaIndex(i_theta + grid.ntheta() / 2);                             \
                                                                                                                               \
                double h1     = 2.0 * grid.radius(0);                                                                          \
                double h2     = grid.radialSpacing(i_r);                                                                       \
                double k1     = grid.angularSpacing(i_theta_M1);                                                               \
                double k2     = grid.angularSpacing(i_theta);                                                                  \
                double coeff1 = 0.5 * (k1 + k2) / h1;                                                                          \
                double coeff2 = 0.5 * (k1 + k2) / h2;                                                                          \
                double coeff3 = 0.5 * (h1 + h2) / k1;                                                                          \
                double coeff4 = 0.5 * (h1 + h2) / k2;                                                                          \
                                                                                                                               \
                const int center_nz_index = getSolverMatrixIndex(i_r, i_theta);                                                \
                const int left_nz_index   = getSolverMatrixIndex(i_r, i_theta_AcrossOrigin);                                   \
                const int right_nz_index  = getSolverMatrixIndex(i_r + 1, i_theta);                                            \
                const int bottom_nz_index = getSolverMatrixIndex(i_r, i_theta_M1);                                             \
                const int top_nz_index    = getSolverMatrixIndex(i_r, i_theta_P1);                                             \
                                                                                                                               \
                int nz_index; /* Current non_zero index in solver_matrix */                                                    \
                                                                                                                               \
                const int center_index = grid.index(i_r, i_theta);                                                             \
                const int left_index   = grid.index(i_r, i_theta_AcrossOrigin);                                                \
                const int right_index  = grid.index(i_r + 1, i_theta);                                                         \
                const int bottom_index = grid.index(i_r, i_theta_M1);                                                          \
                const int top_index    = grid.index(i_r, i_theta_P1);                                                          \
                                                                                                                               \
                /* Fill matrix row of (i,j) */                                                                                 \
                const Stencil& CenterStencil = getStencil(i_r);                                                                \
                                                                                                                               \
                nz_index                          = center_nz_index + CenterStencil[StencilPosition::Center];                      \
                solver_matrix.row_index(nz_index) = center_index;                                                          \
                solver_matrix.col_index(nz_index) = center_index;                                                          \
                solver_matrix.value(nz_index) +=                                                                               \
                    0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */                                  \
                                                                                                                               \
                nz_index                          = center_nz_index + CenterStencil[StencilPosition::Left];                        \
                solver_matrix.row_index(nz_index) = center_index;                                                          \
                solver_matrix.col_index(nz_index) = left_index;                                                            \
                solver_matrix.value(nz_index) += -coeff1 * arr; /* Left */                                                     \
                                                                                                                               \
                nz_index                          = center_nz_index + CenterStencil[StencilPosition::Right];                       \
                solver_matrix.row_index(nz_index) = center_index;                                                          \
                solver_matrix.col_index(nz_index) = right_index;                                                           \
                solver_matrix.value(nz_index) += -coeff2 * arr; /* Right */                                                    \
                                                                                                                               \
                nz_index                          = center_nz_index + CenterStencil[StencilPosition::Bottom];                      \
                solver_matrix.row_index(nz_index) = center_index;                                                          \
                solver_matrix.col_index(nz_index) = bottom_index;                                                          \
                solver_matrix.value(nz_index) += -coeff3 * att; /* Bottom */                                                   \
                                                                                                                               \
                nz_index                          = center_nz_index + CenterStencil[StencilPosition::Top];                         \
                solver_matrix.row_index(nz_index) = center_index;                                                          \
                solver_matrix.col_index(nz_index) = top_index;                                                             \
                solver_matrix.value(nz_index) += -coeff4 * att; /* Top */                                                      \
                                                                                                                               \
                nz_index                          = center_nz_index + CenterStencil[StencilPosition::Center];                      \
                solver_matrix.row_index(nz_index) = center_index;                                                          \
                solver_matrix.col_index(nz_index) = center_index;                                                          \
                /* Center: (Left, Right, Bottom, Top) */                                                                       \
                solver_matrix.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att;                            \
                                                                                                                               \
                /* Fill matrix row of (i-1,j) */                                                                               \
                /* From view the view of the across origin node, the directions are roatated by 180 degrees in the stencil! */ \
                const Stencil& LeftStencil = CenterStencil;                                                                    \
                                                                                                                               \
                nz_index                          = left_nz_index + LeftStencil[StencilPosition::Left];                            \
                solver_matrix.row_index(nz_index) = left_index;                                                            \
                solver_matrix.col_index(nz_index) = center_index;                                                          \
                solver_matrix.value(nz_index) += -coeff1 * arr; /* Right -> Left*/                                             \
                                                                                                                               \
                nz_index                          = left_nz_index + LeftStencil[StencilPosition::Center];                          \
                solver_matrix.row_index(nz_index) = left_index;                                                            \
                solver_matrix.col_index(nz_index) = left_index;                                                            \
                solver_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) -> Center: (Left) */                         \
                                                                                                                               \
                /* Top Right -> Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */                                       \
                /* nz_index = left_nz_index + LeftStencil[StencilPosition::BottomLeft]; */                                         \
                /* solver_matrix.row_index(nz_index) = left_index; */                                                          \
                /* solver_matrix.col_index(nz_index) = top_index; */                                                           \
                /* solver_matrix.value(nz_index) += - 0.25 * art; // Top Right -> Bottom Left*/                                \
                                                                                                                               \
                /* Bottom Right -> Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */                                       \
                /* nz_index = left_nz_index + LeftStencil[StencilPosition::TopLeft]; */                                            \
                /* solver_matrix.row_index(nz_index) = left_index; */                                                          \
                /* solver_matrix.col_index(nz_index) = bottom_index; */                                                        \
                /* solver_matrix.value(nz_index) += 0.25 * art; // Bottom Right -> Top Left */                                 \
                                                                                                                               \
                /* Fill matrix row of (i+1,j) */                                                                               \
                const Stencil& RightStencil = getStencil(i_r + 1);                                                             \
                                                                                                                               \
                nz_index                          = right_nz_index + RightStencil[StencilPosition::Left];                          \
                solver_matrix.row_index(nz_index) = right_index;                                                           \
                solver_matrix.col_index(nz_index) = center_index;                                                          \
                solver_matrix.value(nz_index) += -coeff2 * arr; /* Left */                                                     \
                                                                                                                               \
                nz_index                          = right_nz_index + RightStencil[StencilPosition::Center];                        \
                solver_matrix.row_index(nz_index) = right_index;                                                           \
                solver_matrix.col_index(nz_index) = right_index;                                                           \
                solver_matrix.value(nz_index) += coeff2 * arr; /* Center: (Left) */                                            \
                                                                                                                               \
                nz_index                          = right_nz_index + RightStencil[StencilPosition::TopLeft];                       \
                solver_matrix.row_index(nz_index) = right_index;                                                           \
                solver_matrix.col_index(nz_index) = top_index;                                                             \
                solver_matrix.value(nz_index) += 0.25 * art; /* Top Left */                                                    \
                                                                                                                               \
                nz_index                          = right_nz_index + RightStencil[StencilPosition::BottomLeft];                    \
                solver_matrix.row_index(nz_index) = right_index;                                                           \
                solver_matrix.col_index(nz_index) = bottom_index;                                                          \
                solver_matrix.value(nz_index) += -0.25 * art; /* Bottom Left */                                                \
                                                                                                                               \
                /* Fill matrix row of (i,j-1) */                                                                               \
                const Stencil& BottomStencil = CenterStencil;                                                                  \
                                                                                                                               \
                nz_index                          = bottom_nz_index + BottomStencil[StencilPosition::Top];                         \
                solver_matrix.row_index(nz_index) = bottom_index;                                                          \
                solver_matrix.col_index(nz_index) = center_index;                                                          \
                solver_matrix.value(nz_index) += -coeff3 * att; /* Top */                                                      \
                                                                                                                               \
                nz_index                          = bottom_nz_index + BottomStencil[StencilPosition::Center];                      \
                solver_matrix.row_index(nz_index) = bottom_index;                                                          \
                solver_matrix.col_index(nz_index) = bottom_index;                                                          \
                solver_matrix.value(nz_index) += coeff3 * att; /* Center: (Top) */                                             \
                                                                                                                               \
                nz_index                          = bottom_nz_index + BottomStencil[StencilPosition::TopRight];                    \
                solver_matrix.row_index(nz_index) = bottom_index;                                                          \
                solver_matrix.col_index(nz_index) = right_index;                                                           \
                solver_matrix.value(nz_index) += -0.25 * art; /* Top Right */                                                  \
                                                                                                                               \
                /* REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */                                                                 \
                /* nz_index = bottom_nz_index + BottomStencil[StencilPosition::TopLeft]; */                                        \
                /* solver_matrix.row_index(nz_index) = bottom_index; */                                                    \
                /* solver_matrix.col_index(nz_index) = left_index; */                                                      \
                /* solver_matrix.value(nz_index) += 0.25 * art; // Top Left */                                                 \
                                                                                                                               \
                /* Fill matrix row of (i,j+1) */                                                                               \
                const Stencil& TopStencil = CenterStencil;                                                                     \
                                                                                                                               \
                nz_index                          = top_nz_index + TopStencil[StencilPosition::Bottom];                            \
                solver_matrix.row_index(nz_index) = top_index;                                                             \
                solver_matrix.col_index(nz_index) = center_index;                                                          \
                solver_matrix.value(nz_index) += -coeff4 * att; /* Bottom */                                                   \
                                                                                                                               \
                nz_index                          = top_nz_index + TopStencil[StencilPosition::Center];                            \
                solver_matrix.row_index(nz_index) = top_index;                                                             \
                solver_matrix.col_index(nz_index) = top_index;                                                             \
                solver_matrix.value(nz_index) += coeff4 * att; /* Center: (Bottom) */                                          \
                                                                                                                               \
                nz_index                          = top_nz_index + TopStencil[StencilPosition::BottomRight];                       \
                solver_matrix.row_index(nz_index) = top_index;                                                             \
                solver_matrix.col_index(nz_index) = right_index;                                                           \
                solver_matrix.value(nz_index) += 0.25 * art; /* Bottom Right */                                                \
                                                                                                                               \
                /* REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */                                                                 \
                /* nz_index = top_nz_index + TopStencil[StencilPosition::BottomLeft]; */                                           \
                /* solver_matrix.row_index(nz_index) = top_index; */                                                       \
                /* solver_matrix.col_index(nz_index) = left_index; */                                                      \
                /* solver_matrix.value(nz_index) += - 0.25 * art; // Bottom Left */                                            \
            }                                                                                                                  \
            /* ------------------------------- */                                                                              \
            /* Node next to the inner boundary */                                                                              \
            /* ------------------------------- */                                                                              \
        }                                                                                                                      \
        else if (i_r == 1) {                                                                                                   \
            const double h1     = grid.radialSpacing(i_r - 1);                                                                 \
            const double h2     = grid.radialSpacing(i_r);                                                                     \
            const double k1     = grid.angularSpacing(i_theta - 1);                                                            \
            const double k2     = grid.angularSpacing(i_theta);                                                                \
            const double coeff1 = 0.5 * (k1 + k2) / h1;                                                                        \
            const double coeff2 = 0.5 * (k1 + k2) / h2;                                                                        \
            const double coeff3 = 0.5 * (h1 + h2) / k1;                                                                        \
            const double coeff4 = 0.5 * (h1 + h2) / k2;                                                                        \
                                                                                                                               \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                           \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                           \
                                                                                                                               \
            const int center_nz_index = getSolverMatrixIndex(i_r, i_theta);                                                    \
            const int left_nz_index   = getSolverMatrixIndex(i_r - 1, i_theta);                                                \
            const int right_nz_index  = getSolverMatrixIndex(i_r + 1, i_theta);                                                \
            const int bottom_nz_index = getSolverMatrixIndex(i_r, i_theta_M1);                                                 \
            const int top_nz_index    = getSolverMatrixIndex(i_r, i_theta_P1);                                                 \
                                                                                                                               \
            int nz_index; /* Current non_zero index in solver_matrix */                                                        \
                                                                                                                               \
            const int center_index = grid.index(i_r, i_theta);                                                                 \
            const int left_index   = grid.index(i_r - 1, i_theta);                                                             \
            const int right_index  = grid.index(i_r + 1, i_theta);                                                             \
            const int bottom_index = grid.index(i_r, i_theta_M1);                                                              \
            const int top_index    = grid.index(i_r, i_theta_P1);                                                              \
                                                                                                                               \
            /* Fill matrix row of (i,j) */                                                                                     \
            const Stencil& CenterStencil = getStencil(i_r);                                                                    \
                                                                                                                               \
            nz_index                          = center_nz_index + CenterStencil[StencilPosition::Center];                          \
            solver_matrix.row_index(nz_index) = center_index;                                                              \
            solver_matrix.col_index(nz_index) = center_index;                                                              \
            solver_matrix.value(nz_index) += 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */         \
                                                                                                                               \
            /* REMOVED: Moved to the right hand side to make the matrix symmetric */                                           \
            if (!DirBC_Interior) {                                                                                             \
                nz_index                          = center_nz_index + CenterStencil[StencilPosition::Left];                        \
                solver_matrix.row_index(nz_index) = center_index;                                                          \
                solver_matrix.col_index(nz_index) = left_index;                                                            \
                solver_matrix.value(nz_index) += -coeff1 * arr; /* Left */                                                     \
            }                                                                                                                  \
                                                                                                                               \
            nz_index                          = center_nz_index + CenterStencil[StencilPosition::Right];                           \
            solver_matrix.row_index(nz_index) = center_index;                                                              \
            solver_matrix.col_index(nz_index) = right_index;                                                               \
            solver_matrix.value(nz_index) += -coeff2 * arr; /* Right */                                                        \
                                                                                                                               \
            nz_index                          = center_nz_index + CenterStencil[StencilPosition::Bottom];                          \
            solver_matrix.row_index(nz_index) = center_index;                                                              \
            solver_matrix.col_index(nz_index) = bottom_index;                                                              \
            solver_matrix.value(nz_index) += -coeff3 * att; /* Bottom */                                                       \
                                                                                                                               \
            nz_index                          = center_nz_index + CenterStencil[StencilPosition::Top];                             \
            solver_matrix.row_index(nz_index) = center_index;                                                              \
            solver_matrix.col_index(nz_index) = top_index;                                                                 \
            solver_matrix.value(nz_index) += -coeff4 * att; /* Top */                                                          \
                                                                                                                               \
            nz_index                          = center_nz_index + CenterStencil[StencilPosition::Center];                          \
            solver_matrix.row_index(nz_index) = center_index;                                                              \
            solver_matrix.col_index(nz_index) = center_index;                                                              \
            /* Center: (Left, Right, Bottom, Top) */                                                                           \
            solver_matrix.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att;                                \
                                                                                                                               \
            if (!DirBC_Interior) { /* Don't give to the inner dirichlet boundary! */                                           \
                /* Fill matrix row of (i-1,j) */                                                                               \
                const Stencil& LeftStencil = getStencil(i_r - 1);                                                              \
                                                                                                                               \
                nz_index                          = left_nz_index + LeftStencil[StencilPosition::Right];                           \
                solver_matrix.row_index(nz_index) = left_index;                                                            \
                solver_matrix.col_index(nz_index) = center_index;                                                          \
                solver_matrix.value(nz_index) += -coeff1 * arr; /* Right */                                                    \
                                                                                                                               \
                nz_index                          = left_nz_index + LeftStencil[StencilPosition::Center];                          \
                solver_matrix.row_index(nz_index) = left_index;                                                            \
                solver_matrix.col_index(nz_index) = left_index;                                                            \
                solver_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) */                                           \
                                                                                                                               \
                nz_index                          = left_nz_index + LeftStencil[StencilPosition::TopRight];                        \
                solver_matrix.row_index(nz_index) = left_index;                                                            \
                solver_matrix.col_index(nz_index) = top_index;                                                             \
                solver_matrix.value(nz_index) += -0.25 * art; /* Top Right */                                                  \
                                                                                                                               \
                nz_index                          = left_nz_index + LeftStencil[StencilPosition::BottomRight];                     \
                solver_matrix.row_index(nz_index) = left_index;                                                            \
                solver_matrix.col_index(nz_index) = bottom_index;                                                          \
                solver_matrix.value(nz_index) += 0.25 * art; /* Bottom Right */                                                \
            }                                                                                                                  \
            /* Fill matrix row of (i+1,j) */                                                                                   \
            const Stencil& RightStencil = getStencil(i_r + 1);                                                                 \
                                                                                                                               \
            nz_index                          = right_nz_index + RightStencil[StencilPosition::Left];                              \
            solver_matrix.row_index(nz_index) = right_index;                                                               \
            solver_matrix.col_index(nz_index) = center_index;                                                              \
            solver_matrix.value(nz_index) += -coeff2 * arr; /* Left */                                                         \
                                                                                                                               \
            nz_index                          = right_nz_index + RightStencil[StencilPosition::Center];                            \
            solver_matrix.row_index(nz_index) = right_index;                                                               \
            solver_matrix.col_index(nz_index) = right_index;                                                               \
            solver_matrix.value(nz_index) += coeff2 * arr; /* Center: (Left) */                                                \
                                                                                                                               \
            nz_index                          = right_nz_index + RightStencil[StencilPosition::TopLeft];                           \
            solver_matrix.row_index(nz_index) = right_index;                                                               \
            solver_matrix.col_index(nz_index) = top_index;                                                                 \
            solver_matrix.value(nz_index) += 0.25 * art; /* Top Left */                                                        \
                                                                                                                               \
            nz_index                          = right_nz_index + RightStencil[StencilPosition::BottomLeft];                        \
            solver_matrix.row_index(nz_index) = right_index;                                                               \
            solver_matrix.col_index(nz_index) = bottom_index;                                                              \
            solver_matrix.value(nz_index) += -0.25 * art; /* Bottom Left */                                                    \
                                                                                                                               \
            /* Fill matrix row of (i,j-1) */                                                                                   \
            const Stencil& BottomStencil = CenterStencil;                                                                      \
                                                                                                                               \
            nz_index                          = bottom_nz_index + BottomStencil[StencilPosition::Top];                             \
            solver_matrix.row_index(nz_index) = bottom_index;                                                              \
            solver_matrix.col_index(nz_index) = center_index;                                                              \
            solver_matrix.value(nz_index) += -coeff3 * att; /* Top */                                                          \
                                                                                                                               \
            nz_index                          = bottom_nz_index + BottomStencil[StencilPosition::Center];                          \
            solver_matrix.row_index(nz_index) = bottom_index;                                                              \
            solver_matrix.col_index(nz_index) = bottom_index;                                                              \
            solver_matrix.value(nz_index) += coeff3 * att; /* Center: (Top) */                                                 \
                                                                                                                               \
            nz_index                          = bottom_nz_index + BottomStencil[StencilPosition::TopRight];                        \
            solver_matrix.row_index(nz_index) = bottom_index;                                                              \
            solver_matrix.col_index(nz_index) = right_index;                                                               \
            solver_matrix.value(nz_index) += -0.25 * art; /* Top Right */                                                      \
                                                                                                                               \
            /* REMOVED: Moved to the right hand side to make the matrix symmetric */                                           \
            if (!DirBC_Interior) {                                                                                             \
                nz_index                          = bottom_nz_index + BottomStencil[StencilPosition::TopLeft];                     \
                solver_matrix.row_index(nz_index) = bottom_index;                                                          \
                solver_matrix.col_index(nz_index) = left_index;                                                            \
                solver_matrix.value(nz_index) += 0.25 * art; /* Top Left */                                                    \
            }                                                                                                                  \
                                                                                                                               \
            /* Fill matrix row of (i,j+1) */                                                                                   \
            const Stencil& TopStencil = CenterStencil;                                                                         \
                                                                                                                               \
            nz_index                          = top_nz_index + TopStencil[StencilPosition::Bottom];                                \
            solver_matrix.row_index(nz_index) = top_index;                                                                 \
            solver_matrix.col_index(nz_index) = center_index;                                                              \
            solver_matrix.value(nz_index) += -coeff4 * att; /* Bottom */                                                       \
                                                                                                                               \
            nz_index                          = top_nz_index + TopStencil[StencilPosition::Center];                                \
            solver_matrix.row_index(nz_index) = top_index;                                                                 \
            solver_matrix.col_index(nz_index) = top_index;                                                                 \
            solver_matrix.value(nz_index) += coeff4 * att; /* Center: (Bottom) */                                              \
                                                                                                                               \
            nz_index                          = top_nz_index + TopStencil[StencilPosition::BottomRight];                           \
            solver_matrix.row_index(nz_index) = top_index;                                                                 \
            solver_matrix.col_index(nz_index) = right_index;                                                               \
            solver_matrix.value(nz_index) += 0.25 * art; /* Bottom Right */                                                    \
                                                                                                                               \
            /* REMOVED: Moved to the right hand side to make the matrix symmetric */                                           \
            if (!DirBC_Interior) {                                                                                             \
                nz_index                          = top_nz_index + TopStencil[StencilPosition::BottomLeft];                        \
                solver_matrix.row_index(nz_index) = top_index;                                                             \
                solver_matrix.col_index(nz_index) = left_index;                                                            \
                solver_matrix.value(nz_index) += -0.25 * art; /* Bottom Left */                                                \
            }                                                                                                                  \
                                                                                                                               \
            /* ------------------------------- */                                                                              \
            /* Node next to the outer boundary */                                                                              \
            /* ------------------------------- */                                                                              \
        }                                                                                                                      \
        else if (i_r == grid.nr() - 2) {                                                                                       \
            const double h1     = grid.radialSpacing(i_r - 1);                                                                 \
            const double h2     = grid.radialSpacing(i_r);                                                                     \
            const double k1     = grid.angularSpacing(i_theta - 1);                                                            \
            const double k2     = grid.angularSpacing(i_theta);                                                                \
            const double coeff1 = 0.5 * (k1 + k2) / h1;                                                                        \
            const double coeff2 = 0.5 * (k1 + k2) / h2;                                                                        \
            const double coeff3 = 0.5 * (h1 + h2) / k1;                                                                        \
            const double coeff4 = 0.5 * (h1 + h2) / k2;                                                                        \
                                                                                                                               \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                           \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                           \
                                                                                                                               \
            const int center_nz_index = getSolverMatrixIndex(i_r, i_theta);                                                    \
            const int left_nz_index   = getSolverMatrixIndex(i_r - 1, i_theta);                                                \
            const int right_nz_index  = getSolverMatrixIndex(i_r + 1, i_theta);                                                \
            const int bottom_nz_index = getSolverMatrixIndex(i_r, i_theta_M1);                                                 \
            const int top_nz_index    = getSolverMatrixIndex(i_r, i_theta_P1);                                                 \
                                                                                                                               \
            int nz_index; /* Current non_zero index in solver_matrix */                                                        \
                                                                                                                               \
            const int center_index = grid.index(i_r, i_theta);                                                                 \
            const int left_index   = grid.index(i_r - 1, i_theta);                                                             \
            const int right_index  = grid.index(i_r + 1, i_theta);                                                             \
            const int bottom_index = grid.index(i_r, i_theta_M1);                                                              \
            const int top_index    = grid.index(i_r, i_theta_P1);                                                              \
                                                                                                                               \
            /* Fill matrix row of (i,j) */                                                                                     \
            const Stencil& CenterStencil = getStencil(i_r);                                                                    \
                                                                                                                               \
            nz_index                          = center_nz_index + CenterStencil[StencilPosition::Center];                          \
            solver_matrix.row_index(nz_index) = center_index;                                                              \
            solver_matrix.col_index(nz_index) = center_index;                                                              \
            solver_matrix.value(nz_index) += 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */         \
                                                                                                                               \
            nz_index                          = center_nz_index + CenterStencil[StencilPosition::Left];                            \
            solver_matrix.row_index(nz_index) = center_index;                                                              \
            solver_matrix.col_index(nz_index) = left_index;                                                                \
            solver_matrix.value(nz_index) += -coeff1 * arr; /* Left */                                                         \
                                                                                                                               \
            /* REMOVED: Moved to the right hand side to make the matrix symmetric */                                           \
            /* nz_index = center_nz_index + CenterStencil[StencilPosition::Right]; */                                              \
            /* solver_matrix.row_index(nz_index) = center_index; */                                                        \
            /* solver_matrix.col_index(nz_index) = right_index; */                                                         \
            /* solver_matrix.value(nz_index) += - coeff2 * arr; // Right */                                                    \
                                                                                                                               \
            nz_index                          = center_nz_index + CenterStencil[StencilPosition::Bottom];                          \
            solver_matrix.row_index(nz_index) = center_index;                                                              \
            solver_matrix.col_index(nz_index) = bottom_index;                                                              \
            solver_matrix.value(nz_index) += -coeff3 * att; /* Bottom */                                                       \
                                                                                                                               \
            nz_index                          = center_nz_index + CenterStencil[StencilPosition::Top];                             \
            solver_matrix.row_index(nz_index) = center_index;                                                              \
            solver_matrix.col_index(nz_index) = top_index;                                                                 \
            solver_matrix.value(nz_index) += -coeff4 * att; /* Top */                                                          \
                                                                                                                               \
            nz_index                          = center_nz_index + CenterStencil[StencilPosition::Center];                          \
            solver_matrix.row_index(nz_index) = center_index;                                                              \
            solver_matrix.col_index(nz_index) = center_index;                                                              \
            /* Center: (Left, Right, Bottom, Top) */                                                                           \
            solver_matrix.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att;                                \
                                                                                                                               \
            /* Fill matrix row of (i-1,j) */                                                                                   \
            const Stencil& LeftStencil = getStencil(i_r - 1);                                                                  \
                                                                                                                               \
            nz_index                          = left_nz_index + LeftStencil[StencilPosition::Right];                               \
            solver_matrix.row_index(nz_index) = left_index;                                                                \
            solver_matrix.col_index(nz_index) = center_index;                                                              \
            solver_matrix.value(nz_index) += -coeff1 * arr; /* Right */                                                        \
                                                                                                                               \
            nz_index                          = left_nz_index + LeftStencil[StencilPosition::Center];                              \
            solver_matrix.row_index(nz_index) = left_index;                                                                \
            solver_matrix.col_index(nz_index) = left_index;                                                                \
            solver_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) */                                               \
                                                                                                                               \
            nz_index                          = left_nz_index + LeftStencil[StencilPosition::TopRight];                            \
            solver_matrix.row_index(nz_index) = left_index;                                                                \
            solver_matrix.col_index(nz_index) = top_index;                                                                 \
            solver_matrix.value(nz_index) += -0.25 * art; /* Top Right */                                                      \
                                                                                                                               \
            nz_index                          = left_nz_index + LeftStencil[StencilPosition::BottomRight];                         \
            solver_matrix.row_index(nz_index) = left_index;                                                                \
            solver_matrix.col_index(nz_index) = bottom_index;                                                              \
            solver_matrix.value(nz_index) += 0.25 * art; /* Bottom Right */                                                    \
                                                                                                                               \
            /* Fill matrix row of (i+1,j) */                                                                                   \
            /* Don't give to the outer dirichlet boundary! */                                                                  \
                                                                                                                               \
            /* Fill matrix row of (i,j-1) */                                                                                   \
            const Stencil& BottomStencil = CenterStencil;                                                                      \
                                                                                                                               \
            nz_index                          = bottom_nz_index + BottomStencil[StencilPosition::Top];                             \
            solver_matrix.row_index(nz_index) = bottom_index;                                                              \
            solver_matrix.col_index(nz_index) = center_index;                                                              \
            solver_matrix.value(nz_index) += -coeff3 * att; /* Top */                                                          \
                                                                                                                               \
            nz_index                          = bottom_nz_index + BottomStencil[StencilPosition::Center];                          \
            solver_matrix.row_index(nz_index) = bottom_index;                                                              \
            solver_matrix.col_index(nz_index) = bottom_index;                                                              \
            solver_matrix.value(nz_index) += coeff3 * att; /* Center: (Top) */                                                 \
                                                                                                                               \
            /* REMOVED: Moved to the right hand side to make the matrix symmetric */                                           \
            /* nz_index = bottom_nz_index + BottomStencil[StencilPosition::TopRight]; */                                           \
            /* solver_matrix.row_index(nz_index) = bottom_index; */                                                        \
            /* solver_matrix.col_index(nz_index) = right_index; */                                                         \
            /* solver_matrix.value(nz_index) += - 0.25 * art; // Top Right */                                                  \
                                                                                                                               \
            nz_index                          = bottom_nz_index + BottomStencil[StencilPosition::TopLeft];                         \
            solver_matrix.row_index(nz_index) = bottom_index;                                                              \
            solver_matrix.col_index(nz_index) = left_index;                                                                \
            solver_matrix.value(nz_index) += 0.25 * art; /* Top Left */                                                        \
                                                                                                                               \
            /* Fill matrix row of (i,j+1) */                                                                                   \
            const Stencil& TopStencil = CenterStencil;                                                                         \
                                                                                                                               \
            nz_index                          = top_nz_index + TopStencil[StencilPosition::Bottom];                                \
            solver_matrix.row_index(nz_index) = top_index;                                                                 \
            solver_matrix.col_index(nz_index) = center_index;                                                              \
            solver_matrix.value(nz_index) += -coeff4 * att; /* Bottom */                                                       \
                                                                                                                               \
            nz_index                          = top_nz_index + TopStencil[StencilPosition::Center];                                \
            solver_matrix.row_index(nz_index) = top_index;                                                                 \
            solver_matrix.col_index(nz_index) = top_index;                                                                 \
            solver_matrix.value(nz_index) += coeff4 * att; /* Center: (Bottom) */                                              \
                                                                                                                               \
            /* REMOVED: Moved to the right hand side to make the matrix symmetric */                                           \
            /* nz_index = top_nz_index + TopStencil[StencilPosition::BottomRight]; */                                              \
            /* solver_matrix.row_index(nz_index) = top_index; */                                                           \
            /* solver_matrix.col_index(nz_index) = right_index; */                                                         \
            /* solver_matrix.value(nz_index) += 0.25 * art; // Bottom Right */                                                 \
                                                                                                                               \
            nz_index                          = top_nz_index + TopStencil[StencilPosition::BottomLeft];                            \
            solver_matrix.row_index(nz_index) = top_index;                                                                 \
            solver_matrix.col_index(nz_index) = left_index;                                                                \
            solver_matrix.value(nz_index) += -0.25 * art; /* Bottom Left */                                                    \
                                                                                                                               \
            /* ------------------------------------ */                                                                         \
            /* Node on the outer dirichlet boundary */                                                                         \
            /* ------------------------------------ */                                                                         \
        }                                                                                                                      \
        else if (i_r == grid.nr() - 1) {                                                                                       \
            double h1     = grid.radialSpacing(i_r - 1);                                                                       \
            double k1     = grid.angularSpacing(i_theta - 1);                                                                  \
            double k2     = grid.angularSpacing(i_theta);                                                                      \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                              \
                                                                                                                               \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                           \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                           \
                                                                                                                               \
            int nz_index; /* Current non_zero index in solver_matrix */                                                        \
                                                                                                                               \
            const int center_nz_index = getSolverMatrixIndex(i_r, i_theta);                                                    \
            const int left_nz_index   = getSolverMatrixIndex(i_r - 1, i_theta);                                                \
                                                                                                                               \
            const int center_index = grid.index(i_r, i_theta);                                                                 \
            const int left_index   = grid.index(i_r - 1, i_theta);                                                             \
            const int bottom_index = grid.index(i_r, i_theta_M1);                                                              \
            const int top_index    = grid.index(i_r, i_theta_P1);                                                              \
                                                                                                                               \
            /* Fill matrix row of (i,j) */                                                                                     \
            const Stencil& CenterStencil = getStencil(i_r);                                                                    \
                                                                                                                               \
            nz_index                          = center_nz_index + CenterStencil[StencilPosition::Center];                          \
            solver_matrix.row_index(nz_index) = center_index;                                                              \
            solver_matrix.col_index(nz_index) = center_index;                                                              \
            solver_matrix.value(nz_index) += 1.0;                                                                              \
                                                                                                                               \
            /* Give value to the interior nodes! */                                                                            \
            /* Fill matrix row of (i-1,j) */                                                                                   \
            const Stencil& LeftStencil = getStencil(i_r - 1);                                                                  \
                                                                                                                               \
            /* REMOVED: Moved to the right hand side to make the matrix symmetric */                                           \
            /* nz_index = left_nz_index + LeftStencil[StencilPosition::Right]; */                                                  \
            /* solver_matrix.row_index(nz_index) = left_index; */                                                          \
            /* solver_matrix.col_index(nz_index) = center_index; */                                                        \
            /* solver_matrix.value(nz_index) += - coeff1 * arr; // Right */                                                    \
                                                                                                                               \
            nz_index                          = left_nz_index + LeftStencil[StencilPosition::Center];                              \
            solver_matrix.row_index(nz_index) = left_index;                                                                \
            solver_matrix.col_index(nz_index) = left_index;                                                                \
            solver_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) */                                               \
                                                                                                                               \
            /* REMOVED: Moved to the right hand side to make the matrix symmetric */                                           \
            /* nz_index = left_nz_index + LeftStencil[StencilPosition::TopRight]; */                                               \
            /* solver_matrix.row_index(nz_index) = left_index; */                                                          \
            /* solver_matrix.col_index(nz_index) = top_index; */                                                           \
            /* solver_matrix.value(nz_index) += - 0.25 * art; // Top Right */                                                  \
                                                                                                                               \
            /* REMOVED: Moved to the right hand side to make the matrix symmetric */                                           \
            /* nz_index = left_nz_index + LeftStencil[StencilPosition::BottomRight]; */                                            \
            /* solver_matrix.row_index(nz_index) = left_index; */                                                          \
            /* solver_matrix.col_index(nz_index) = bottom_index; */                                                        \
            /* solver_matrix.value(nz_index) += 0.25 * art; // Bottom Right */                                                 \
        }                                                                                                                      \
    } while (0)

void DirectSolverGive::buildSolverMatrixCircleSection(const int i_r, SparseMatrixCOO<double>& solver_matrix)
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
            compute_jacobian_elements(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art,
                                      detDF);
        }

        // Build solver matrix at the current node
        NODE_BUILD_SOLVER_MATRIX_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, grid_, DirBC_Interior_,
                                      solver_matrix, arr, att, art, detDF, coeff_beta);
    }
}

void DirectSolverGive::buildSolverMatrixRadialSection(const int i_theta, SparseMatrixCOO<double>& solver_matrix)
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
            compute_jacobian_elements(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art,
                                      detDF);
        }

        // Build solver matrix at the current node
        NODE_BUILD_SOLVER_MATRIX_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, grid_, DirBC_Interior_,
                                      solver_matrix, arr, att, art, detDF, coeff_beta);
    }
}

// clang-format off

/* ------------------------------------------------------------------------ */
/* If the indexing is not smoother-based, please adjust the access patterns */
SparseMatrixCOO<double> DirectSolverGive::buildSolverMatrix()
{
    omp_set_num_threads(num_omp_threads_);

    const int n   = grid_.numberOfNodes();
    const int nnz = getNonZeroCountSolverMatrix();

    // Although the matrix is symmetric, we need to store all its entries, so we disable the symmetry.
    SparseMatrixCOO<double> solver_matrix(n, n, nnz);
    solver_matrix.is_symmetric(false);

    #pragma omp parallel for if (nnz > 10'000)
    for (int i = 0; i < nnz; i++) {
        solver_matrix.value(i) = 0.0;
    }

    if (omp_get_max_threads() == 1) {
        /* Single-threaded execution */
        for (int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
            buildSolverMatrixCircleSection(i_r, solver_matrix);
        }
        for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
            buildSolverMatrixRadialSection(i_theta, solver_matrix);
        }
    }
    else {
        /*  Multi-threaded execution: For Loops */
        const int num_circle_tasks        = grid_.numberSmootherCircles();
        const int additional_radial_tasks = grid_.ntheta() % 3;
        const int num_radial_tasks        = grid_.ntheta() - additional_radial_tasks;

        #pragma omp parallel
        {
            #pragma omp for
            for (int circle_task = 0; circle_task < num_circle_tasks; circle_task += 3) {
                int i_r = grid_.numberSmootherCircles() - circle_task - 1;
                buildSolverMatrixCircleSection(i_r, solver_matrix);
            }
            #pragma omp for
            for (int circle_task = 1; circle_task < num_circle_tasks; circle_task += 3) {
                int i_r = grid_.numberSmootherCircles() - circle_task - 1;
                buildSolverMatrixCircleSection(i_r, solver_matrix);
            }
            #pragma omp for nowait
            for (int circle_task = 2; circle_task < num_circle_tasks; circle_task += 3) {
                int i_r = grid_.numberSmootherCircles() - circle_task - 1;
                buildSolverMatrixCircleSection(i_r, solver_matrix);
            }

            #pragma omp for
            for (int radial_task = 0; radial_task < num_radial_tasks; radial_task += 3) {
                if (radial_task > 0) {
                    int i_theta = radial_task + additional_radial_tasks;
                    buildSolverMatrixRadialSection(i_theta, solver_matrix);
                }
                else {
                    if (additional_radial_tasks == 0) {
                        buildSolverMatrixRadialSection(0, solver_matrix);
                    }
                    else if (additional_radial_tasks >= 1) {
                        buildSolverMatrixRadialSection(0, solver_matrix);
                        buildSolverMatrixRadialSection(1, solver_matrix);
                    }
                }
            }
            #pragma omp for
            for (int radial_task = 1; radial_task < num_radial_tasks; radial_task += 3) {
                if (radial_task > 1) {
                    int i_theta = radial_task + additional_radial_tasks;
                    buildSolverMatrixRadialSection(i_theta, solver_matrix);
                }
                else {
                    if (additional_radial_tasks == 0) {
                        buildSolverMatrixRadialSection(1, solver_matrix);
                    }
                    else if (additional_radial_tasks == 1) {
                        buildSolverMatrixRadialSection(2, solver_matrix);
                    }
                    else if (additional_radial_tasks == 2) {
                        buildSolverMatrixRadialSection(2, solver_matrix);
                        buildSolverMatrixRadialSection(3, solver_matrix);
                    }
                }
            }
            #pragma omp for
            for (int radial_task = 2; radial_task < num_radial_tasks; radial_task += 3) {
                int i_theta = radial_task + additional_radial_tasks;
                buildSolverMatrixRadialSection(i_theta, solver_matrix);
            }
        }
    }

    /* Mumps: In the case of symmetric matrices, only half of the matrix should be provided. */
    /* Speeds up factorization time. */
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