#include "../../../include/DirectSolver/DirectSolverTake/directSolverTake.h"


#define NODE_BUILD_SOLVER_MATRIX_TAKE(i_r, i_theta, \
    grid, DirBC_Interior, solver_matrix, \
    arr, att, art, detDF, coeff_beta) \
do { \
    /* -------------------- */ \
    /* Node in the interior */ \
    /* -------------------- */ \
    if (i_r > 1 && i_r < grid.nr() - 2) { \
        const int i_theta_M1 = grid.wrapThetaIndex(i_theta-1); \
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta+1); \
        \
        const double h1 = grid.radialSpacing(i_r-1); \
        const double h2 = grid.radialSpacing(i_r); \
        const double k1 = grid.angularSpacing(i_theta_M1); \
        const double k2 = grid.angularSpacing(i_theta); \
        const double coeff1 = 0.5*(k1+k2)/h1; \
        const double coeff2 = 0.5*(k1+k2)/h2; \
        const double coeff3 = 0.5*(h1+h2)/k1; \
        const double coeff4 = 0.5*(h1+h2)/k2; \
        \
        const int center_nz_index = getSolverMatrixIndex(i_r, i_theta); \
        \
        int nz_index; /* Current non_zero index in solver_matrix */ \
        \
        const int center_index = grid.index(i_r, i_theta); \
        const int left_index = grid.index(i_r-1, i_theta); \
        const int right_index = grid.index(i_r+1, i_theta); \
        const int bottom_index = grid.index(i_r, i_theta_M1); \
        const int top_index = grid.index(i_r, i_theta_P1); \
        const int bottom_left_index = grid.index(i_r-1, i_theta_M1); \
        const int bottom_right_index = grid.index(i_r+1, i_theta_M1); \
        const int top_left_index = grid.index(i_r-1, i_theta_P1); \
        const int top_right_index = grid.index(i_r+1, i_theta_P1); \
        \
        const double left_value = - coeff1 * (arr[center_index] + arr[left_index]); /* Left */ \
        const double right_value = - coeff2 * (arr[center_index] + arr[right_index]); /* Right */ \
        const double bottom_value = - coeff3 * (att[center_index] + att[bottom_index]); /* Bottom */ \
        const double top_value = - coeff4 * (att[center_index] + att[top_index]); /* Top */ \
        \
        const double center_value = ( \
            + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center_index]) /* beta_{i,j} */ \
            - left_value /* Center: (Left) */ \
            - right_value /* Center: (Right) */ \
            - bottom_value /* Center: (Bottom) */ \
            - top_value /* Center: (Top) */ \
        ); \
        \
        const double bottom_left_value = - 0.25 * (art[left_index] + art[bottom_index]); /* Bottom Left */ \
        const double bottom_right_value = + 0.25 * (art[right_index] + art[bottom_index]); /* Bottom Right */ \
        const double top_left_value = + 0.25 * (art[left_index] + art[top_index]); /* Top Left */ \
        const double top_right_value = - 0.25 * (art[right_index] + art[top_index]); /* Top Right */ \
        \
        /* Fill matrix row of (i,j) */ \
        const Stencil& CenterStencil = getStencil(i_r); \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        solver_matrix.value(nz_index) = center_value; \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Left]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = left_index + 1; \
        solver_matrix.value(nz_index) = left_value; \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Right]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = right_index + 1; \
        solver_matrix.value(nz_index) = right_value; \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Bottom]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = bottom_index + 1; \
        solver_matrix.value(nz_index) = bottom_value; \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Top]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = top_index + 1; \
        solver_matrix.value(nz_index) = top_value; \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::BottomLeft]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = bottom_left_index + 1; \
        solver_matrix.value(nz_index) = bottom_left_value; \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::BottomRight]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = bottom_right_index + 1; \
        solver_matrix.value(nz_index) = bottom_right_value; \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::TopLeft]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = top_left_index + 1; \
        solver_matrix.value(nz_index) = top_left_value; \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::TopRight]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = top_right_index + 1; \
        solver_matrix.value(nz_index) = top_right_value; \
        \
    /* -------------------------- */ \
    /* Node on the inner boundary */ \
    /* -------------------------- */ \
    } else if (i_r == 0) { \
        /* ------------------------------------------------ */ \
        /* Case 1: Dirichlet boundary on the inner boundary */ \
        /* ------------------------------------------------ */ \
        if(DirBC_Interior){ \
            const int center_nz_index = getSolverMatrixIndex(i_r, i_theta); \
            int nz_index; /* Current non_zero index in solver_matrix */ \
            const int center_index = grid.index(i_r, i_theta); \
            \
            /* Fill matrix row of (i,j) */ \
            const Stencil& CenterStencil = getStencil(i_r); \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
            solver_matrix.row_index(nz_index) = center_index + 1; \
            solver_matrix.col_index(nz_index) = center_index + 1; \
            solver_matrix.value(nz_index) = 1.0; \
            \
        } else{ \
            /* ------------------------------------------------------------- */ \
            /* Case 2: Across origin discretization on the interior boundary */ \
            /* ------------------------------------------------------------- */ \
            /* h1 gets replaced with 2 * R0. */ \
            /* (i_r-1,i_theta) gets replaced with (i_r, i_theta + grid.ntheta()/2). */ \
            /* Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil. */ \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta-1); \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta+1); \
            \
            assert(grid_.ntheta() % 2 == 0); \
            const int i_theta_AcrossOrigin = grid.wrapThetaIndex(i_theta + grid.ntheta()/2); \
            \
            double h1 = 2.0 * grid.radius(0); \
            double h2 = grid.radialSpacing(i_r); \
            double k1 = grid.angularSpacing(i_theta_M1); \
            double k2 = grid.angularSpacing(i_theta); \
            \
            const double coeff1 = 0.5*(k1+k2)/h1; \
            const double coeff2 = 0.5*(k1+k2)/h2; \
            const double coeff3 = 0.5*(h1+h2)/k1; \
            const double coeff4 = 0.5*(h1+h2)/k2; \
            \
            const int center_nz_index = getSolverMatrixIndex(i_r, i_theta); \
            \
            int nz_index; /* Current non_zero index in solver_matrix */ \
            \
            const int center_index = grid.index(i_r, i_theta); \
            const int left_index = grid.index(i_r, i_theta_AcrossOrigin); \
            const int right_index = grid.index(i_r+1, i_theta); \
            const int bottom_index = grid.index(i_r, i_theta_M1); \
            const int top_index = grid.index(i_r, i_theta_P1); \
            const int bottom_right_index = grid.index(i_r+1, i_theta_M1); \
            const int top_right_index = grid.index(i_r+1, i_theta_P1); \
            \
            const double left_value = - coeff1 * (arr[center_index] + arr[left_index]); /* Left */ \
            const double right_value = - coeff2 * (arr[center_index] + arr[right_index]); /* Right */ \
            const double bottom_value = - coeff3 * (att[center_index] + att[bottom_index]); /* Bottom */ \
            const double top_value = - coeff4 * (att[center_index] + att[top_index]); /* Top */ \
            \
            const double center_value = ( \
                + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center_index]) /* beta_{i,j} */ \
                - left_value /* Center: (Left) */ \
                - right_value /* Center: (Right) */ \
                - bottom_value /* Center: (Bottom) */ \
                - top_value /* Center: (Top) */ \
            ); \
            \
            const double bottom_right_value = + 0.25 * (art[right_index] + art[bottom_index]); /* Bottom Right */ \
            const double top_right_value = - 0.25 * (art[right_index] + art[top_index]); /* Top Right */ \
            \
            /* Fill matrix row of (i,j) */ \
            const Stencil& CenterStencil = getStencil(i_r); \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
            solver_matrix.row_index(nz_index) = center_index + 1; \
            solver_matrix.col_index(nz_index) = center_index + 1; \
            solver_matrix.value(nz_index) = center_value; \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Left]; \
            solver_matrix.row_index(nz_index) = center_index + 1; \
            solver_matrix.col_index(nz_index) = left_index + 1; \
            solver_matrix.value(nz_index) = left_value; \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Right]; \
            solver_matrix.row_index(nz_index) = center_index + 1; \
            solver_matrix.col_index(nz_index) = right_index + 1; \
            solver_matrix.value(nz_index) = right_value; \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Bottom]; \
            solver_matrix.row_index(nz_index) = center_index + 1; \
            solver_matrix.col_index(nz_index) = bottom_index + 1; \
            solver_matrix.value(nz_index) = bottom_value; \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Top]; \
            solver_matrix.row_index(nz_index) = center_index + 1; \
            solver_matrix.col_index(nz_index) = top_index + 1; \
            solver_matrix.value(nz_index) = top_value; \
            \
            /* REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = center_nz_index + CenterStencil[StencilType::BottomLeft]; */ \
            /* solver_matrix.row_index(nz_index) = center_index + 1; */ \
            /* solver_matrix.col_index(nz_index) = bottom_left_index + 1; */ \
            /* solver_matrix.value(nz_index) = bottom_left_value; */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::BottomRight]; \
            solver_matrix.row_index(nz_index) = center_index + 1; \
            solver_matrix.col_index(nz_index) = bottom_right_index + 1; \
            solver_matrix.value(nz_index) = bottom_right_value; \
            \
            /* REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = center_nz_index + CenterStencil[StencilType::TopLeft]; */ \
            /* solver_matrix.row_index(nz_index) = center_index + 1; */ \
            /* solver_matrix.col_index(nz_index) = top_left_index + 1; */ \
            /* solver_matrix.value(nz_index) = top_left_value; */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::TopRight]; \
            solver_matrix.row_index(nz_index) = center_index + 1; \
            solver_matrix.col_index(nz_index) = top_right_index + 1; \
            solver_matrix.value(nz_index) = top_right_value; \
            \
        } \
    /* ------------------------------- */ \
    /* Node next to the inner boundary */ \
    /* ------------------------------- */ \
    } else if (i_r == 1) { \
        const int i_theta_M1 = grid.wrapThetaIndex(i_theta-1); \
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta+1); \
        \
        const double h1 = grid.radialSpacing(i_r-1); \
        const double h2 = grid.radialSpacing(i_r); \
        const double k1 = grid.angularSpacing(i_theta_M1); \
        const double k2 = grid.angularSpacing(i_theta); \
        const double coeff1 = 0.5*(k1+k2)/h1; \
        const double coeff2 = 0.5*(k1+k2)/h2; \
        const double coeff3 = 0.5*(h1+h2)/k1; \
        const double coeff4 = 0.5*(h1+h2)/k2; \
        \
        const int center_nz_index = getSolverMatrixIndex(i_r, i_theta); \
        \
        int nz_index; /* Current non_zero index in solver_matrix */ \
        \
        const int center_index = grid.index(i_r, i_theta); \
        const int left_index = grid.index(i_r-1, i_theta); \
        const int right_index = grid.index(i_r+1, i_theta); \
        const int bottom_index = grid.index(i_r, i_theta_M1); \
        const int top_index = grid.index(i_r, i_theta_P1); \
        const int bottom_left_index = grid.index(i_r-1, i_theta_M1); \
        const int bottom_right_index = grid.index(i_r+1, i_theta_M1); \
        const int top_left_index = grid.index(i_r-1, i_theta_P1); \
        const int top_right_index = grid.index(i_r+1, i_theta_P1); \
        \
        const double left_value = - coeff1 * (arr[center_index] + arr[left_index]); /* Left */ \
        const double right_value = - coeff2 * (arr[center_index] + arr[right_index]); /* Right */ \
        const double bottom_value = - coeff3 * (att[center_index] + att[bottom_index]); /* Bottom */ \
        const double top_value = - coeff4 * (att[center_index] + att[top_index]); /* Top */ \
        \
        const double center_value = ( \
            + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center_index]) /* beta_{i,j} */ \
            - left_value /* Center: (Left) */ \
            - right_value /* Center: (Right) */ \
            - bottom_value /* Center: (Bottom) */ \
            - top_value /* Center: (Top) */ \
        ); \
        \
        const double bottom_left_value = - 0.25 * (art[left_index] + art[bottom_index]); /* Bottom Left */ \
        const double bottom_right_value = + 0.25 * (art[right_index] + art[bottom_index]); /* Bottom Right */ \
        const double top_left_value = + 0.25 * (art[left_index] + art[top_index]); /* Top Left */ \
        const double top_right_value = - 0.25 * (art[right_index] + art[top_index]); /* Top Right */ \
        \
        /* Fill matrix row of (i,j) */ \
        const Stencil& CenterStencil = getStencil(i_r); \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        solver_matrix.value(nz_index) = center_value; \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        if(!DirBC_Interior) { \
            nz_index = center_nz_index + CenterStencil[StencilType::Left]; \
            solver_matrix.row_index(nz_index) = center_index + 1; \
            solver_matrix.col_index(nz_index) = left_index + 1; \
            solver_matrix.value(nz_index) = left_value; \
        } \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Right]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = right_index + 1; \
        solver_matrix.value(nz_index) = right_value; \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Bottom]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = bottom_index + 1; \
        solver_matrix.value(nz_index) = bottom_value; \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Top]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = top_index + 1; \
        solver_matrix.value(nz_index) = top_value; \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        if(!DirBC_Interior) { \
            nz_index = center_nz_index + CenterStencil[StencilType::BottomLeft]; \
            solver_matrix.row_index(nz_index) = center_index + 1; \
            solver_matrix.col_index(nz_index) = bottom_left_index + 1; \
            solver_matrix.value(nz_index) = bottom_left_value; \
        } \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::BottomRight]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = bottom_right_index + 1; \
        solver_matrix.value(nz_index) = bottom_right_value; \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        if(!DirBC_Interior) { \
            nz_index = center_nz_index + CenterStencil[StencilType::TopLeft]; \
            solver_matrix.row_index(nz_index) = center_index + 1; \
            solver_matrix.col_index(nz_index) = top_left_index + 1; \
            solver_matrix.value(nz_index) = top_left_value; \
        } \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::TopRight]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = top_right_index + 1; \
        solver_matrix.value(nz_index) = top_right_value; \
        \
    /* ------------------------------- */ \
    /* Node next to the outer boundary */ \
    /* ------------------------------- */ \
    } else if (i_r == grid.nr() - 2) { \
        const int i_theta_M1 = grid.wrapThetaIndex(i_theta-1); \
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta+1); \
        \
        const double h1 = grid.radialSpacing(i_r-1); \
        const double h2 = grid.radialSpacing(i_r); \
        const double k1 = grid.angularSpacing(i_theta_M1); \
        const double k2 = grid.angularSpacing(i_theta); \
        \
        const double coeff1 = 0.5*(k1+k2)/h1; \
        const double coeff2 = 0.5*(k1+k2)/h2; \
        const double coeff3 = 0.5*(h1+h2)/k1; \
        const double coeff4 = 0.5*(h1+h2)/k2; \
        \
        const int center_nz_index = getSolverMatrixIndex(i_r, i_theta); \
        \
        int nz_index; /* Current non_zero index in solver_matrix */ \
        \
        const int center_index = grid.index(i_r, i_theta); \
        const int left_index = grid.index(i_r-1, i_theta); \
        const int right_index = grid.index(i_r+1, i_theta); \
        const int bottom_index = grid.index(i_r, i_theta_M1); \
        const int top_index = grid.index(i_r, i_theta_P1); \
        const int bottom_left_index = grid.index(i_r-1, i_theta_M1); \
        const int bottom_right_index = grid.index(i_r+1, i_theta_M1); \
        const int top_left_index = grid.index(i_r-1, i_theta_P1); \
        const int top_right_index = grid.index(i_r+1, i_theta_P1); \
        \
        const double left_value = - coeff1 * (arr[center_index] + arr[left_index]); /* Left */ \
        const double right_value = - coeff2 * (arr[center_index] + arr[right_index]); /* Right */ \
        const double bottom_value = - coeff3 * (att[center_index] + att[bottom_index]); /* Bottom */ \
        const double top_value = - coeff4 * (att[center_index] + att[top_index]); /* Top */ \
        \
        const double center_value = ( \
            + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center_index]) /* beta_{i,j} */ \
            - left_value /* Center: (Left) */ \
            - right_value /* Center: (Right) */ \
            - bottom_value /* Center: (Bottom) */ \
            - top_value /* Center: (Top) */ \
        ); \
        \
        const double bottom_left_value = - 0.25 * (art[left_index] + art[bottom_index]); /* Bottom Left */ \
        const double bottom_right_value = + 0.25 * (art[right_index] + art[bottom_index]); /* Bottom Right */ \
        const double top_left_value = + 0.25 * (art[left_index] + art[top_index]); /* Top Left */ \
        const double top_right_value = - 0.25 * (art[right_index] + art[top_index]); /* Top Right */ \
        \
        /* Fill matrix row of (i,j) */ \
        const Stencil& CenterStencil = getStencil(i_r); \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        solver_matrix.value(nz_index) = center_value; \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Left]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = left_index + 1; \
        solver_matrix.value(nz_index) = left_value; \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        /* nz_index = center_nz_index + CenterStencil[StencilType::Right]; */ \
        /* solver_matrix.row_index(nz_index) = center_index + 1; */ \
        /* solver_matrix.col_index(nz_index) = right_index + 1; */ \
        /* solver_matrix.value(nz_index) = right_value; */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Bottom]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = bottom_index + 1; \
        solver_matrix.value(nz_index) = bottom_value; \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Top]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = top_index + 1; \
        solver_matrix.value(nz_index) = top_value; \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::BottomLeft]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = bottom_left_index + 1; \
        solver_matrix.value(nz_index) = bottom_left_value; \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        /* nz_index = center_nz_index + CenterStencil[StencilType::BottomRight]; */ \
        /* solver_matrix.row_index(nz_index) = center_index + 1; */ \
        /* solver_matrix.col_index(nz_index) = bottom_right_index + 1; */ \
        /* solver_matrix.value(nz_index) = bottom_right_value; */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::TopLeft]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = top_left_index + 1; \
        solver_matrix.value(nz_index) = top_left_value; \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        /* nz_index = center_nz_index + CenterStencil[StencilType::TopRight]; */ \
        /* solver_matrix.row_index(nz_index) = center_index + 1; */ \
        /* solver_matrix.col_index(nz_index) = top_right_index + 1; */ \
        /* solver_matrix.value(nz_index) = top_right_value; */ \
        \
    /* ------------------------------------ */ \
    /* Node on the outer dirichlet boundary */ \
    /* ------------------------------------ */ \
    } else if (i_r == grid.nr() - 1) { \
        int nz_index; /* Current non_zero index in solver_matrix */ \
        const int center_nz_index = getSolverMatrixIndex(i_r, i_theta); \
        const int center_index = grid.index(i_r, i_theta); \
        \
        /* Fill matrix row of (i,j) */ \
        const Stencil& CenterStencil = getStencil(i_r); \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        solver_matrix.value(nz_index) = 1.0; \
        \
    } \
} while(0)


void DirectSolverTake::buildSolverMatrixCircleSection(const int i_r, SparseMatrix<double>& solver_matrix){
    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    const auto& arr = level_cache_.arr();
    const auto& att = level_cache_.att();
    const auto& art = level_cache_.art();
    const auto& detDF = level_cache_.detDF();
    const auto& coeff_beta = level_cache_.coeff_beta();

    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++){
        // Build solver matrix at the current node
        NODE_BUILD_SOLVER_MATRIX_TAKE(i_r, i_theta,
            grid_, DirBC_Interior_, solver_matrix,
            arr, att, art, detDF, coeff_beta);
    }
}

void DirectSolverTake::buildSolverMatrixRadialSection(const int i_theta, SparseMatrix<double>& solver_matrix){
    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    const auto& arr = level_cache_.arr();
    const auto& att = level_cache_.att();
    const auto& art = level_cache_.art();
    const auto& detDF = level_cache_.detDF();
    const auto& coeff_beta = level_cache_.coeff_beta();

    for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++){
        // Build solver matrix at the current node
        NODE_BUILD_SOLVER_MATRIX_TAKE(i_r, i_theta,
            grid_, DirBC_Interior_, solver_matrix,
            arr, att, art, detDF, coeff_beta);
    }
}


/* ------------------------------------------------------------------------ */
/* If the indexing is not smoother-based, please adjust the access patterns */
SparseMatrix<double> DirectSolverTake::buildSolverMatrix()
{
    omp_set_num_threads(num_omp_threads_);

    const int n = grid_.numberOfNodes();
    const int nnz = getNonZeroCountSolverMatrix();

    // Although the matrix is symmetric, we need to store all its entries, so we disable the symmetry.
    SparseMatrix<double> solver_matrix(n, n, nnz);
    solver_matrix.is_symmetric(false);

    if(omp_get_max_threads() == 1) {
        /* Single-threaded execution */
        for(int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
            buildSolverMatrixCircleSection(i_r, solver_matrix);
        }
        for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
            buildSolverMatrixRadialSection(i_theta, solver_matrix);
        }
    }
    else{
        /* Multi-threaded execution */
        #pragma omp parallel
        {
            /* Circle Section */
            #pragma omp for nowait
            for(int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
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

    if(!construct_symmetric){
        return solver_matrix;
    }

    /* Only store the upper tridiagonal entries of the symmetric solver_matrix */
    const int symmetric_nnz = nnz - (nnz - n) / 2;
    SparseMatrix<double> symmetric_solver_matrix(n, n, symmetric_nnz);
    symmetric_solver_matrix.is_symmetric(true);

    int current_nz = 0;
    for (int nz_index = 0; nz_index < nnz; nz_index++) {
        const int row = solver_matrix.row_index(nz_index);
        const int col = solver_matrix.col_index(nz_index);
        if (row <= col) {
            symmetric_solver_matrix.row_index(current_nz) = row;
            symmetric_solver_matrix.col_index(current_nz) = col;
            symmetric_solver_matrix.value(current_nz) = std::move(solver_matrix.value(nz_index));
            current_nz++;
        }
    }

    return symmetric_solver_matrix;
}