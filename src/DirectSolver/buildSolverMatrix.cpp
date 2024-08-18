#include "../../include/DirectSolver/directSolver.h"


#define COMPUTE_JACOBIAN_ELEMENTS(domain_geometry, r, theta, sin_theta, cos_theta, coeff_alpha, \
    arr, att, art, detDF) \
do { \
    /* Calculate the elements of the Jacobian matrix for the transformation mapping */ \
    /* The Jacobian matrix is: */ \
    /* [Jrr, Jrt] */ \
    /* [Jtr, Jtt] */ \
    const double Jrr = domain_geometry.dFx_dr(r, theta, sin_theta, cos_theta); \
    const double Jtr = domain_geometry.dFy_dr(r, theta, sin_theta, cos_theta); \
    const double Jrt = domain_geometry.dFx_dt(r, theta, sin_theta, cos_theta); \
    const double Jtt = domain_geometry.dFy_dt(r, theta, sin_theta, cos_theta); \
    /* Compute the determinant of the Jacobian matrix */ \
    detDF = Jrr * Jtt - Jrt * Jtr; \
    /* Compute the elements of the symmetric matrix: */ \
    /* 0.5 * alpha * DF^{-1} * DF^{-T} * |det(DF)| */ \
    /* which is represented by: */ \
    /* [arr, 0.5*art] */ \
    /* [0.5*atr, att] */ \
    arr = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / fabs(detDF); \
    att = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / fabs(detDF); \
    art = (- Jtt * Jtr - Jrt * Jrr) * coeff_alpha / fabs(detDF); \
    /* Note that the inverse Jacobian matrix DF^{-1} is: */ \
    /* 1.0 / det(DF) *   */ \
    /* [Jtt, -Jrt] */ \
    /* [-Jtr, Jrr] */ \
} while(0) \


#define NODE_BUILD_SOLVER_MATRIX_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, \
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
        const int left_nz_index = getSolverMatrixIndex(i_r-1, i_theta); \
        const int right_nz_index = getSolverMatrixIndex(i_r+1, i_theta); \
        const int bottom_nz_index = getSolverMatrixIndex(i_r, i_theta_M1); \
        const int top_nz_index = getSolverMatrixIndex(i_r, i_theta_P1); \
        \
        int nz_index; /* Current non_zero index in solver_matrix */ \
        \
        const int center_index = grid.index(i_r, i_theta); \
        const int left_index = grid.index(i_r-1, i_theta); \
        const int right_index = grid.index(i_r+1, i_theta); \
        const int bottom_index = grid.index(i_r, i_theta_M1); \
        const int top_index = grid.index(i_r, i_theta_P1); \
        \
        /* Fill matrix row of (i,j) */ \
        const Stencil& CenterStencil = getStencil(i_r); \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        solver_matrix.value(nz_index) += 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Left]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = left_index + 1; \
        solver_matrix.value(nz_index) += - coeff1 * arr; /* Left */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Right]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = right_index + 1; \
        solver_matrix.value(nz_index) += - coeff2 * arr; /* Right */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Bottom]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = bottom_index + 1; \
        solver_matrix.value(nz_index) += - coeff3 * att; /* Bottom */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Top]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = top_index + 1; \
        solver_matrix.value(nz_index) += - coeff4 * att; /* Top */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        /* Center: (Left, Right, Bottom, Top) */ \
        solver_matrix.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; \
        \
        /* Fill matrix row of (i-1,j) */ \
        const Stencil& LeftStencil = getStencil(i_r-1); \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::Right]; \
        solver_matrix.row_index(nz_index) = left_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        solver_matrix.value(nz_index) += - coeff1 * arr; /* Right */ \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = left_index + 1; \
        solver_matrix.col_index(nz_index) = left_index + 1; \
        solver_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) */ \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::TopRight]; \
        solver_matrix.row_index(nz_index) = left_index + 1; \
        solver_matrix.col_index(nz_index) = top_index + 1; \
        solver_matrix.value(nz_index) += - 0.25 * art; /* Top Right */ \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::BottomRight]; \
        solver_matrix.row_index(nz_index) = left_index + 1; \
        solver_matrix.col_index(nz_index) = bottom_index + 1; \
        solver_matrix.value(nz_index) += 0.25 * art; /* Bottom Right */ \
        \
        /* Fill matrix row of (i+1,j) */ \
        const Stencil& RightStencil = getStencil(i_r+1); \
        \
        nz_index = right_nz_index + RightStencil[StencilType::Left]; \
        solver_matrix.row_index(nz_index) = right_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        solver_matrix.value(nz_index) += - coeff2 * arr; /* Left */ \
        \
        nz_index = right_nz_index + RightStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = right_index + 1; \
        solver_matrix.col_index(nz_index) = right_index + 1; \
        solver_matrix.value(nz_index) += coeff2 * arr; /* Center: (Left) */ \
        \
        nz_index = right_nz_index + RightStencil[StencilType::TopLeft]; \
        solver_matrix.row_index(nz_index) = right_index + 1; \
        solver_matrix.col_index(nz_index) = top_index + 1; \
        solver_matrix.value(nz_index) += 0.25 * art; /* Top Left */ \
        \
        nz_index = right_nz_index + RightStencil[StencilType::BottomLeft]; \
        solver_matrix.row_index(nz_index) = right_index + 1; \
        solver_matrix.col_index(nz_index) = bottom_index + 1; \
        solver_matrix.value(nz_index) += - 0.25 * art; /* Bottom Left */ \
        \
        /* Fill matrix row of (i,j-1) */ \
        const Stencil& BottomStencil = CenterStencil; \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::Top]; \
        solver_matrix.row_index(nz_index) = bottom_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        solver_matrix.value(nz_index) += - coeff3 * att; /* Top */ \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = bottom_index + 1; \
        solver_matrix.col_index(nz_index) = bottom_index + 1; \
        solver_matrix.value(nz_index) += coeff3 * att; /* Center: (Top) */ \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::TopRight]; \
        solver_matrix.row_index(nz_index) = bottom_index + 1; \
        solver_matrix.col_index(nz_index) = right_index + 1; \
        solver_matrix.value(nz_index) += - 0.25 * art; /* Top Right */ \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::TopLeft]; \
        solver_matrix.row_index(nz_index) = bottom_index + 1; \
        solver_matrix.col_index(nz_index) = left_index + 1; \
        solver_matrix.value(nz_index) += 0.25 * art; /* Top Left */ \
        \
        /* Fill matrix row of (i,j+1) */ \
        const Stencil& TopStencil = CenterStencil; \
        \
        nz_index = top_nz_index + TopStencil[StencilType::Bottom]; \
        solver_matrix.row_index(nz_index) = top_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        solver_matrix.value(nz_index) += - coeff4 * att; /* Bottom */ \
        \
        nz_index = top_nz_index + TopStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = top_index + 1; \
        solver_matrix.col_index(nz_index) = top_index + 1; \
        solver_matrix.value(nz_index) += coeff4 * att; /* Center: (Bottom) */ \
        \
        nz_index = top_nz_index + TopStencil[StencilType::BottomRight]; \
        solver_matrix.row_index(nz_index) = top_index + 1; \
        solver_matrix.col_index(nz_index) = right_index + 1; \
        solver_matrix.value(nz_index) += 0.25 * art; /* Bottom Right */ \
        \
        nz_index = top_nz_index + TopStencil[StencilType::BottomLeft]; \
        solver_matrix.row_index(nz_index) = top_index + 1; \
        solver_matrix.col_index(nz_index) = left_index + 1; \
        solver_matrix.value(nz_index) += - 0.25 * art; /* Bottom Left */ \
        \
    /* -------------------------- */ \
    /* Node on the inner boundary */ \
    /* -------------------------- */ \
    } else if (i_r == 0) { \
        /* ------------------------------------------------ */ \
        /* Case 1: Dirichlet boundary on the inner boundary */ \
        /* ------------------------------------------------ */ \
        if(DirBC_Interior){ \
            const double h2 = grid.radialSpacing(i_r); \
            const double k1 = grid.angularSpacing(i_theta-1); \
            const double k2 = grid.angularSpacing(i_theta); \
            const double coeff2 = 0.5*(k1+k2)/h2; \
            \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta-1); \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta+1); \
            \
            const int center_nz_index = getSolverMatrixIndex(i_r, i_theta); \
            const int right_nz_index = getSolverMatrixIndex(i_r+1, i_theta); \
            \
            int nz_index; /* Current non_zero index in solver_matrix */ \
            \
            const int center_index = grid.index(i_r, i_theta); \
            const int right_index = grid.index(i_r+1, i_theta); \
            const int bottom_index = grid.index(i_r, i_theta_M1); \
            const int top_index = grid.index(i_r, i_theta_P1); \
            \
            /* Fill matrix row of (i,j) */ \
            const Stencil& CenterStencil = getStencil(i_r); \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
            solver_matrix.row_index(nz_index) = center_index + 1; \
            solver_matrix.col_index(nz_index) = center_index + 1; \
            solver_matrix.value(nz_index) += 1.0; \
            \
            /* Fill matrix row of (i+1,j) */ \
            const Stencil& RightStencil = getStencil(i_r+1); \
            \
            /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
            /* nz_index = right_nz_index + RightStencil[StencilType::Left]; */ \
            /* solver_matrix.row_index(nz_index) = right_index + 1; */ \
            /* solver_matrix.col_index(nz_index) = center_index + 1; */ \
            /* solver_matrix.value(nz_index) += - coeff2 * arr; // Left */ \
            \
            nz_index = right_nz_index + RightStencil[StencilType::Center]; \
            solver_matrix.row_index(nz_index) = right_index + 1; \
            solver_matrix.col_index(nz_index) = right_index + 1; \
            solver_matrix.value(nz_index) += coeff2 * arr; /* Center: (Left) */ \
            \
            /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
            /* nz_index = right_nz_index + RightStencil[StencilType::TopLeft]; */ \
            /* solver_matrix.row_index(nz_index) = right_index + 1; */ \
            /* solver_matrix.col_index(nz_index) = top_index + 1; */ \
            /* solver_matrix.value(nz_index) += 0.25 * art; // Top Left */ \
            \
            /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
            /* nz_index = right_nz_index + RightStencil[StencilType::BottomLeft]; */ \
            /* solver_matrix.row_index(nz_index) = right_index + 1; */ \
            /* solver_matrix.col_index(nz_index) = bottom_index + 1; */ \
            /* solver_matrix.value(nz_index) += - 0.25 * art; // Bottom Left */ \
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
            double coeff1 = 0.5*(k1+k2)/h1; \
            double coeff2 = 0.5*(k1+k2)/h2; \
            double coeff3 = 0.5*(h1+h2)/k1; \
            double coeff4 = 0.5*(h1+h2)/k2; \
            \
            const int center_nz_index = getSolverMatrixIndex(i_r, i_theta); \
            const int left_nz_index = getSolverMatrixIndex(i_r, i_theta_AcrossOrigin); \
            const int right_nz_index = getSolverMatrixIndex(i_r+1, i_theta); \
            const int bottom_nz_index = getSolverMatrixIndex(i_r, i_theta_M1); \
            const int top_nz_index = getSolverMatrixIndex(i_r, i_theta_P1); \
            \
            int nz_index; /* Current non_zero index in solver_matrix */ \
            \
            const int center_index = grid.index(i_r, i_theta); \
            const int left_index = grid.index(i_r, i_theta_AcrossOrigin); \
            const int right_index = grid.index(i_r+1, i_theta); \
            const int bottom_index = grid.index(i_r, i_theta_M1); \
            const int top_index = grid.index(i_r, i_theta_P1); \
            \
            /* Fill matrix row of (i,j) */ \
            const Stencil& CenterStencil = getStencil(i_r); \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
            solver_matrix.row_index(nz_index) = center_index + 1; \
            solver_matrix.col_index(nz_index) = center_index + 1; \
            solver_matrix.value(nz_index) += 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Left]; \
            solver_matrix.row_index(nz_index) = center_index + 1; \
            solver_matrix.col_index(nz_index) = left_index + 1; \
            solver_matrix.value(nz_index) += - coeff1 * arr; /* Left */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Right]; \
            solver_matrix.row_index(nz_index) = center_index + 1; \
            solver_matrix.col_index(nz_index) = right_index + 1; \
            solver_matrix.value(nz_index) += - coeff2 * arr; /* Right */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Bottom]; \
            solver_matrix.row_index(nz_index) = center_index + 1; \
            solver_matrix.col_index(nz_index) = bottom_index + 1; \
            solver_matrix.value(nz_index) += - coeff3 * att; /* Bottom */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Top]; \
            solver_matrix.row_index(nz_index) = center_index + 1; \
            solver_matrix.col_index(nz_index) = top_index + 1; \
            solver_matrix.value(nz_index) += - coeff4 * att; /* Top */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
            solver_matrix.row_index(nz_index) = center_index + 1; \
            solver_matrix.col_index(nz_index) = center_index + 1; \
            /* Center: (Left, Right, Bottom, Top) */ \
            solver_matrix.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; \
            \
            /* Fill matrix row of (i-1,j) */ \
            /* From view the view of the across origin node, the directions are roatated by 180 degrees in the stencil! */ \
            const Stencil& LeftStencil = CenterStencil; \
            \
            nz_index = left_nz_index + LeftStencil[StencilType::Left]; \
            solver_matrix.row_index(nz_index) = left_index + 1; \
            solver_matrix.col_index(nz_index) = center_index + 1; \
            solver_matrix.value(nz_index) += - coeff1 * arr; /* Right -> Left*/ \
            \
            nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
            solver_matrix.row_index(nz_index) = left_index + 1; \
            solver_matrix.col_index(nz_index) = left_index + 1; \
            solver_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) -> Center: (Left) */ \
            \
            /* Top Right -> Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = left_nz_index + LeftStencil[StencilType::BottomLeft]; */  \
            /* solver_matrix.row_index(nz_index) = left_index; */  \
            /* solver_matrix.col_index(nz_index) = top_index; */  \
            /* solver_matrix.value(nz_index) += - 0.25 * art; // Top Right -> Bottom Left*/ \
            \
            /* Bottom Right -> Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = left_nz_index + LeftStencil[StencilType::TopLeft]; */ \
            /* solver_matrix.row_index(nz_index) = left_index; */  \
            /* solver_matrix.col_index(nz_index) = bottom_index; */  \
            /* solver_matrix.value(nz_index) += 0.25 * art; // Bottom Right -> Top Left */ \
            \
            /* Fill matrix row of (i+1,j) */ \
            const Stencil& RightStencil = getStencil(i_r+1); \
            \
            nz_index = right_nz_index + RightStencil[StencilType::Left]; \
            solver_matrix.row_index(nz_index) = right_index + 1; \
            solver_matrix.col_index(nz_index) = center_index + 1; \
            solver_matrix.value(nz_index) += - coeff2 * arr; /* Left */ \
            \
            nz_index = right_nz_index + RightStencil[StencilType::Center]; \
            solver_matrix.row_index(nz_index) = right_index + 1; \
            solver_matrix.col_index(nz_index) = right_index + 1; \
            solver_matrix.value(nz_index) += coeff2 * arr; /* Center: (Left) */ \
            \
            nz_index = right_nz_index + RightStencil[StencilType::TopLeft]; \
            solver_matrix.row_index(nz_index) = right_index + 1; \
            solver_matrix.col_index(nz_index) = top_index + 1; \
            solver_matrix.value(nz_index) += 0.25 * art; /* Top Left */ \
            \
            nz_index = right_nz_index + RightStencil[StencilType::BottomLeft]; \
            solver_matrix.row_index(nz_index) = right_index + 1; \
            solver_matrix.col_index(nz_index) = bottom_index + 1; \
            solver_matrix.value(nz_index) += - 0.25 * art; /* Bottom Left */ \
            \
            /* Fill matrix row of (i,j-1) */ \
            const Stencil& BottomStencil = CenterStencil; \
            \
            nz_index = bottom_nz_index + BottomStencil[StencilType::Top]; \
            solver_matrix.row_index(nz_index) = bottom_index + 1; \
            solver_matrix.col_index(nz_index) = center_index + 1; \
            solver_matrix.value(nz_index) += - coeff3 * att; /* Top */ \
            \
            nz_index = bottom_nz_index + BottomStencil[StencilType::Center]; \
            solver_matrix.row_index(nz_index) = bottom_index + 1; \
            solver_matrix.col_index(nz_index) = bottom_index + 1; \
            solver_matrix.value(nz_index) += coeff3 * att; /* Center: (Top) */ \
            \
            nz_index = bottom_nz_index + BottomStencil[StencilType::TopRight]; \
            solver_matrix.row_index(nz_index) = bottom_index + 1; \
            solver_matrix.col_index(nz_index) = right_index + 1; \
            solver_matrix.value(nz_index) += - 0.25 * art; /* Top Right */ \
            \
            /* REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = bottom_nz_index + BottomStencil[StencilType::TopLeft]; */ \
            /* solver_matrix.row_index(nz_index) = bottom_index + 1; */ \
            /* solver_matrix.col_index(nz_index) = left_index + 1; */ \
            /* solver_matrix.value(nz_index) += 0.25 * art; // Top Left */ \
            \
            /* Fill matrix row of (i,j+1) */ \
            const Stencil& TopStencil = CenterStencil; \
            \
            nz_index = top_nz_index + TopStencil[StencilType::Bottom]; \
            solver_matrix.row_index(nz_index) = top_index + 1; \
            solver_matrix.col_index(nz_index) = center_index + 1; \
            solver_matrix.value(nz_index) += - coeff4 * att; /* Bottom */ \
            \
            nz_index = top_nz_index + TopStencil[StencilType::Center]; \
            solver_matrix.row_index(nz_index) = top_index + 1; \
            solver_matrix.col_index(nz_index) = top_index + 1; \
            solver_matrix.value(nz_index) += coeff4 * att; /* Center: (Bottom) */ \
            \
            nz_index = top_nz_index + TopStencil[StencilType::BottomRight]; \
            solver_matrix.row_index(nz_index) = top_index + 1; \
            solver_matrix.col_index(nz_index) = right_index + 1; \
            solver_matrix.value(nz_index) += 0.25 * art; /* Bottom Right */ \
            \
            /* REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = top_nz_index + TopStencil[StencilType::BottomLeft]; */ \
            /* solver_matrix.row_index(nz_index) = top_index + 1; */ \
            /* solver_matrix.col_index(nz_index) = left_index + 1; */ \
            /* solver_matrix.value(nz_index) += - 0.25 * art; // Bottom Left */ \
            \
        } \
    /* ------------------------------- */ \
    /* Node next to the inner boundary */ \
    /* ------------------------------- */ \
    } else if (i_r == 1) { \
        const double h1 = grid.radialSpacing(i_r-1); \
        const double h2 = grid.radialSpacing(i_r); \
        const double k1 = grid.angularSpacing(i_theta-1); \
        const double k2 = grid.angularSpacing(i_theta); \
        const double coeff1 = 0.5*(k1+k2)/h1; \
        const double coeff2 = 0.5*(k1+k2)/h2; \
        const double coeff3 = 0.5*(h1+h2)/k1; \
        const double coeff4 = 0.5*(h1+h2)/k2; \
        \
        const int i_theta_M1 = grid.wrapThetaIndex(i_theta-1); \
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta+1); \
        \
        const int center_nz_index = getSolverMatrixIndex(i_r, i_theta); \
        const int left_nz_index = getSolverMatrixIndex(i_r-1, i_theta); \
        const int right_nz_index = getSolverMatrixIndex(i_r+1, i_theta); \
        const int bottom_nz_index = getSolverMatrixIndex(i_r, i_theta_M1); \
        const int top_nz_index = getSolverMatrixIndex(i_r, i_theta_P1); \
        \
        int nz_index; /* Current non_zero index in solver_matrix */ \
        \
        const int center_index = grid.index(i_r, i_theta); \
        const int left_index = grid.index(i_r-1, i_theta); \
        const int right_index = grid.index(i_r+1, i_theta); \
        const int bottom_index = grid.index(i_r, i_theta_M1); \
        const int top_index = grid.index(i_r, i_theta_P1); \
        \
        /* Fill matrix row of (i,j) */ \
        const Stencil& CenterStencil = getStencil(i_r); \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        solver_matrix.value(nz_index) += 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */ \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        if(!DirBC_Interior) { \
            nz_index = center_nz_index + CenterStencil[StencilType::Left]; \
            solver_matrix.row_index(nz_index) = center_index + 1; \
            solver_matrix.col_index(nz_index) = left_index + 1; \
            solver_matrix.value(nz_index) += - coeff1 * arr; /* Left */ \
        } \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Right]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = right_index + 1; \
        solver_matrix.value(nz_index) += - coeff2 * arr; /* Right */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Bottom]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = bottom_index + 1; \
        solver_matrix.value(nz_index) += - coeff3 * att; /* Bottom */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Top]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = top_index + 1; \
        solver_matrix.value(nz_index) += - coeff4 * att; /* Top */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        /* Center: (Left, Right, Bottom, Top) */ \
        solver_matrix.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; \
        \
        if(!DirBC_Interior){ /* Don't give to the inner dirichlet boundary! */ \
            /* Fill matrix row of (i-1,j) */ \
            const Stencil& LeftStencil = getStencil(i_r-1); \
            \
            nz_index = left_nz_index + LeftStencil[StencilType::Right]; \
            solver_matrix.row_index(nz_index) = left_index + 1; \
            solver_matrix.col_index(nz_index) = center_index + 1; \
            solver_matrix.value(nz_index) += - coeff1 * arr; /* Right */ \
            \
            nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
            solver_matrix.row_index(nz_index) = left_index + 1; \
            solver_matrix.col_index(nz_index) = left_index + 1; \
            solver_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) */ \
            \
            nz_index = left_nz_index + LeftStencil[StencilType::TopRight]; \
            solver_matrix.row_index(nz_index) = left_index + 1; \
            solver_matrix.col_index(nz_index) = top_index + 1; \
            solver_matrix.value(nz_index) += - 0.25 * art; /* Top Right */ \
            \
            nz_index = left_nz_index + LeftStencil[StencilType::BottomRight]; \
            solver_matrix.row_index(nz_index) = left_index + 1; \
            solver_matrix.col_index(nz_index) = bottom_index + 1; \
            solver_matrix.value(nz_index) += 0.25 * art; /* Bottom Right */ \
            \
        } \
        /* Fill matrix row of (i+1,j) */ \
        const Stencil& RightStencil = getStencil(i_r+1); \
        \
        nz_index = right_nz_index + RightStencil[StencilType::Left]; \
        solver_matrix.row_index(nz_index) = right_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        solver_matrix.value(nz_index) += - coeff2 * arr; /* Left */ \
        \
        nz_index = right_nz_index + RightStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = right_index + 1; \
        solver_matrix.col_index(nz_index) = right_index + 1; \
        solver_matrix.value(nz_index) += coeff2 * arr; /* Center: (Left) */ \
        \
        nz_index = right_nz_index + RightStencil[StencilType::TopLeft]; \
        solver_matrix.row_index(nz_index) = right_index + 1; \
        solver_matrix.col_index(nz_index) = top_index + 1; \
        solver_matrix.value(nz_index) += 0.25 * art; /* Top Left */ \
        \
        nz_index = right_nz_index + RightStencil[StencilType::BottomLeft]; \
        solver_matrix.row_index(nz_index) = right_index + 1; \
        solver_matrix.col_index(nz_index) = bottom_index + 1; \
        solver_matrix.value(nz_index) += - 0.25 * art; /* Bottom Left */ \
        \
        /* Fill matrix row of (i,j-1) */ \
        const Stencil& BottomStencil = CenterStencil; \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::Top]; \
        solver_matrix.row_index(nz_index) = bottom_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        solver_matrix.value(nz_index) += - coeff3 * att; /* Top */ \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = bottom_index + 1; \
        solver_matrix.col_index(nz_index) = bottom_index + 1; \
        solver_matrix.value(nz_index) += coeff3 * att; /* Center: (Top) */ \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::TopRight]; \
        solver_matrix.row_index(nz_index) = bottom_index + 1; \
        solver_matrix.col_index(nz_index) = right_index + 1; \
        solver_matrix.value(nz_index) += - 0.25 * art; /* Top Right */ \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        if(!DirBC_Interior) { \
            nz_index = bottom_nz_index + BottomStencil[StencilType::TopLeft]; \
            solver_matrix.row_index(nz_index) = bottom_index + 1; \
            solver_matrix.col_index(nz_index) = left_index + 1; \
            solver_matrix.value(nz_index) += 0.25 * art; /* Top Left */ \
        } \
        \
        /* Fill matrix row of (i,j+1) */ \
        const Stencil& TopStencil = CenterStencil; \
        \
        nz_index = top_nz_index + TopStencil[StencilType::Bottom]; \
        solver_matrix.row_index(nz_index) = top_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        solver_matrix.value(nz_index) += - coeff4 * att; /* Bottom */ \
        \
        nz_index = top_nz_index + TopStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = top_index + 1; \
        solver_matrix.col_index(nz_index) = top_index + 1; \
        solver_matrix.value(nz_index) += coeff4 * att; /* Center: (Bottom) */ \
        \
        nz_index = top_nz_index + TopStencil[StencilType::BottomRight]; \
        solver_matrix.row_index(nz_index) = top_index + 1; \
        solver_matrix.col_index(nz_index) = right_index + 1; \
        solver_matrix.value(nz_index) += 0.25 * art; /* Bottom Right */ \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        if(!DirBC_Interior) { \
            nz_index = top_nz_index + TopStencil[StencilType::BottomLeft]; \
            solver_matrix.row_index(nz_index) = top_index + 1; \
            solver_matrix.col_index(nz_index) = left_index + 1; \
            solver_matrix.value(nz_index) += - 0.25 * art; /* Bottom Left */ \
        } \
        \
    /* ------------------------------- */ \
    /* Node next to the outer boundary */ \
    /* ------------------------------- */ \
    } else if (i_r == grid.nr() - 2) { \
        const double h1 = grid.radialSpacing(i_r-1); \
        const double h2 = grid.radialSpacing(i_r); \
        const double k1 = grid.angularSpacing(i_theta-1); \
        const double k2 = grid.angularSpacing(i_theta); \
        const double coeff1 = 0.5*(k1+k2)/h1; \
        const double coeff2 = 0.5*(k1+k2)/h2; \
        const double coeff3 = 0.5*(h1+h2)/k1; \
        const double coeff4 = 0.5*(h1+h2)/k2; \
        \
        const int i_theta_M1 = grid.wrapThetaIndex(i_theta-1); \
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta+1); \
        \
        const int center_nz_index = getSolverMatrixIndex(i_r, i_theta); \
        const int left_nz_index = getSolverMatrixIndex(i_r-1, i_theta); \
        const int right_nz_index = getSolverMatrixIndex(i_r+1, i_theta); \
        const int bottom_nz_index = getSolverMatrixIndex(i_r, i_theta_M1); \
        const int top_nz_index = getSolverMatrixIndex(i_r, i_theta_P1); \
        \
        int nz_index; /* Current non_zero index in solver_matrix */ \
        \
        const int center_index = grid.index(i_r, i_theta); \
        const int left_index = grid.index(i_r-1, i_theta); \
        const int right_index = grid.index(i_r+1, i_theta); \
        const int bottom_index = grid.index(i_r, i_theta_M1); \
        const int top_index = grid.index(i_r, i_theta_P1); \
        \
        /* Fill matrix row of (i,j) */ \
        const Stencil& CenterStencil = getStencil(i_r); \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        solver_matrix.value(nz_index) += 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Left]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = left_index + 1; \
        solver_matrix.value(nz_index) += - coeff1 * arr; /* Left */ \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        /* nz_index = center_nz_index + CenterStencil[StencilType::Right]; */ \
        /* solver_matrix.row_index(nz_index) = center_index + 1; */ \
        /* solver_matrix.col_index(nz_index) = right_index + 1; */ \
        /* solver_matrix.value(nz_index) += - coeff2 * arr; // Right */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Bottom]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = bottom_index + 1; \
        solver_matrix.value(nz_index) += - coeff3 * att; /* Bottom */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Top]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = top_index + 1; \
        solver_matrix.value(nz_index) += - coeff4 * att; /* Top */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        /* Center: (Left, Right, Bottom, Top) */ \
        solver_matrix.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; \
        \
        /* Fill matrix row of (i-1,j) */ \
        const Stencil& LeftStencil = getStencil(i_r-1); \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::Right]; \
        solver_matrix.row_index(nz_index) = left_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        solver_matrix.value(nz_index) += - coeff1 * arr; /* Right */ \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = left_index + 1; \
        solver_matrix.col_index(nz_index) = left_index + 1; \
        solver_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) */ \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::TopRight]; \
        solver_matrix.row_index(nz_index) = left_index + 1; \
        solver_matrix.col_index(nz_index) = top_index + 1; \
        solver_matrix.value(nz_index) += - 0.25 * art; /* Top Right */ \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::BottomRight]; \
        solver_matrix.row_index(nz_index) = left_index + 1; \
        solver_matrix.col_index(nz_index) = bottom_index + 1; \
        solver_matrix.value(nz_index) += 0.25 * art; /* Bottom Right */ \
        \
        /* Fill matrix row of (i+1,j) */ \
        /* Don't give to the outer dirichlet boundary! */ \
        \
        /* Fill matrix row of (i,j-1) */ \
        const Stencil& BottomStencil = CenterStencil; \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::Top]; \
        solver_matrix.row_index(nz_index) = bottom_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        solver_matrix.value(nz_index) += - coeff3 * att; /* Top */ \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = bottom_index + 1; \
        solver_matrix.col_index(nz_index) = bottom_index + 1; \
        solver_matrix.value(nz_index) += coeff3 * att; /* Center: (Top) */ \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        /* nz_index = bottom_nz_index + BottomStencil[StencilType::TopRight]; */ \
        /* solver_matrix.row_index(nz_index) = bottom_index + 1; */ \
        /* solver_matrix.col_index(nz_index) = right_index + 1; */ \
        /* solver_matrix.value(nz_index) += - 0.25 * art; // Top Right */ \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::TopLeft]; \
        solver_matrix.row_index(nz_index) = bottom_index + 1; \
        solver_matrix.col_index(nz_index) = left_index + 1; \
        solver_matrix.value(nz_index) += 0.25 * art; /* Top Left */ \
        \
        /* Fill matrix row of (i,j+1) */ \
        const Stencil& TopStencil = CenterStencil; \
        \
        nz_index = top_nz_index + TopStencil[StencilType::Bottom]; \
        solver_matrix.row_index(nz_index) = top_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        solver_matrix.value(nz_index) += - coeff4 * att; /* Bottom */ \
        \
        nz_index = top_nz_index + TopStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = top_index + 1; \
        solver_matrix.col_index(nz_index) = top_index + 1; \
        solver_matrix.value(nz_index) += coeff4 * att; /* Center: (Bottom) */ \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        /* nz_index = top_nz_index + TopStencil[StencilType::BottomRight]; */ \
        /* solver_matrix.row_index(nz_index) = top_index + 1; */ \
        /* solver_matrix.col_index(nz_index) = right_index + 1; */ \
        /* solver_matrix.value(nz_index) += 0.25 * art; // Bottom Right */ \
        \
        nz_index = top_nz_index + TopStencil[StencilType::BottomLeft]; \
        solver_matrix.row_index(nz_index) = top_index + 1; \
        solver_matrix.col_index(nz_index) = left_index + 1; \
        solver_matrix.value(nz_index) += - 0.25 * art; /* Bottom Left */ \
        \
    /* ------------------------------------ */ \
    /* Node on the outer dirichlet boundary */ \
    /* ------------------------------------ */ \
    } else if (i_r == grid.nr() - 1) { \
        double h1 = grid.radialSpacing(i_r-1); \
        double k1 = grid.angularSpacing(i_theta-1); \
        double k2 = grid.angularSpacing(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        \
        const int i_theta_M1 = grid.wrapThetaIndex(i_theta-1); \
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta+1); \
        \
        int nz_index; /* Current non_zero index in solver_matrix */ \
        \
        const int center_nz_index = getSolverMatrixIndex(i_r, i_theta); \
        const int left_nz_index = getSolverMatrixIndex(i_r-1, i_theta); \
        \
        const int center_index = grid.index(i_r, i_theta); \
        const int left_index = grid.index(i_r-1, i_theta); \
        const int bottom_index = grid.index(i_r, i_theta_M1); \
        const int top_index = grid.index(i_r, i_theta_P1); \
        \
        /* Fill matrix row of (i,j) */ \
        const Stencil& CenterStencil = getStencil(i_r); \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = center_index + 1; \
        solver_matrix.col_index(nz_index) = center_index + 1; \
        solver_matrix.value(nz_index) += 1.0; \
        \
        /* Give value to the interior nodes! */ \
        /* Fill matrix row of (i-1,j) */ \
        const Stencil& LeftStencil = getStencil(i_r-1); \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        /* nz_index = left_nz_index + LeftStencil[StencilType::Right]; */ \
        /* solver_matrix.row_index(nz_index) = left_index + 1; */ \
        /* solver_matrix.col_index(nz_index) = center_index + 1; */ \
        /* solver_matrix.value(nz_index) += - coeff1 * arr; // Right */ \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
        solver_matrix.row_index(nz_index) = left_index + 1; \
        solver_matrix.col_index(nz_index) = left_index + 1; \
        solver_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) */ \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        /* nz_index = left_nz_index + LeftStencil[StencilType::TopRight]; */ \
        /* solver_matrix.row_index(nz_index) = left_index + 1; */ \
        /* solver_matrix.col_index(nz_index) = top_index + 1; */ \
        /* solver_matrix.value(nz_index) += - 0.25 * art; // Top Right */ \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        /* nz_index = left_nz_index + LeftStencil[StencilType::BottomRight]; */ \
        /* solver_matrix.row_index(nz_index) = left_index + 1; */ \
        /* solver_matrix.col_index(nz_index) = bottom_index + 1; */ \
        /* solver_matrix.value(nz_index) += 0.25 * art; // Bottom Right */ \
        \
    } \
} while(0)


void DirectSolver::buildSolverMatrixCircleSection(const int i_r, SparseMatrix<double>& solver_matrix){
    const double r = grid_.radius(i_r);
    const double coeff_alpha = coeff_alpha_cache_[i_r];
    const double coeff_beta = coeff_beta_cache_[i_r];
    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++){
        const double theta = grid_.theta(i_theta);
        const double sin_theta = sin_theta_cache_[i_theta];
        const double cos_theta = cos_theta_cache_[i_theta];
        // Compute arr, att, art, detDF value at the current node 
        double arr, att, art, detDF;
        COMPUTE_JACOBIAN_ELEMENTS(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha,
            arr, att, art, detDF);
        // Build solver matrix at the current node
        NODE_BUILD_SOLVER_MATRIX_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta,
            grid_, DirBC_Interior_, solver_matrix,
            arr, att, art, detDF, coeff_beta);
    }
}

void DirectSolver::buildSolverMatrixRadialSection(const int i_theta, SparseMatrix<double>& solver_matrix){
    const double theta = grid_.theta(i_theta);
    const double sin_theta = sin_theta_cache_[i_theta];
    const double cos_theta = cos_theta_cache_[i_theta];
    for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++){
        const double r = grid_.radius(i_r);
        const double coeff_alpha = coeff_alpha_cache_[i_r];
        const double coeff_beta = coeff_beta_cache_[i_r];
        // Compute arr, att, art, detDF value at the current node 
        double arr, att, art, detDF;
        COMPUTE_JACOBIAN_ELEMENTS(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha,
            arr, att, art, detDF);
        // Build solver matrix at the current node
        NODE_BUILD_SOLVER_MATRIX_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta,
            grid_, DirBC_Interior_, solver_matrix,
            arr, att, art, detDF, coeff_beta);
    }
}


/* ------------------------------------------------------------------------ */
/* If the indexing is not smoother-based, please adjust the access patterns */
SparseMatrix<double> DirectSolver::buildSolverMatrix()
{
    omp_set_num_threads(num_omp_threads_);

    const int n = grid_.numberOfNodes();
    const int nnz = getNonZeroCountSolverMatrix();

    // Although the matrix is symmetric, we need to store all its entries, so we disable the symmetry.
    SparseMatrix<double> solver_matrix(n, n, nnz);
    solver_matrix.is_symmetric(false);

    #pragma omp parallel for if(nnz > 100'000)
    for (int i = 0; i < nnz; i++){
        solver_matrix.value(i) = 0.0;
    }

    if(num_omp_threads_ == 1) {
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
        const int num_circle_tasks = grid_.numberSmootherCircles();
        const int additional_radial_tasks = grid_.ntheta() % 3;
        const int num_radial_tasks = grid_.ntheta() - additional_radial_tasks;

        assert(num_circle_tasks >= 2);
        assert(num_radial_tasks >= 3 && num_radial_tasks % 3 == 0);

        /* Make sure to deallocate at the end */
        const int boundary_margin = 2; // Additional space to ensure safe access
        int* circle_dep = new int[num_circle_tasks + boundary_margin];
        int* radial_dep = new int[num_radial_tasks];

        #pragma omp parallel
        {
            #pragma omp single
            {
                /* ------------ */
                /* Circle Tasks */
                /* ------------ */

                /* Mod 0 Circles */
                for(int circle_task = 0; circle_task < num_circle_tasks; circle_task += 3) {
                    #pragma omp task \
                        depend(out: circle_dep[circle_task])
                    {
                        const int i_r = grid_.numberSmootherCircles() - circle_task - 1;    
                        buildSolverMatrixCircleSection(i_r, solver_matrix);
                    }
                }
                /* Mod 2 Circles */
                for(int circle_task = 1; circle_task < num_circle_tasks; circle_task += 3) {
                    #pragma omp task \
                        depend(out: circle_dep[circle_task]) \
                        depend(in: circle_dep[circle_task-1], circle_dep[circle_task+2])   
                    {
                        const int i_r = grid_.numberSmootherCircles() - circle_task - 1;    
                        buildSolverMatrixCircleSection(i_r, solver_matrix);
                    }
                    
                }
                /* Mod 2 Circles */
                for(int circle_task = 2; circle_task < num_circle_tasks; circle_task += 3) {
                        #pragma omp task \
                            depend(out: circle_dep[circle_task]) \
                            depend(in: circle_dep[circle_task-1], circle_dep[circle_task+2])   
                        {
                            const int i_r = grid_.numberSmootherCircles() - circle_task - 1;    
                            buildSolverMatrixCircleSection(i_r, solver_matrix);
                        }
                }

                /* ------------ */
                /* Radial Tasks */
                /* ------------ */

                /* Mod 0 Radials */
                for(int radial_task = 0; radial_task < num_radial_tasks; radial_task += 3) {
                    #pragma omp task \
                        depend(out: radial_dep[radial_task]) \
                        depend(in: circle_dep[1])
                    {
                        if(radial_task > 2){
                            const int i_theta = radial_task + additional_radial_tasks;    
                            buildSolverMatrixRadialSection(i_theta, solver_matrix);
                        } else{
                            if(additional_radial_tasks == 0){
                                buildSolverMatrixRadialSection(0, solver_matrix);
                            } 
                            else if(additional_radial_tasks >= 1){
                                buildSolverMatrixRadialSection(0, solver_matrix);
                                buildSolverMatrixRadialSection(1, solver_matrix);
                            }
                        }
                    }
                }
                /* Mod 1 Radials */
                for(int radial_task = 1; radial_task < num_radial_tasks; radial_task += 3) {
                    #pragma omp task \
                        depend(out: radial_dep[radial_task]) \
                        depend(in: radial_dep[radial_task-1], radial_dep[(radial_task+2) % num_radial_tasks])   
                    {
                        if(radial_task > 2){
                            int i_theta = radial_task + additional_radial_tasks;    
                            buildSolverMatrixRadialSection(i_theta, solver_matrix);
                        } else {
                            if(additional_radial_tasks == 0){
                                buildSolverMatrixRadialSection(1, solver_matrix);
                            } 
                            else if(additional_radial_tasks == 1){
                                buildSolverMatrixRadialSection(2, solver_matrix);
                            }
                            else if(additional_radial_tasks == 2){
                                buildSolverMatrixRadialSection(2, solver_matrix);
                                buildSolverMatrixRadialSection(3, solver_matrix);
                            }
                        }
                    }
                }
                /* Mod 2 Radials */
                for(int radial_task = 2; radial_task < num_radial_tasks; radial_task += 3) {
                    #pragma omp task \
                        depend(out: radial_dep[radial_task]) \
                        depend(in: radial_dep[radial_task-1], radial_dep[(radial_task+2) % num_radial_tasks])   
                    {
                        int i_theta = radial_task + additional_radial_tasks;    
                        buildSolverMatrixRadialSection(i_theta, solver_matrix);
                    }
                }
            }
        }
        delete[] circle_dep;
        delete[] radial_dep;
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