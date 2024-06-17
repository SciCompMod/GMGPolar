#include "../../include/CoarseSolver/coarseSolver.h"


#define ARR_ATT_ART(domain_geometry, r, theta, sin_theta, cos_theta, coeff_alpha, \
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



#define NODE_BUILD_A_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, \
    grid, DirBC_Interior, \
    matrixA, nz_index, \
    center_nz_index, left_nz_index, right_nz_index, bottom_nz_index, top_nz_index, \
    center_index, left_index, right_index, bottom_index, top_index, \
    arr, att, art, coeff_beta, detDF) \
do { \
    /* -------------------- */ \
    /* Node in the interior */ \
    /* -------------------- */ \
    if (i_r > 1 && i_r < grid.nr() - 2) { \
        double h1 = grid.r_dist(i_r-1); \
        double h2 = grid.r_dist(i_r); \
        double k1 = grid.theta_dist(i_theta-1); \
        double k2 = grid.theta_dist(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        \
        center_nz_index = ptr_nz_index_matrixA(i_r, i_theta); \
        left_nz_index = ptr_nz_index_matrixA(i_r-1, i_theta); \
        right_nz_index = ptr_nz_index_matrixA(i_r+1, i_theta); \
        bottom_nz_index = ptr_nz_index_matrixA(i_r, grid.wrap_theta_index(i_theta-1)); \
        top_nz_index = ptr_nz_index_matrixA(i_r, grid.wrap_theta_index(i_theta+1)); \
        \
        center_index = grid.index(i_r,i_theta); \
        left_index = grid.index(i_r-1,i_theta); \
        right_index = grid.index(i_r+1,i_theta); \
        bottom_index = grid.index(i_r,i_theta-1); \
        top_index = grid.index(i_r,i_theta+1); \
        \
        /* Fill matrix row of (i,j) */ \
        const Stencil& CenterStencil = get_stencil(i_r); \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = center_index + 1; \
        matrixA.col_index(nz_index) = center_index + 1; \
        matrixA.value(nz_index) += 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Left]; \
        matrixA.row_index(nz_index) = center_index + 1; \
        matrixA.col_index(nz_index) = left_index + 1; \
        matrixA.value(nz_index) += - coeff1 * arr; /* Left */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Right]; \
        matrixA.row_index(nz_index) = center_index + 1; \
        matrixA.col_index(nz_index) = right_index + 1; \
        matrixA.value(nz_index) += - coeff2 * arr; /* Right */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Bottom]; \
        matrixA.row_index(nz_index) = center_index + 1; \
        matrixA.col_index(nz_index) = bottom_index + 1; \
        matrixA.value(nz_index) += - coeff3 * att; /* Bottom */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Top]; \
        matrixA.row_index(nz_index) = center_index + 1; \
        matrixA.col_index(nz_index) = top_index + 1; \
        matrixA.value(nz_index) += - coeff4 * att; /* Top */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = center_index + 1; \
        matrixA.col_index(nz_index) = center_index + 1; \
        /* Center: (Left, Right, Bottom, Top) */ \
        matrixA.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; \
        \
        /* Fill matrix row of (i-1,j) */ \
        const Stencil& LeftStencil = get_stencil(i_r-1); \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::Right]; \
        matrixA.row_index(nz_index) = left_index + 1; \
        matrixA.col_index(nz_index) = center_index + 1; \
        matrixA.value(nz_index) += - coeff1 * arr; /* Right */ \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = left_index + 1; \
        matrixA.col_index(nz_index) = left_index + 1; \
        matrixA.value(nz_index) += coeff1 * arr; /* Center: (Right) */ \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::TopRight]; \
        matrixA.row_index(nz_index) = left_index + 1; \
        matrixA.col_index(nz_index) = top_index + 1; \
        matrixA.value(nz_index) += - 0.25 * art; /* Top Right */ \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::BottomRight]; \
        matrixA.row_index(nz_index) = left_index + 1; \
        matrixA.col_index(nz_index) = bottom_index + 1; \
        matrixA.value(nz_index) += 0.25 * art; /* Bottom Right */ \
        \
        /* Fill matrix row of (i+1,j) */ \
        const Stencil& RightStencil = get_stencil(i_r+1); \
        \
        nz_index = right_nz_index + RightStencil[StencilType::Left]; \
        matrixA.row_index(nz_index) = right_index + 1; \
        matrixA.col_index(nz_index) = center_index + 1; \
        matrixA.value(nz_index) += - coeff2 * arr; /* Left */ \
        \
        nz_index = right_nz_index + RightStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = right_index + 1; \
        matrixA.col_index(nz_index) = right_index + 1; \
        matrixA.value(nz_index) += coeff2 * arr; /* Center: (Left) */ \
        \
        nz_index = right_nz_index + RightStencil[StencilType::TopLeft]; \
        matrixA.row_index(nz_index) = right_index + 1; \
        matrixA.col_index(nz_index) = top_index + 1; \
        matrixA.value(nz_index) += 0.25 * art; /* Top Left */ \
        \
        nz_index = right_nz_index + RightStencil[StencilType::BottomLeft]; \
        matrixA.row_index(nz_index) = right_index + 1; \
        matrixA.col_index(nz_index) = bottom_index + 1; \
        matrixA.value(nz_index) += - 0.25 * art; /* Bottom Left */ \
        \
        /* Fill matrix row of (i,j-1) */ \
        const Stencil& BottomStencil = CenterStencil; \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::Top]; \
        matrixA.row_index(nz_index) = bottom_index + 1; \
        matrixA.col_index(nz_index) = center_index + 1; \
        matrixA.value(nz_index) += - coeff3 * att; /* Top */ \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = bottom_index + 1; \
        matrixA.col_index(nz_index) = bottom_index + 1; \
        matrixA.value(nz_index) += coeff3 * att; /* Center: (Top) */ \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::TopRight]; \
        matrixA.row_index(nz_index) = bottom_index + 1; \
        matrixA.col_index(nz_index) = right_index + 1; \
        matrixA.value(nz_index) += - 0.25 * art; /* Top Right */ \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::TopLeft]; \
        matrixA.row_index(nz_index) = bottom_index + 1; \
        matrixA.col_index(nz_index) = left_index + 1; \
        matrixA.value(nz_index) += 0.25 * art; /* Top Left */ \
        \
        /* Fill matrix row of (i,j+1) */ \
        const Stencil& TopStencil = CenterStencil; \
        \
        nz_index = top_nz_index + TopStencil[StencilType::Bottom]; \
        matrixA.row_index(nz_index) = top_index + 1; \
        matrixA.col_index(nz_index) = center_index + 1; \
        matrixA.value(nz_index) += - coeff4 * att; /* Bottom */ \
        \
        nz_index = top_nz_index + TopStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = top_index + 1; \
        matrixA.col_index(nz_index) = top_index + 1; \
        matrixA.value(nz_index) += coeff4 * att; /* Center: (Bottom) */ \
        \
        nz_index = top_nz_index + TopStencil[StencilType::BottomRight]; \
        matrixA.row_index(nz_index) = top_index + 1; \
        matrixA.col_index(nz_index) = right_index + 1; \
        matrixA.value(nz_index) += 0.25 * art; /* Bottom Right */ \
        \
        nz_index = top_nz_index + TopStencil[StencilType::BottomLeft]; \
        matrixA.row_index(nz_index) = top_index + 1; \
        matrixA.col_index(nz_index) = left_index + 1; \
        matrixA.value(nz_index) += - 0.25 * art; /* Bottom Left */ \
        \
    /* -------------------------- */ \
    /* Node on the inner boundary */ \
    /* -------------------------- */ \
    } else if (i_r == 0) { \
        /* ------------------------------------------------ */ \
        /* Case 1: Dirichlet boundary on the inner boundary */ \
        /* ------------------------------------------------ */ \
        if(DirBC_Interior){ \
            double h2 = grid.r_dist(i_r); \
            double k1 = grid.theta_dist(i_theta-1); \
            double k2 = grid.theta_dist(i_theta); \
            double coeff2 = 0.5*(k1+k2)/h2; \
            \
            center_nz_index = ptr_nz_index_matrixA(i_r, i_theta); \
            right_nz_index = ptr_nz_index_matrixA(i_r+1, i_theta); \
            \
            center_index = grid.index(i_r,i_theta); \
            right_index = grid.index(i_r+1,i_theta); \
            bottom_index = grid.index(i_r,i_theta-1); \
            top_index = grid.index(i_r,i_theta+1); \
            \
            /* Fill matrix row of (i,j) */ \
            const Stencil& CenterStencil = get_stencil(i_r); \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
            matrixA.row_index(nz_index) = center_index + 1; \
            matrixA.col_index(nz_index) = center_index + 1; \
            matrixA.value(nz_index) += 1.0; \
            \
            /* Give values to the interior nodes! */ \
            /* Fill matrix row of (i+1,j) */ \
            const Stencil& RightStencil = get_stencil(i_r+1); \
            \
            /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
            /* nz_index = right_nz_index + RightStencil[StencilType::Left]; */ \
            /* matrixA.row_index(nz_index) = right_index + 1; */ \
            /* matrixA.col_index(nz_index) = center_index + 1; */ \
            /* matrixA.value(nz_index) += - coeff2 * arr; // Left */ \
            \
            nz_index = right_nz_index + RightStencil[StencilType::Center]; \
            matrixA.row_index(nz_index) = right_index + 1; \
            matrixA.col_index(nz_index) = right_index + 1; \
            matrixA.value(nz_index) += coeff2 * arr; /* Center: (Left) */ \
            \
            /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
            /* nz_index = right_nz_index + RightStencil[StencilType::TopLeft]; */ \
            /* matrixA.row_index(nz_index) = right_index + 1; */ \
            /* matrixA.col_index(nz_index) = top_index + 1; */ \
            /* matrixA.value(nz_index) += 0.25 * art; // Top Left */ \
            \
            /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
            /* nz_index = right_nz_index + RightStencil[StencilType::BottomLeft]; */ \
            /* matrixA.row_index(nz_index) = right_index + 1; */ \
            /* matrixA.col_index(nz_index) = bottom_index + 1; */ \
            /* matrixA.value(nz_index) += - 0.25 * art; // Bottom Left */ \
            \
        } else{ \
            /* ------------------------------------------------------------- */ \
            /* Case 2: Across origin discretization on the interior boundary */ \
            /* ------------------------------------------------------------- */ \
            /* h1 gets replaced with 2 * R0. */ \
            /* (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()>>1)). */ \
            /* Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil. */ \
            double h1 = 2 * grid.radius(0); \
            double h2 = grid.r_dist(i_r); \
            double k1 = grid.theta_dist(i_theta-1); \
            double k2 = grid.theta_dist(i_theta); \
            double coeff1 = 0.5*(k1+k2)/h1; \
            double coeff2 = 0.5*(k1+k2)/h2; \
            double coeff3 = 0.5*(h1+h2)/k1; \
            double coeff4 = 0.5*(h1+h2)/k2; \
            \
            center_nz_index = ptr_nz_index_matrixA(i_r, i_theta); \
            left_nz_index = ptr_nz_index_matrixA(i_r, grid.wrap_theta_index(i_theta + (grid.ntheta()>>1))); \
            right_nz_index = ptr_nz_index_matrixA(i_r+1, i_theta); \
            bottom_nz_index = ptr_nz_index_matrixA(i_r, grid.wrap_theta_index(i_theta-1)); \
            top_nz_index = ptr_nz_index_matrixA(i_r, grid.wrap_theta_index(i_theta+1)); \
            \
            center_index = grid.index(i_r,i_theta); \
            left_index = grid.index(i_r, i_theta + (grid.ntheta()>>1)); \
            right_index = grid.index(i_r+1,i_theta); \
            bottom_index = grid.index(i_r,i_theta-1); \
            top_index = grid.index(i_r,i_theta+1); \
            \
            /* Fill matrix row of (i,j) */ \
            const Stencil& CenterStencil = get_stencil(i_r); \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
            matrixA.row_index(nz_index) = center_index + 1; \
            matrixA.col_index(nz_index) = center_index + 1; \
            matrixA.value(nz_index) += 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Left]; \
            matrixA.row_index(nz_index) = center_index + 1; \
            matrixA.col_index(nz_index) = left_index + 1; \
            matrixA.value(nz_index) += - coeff1 * arr; /* Left */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Right]; \
            matrixA.row_index(nz_index) = center_index + 1; \
            matrixA.col_index(nz_index) = right_index + 1; \
            matrixA.value(nz_index) += - coeff2 * arr; /* Right */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Bottom]; \
            matrixA.row_index(nz_index) = center_index + 1; \
            matrixA.col_index(nz_index) = bottom_index + 1; \
            matrixA.value(nz_index) += - coeff3 * att; /* Bottom */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Top]; \
            matrixA.row_index(nz_index) = center_index + 1; \
            matrixA.col_index(nz_index) = top_index + 1; \
            matrixA.value(nz_index) += - coeff4 * att; /* Top */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
            matrixA.row_index(nz_index) = center_index + 1; \
            matrixA.col_index(nz_index) = center_index + 1; \
            /* Center: (Left, Right, Bottom, Top) */ \
            matrixA.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; \
            \
            /* Fill matrix row of (i-1,j) */ \
            /* From view the view of the across origin node, the directions are roatated by 180 degrees in the stencil! */ \
            const Stencil& LeftStencil = CenterStencil; \
            \
            nz_index = left_nz_index + LeftStencil[StencilType::Left]; \
            matrixA.row_index(nz_index) = left_index + 1; \
            matrixA.col_index(nz_index) = center_index + 1; \
            matrixA.value(nz_index) += - coeff1 * arr; /* Right -> Left*/ \
            \
            nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
            matrixA.row_index(nz_index) = left_index + 1; \
            matrixA.col_index(nz_index) = left_index + 1; \
            matrixA.value(nz_index) += coeff1 * arr; /* Center: (Right) -> Center: (Left) */ \
            \
            /* Top Right -> Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = left_nz_index + LeftStencil[StencilType::BottomLeft]; */  \
            /* matrixA.row_index(nz_index) = left_index; */  \
            /* matrixA.col_index(nz_index) = top_index; */  \
            /* matrixA.value(nz_index) += - 0.25 * art; // Top Right -> Bottom Left*/ \
            \
            /* Bottom Right -> Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = left_nz_index + LeftStencil[StencilType::TopLeft]; */ \
            /* matrixA.row_index(nz_index) = left_index; */  \
            /* matrixA.col_index(nz_index) = bottom_index; */  \
            /* matrixA.value(nz_index) += 0.25 * art; // Bottom Right -> Top Left */ \
            \
            /* Fill matrix row of (i+1,j) */ \
            const Stencil& RightStencil = get_stencil(i_r+1); \
            \
            nz_index = right_nz_index + RightStencil[StencilType::Left]; \
            matrixA.row_index(nz_index) = right_index + 1; \
            matrixA.col_index(nz_index) = center_index + 1; \
            matrixA.value(nz_index) += - coeff2 * arr; /* Left */ \
            \
            nz_index = right_nz_index + RightStencil[StencilType::Center]; \
            matrixA.row_index(nz_index) = right_index + 1; \
            matrixA.col_index(nz_index) = right_index + 1; \
            matrixA.value(nz_index) += coeff2 * arr; /* Center: (Left) */ \
            \
            nz_index = right_nz_index + RightStencil[StencilType::TopLeft]; \
            matrixA.row_index(nz_index) = right_index + 1; \
            matrixA.col_index(nz_index) = top_index + 1; \
            matrixA.value(nz_index) += 0.25 * art; /* Top Left */ \
            \
            nz_index = right_nz_index + RightStencil[StencilType::BottomLeft]; \
            matrixA.row_index(nz_index) = right_index + 1; \
            matrixA.col_index(nz_index) = bottom_index + 1; \
            matrixA.value(nz_index) += - 0.25 * art; /* Bottom Left */ \
            \
            /* Fill matrix row of (i,j-1) */ \
            const Stencil& BottomStencil = CenterStencil; \
            \
            nz_index = bottom_nz_index + BottomStencil[StencilType::Top]; \
            matrixA.row_index(nz_index) = bottom_index + 1; \
            matrixA.col_index(nz_index) = center_index + 1; \
            matrixA.value(nz_index) += - coeff3 * att; /* Top */ \
            \
            nz_index = bottom_nz_index + BottomStencil[StencilType::Center]; \
            matrixA.row_index(nz_index) = bottom_index + 1; \
            matrixA.col_index(nz_index) = bottom_index + 1; \
            matrixA.value(nz_index) += coeff3 * att; /* Center: (Top) */ \
            \
            nz_index = bottom_nz_index + BottomStencil[StencilType::TopRight]; \
            matrixA.row_index(nz_index) = bottom_index + 1; \
            matrixA.col_index(nz_index) = right_index + 1; \
            matrixA.value(nz_index) += - 0.25 * art; /* Top Right */ \
            \
            /* REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = bottom_nz_index + BottomStencil[StencilType::TopLeft]; */ \
            /* matrixA.row_index(nz_index) = bottom_index + 1; */ \
            /* matrixA.col_index(nz_index) = left_index + 1; */ \
            /* matrixA.value(nz_index) += 0.25 * art; // Top Left */ \
            \
            /* Fill matrix row of (i,j+1) */ \
            const Stencil& TopStencil = CenterStencil; \
            \
            nz_index = top_nz_index + TopStencil[StencilType::Bottom]; \
            matrixA.row_index(nz_index) = top_index + 1; \
            matrixA.col_index(nz_index) = center_index + 1; \
            matrixA.value(nz_index) += - coeff4 * att; /* Bottom */ \
            \
            nz_index = top_nz_index + TopStencil[StencilType::Center]; \
            matrixA.row_index(nz_index) = top_index + 1; \
            matrixA.col_index(nz_index) = top_index + 1; \
            matrixA.value(nz_index) += coeff4 * att; /* Center: (Bottom) */ \
            \
            nz_index = top_nz_index + TopStencil[StencilType::BottomRight]; \
            matrixA.row_index(nz_index) = top_index + 1; \
            matrixA.col_index(nz_index) = right_index + 1; \
            matrixA.value(nz_index) += 0.25 * art; /* Bottom Right */ \
            \
            /* REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = top_nz_index + TopStencil[StencilType::BottomLeft]; */ \
            /* matrixA.row_index(nz_index) = top_index + 1; */ \
            /* matrixA.col_index(nz_index) = left_index + 1; */ \
            /* matrixA.value(nz_index) += - 0.25 * art; // Bottom Left */ \
            \
        } \
    /* ------------------------------- */ \
    /* Node next to the inner boundary */ \
    /* ------------------------------- */ \
    } else if (i_r == 1) { \
        double h1 = grid.r_dist(i_r-1); \
        double h2 = grid.r_dist(i_r); \
        double k1 = grid.theta_dist(i_theta-1); \
        double k2 = grid.theta_dist(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        \
        center_nz_index = ptr_nz_index_matrixA(i_r, i_theta); \
        left_nz_index = ptr_nz_index_matrixA(i_r-1, i_theta); \
        right_nz_index = ptr_nz_index_matrixA(i_r+1, i_theta); \
        bottom_nz_index = ptr_nz_index_matrixA(i_r, grid.wrap_theta_index(i_theta-1)); \
        top_nz_index = ptr_nz_index_matrixA(i_r, grid.wrap_theta_index(i_theta+1)); \
        \
        center_index = grid.index(i_r,i_theta); \
        left_index = grid.index(i_r-1,i_theta); \
        right_index = grid.index(i_r+1,i_theta); \
        bottom_index = grid.index(i_r,i_theta-1); \
        top_index = grid.index(i_r,i_theta+1); \
        \
        /* Fill matrix row of (i,j) */ \
        const Stencil& CenterStencil = get_stencil(i_r); \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = center_index + 1; \
        matrixA.col_index(nz_index) = center_index + 1; \
        matrixA.value(nz_index) += 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */ \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        if(!DirBC_Interior) { \
            nz_index = center_nz_index + CenterStencil[StencilType::Left]; \
            matrixA.row_index(nz_index) = center_index + 1; \
            matrixA.col_index(nz_index) = left_index + 1; \
            matrixA.value(nz_index) += - coeff1 * arr; /* Left */ \
        } \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Right]; \
        matrixA.row_index(nz_index) = center_index + 1; \
        matrixA.col_index(nz_index) = right_index + 1; \
        matrixA.value(nz_index) += - coeff2 * arr; /* Right */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Bottom]; \
        matrixA.row_index(nz_index) = center_index + 1; \
        matrixA.col_index(nz_index) = bottom_index + 1; \
        matrixA.value(nz_index) += - coeff3 * att; /* Bottom */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Top]; \
        matrixA.row_index(nz_index) = center_index + 1; \
        matrixA.col_index(nz_index) = top_index + 1; \
        matrixA.value(nz_index) += - coeff4 * att; /* Top */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = center_index + 1; \
        matrixA.col_index(nz_index) = center_index + 1; \
        /* Center: (Left, Right, Bottom, Top) */ \
        matrixA.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; \
        \
        if(!DirBC_Interior){ /* Don't give to the inner dirichlet boundary! */ \
            /* Fill matrix row of (i-1,j) */ \
            const Stencil& LeftStencil = get_stencil(i_r-1); \
            \
            nz_index = left_nz_index + LeftStencil[StencilType::Right]; \
            matrixA.row_index(nz_index) = left_index + 1; \
            matrixA.col_index(nz_index) = center_index + 1; \
            matrixA.value(nz_index) += - coeff1 * arr; /* Right */ \
            \
            nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
            matrixA.row_index(nz_index) = left_index + 1; \
            matrixA.col_index(nz_index) = left_index + 1; \
            matrixA.value(nz_index) += coeff1 * arr; /* Center: (Right) */ \
            \
            nz_index = left_nz_index + LeftStencil[StencilType::TopRight]; \
            matrixA.row_index(nz_index) = left_index + 1; \
            matrixA.col_index(nz_index) = top_index + 1; \
            matrixA.value(nz_index) += - 0.25 * art; /* Top Right */ \
            \
            nz_index = left_nz_index + LeftStencil[StencilType::BottomRight]; \
            matrixA.row_index(nz_index) = left_index + 1; \
            matrixA.col_index(nz_index) = bottom_index + 1; \
            matrixA.value(nz_index) += 0.25 * art; /* Bottom Right */ \
            \
        } \
        /* Fill matrix row of (i+1,j) */ \
        const Stencil& RightStencil = get_stencil(i_r+1); \
        \
        nz_index = right_nz_index + RightStencil[StencilType::Left]; \
        matrixA.row_index(nz_index) = right_index + 1; \
        matrixA.col_index(nz_index) = center_index + 1; \
        matrixA.value(nz_index) += - coeff2 * arr; /* Left */ \
        \
        nz_index = right_nz_index + RightStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = right_index + 1; \
        matrixA.col_index(nz_index) = right_index + 1; \
        matrixA.value(nz_index) += coeff2 * arr; /* Center: (Left) */ \
        \
        nz_index = right_nz_index + RightStencil[StencilType::TopLeft]; \
        matrixA.row_index(nz_index) = right_index + 1; \
        matrixA.col_index(nz_index) = top_index + 1; \
        matrixA.value(nz_index) += 0.25 * art; /* Top Left */ \
        \
        nz_index = right_nz_index + RightStencil[StencilType::BottomLeft]; \
        matrixA.row_index(nz_index) = right_index + 1; \
        matrixA.col_index(nz_index) = bottom_index + 1; \
        matrixA.value(nz_index) += - 0.25 * art; /* Bottom Left */ \
        \
        /* Fill matrix row of (i,j-1) */ \
        const Stencil& BottomStencil = CenterStencil; \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::Top]; \
        matrixA.row_index(nz_index) = bottom_index + 1; \
        matrixA.col_index(nz_index) = center_index + 1; \
        matrixA.value(nz_index) += - coeff3 * att; /* Top */ \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = bottom_index + 1; \
        matrixA.col_index(nz_index) = bottom_index + 1; \
        matrixA.value(nz_index) += coeff3 * att; /* Center: (Top) */ \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::TopRight]; \
        matrixA.row_index(nz_index) = bottom_index + 1; \
        matrixA.col_index(nz_index) = right_index + 1; \
        matrixA.value(nz_index) += - 0.25 * art; /* Top Right */ \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        if(!DirBC_Interior) { \
            nz_index = bottom_nz_index + BottomStencil[StencilType::TopLeft]; \
            matrixA.row_index(nz_index) = bottom_index + 1; \
            matrixA.col_index(nz_index) = left_index + 1; \
            matrixA.value(nz_index) += 0.25 * art; /* Top Left */ \
        } \
        \
        /* Fill matrix row of (i,j+1) */ \
        const Stencil& TopStencil = CenterStencil; \
        \
        nz_index = top_nz_index + TopStencil[StencilType::Bottom]; \
        matrixA.row_index(nz_index) = top_index + 1; \
        matrixA.col_index(nz_index) = center_index + 1; \
        matrixA.value(nz_index) += - coeff4 * att; /* Bottom */ \
        \
        nz_index = top_nz_index + TopStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = top_index + 1; \
        matrixA.col_index(nz_index) = top_index + 1; \
        matrixA.value(nz_index) += coeff4 * att; /* Center: (Bottom) */ \
        \
        nz_index = top_nz_index + TopStencil[StencilType::BottomRight]; \
        matrixA.row_index(nz_index) = top_index + 1; \
        matrixA.col_index(nz_index) = right_index + 1; \
        matrixA.value(nz_index) += 0.25 * art; /* Bottom Right */ \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        if(!DirBC_Interior) { \
            nz_index = top_nz_index + TopStencil[StencilType::BottomLeft]; \
            matrixA.row_index(nz_index) = top_index + 1; \
            matrixA.col_index(nz_index) = left_index + 1; \
            matrixA.value(nz_index) += - 0.25 * art; /* Bottom Left */ \
        } \
        \
    /* ------------------------------- */ \
    /* Node next to the outer boundary */ \
    /* ------------------------------- */ \
    } else if (i_r == grid.nr() - 2) { \
        double h1 = grid.r_dist(i_r-1); \
        double h2 = grid.r_dist(i_r); \
        double k1 = grid.theta_dist(i_theta-1); \
        double k2 = grid.theta_dist(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        \
        center_nz_index = ptr_nz_index_matrixA(i_r, i_theta); \
        left_nz_index = ptr_nz_index_matrixA(i_r-1, i_theta); \
        right_nz_index = ptr_nz_index_matrixA(i_r+1, i_theta); \
        bottom_nz_index = ptr_nz_index_matrixA(i_r, grid.wrap_theta_index(i_theta-1)); \
        top_nz_index = ptr_nz_index_matrixA(i_r, grid.wrap_theta_index(i_theta+1)); \
        \
        center_index = grid.index(i_r,i_theta); \
        left_index = grid.index(i_r-1,i_theta); \
        right_index = grid.index(i_r+1,i_theta); \
        bottom_index = grid.index(i_r,i_theta-1); \
        top_index = grid.index(i_r,i_theta+1); \
        \
        /* Fill matrix row of (i,j) */ \
        const Stencil& CenterStencil = get_stencil(i_r); \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = center_index + 1; \
        matrixA.col_index(nz_index) = center_index + 1; \
        matrixA.value(nz_index) += 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Left]; \
        matrixA.row_index(nz_index) = center_index + 1; \
        matrixA.col_index(nz_index) = left_index + 1; \
        matrixA.value(nz_index) += - coeff1 * arr; /* Left */ \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        /* nz_index = center_nz_index + CenterStencil[StencilType::Right]; */ \
        /* matrixA.row_index(nz_index) = center_index + 1; */ \
        /* matrixA.col_index(nz_index) = right_index + 1; */ \
        /* matrixA.value(nz_index) += - coeff2 * arr; // Right */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Bottom]; \
        matrixA.row_index(nz_index) = center_index + 1; \
        matrixA.col_index(nz_index) = bottom_index + 1; \
        matrixA.value(nz_index) += - coeff3 * att; /* Bottom */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Top]; \
        matrixA.row_index(nz_index) = center_index + 1; \
        matrixA.col_index(nz_index) = top_index + 1; \
        matrixA.value(nz_index) += - coeff4 * att; /* Top */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = center_index + 1; \
        matrixA.col_index(nz_index) = center_index + 1; \
        /* Center: (Left, Right, Bottom, Top) */ \
        matrixA.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; \
        \
        /* Fill matrix row of (i-1,j) */ \
        const Stencil& LeftStencil = get_stencil(i_r-1); \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::Right]; \
        matrixA.row_index(nz_index) = left_index + 1; \
        matrixA.col_index(nz_index) = center_index + 1; \
        matrixA.value(nz_index) += - coeff1 * arr; /* Right */ \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = left_index + 1; \
        matrixA.col_index(nz_index) = left_index + 1; \
        matrixA.value(nz_index) += coeff1 * arr; /* Center: (Right) */ \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::TopRight]; \
        matrixA.row_index(nz_index) = left_index + 1; \
        matrixA.col_index(nz_index) = top_index + 1; \
        matrixA.value(nz_index) += - 0.25 * art; /* Top Right */ \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::BottomRight]; \
        matrixA.row_index(nz_index) = left_index + 1; \
        matrixA.col_index(nz_index) = bottom_index + 1; \
        matrixA.value(nz_index) += 0.25 * art; /* Bottom Right */ \
        \
        /* Don't give to the outer dirichlet boundary! */ \
        /* Fill matrix row of (i+1,j) */ \
        /* const Stencil& RightStencil = get_stencil(i_r+1); */ \
        /* nz_index = right_nz_index + RightStencil[StencilType::Left]; */ \
        /* matrixA.row_index(nz_index) = right_index + 1; */ \
        /* matrixA.col_index(nz_index) = center_index + 1; */ \
        /* matrixA.value(nz_index) += - coeff2 * arr; // Left */ \
        /* nz_index = right_nz_index + RightStencil[StencilType::Center]; */ \
        /* matrixA.row_index(nz_index) = right_index + 1; */ \
        /* matrixA.col_index(nz_index) = right_index + 1; */ \
        /* matrixA.value(nz_index) += coeff2 * arr; // Center: (Left) */ \
        /* nz_index = right_nz_index + RightStencil[StencilType::TopLeft]; */ \
        /* matrixA.row_index(nz_index) = right_index + 1; */ \
        /* matrixA.col_index(nz_index) = top_index + 1; */ \
        /* matrixA.value(nz_index) += 0.25 * art; // Top Left */ \
        /* nz_index = right_nz_index + RightStencil[StencilType::BottomLeft]; */ \
        /* matrixA.row_index(nz_index) = right_index + 1; */ \
        /* matrixA.col_index(nz_index) = bottom_index + 1; */ \
        /* matrixA.value(nz_index) += - 0.25 * art; // Bottom Left */ \
        \
        /* Fill matrix row of (i,j-1) */ \
        const Stencil& BottomStencil = CenterStencil; \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::Top]; \
        matrixA.row_index(nz_index) = bottom_index + 1; \
        matrixA.col_index(nz_index) = center_index + 1; \
        matrixA.value(nz_index) += - coeff3 * att; /* Top */ \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = bottom_index + 1; \
        matrixA.col_index(nz_index) = bottom_index + 1; \
        matrixA.value(nz_index) += coeff3 * att; /* Center: (Top) */ \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        /* nz_index = bottom_nz_index + BottomStencil[StencilType::TopRight]; */ \
        /* matrixA.row_index(nz_index) = bottom_index + 1; */ \
        /* matrixA.col_index(nz_index) = right_index + 1; */ \
        /* matrixA.value(nz_index) += - 0.25 * art; // Top Right */ \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::TopLeft]; \
        matrixA.row_index(nz_index) = bottom_index + 1; \
        matrixA.col_index(nz_index) = left_index + 1; \
        matrixA.value(nz_index) += 0.25 * art; /* Top Left */ \
        \
        /* Fill matrix row of (i,j+1) */ \
        const Stencil& TopStencil = CenterStencil; \
        \
        nz_index = top_nz_index + TopStencil[StencilType::Bottom]; \
        matrixA.row_index(nz_index) = top_index + 1; \
        matrixA.col_index(nz_index) = center_index + 1; \
        matrixA.value(nz_index) += - coeff4 * att; /* Bottom */ \
        \
        nz_index = top_nz_index + TopStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = top_index + 1; \
        matrixA.col_index(nz_index) = top_index + 1; \
        matrixA.value(nz_index) += coeff4 * att; /* Center: (Bottom) */ \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        /* nz_index = top_nz_index + TopStencil[StencilType::BottomRight]; */ \
        /* matrixA.row_index(nz_index) = top_index + 1; */ \
        /* matrixA.col_index(nz_index) = right_index + 1; */ \
        /* matrixA.value(nz_index) += 0.25 * art; // Bottom Right */ \
        \
        nz_index = top_nz_index + TopStencil[StencilType::BottomLeft]; \
        matrixA.row_index(nz_index) = top_index + 1; \
        matrixA.col_index(nz_index) = left_index + 1; \
        matrixA.value(nz_index) += - 0.25 * art; /* Bottom Left */ \
        \
    /* ------------------------------------ */ \
    /* Node on the outer dirichlet boundary */ \
    /* ------------------------------------ */ \
    } else if (i_r == grid.nr() - 1) { \
        double h1 = grid.r_dist(i_r-1); \
        double k1 = grid.theta_dist(i_theta-1); \
        double k2 = grid.theta_dist(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        \
        center_nz_index = ptr_nz_index_matrixA(i_r, i_theta); \
        left_nz_index = ptr_nz_index_matrixA(i_r-1, i_theta); \
        \
        center_index = grid.index(i_r,i_theta); \
        left_index = grid.index(i_r-1,i_theta); \
        bottom_index = grid.index(i_r,i_theta-1); \
        top_index = grid.index(i_r,i_theta+1); \
        \
        /* Fill matrix row of (i,j) */ \
        const Stencil& CenterStencil = get_stencil(i_r); \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = center_index + 1; \
        matrixA.col_index(nz_index) = center_index + 1; \
        matrixA.value(nz_index) += 1.0; \
        \
        /* Give value to the interior nodes! */ \
        /* Fill matrix row of (i-1,j) */ \
        const Stencil& LeftStencil = get_stencil(i_r-1); \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        /* nz_index = left_nz_index + LeftStencil[StencilType::Right]; */ \
        /* matrixA.row_index(nz_index) = left_index + 1; */ \
        /* matrixA.col_index(nz_index) = center_index + 1; */ \
        /* matrixA.value(nz_index) += - coeff1 * arr; // Right */ \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = left_index + 1; \
        matrixA.col_index(nz_index) = left_index + 1; \
        matrixA.value(nz_index) += coeff1 * arr; /* Center: (Right) */ \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        /* nz_index = left_nz_index + LeftStencil[StencilType::TopRight]; */ \
        /* matrixA.row_index(nz_index) = left_index + 1; */ \
        /* matrixA.col_index(nz_index) = top_index + 1; */ \
        /* matrixA.value(nz_index) += - 0.25 * art; // Top Right */ \
        \
        /* REMOVED: Moved to the right hand side to make the matrix symmetric */ \
        /* nz_index = left_nz_index + LeftStencil[StencilType::BottomRight]; */ \
        /* matrixA.row_index(nz_index) = left_index + 1; */ \
        /* matrixA.col_index(nz_index) = bottom_index + 1; */ \
        /* matrixA.value(nz_index) += 0.25 * art; // Bottom Right */ \
        \
    } \
} while(0)



#define CIRCLE_SECTION_BUILD_A_GIVE(i_r) \
do { \
    r = grid_.radius(i_r); \
    coeff_alpha = system_parameters_.alpha(r); \
    coeff_beta = system_parameters_.beta(r); \
    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++){ \
        theta = grid_.theta(i_theta); \
        sin_theta = sin_theta_[i_theta]; \
        cos_theta = cos_theta_[i_theta]; \
        /* Get arr, att, art, detDF value at the current node */ \
        ARR_ATT_ART(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, \
            arr, att, art, detDF); \
        \
        NODE_BUILD_A_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, \
            grid_, DirBC_Interior_, \
            matrixA, nz_index, \
            center_nz_index, left_nz_index, right_nz_index, bottom_nz_index, top_nz_index, \
            center_index, left_index, right_index, bottom_index, top_index, \
            arr, att, art, coeff_beta, detDF); \
    } \
} while(0)


#define RADIAL_SECTION_NODE_BUILD_A_GIVE(i_theta) \
do { \
    theta = grid_.theta(i_theta); \
    sin_theta = sin_theta_[i_theta]; \
    cos_theta = cos_theta_[i_theta]; \
    for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++){ \
        r = grid_.radius(i_r); \
        coeff_alpha = system_parameters_.alpha(r); \
        coeff_beta = system_parameters_.beta(r); \
        /* Get arr, att, art, detDF value at the current node */ \
        ARR_ATT_ART(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, \
            arr, att, art, detDF); \
        \
        NODE_BUILD_A_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, \
            grid_, DirBC_Interior_, \
            matrixA, nz_index, \
            center_nz_index, left_nz_index, right_nz_index, bottom_nz_index, top_nz_index, \
            center_index, left_index, right_index, bottom_index, top_index, \
            arr, att, art, coeff_beta, detDF); \
    } \
} while(0)


/* ---------------------------------- */
/* Build symmetric matrix A           */
/* -> shift non-symmetric part to rhs */
/* ---------------------------------- */

void CoarseSolver::buildMatrixA(SparseMatrix<double>& symetric_matrixA)
{
    omp_set_num_threads(maxOpenMPThreads_);

    const int n = grid_.number_of_nodes();
    const int matrixA_nnz = nnz_matrixA();

    SparseMatrix<double> matrixA(n, n, matrixA_nnz);

    #pragma omp parallel for
    for (int i = 0; i < matrixA_nnz; i++){
        matrixA.value(i) = 0.0;
    }

    const int numCircleTasks = grid_.numberSmootherCircles();
    const int additionalRadialTasks = grid_.ntheta() % 3;
    const int numRadialTasks = grid_.ntheta() - additionalRadialTasks;

    assert(numCircleTasks >= 2);
    assert(numRadialTasks >= 3 && numRadialTasks % 3 == 0);

    /* Make sure to deallocate at the end */
    int* dep = new int[numCircleTasks + numRadialTasks];

    omp_set_num_threads(openMPTaskThreads_);
    #pragma omp parallel num_threads(openMPTaskThreads_) /* Outside variable are shared by default */
    {
        /* Define thread-local variables */
        double r, theta;
        double sin_theta, cos_theta;
        double arr, att, art;
        double coeff_alpha, coeff_beta;
        double detDF;

        int center_nz_index;
        int left_nz_index, right_nz_index;
        int bottom_nz_index, top_nz_index;

        int center_index;
        int left_index, right_index;
        int bottom_index, top_index;

        int nz_index;
        #pragma omp single
        {
            /* ------------ */
            /* Circle Tasks */
            /* ------------ */

            /* Mod 0 Circles */
            for(int circle_task = 0; circle_task < numCircleTasks; circle_task += 3) {
                #pragma omp task \
                    depend(out: dep[circle_task])
                {
                    int i_r = grid_.numberSmootherCircles() - circle_task - 1;    
                    CIRCLE_SECTION_BUILD_A_GIVE(i_r);
                }
            }
            /* Mod 2 Circles */
            for(int circle_task = 1; circle_task < numCircleTasks; circle_task += 3) {
                #pragma omp task \
                    depend(out: dep[circle_task]) \
                    depend(in: dep[circle_task-1], dep[circle_task+2])   
                {
                    int i_r = grid_.numberSmootherCircles() - circle_task - 1;    
                    CIRCLE_SECTION_BUILD_A_GIVE(i_r);
                }
                
            }
            /* Mod 2 Circles */
            for(int circle_task = 2; circle_task < numCircleTasks; circle_task += 3) {
                    #pragma omp task \
                        depend(out: dep[circle_task]) \
                        depend(in: dep[circle_task-1], dep[circle_task+2])   
                    {
                        int i_r = grid_.numberSmootherCircles() - circle_task - 1;    
                        CIRCLE_SECTION_BUILD_A_GIVE(i_r);
                    }
            }

            /* ------------ */
            /* Radial Tasks */
            /* ------------ */

            /* Mod 0 Radials */
            for(int radial_task = 0; radial_task < numRadialTasks; radial_task += 3) {
                #pragma omp task \
                    depend(out: dep[numCircleTasks+radial_task]) \
                    depend(in: dep[1])   
                {
                    if(radial_task > 0){
                        int i_theta = radial_task + additionalRadialTasks;    
                        RADIAL_SECTION_NODE_BUILD_A_GIVE(i_theta);
                    } else{
                        if(additionalRadialTasks == 0){
                            RADIAL_SECTION_NODE_BUILD_A_GIVE(0);
                        } 
                        else if(additionalRadialTasks >= 1){
                            RADIAL_SECTION_NODE_BUILD_A_GIVE(0);
                            RADIAL_SECTION_NODE_BUILD_A_GIVE(1);
                        }
                    }
                }
            }
            /* Mod 1 Radials */
            for(int radial_task = 1; radial_task < numRadialTasks; radial_task += 3) {
                #pragma omp task \
                    depend(out: dep[numCircleTasks+radial_task]) \
                    depend(in: dep[1], dep[numCircleTasks+radial_task-1], dep[numCircleTasks+(radial_task+2) % numRadialTasks])   
                {
                    if(radial_task > 1){
                        int i_theta = radial_task + additionalRadialTasks;    
                        RADIAL_SECTION_NODE_BUILD_A_GIVE(i_theta);
                    } else {
                        if(additionalRadialTasks == 0){
                            RADIAL_SECTION_NODE_BUILD_A_GIVE(1);
                        } 
                        else if(additionalRadialTasks == 1){
                            RADIAL_SECTION_NODE_BUILD_A_GIVE(2);
                        }
                        else if(additionalRadialTasks == 2){
                            RADIAL_SECTION_NODE_BUILD_A_GIVE(2);
                            RADIAL_SECTION_NODE_BUILD_A_GIVE(3);
                        }
                    }
                }
            }
            /* Mod 2 Radials */
            for(int radial_task = 2; radial_task < numRadialTasks; radial_task += 3) {
                #pragma omp task \
                    depend(out: dep[numCircleTasks+radial_task]) \
                    depend(in:  dep[1], dep[numCircleTasks+radial_task-1], dep[numCircleTasks+(radial_task+2) % numRadialTasks])   
                {
                    int i_theta = radial_task + additionalRadialTasks;    
                    RADIAL_SECTION_NODE_BUILD_A_GIVE(i_theta);
                }
            }
        }
    }

    omp_set_num_threads(maxOpenMPThreads_);
    delete[] dep;

    const int symmetric_matrixA_nnz = matrixA_nnz - (matrixA_nnz - n) / 2;
    symetric_matrixA = SparseMatrix<double> (n, n, symmetric_matrixA_nnz);
    symetric_matrixA.is_symmetric(true);

    int current_nz = 0;
    for (int nz_index = 0; nz_index < matrixA.non_zero_size(); nz_index++) {
        int current_matrixA_row = matrixA.row_index(nz_index);
        int current_matrixA_col = matrixA.col_index(nz_index);
        if (current_matrixA_row <= current_matrixA_col) {
            symetric_matrixA.row_index(current_nz) = current_matrixA_row;
            symetric_matrixA.col_index(current_nz) = current_matrixA_col;
            symetric_matrixA.value(current_nz) = std::move(matrixA.value(nz_index));
            current_nz++;
        }
    }
}
