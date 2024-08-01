#include "../../include/Smoother/smoother.h"

/* arr, att, art and detDF get computed */
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



#define NODE_BUILD_SMOOTHER_GIVE(i_r, i_theta, grid, DirBC_Interior, \
    nonSymmetric_inner_boundary_circle_Asc_matrix, \
    circle_symmetric_cyclic_tridiagonal_solver, \
    radial_symmetric_tridiagonal_solver \
) \
do { \
    assert(i_r >= 0 && i_r < grid.nr()); \
    assert(i_theta >= 0 && i_theta < grid.ntheta()); \
    \
    const int numberSmootherCircles = grid.numberSmootherCircles(); \
    const int lengthSmootherRadial = grid.lengthSmootherRadial(); \
    \
    assert(numberSmootherCircles >= 2); \
    assert(lengthSmootherRadial >= 3); \
    \
    int row, column; \
    double value; \
    /* ------------------------------------------ */ \
    /* Node in the interior of the Circle Section */ \
    /* ------------------------------------------ */ \
    if (i_r > 0 && i_r < numberSmootherCircles) { \
        const double h1 = grid.r_dist(i_r-1); \
        const double h2 = grid.r_dist(i_r); \
        const double k1 = grid.theta_dist(i_theta-1); \
        const double k2 = grid.theta_dist(i_theta); \
        const double coeff1 = 0.5*(k1+k2)/h1; \
        const double coeff2 = 0.5*(k1+k2)/h2; \
        const double coeff3 = 0.5*(h1+h2)/k1; \
        const double coeff4 = 0.5*(h1+h2)/k2; \
        \
        const int i_theta_M1 = grid.wrap_theta_index(i_theta-1); \
        const int i_theta_P1 = grid.wrap_theta_index(i_theta+1); \
        \
        const int center_index = i_theta; \
        const int left_index = i_theta; \
        const int right_index = (i_r + 1 == numberSmootherCircles) ? 0 : i_theta; \
        const int bottom_index = i_theta_M1; \
        const int top_index = i_theta_P1; \
        \
        const int circle_tridiagonal_solver_index = i_r - 1; \
        auto& left_matrix = circle_symmetric_cyclic_tridiagonal_solver[circle_tridiagonal_solver_index-1]; \
        auto& center_matrix = circle_symmetric_cyclic_tridiagonal_solver[circle_tridiagonal_solver_index]; \
        auto& right_matrix = (i_r + 1 == numberSmootherCircles) ? \
            radial_symmetric_tridiagonal_solver[i_theta] : \
            circle_symmetric_cyclic_tridiagonal_solver[circle_tridiagonal_solver_index+1]; \
        \
        /* Fill matrix row of (i,j) */ \
        row = center_index; \
        column = center_index; \
        value = 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        row = center_index; \
        column = bottom_index; \
        value = - coeff3 * att; /* Bottom */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        row = center_index; \
        column = top_index; \
        value = - coeff4 * att; /* Top */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        row = center_index; \
        column = center_index; \
        value = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        /* Fill matrix row of (i,j-1) */ \
        row = bottom_index; \
        column = center_index; \
        value = - coeff3 * att; /* Top */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        row = bottom_index; \
        column = bottom_index; \
        value = coeff3 * att; /* Center: (Top) */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        /* Fill matrix row of (i,j+1) */ \
        row = top_index; \
        column = center_index; \
        value = - coeff4 * att; /* Bottom */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        row = top_index; \
        column = top_index; \
        value = coeff4 * att; /* Center: (Bottom) */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        /* Fill matrix row of (i-1,j) */ \
        if(!DirBC_Interior && i_r == 1) { \
            const Stencil& LeftStencil = get_stencil(i_r-1); \
            int left_nz_index = ptr_nz_index_circle_Asc(i_r-1, i_theta); \
            int nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
            nonSymmetric_inner_boundary_circle_Asc_matrix.row_index(nz_index) = left_index + 1; \
            nonSymmetric_inner_boundary_circle_Asc_matrix.col_index(nz_index) = left_index + 1; \
            nonSymmetric_inner_boundary_circle_Asc_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) */ \
        } \
        if(i_r > 1) { \
            row = left_index; \
            column = left_index; \
            value = coeff1 * arr; /* Center: (Right) */ \
            if(row == column) left_matrix.main_diagonal(row) += value; \
            else if(row == column - 1) left_matrix.sub_diagonal(row) += value; \
            else if(row == 0 && column == left_matrix.columns()-1) left_matrix.cyclic_corner_element() += value; \
        } \
        \
        /* Fill matrix row of (i+1,j) */ \
        row = right_index; \
        column = right_index; \
        value = coeff2 * arr; /* Center: (Left) */ \
        if(row == column) right_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) right_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == right_matrix.columns()-1) right_matrix.cyclic_corner_element() += value; \
    } \
    /* ------------------------------------------ */ \
    /* Node in the interior of the Radial Section */ \
    /* ------------------------------------------ */ \
    else if(i_r > numberSmootherCircles && i_r < grid.nr() - 2){ \
        const double h1 = grid.r_dist(i_r-1); \
        const double h2 = grid.r_dist(i_r); \
        const double k1 = grid.theta_dist(i_theta-1); \
        const double k2 = grid.theta_dist(i_theta); \
        const double coeff1 = 0.5*(k1+k2)/h1; \
        const double coeff2 = 0.5*(k1+k2)/h2; \
        const double coeff3 = 0.5*(h1+h2)/k1; \
        const double coeff4 = 0.5*(h1+h2)/k2; \
        \
        const int i_theta_M1 = grid.wrap_theta_index(i_theta-1); \
        const int i_theta_P1 = grid.wrap_theta_index(i_theta+1); \
        \
        auto& bottom_matrix = radial_symmetric_tridiagonal_solver[i_theta_M1]; \
        auto& center_matrix = radial_symmetric_tridiagonal_solver[i_theta]; \
        auto& top_matrix = radial_symmetric_tridiagonal_solver[i_theta_P1]; \
        \
        const int center_index = i_r - numberSmootherCircles; \
        const int left_index = i_r - numberSmootherCircles - 1; \
        const int right_index = i_r - numberSmootherCircles + 1; \
        const int bottom_index = i_r - numberSmootherCircles; \
        const int top_index = i_r - numberSmootherCircles; \
        \
        /* Fill matrix row of (i,j) */ \
        row = center_index; \
        column = center_index; \
        value = 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        row = center_index; \
        column = left_index; \
        value = - coeff1 * arr; /* Left */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        row = center_index; \
        column = right_index; \
        value = - coeff2 * arr; /* Right */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        row = center_index; \
        column = center_index; \
        value = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        /* Fill matrix row of (i-1,j) */ \
        row = left_index; \
        column = center_index; \
        value = - coeff1 * arr; /* Right */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        row = left_index; \
        column = left_index; \
        value = coeff1 * arr; /* Center: (Right) */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        /* Fill matrix row of (i+1,j) */ \
        row = right_index; \
        column = center_index; \
        value = - coeff2 * arr; /* Left */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        row = right_index; \
        column = right_index; \
        value = coeff2 * arr; /* Center: (Left) */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        /* Fill matrix row of (i,j-1) */ \
        row = bottom_index; \
        column = bottom_index; \
        value = coeff3 * att; /* Center: (Top) */ \
        if(row == column) bottom_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) bottom_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == bottom_matrix.columns()-1) bottom_matrix.cyclic_corner_element() += value; \
        \
        /* Fill matrix row of (i,j+1) */ \
        row = top_index; \
        column = top_index; \
        value = coeff4 * att; /* Center: (Bottom) */ \
        if(row == column) top_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) top_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == top_matrix.columns()-1) top_matrix.cyclic_corner_element() += value; \
    } \
    /* ------------------------------------------ */ \
    /* Circle Section: Node in the inner boundary */ \
    /* ------------------------------------------ */ \
    else if(i_r == 0){ \
        /* ------------------------------------------------ */ \
        /* Case 1: Dirichlet boundary on the inner boundary */ \
        /* ------------------------------------------------ */ \
        if(DirBC_Interior){ \
            /* Fill result(i,j) */ \
            const double h2 = grid.r_dist(i_r); \
            const double k1 = grid.theta_dist(i_theta-1); \
            const double k2 = grid.theta_dist(i_theta); \
            const double coeff2 = 0.5*(k1+k2)/h2; \
            \
            const int i_theta_M1 = grid.wrap_theta_index(i_theta-1); \
            const int i_theta_P1 = grid.wrap_theta_index(i_theta+1); \
            \
            const int circle_tridiagonal_solver_index = i_r - 1; \
            auto& center_matrix = nonSymmetric_inner_boundary_circle_Asc_matrix; \
            auto& right_matrix = circle_symmetric_cyclic_tridiagonal_solver[circle_tridiagonal_solver_index+1]; \
            \
            const int center_index = i_theta; \
            const int right_index = i_theta; \
            const int bottom_index = i_theta_M1; \
            const int top_index = i_theta_P1; \
            \
            /* Fill matrix row of (i,j) */ \
            const Stencil& CenterStencil = get_stencil(i_r); \
            int center_nz_index = ptr_nz_index_circle_Asc(i_r, i_theta); \
            int nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
            center_matrix.row_index(nz_index) = center_index + 1; \
            center_matrix.col_index(nz_index) = center_index + 1; \
            center_matrix.value(nz_index) += 1.0; \
            /* Fill matrix row of (i+1,j) */ \
            row = right_index; \
            column = right_index; \
            value = coeff2 * arr; /* Center: (Left) */ \
            if(row == column) right_matrix.main_diagonal(row) += value; \
            else if(row == column - 1) right_matrix.sub_diagonal(row) += value; \
            else if(row == 0 && column == right_matrix.columns()-1) right_matrix.cyclic_corner_element() += value; \
        } else{ \
            /* ------------------------------------------------------------- */ \
            /* Case 2: Across origin discretization on the interior boundary */ \
            /* ------------------------------------------------------------- */ \
            /* h1 gets replaced with 2 * R0. */ \
            /* (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()>>1)). */ \
            /* Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil. */ \
            const double h1 = 2 * grid.radius(0); \
            const double h2 = grid.r_dist(i_r); \
            const double k1 = grid.theta_dist(i_theta-1); \
            const double k2 = grid.theta_dist(i_theta); \
            const double coeff1 = 0.5*(k1+k2)/h1; \
            const double coeff2 = 0.5*(k1+k2)/h2; \
            const double coeff3 = 0.5*(h1+h2)/k1; \
            const double coeff4 = 0.5*(h1+h2)/k2; \
            \
            const int circle_tridiagonal_solver_index = i_r - 1; \
            auto& center_matrix = nonSymmetric_inner_boundary_circle_Asc_matrix; \
            auto& right_matrix = circle_symmetric_cyclic_tridiagonal_solver[circle_tridiagonal_solver_index+1]; \
            auto& left_matrix = nonSymmetric_inner_boundary_circle_Asc_matrix; \
            \
            const int i_theta_M1 = grid.wrap_theta_index(i_theta-1); \
            const int i_theta_P1 = grid.wrap_theta_index(i_theta+1); \
            const int i_theta_AcrossOrigin = grid.wrap_theta_index(i_theta + (grid.ntheta()>>1)); \
            \
            const int center_index = i_theta; \
            const int left_index = i_theta_AcrossOrigin; \
            const int right_index = i_theta; \
            const int bottom_index = i_theta_M1; \
            const int top_index = i_theta_P1; \
            \
            const int center_nz_index = ptr_nz_index_circle_Asc(i_r, i_theta); \
            const int bottom_nz_index = ptr_nz_index_circle_Asc(i_r, i_theta_M1); \
            const int top_nz_index = ptr_nz_index_circle_Asc(i_r, i_theta_P1); \
            const int left_nz_index = ptr_nz_index_circle_Asc(i_r, i_theta_AcrossOrigin); \
            \
            int nz_index; \
            /* Fill matrix row of (i,j) */ \
            const Stencil& CenterStencil = get_stencil(i_r); \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
            center_matrix.row_index(nz_index) = center_index + 1; \
            center_matrix.col_index(nz_index) = center_index + 1; \
            center_matrix.value(nz_index) += 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Left]; \
            center_matrix.row_index(nz_index) = center_index + 1; \
            center_matrix.col_index(nz_index) = left_index + 1; \
            center_matrix.value(nz_index) += - coeff1 * arr; /* Left */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Bottom]; \
            center_matrix.row_index(nz_index) = center_index + 1; \
            center_matrix.col_index(nz_index) = bottom_index + 1; \
            center_matrix.value(nz_index) += - coeff3 * att; /* Bottom */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Top]; \
            center_matrix.row_index(nz_index) = center_index + 1; \
            center_matrix.col_index(nz_index) = top_index + 1; \
            center_matrix.value(nz_index) += - coeff4 * att; /* Top */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
            center_matrix.row_index(nz_index) = center_index + 1; \
            center_matrix.col_index(nz_index) = center_index + 1; \
            /* Center: (Left, Right, Bottom, Top) */ \
            center_matrix.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; \
            \
            /* Fill matrix row of (i-1,j) */ \
            /* From view the view of the across origin node, the directions are roatated by 180 degrees in the stencil! */ \
            const Stencil& LeftStencil = CenterStencil; \
            \
            nz_index = left_nz_index + LeftStencil[StencilType::Left]; \
            left_matrix.row_index(nz_index) = left_index + 1; \
            left_matrix.col_index(nz_index) = center_index + 1; \
            left_matrix.value(nz_index) += - coeff1 * arr; /* Right -> Left*/ \
            \
            nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
            left_matrix.row_index(nz_index) = left_index + 1; \
            left_matrix.col_index(nz_index) = left_index + 1; \
            left_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) -> Center: (Left) */ \
            \
            /* Top Right -> Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = left_nz_index + LeftStencil[StencilType::BottomLeft]; */  \
            /* left_matrix.row_index(nz_index) = left_index; */  \
            /* left_matrix.col_index(nz_index) = top_index; */  \
            /* left_matrix.value(nz_index) += - 0.25 * art; // Top Right -> Bottom Left*/ \
            \
            /* Bottom Right -> Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = left_nz_index + LeftStencil[StencilType::TopLeft]; */ \
            /* left_matrix.row_index(nz_index) = left_index; */  \
            /* left_matrix.col_index(nz_index) = bottom_index; */  \
            /* lefc_matrix.value(nz_index) += 0.25 * art; // Bottom Right -> Top Left */ \
            \
            /* Fill matrix row of (i+1,j) */ \
            row = right_index; \
            column = right_index; \
            value = coeff2 * arr; /* Center: (Left) */ \
            if(row == column) right_matrix.main_diagonal(row) += value; \
            else if(row == column - 1) right_matrix.sub_diagonal(row) += value; \
            else if(row == 0 && column == right_matrix.columns()-1) right_matrix.cyclic_corner_element() += value; \
            \
            /* Fill matrix row of (i,j-1) */ \
            const Stencil& BottomStencil = CenterStencil; \
            \
            nz_index = bottom_nz_index + BottomStencil[StencilType::Top]; \
            center_matrix.row_index(nz_index) = bottom_index + 1; \
            center_matrix.col_index(nz_index) = center_index + 1; \
            center_matrix.value(nz_index) += - coeff3 * att; /* Top */ \
            \
            nz_index = bottom_nz_index + BottomStencil[StencilType::Center]; \
            center_matrix.row_index(nz_index) = bottom_index + 1; \
            center_matrix.col_index(nz_index) = bottom_index + 1; \
            center_matrix.value(nz_index) += coeff3 * att; /* Center: (Top) */ \
            \
            /* REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = bottom_nz_index + BottomStencil[StencilType::TopLeft]; */ \
            /* center_matrix.row_index(nz_index) = bottom_index + 1; */ \
            /* center_matrix.col_index(nz_index) = left_index + 1; */ \
            /* center_matrix.value(nz_index) += 0.25 * art; // Top Left */ \
            \
            /* Fill matrix row of (i,j+1) */ \
            const Stencil& TopStencil = CenterStencil; \
            \
            nz_index = top_nz_index + TopStencil[StencilType::Bottom]; \
            center_matrix.row_index(nz_index) = top_index + 1; \
            center_matrix.col_index(nz_index) = center_index + 1; \
            center_matrix.value(nz_index) += - coeff4 * att; /* Bottom */ \
            \
            nz_index = top_nz_index + TopStencil[StencilType::Center]; \
            center_matrix.row_index(nz_index) = top_index + 1; \
            center_matrix.col_index(nz_index) = top_index + 1; \
            center_matrix.value(nz_index) += coeff4 * att; /* Center: (Bottom) */ \
            \
            /* REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = top_nz_index + TopStencil[StencilType::BottomLeft]; */ \
            /* center_matrix.row_index(nz_index) = top_index + 1; */ \
            /* center_matrix.col_index(nz_index) = left_index + 1; */ \
            /* center_matrix.value(nz_index) += - 0.25 * art; // Bottom Left */ \
        } \
    } \
    /* --------------------------------------------- */ \
    /* Radial Section: Node next to circular section */ \
    /* --------------------------------------------- */ \
    else if(i_r == numberSmootherCircles){ \
        double h1 = grid.r_dist(i_r-1); \
        double h2 = grid.r_dist(i_r); \
        double k1 = grid.theta_dist(i_theta-1); \
        double k2 = grid.theta_dist(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        \
        const int i_theta_M1 = grid.wrap_theta_index(i_theta-1); \
        const int i_theta_P1 = grid.wrap_theta_index(i_theta+1); \
        \
        const int circle_tridiagonal_solver_index = i_r - 1; \
        auto& bottom_matrix = radial_symmetric_tridiagonal_solver[i_theta_M1]; \
        auto& center_matrix = radial_symmetric_tridiagonal_solver[i_theta]; \
        auto& top_matrix = radial_symmetric_tridiagonal_solver[i_theta_P1]; \
        auto& left_matrix = circle_symmetric_cyclic_tridiagonal_solver[circle_tridiagonal_solver_index-1]; \
        \
        const int center_index = i_r - numberSmootherCircles; \
        const int left_index = i_theta; \
        const int right_index = i_r - numberSmootherCircles + 1; \
        const int bottom_index = i_r - numberSmootherCircles; \
        const int top_index = i_r - numberSmootherCircles; \
        \
        /* Fill matrix row of (i,j) */ \
        row = center_index; \
        column = center_index; \
        value = 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        row = center_index; \
        column = right_index; \
        value = - coeff2 * arr; /* Right */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        row = center_index; \
        column = center_index; \
        value = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        /* Fill matrix row of (i-1,j) */ \
        row = left_index; \
        column = left_index; \
        value = coeff1 * arr; /* Center: (Right) */ \
        if(row == column) left_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) left_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == left_matrix.columns()-1) left_matrix.cyclic_corner_element() += value; \
        \
        /* Fill matrix row of (i+1,j) */ \
        row = right_index; \
        column = center_index; \
        value = - coeff2 * arr; /* Left */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        row = right_index; \
        column = right_index; \
        value = coeff2 * arr; /* Center: (Left) */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        /* Fill matrix row of (i,j-1) */ \
        row = bottom_index; \
        column = bottom_index; \
        value = coeff3 * att; /* Center: (Top) */ \
        if(row == column) bottom_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) bottom_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == bottom_matrix.columns()-1) bottom_matrix.cyclic_corner_element() += value; \
        \
        /* Fill matrix row of (i,j+1) */ \
        row = top_index; \
        column = top_index; \
        value = coeff4 * att; /* Center: (Bottom) */ \
        if(row == column) top_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) top_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == top_matrix.columns()-1) top_matrix.cyclic_corner_element() += value; \
    } \
    /* ------------------------------------------- */ \
    /* Radial Section: Node next to outer boundary */ \
    /* ------------------------------------------- */ \
    else if(i_r == grid.nr() - 2){ \
        const double h1 = grid.r_dist(i_r-1); \
        const double h2 = grid.r_dist(i_r); \
        const double k1 = grid.theta_dist(i_theta-1); \
        const double k2 = grid.theta_dist(i_theta); \
        const double coeff1 = 0.5*(k1+k2)/h1; \
        const double coeff2 = 0.5*(k1+k2)/h2; \
        const double coeff3 = 0.5*(h1+h2)/k1; \
        const double coeff4 = 0.5*(h1+h2)/k2; \
        \
        const int i_theta_M1 = grid.wrap_theta_index(i_theta-1); \
        const int i_theta_P1 = grid.wrap_theta_index(i_theta+1); \
        \
        auto& bottom_matrix = radial_symmetric_tridiagonal_solver[i_theta_M1]; \
        auto& center_matrix = radial_symmetric_tridiagonal_solver[i_theta]; \
        auto& top_matrix = radial_symmetric_tridiagonal_solver[i_theta_P1]; \
        \
        const int center_index = i_r - numberSmootherCircles; \
        const int left_index = i_r - numberSmootherCircles - 1; \
        const int right_index = i_r - numberSmootherCircles + 1; \
        const int bottom_index = i_r - numberSmootherCircles; \
        const int top_index = i_r - numberSmootherCircles; \
        \
        /* ---------------------------- */ \
        /* Give values to center matrix */ \
        /* ---------------------------- */ \
        /* Fill matrix row of (i,j) */ \
        row = center_index; \
        column = center_index; \
        value = 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        row = center_index; \
        column = left_index; \
        value = - coeff1 * arr; /* Left */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        row = center_index; \
        column = center_index; \
        value = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        /* Fill matrix row of (i-1,j) */ \
        row = left_index; \
        column = center_index; \
        value = - coeff1 * arr; /* Right */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        row = left_index; \
        column = left_index; \
        value = coeff1 * arr; /* Center: (Right) */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        /* Fill matrix row of (i,j-1) */ \
        row = bottom_index; \
        column = bottom_index; \
        value = coeff3 * att; /* Center: (Top) */ \
        if(row == column) bottom_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) bottom_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == bottom_matrix.columns()-1) bottom_matrix.cyclic_corner_element() += value; \
        \
        /* Fill matrix row of (i,j+1) */ \
        row = top_index; \
        column = top_index; \
        value = coeff4 * att; /* Center: (Bottom) */ \
        if(row == column) top_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) top_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == top_matrix.columns()-1) top_matrix.cyclic_corner_element() += value; \
    } \
    /* ------------------------------------------ */ \
    /* Radial Section: Node on the outer boundary */ \
    /* ------------------------------------------ */ \
    else if(i_r == grid.nr() - 1){ \
        double h1 = grid.r_dist(i_r-1); \
        double k1 = grid.theta_dist(i_theta-1); \
        double k2 = grid.theta_dist(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        \
        auto& center_matrix = radial_symmetric_tridiagonal_solver[i_theta]; \
        \
        const int center_index = i_r - numberSmootherCircles; \
        const int left_index = i_r - numberSmootherCircles - 1; \
        \
        /* Fill matrix row of (i,j) */ \
        row = center_index; \
        column = center_index; \
        value = 1.0; \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
        \
        /* Fill matrix row of (i-1,j) */ \
        row = left_index; \
        column = left_index; \
        value = coeff1 * arr; /* Center: (Right) */ \
        if(row == column) center_matrix.main_diagonal(row) += value; \
        else if(row == column - 1) center_matrix.sub_diagonal(row) += value; \
        else if(row == 0 && column == center_matrix.columns()-1) center_matrix.cyclic_corner_element() += value; \
    } \
} while(0)



#define CIRCLE_SECTION_BUILD_SMOOTHER_GIVE(i_r) \
do { \
    const double r = grid_.radius(i_r); \
    const double coeff_alpha = system_parameters_.alpha(r); \
    const double coeff_beta = system_parameters_.beta(r); \
    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++){ \
        const double theta = grid_.theta(i_theta); \
        const double sin_theta = sin_theta_[i_theta]; \
        const double cos_theta = cos_theta_[i_theta]; \
        /* Get arr, att, art, detDF value at the current node */ \
        double arr, att, art, detDF; \
        ARR_ATT_ART(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDF); \
        \
        NODE_BUILD_SMOOTHER_GIVE(i_r, i_theta, grid_, DirBC_Interior_, \
            nonSymmetric_inner_boundary_circle_Asc_matrix, \
            circle_symmetric_cyclic_tridiagonal_solver_, \
            radial_symmetric_tridiagonal_solver_ \
        ); \
    } \
} while(0)



#define RADIAL_SECTION_BUILD_SMOOTHER_GIVE(i_theta) \
do { \
    const double theta = grid_.theta(i_theta); \
    const double sin_theta = sin_theta_[i_theta]; \
    const double cos_theta = cos_theta_[i_theta]; \
    for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++){ \
        const double r = grid_.radius(i_r); \
        const double coeff_alpha = system_parameters_.alpha(r); \
        const double coeff_beta = system_parameters_.beta(r); \
        /* Get arr, att, art, detDF value at the current node */ \
        double arr, att, art, detDF; \
        ARR_ATT_ART(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDF); \
        \
        NODE_BUILD_SMOOTHER_GIVE(i_r, i_theta, grid_, DirBC_Interior_, \
            nonSymmetric_inner_boundary_circle_Asc_matrix, \
            circle_symmetric_cyclic_tridiagonal_solver_, \
            radial_symmetric_tridiagonal_solver_ \
        ); \
    } \
} while(0)



void Smoother::build_Asc_matrices()
{
    omp_set_num_threads(maxOpenMPThreads_);

    /* -------------------------------------- */
    /* Part 1: Allocate Asc Smoother matrices */
    /* -------------------------------------- */

    const int numberSmootherCircles = grid_.numberSmootherCircles();
    const int lengthSmootherRadial = grid_.lengthSmootherRadial();

    SparseMatrix<double> nonSymmetric_inner_boundary_circle_Asc_matrix;

    const int num_circle_nodes = grid_.ntheta();
    circle_symmetric_cyclic_tridiagonal_solver_.resize(numberSmootherCircles-1);

    const int num_radial_nodes = lengthSmootherRadial;
    radial_symmetric_tridiagonal_solver_.resize(grid_.ntheta());

    #pragma omp parallel num_threads(maxOpenMPThreads_) if(grid_.number_of_nodes() > 100'000)
    {
        // ---------------- //
        // Circular Section //
        #pragma omp for nowait
        for (int circle_Asc_index = 0; circle_Asc_index < numberSmootherCircles; circle_Asc_index++){
            
            /* Inner boundary circle */
            if(circle_Asc_index == 0){
                const int nnz = nnz_circle_Asc(circle_Asc_index);
                nonSymmetric_inner_boundary_circle_Asc_matrix = SparseMatrix<double>(num_circle_nodes, num_circle_nodes, nnz);
                for (int i = 0; i < nnz; i++){
                    nonSymmetric_inner_boundary_circle_Asc_matrix.value(i) = 0.0;
                }
            }

            /* Interior Circle Section */
            else{
                const int circle_tridiagonal_solver_index = circle_Asc_index - 1;
                auto& solverMatrix = circle_symmetric_cyclic_tridiagonal_solver_[circle_tridiagonal_solver_index];

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
        for (int radial_Asc_index = 0; radial_Asc_index < grid_.ntheta(); radial_Asc_index++){
            auto& solverMatrix = radial_symmetric_tridiagonal_solver_[radial_Asc_index];

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
                    CIRCLE_SECTION_BUILD_SMOOTHER_GIVE(i_r);
                }
            }
            /* Mod 2 Circles */
            for(int circle_task = 1; circle_task < numCircleTasks; circle_task += 3) {
                #pragma omp task \
                    depend(out: dep[circle_task]) \
                    depend(in: dep[circle_task-1], dep[circle_task+2])   
                {
                    int i_r = grid_.numberSmootherCircles() - circle_task - 1;    
                    CIRCLE_SECTION_BUILD_SMOOTHER_GIVE(i_r);
                }
                
            }
            /* Mod 2 Circles */
            for(int circle_task = 2; circle_task < numCircleTasks; circle_task += 3) {
                #pragma omp task \
                    depend(out: dep[circle_task]) \
                    depend(in: dep[circle_task-1], dep[circle_task+2])   
                {
                    int i_r = grid_.numberSmootherCircles() - circle_task - 1;    
                    CIRCLE_SECTION_BUILD_SMOOTHER_GIVE(i_r);
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
                        RADIAL_SECTION_BUILD_SMOOTHER_GIVE(i_theta);
                    } else{
                        if(additionalRadialTasks == 0){
                            RADIAL_SECTION_BUILD_SMOOTHER_GIVE(0);
                        } 
                        else if(additionalRadialTasks >= 1){
                            RADIAL_SECTION_BUILD_SMOOTHER_GIVE(0);
                            RADIAL_SECTION_BUILD_SMOOTHER_GIVE(1);
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
                        RADIAL_SECTION_BUILD_SMOOTHER_GIVE(i_theta);
                    } else {
                        if(additionalRadialTasks == 0){
                            RADIAL_SECTION_BUILD_SMOOTHER_GIVE(1);
                        } 
                        else if(additionalRadialTasks == 1){
                            RADIAL_SECTION_BUILD_SMOOTHER_GIVE(2);
                        }
                        else if(additionalRadialTasks == 2){
                            RADIAL_SECTION_BUILD_SMOOTHER_GIVE(2);
                            RADIAL_SECTION_BUILD_SMOOTHER_GIVE(3);
                        }
                    }
                }
            }
            /* Mod 2 Radials */
            for(int radial_task = 2; radial_task < numRadialTasks; radial_task += 3) {
                #pragma omp task \
                    depend(out: dep[numCircleTasks+radial_task]) \
                    depend(in: dep[1], dep[numCircleTasks+radial_task-1], dep[numCircleTasks+(radial_task+2) % numRadialTasks])   
                {
                    int i_theta = radial_task + additionalRadialTasks;    
                    RADIAL_SECTION_BUILD_SMOOTHER_GIVE(i_theta);
                }
            }
        }
    }

    omp_set_num_threads(maxOpenMPThreads_);
    delete[] dep;

    /* ----------------------------------------------------------------------------------- */
    /* Part 3: Convert nonSymmetric_inner_boundary_circle_Asc_matrix to a symmetric matrix */
    /* ----------------------------------------------------------------------------------- */

    auto& unsymmetric_Asc_matrix = nonSymmetric_inner_boundary_circle_Asc_matrix;
    auto& symmetric_Asc_matrix = inner_boundary_circle_Asc_matrix_;

    const int nnz = unsymmetric_Asc_matrix.non_zero_size();
    const int numRows = unsymmetric_Asc_matrix.rows();
    const int numColumns = unsymmetric_Asc_matrix.columns();
    const int symmetric_nnz = nnz - (nnz - numRows) / 2;

    symmetric_Asc_matrix = SparseMatrix<double>(numRows, numColumns, symmetric_nnz);
    symmetric_Asc_matrix.is_symmetric(true);

    int current_nz = 0;
    for (int nz_index = 0; nz_index < unsymmetric_Asc_matrix.non_zero_size(); nz_index++) {
        int current_row = unsymmetric_Asc_matrix.row_index(nz_index);
        int current_col = unsymmetric_Asc_matrix.col_index(nz_index);
        if (current_row <= current_col) {
            symmetric_Asc_matrix.row_index(current_nz) = current_row;
            symmetric_Asc_matrix.col_index(current_nz) = current_col;
            symmetric_Asc_matrix.value(current_nz) = std::move(unsymmetric_Asc_matrix.value(nz_index));
            current_nz++;
        }
    }
}