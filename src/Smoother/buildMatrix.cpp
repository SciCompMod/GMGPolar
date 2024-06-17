#include "../../include/Smoother/smoother.h"


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



#define NODE_BUILD_SMOOTHER_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, \
    grid, DirBC_Interior, \
    circle_Asc_matrix, radial_Asc_matrix, nz_index, \
    center_nz_index, left_nz_index, right_nz_index, bottom_nz_index, top_nz_index, \
    center_index, left_index, right_index, bottom_index, top_index, \
    arr, att, art, coeff_beta, detDF) \
do { \
    assert(i_r >= 0 && i_r < grid.nr()); \
    assert(i_theta >= 0 && i_theta < grid.ntheta()); \
    \
    const int numberSmootherCircles = grid.numberSmootherCircles(); \
    const int lengthSmootherRadial = grid.lengthSmootherRadial(); \
    assert(numberSmootherCircles >= 2); \
    assert(lengthSmootherRadial >= 3); \
    /* ------------------------------------------ */ \
    /* Node in the interior of the Circle Section */ \
    /* ------------------------------------------ */ \
    if (i_r > 0 && i_r < numberSmootherCircles) { \
        double h1 = grid.r_dist(i_r-1); \
        double h2 = grid.r_dist(i_r); \
        double k1 = grid.theta_dist(i_theta-1); \
        double k2 = grid.theta_dist(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        \
        auto& left_Asc_matrix = circle_Asc_matrix[i_r-1]; \
        auto& center_Asc_matrix = circle_Asc_matrix[i_r]; \
        auto& right_Asc_matrix = (i_r + 1 == numberSmootherCircles) ? \
            radial_Asc_matrix[i_theta] : \
            circle_Asc_matrix[i_r+1]; \
        \
        center_nz_index = ptr_nz_index_circle_Asc(i_r, i_theta); \
        left_nz_index = ptr_nz_index_circle_Asc(i_r-1, i_theta); \
        right_nz_index = (i_r + 1 == numberSmootherCircles) ? \
            ptr_nz_index_radial_Asc(i_r+1, i_theta) : \
            ptr_nz_index_circle_Asc(i_r+1, i_theta); \
        bottom_nz_index = ptr_nz_index_circle_Asc(i_r, grid.wrap_theta_index(i_theta-1)); \
        top_nz_index = ptr_nz_index_circle_Asc(i_r, grid.wrap_theta_index(i_theta+1)); \
        \
        center_index = i_theta; \
        left_index = i_theta; \
        right_index = (i_r + 1 == numberSmootherCircles) ? 0 : i_theta; \
        bottom_index = grid.wrap_theta_index(i_theta-1); \
        top_index = grid.wrap_theta_index(i_theta+1); \
        \
        /* Fill matrix row of (i,j) */ \
        const Stencil& CenterStencil = get_stencil(i_r); \
        \
        assert(CenterStencil[StencilType::Center] != -1); \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        center_Asc_matrix.row_index(nz_index) = center_index + 1; \
        center_Asc_matrix.col_index(nz_index) = center_index + 1; \
        center_Asc_matrix.value(nz_index) += 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */ \
        \
        assert(CenterStencil[StencilType::Bottom] != -1); \
        nz_index = center_nz_index + CenterStencil[StencilType::Bottom]; \
        center_Asc_matrix.row_index(nz_index) = center_index + 1; \
        center_Asc_matrix.col_index(nz_index) = bottom_index + 1; \
        center_Asc_matrix.value(nz_index) += - coeff3 * att; /* Bottom */ \
        \
        assert(CenterStencil[StencilType::Top] != -1); \
        nz_index = center_nz_index + CenterStencil[StencilType::Top]; \
        center_Asc_matrix.row_index(nz_index) = center_index + 1; \
        center_Asc_matrix.col_index(nz_index) = top_index + 1; \
        center_Asc_matrix.value(nz_index) += - coeff4 * att; /* Top */ \
        \
        assert(CenterStencil[StencilType::Center] != -1); \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        center_Asc_matrix.row_index(nz_index) = center_index + 1; \
        center_Asc_matrix.col_index(nz_index) = center_index + 1; \
        /* Center: (Left, Right, Bottom, Top) */ \
        center_Asc_matrix.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; \
        \
        /* Fill matrix row of (i-1,j) */ \
        const Stencil& LeftStencil = get_stencil(i_r-1); \
        \
        if(!DirBC_Interior || i_r > 1) { \
            assert(LeftStencil[StencilType::Center] != -1); \
            nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
            left_Asc_matrix.row_index(nz_index) = left_index + 1; \
            left_Asc_matrix.col_index(nz_index) = left_index + 1; \
            left_Asc_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) */ \
        } \
        \
        /* Fill matrix row of (i+1,j) */ \
        const Stencil& RightStencil = get_stencil(i_r+1); \
        \
        assert(RightStencil[StencilType::Center] != -1); \
        nz_index = right_nz_index + RightStencil[StencilType::Center]; \
        right_Asc_matrix.row_index(nz_index) = right_index + 1; \
        right_Asc_matrix.col_index(nz_index) = right_index + 1; \
        right_Asc_matrix.value(nz_index) += coeff2 * arr; /* Center: (Left) */ \
        \
        /* Fill matrix row of (i,j-1) */ \
        const Stencil& BottomStencil = CenterStencil; \
        \
        assert(BottomStencil[StencilType::Top] != -1); \
        nz_index = bottom_nz_index + BottomStencil[StencilType::Top]; \
        center_Asc_matrix.row_index(nz_index) = bottom_index + 1; \
        center_Asc_matrix.col_index(nz_index) = center_index + 1; \
        center_Asc_matrix.value(nz_index) += - coeff3 * att; /* Top */ \
        \
        assert(BottomStencil[StencilType::Center] != -1); \
        nz_index = bottom_nz_index + BottomStencil[StencilType::Center]; \
        center_Asc_matrix.row_index(nz_index) = bottom_index + 1; \
        center_Asc_matrix.col_index(nz_index) = bottom_index + 1; \
        center_Asc_matrix.value(nz_index) += coeff3 * att; /* Center: (Top) */ \
        \
        /* Fill matrix row of (i,j+1) */ \
        const Stencil& TopStencil = CenterStencil; \
        \
        assert(TopStencil[StencilType::Bottom] != -1); \
        nz_index = top_nz_index + TopStencil[StencilType::Bottom]; \
        center_Asc_matrix.row_index(nz_index) = top_index + 1; \
        center_Asc_matrix.col_index(nz_index) = center_index + 1; \
        center_Asc_matrix.value(nz_index) += - coeff4 * att; /* Bottom */ \
        \
        assert(TopStencil[StencilType::Center] != -1); \
        nz_index = top_nz_index + TopStencil[StencilType::Center]; \
        center_Asc_matrix.row_index(nz_index) = top_index + 1; \
        center_Asc_matrix.col_index(nz_index) = top_index + 1; \
        center_Asc_matrix.value(nz_index) += coeff4 * att; /* Center: (Bottom) */ \
    } \
    /* ------------------------------------------ */ \
    /* Node in the interior of the Radial Section */ \
    /* ------------------------------------------ */ \
    else if(i_r > numberSmootherCircles && i_r < grid.nr() - 2){ \
        double h1 = grid.r_dist(i_r-1); \
        double h2 = grid.r_dist(i_r); \
        double k1 = grid.theta_dist(i_theta-1); \
        double k2 = grid.theta_dist(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        \
        auto& bottom_Asc_matrix = radial_Asc_matrix[grid.wrap_theta_index(i_theta-1)]; \
        auto& center_Asc_matrix = radial_Asc_matrix[i_theta]; \
        auto& top_Asc_matrix = radial_Asc_matrix[grid.wrap_theta_index(i_theta+1)]; \
        \
        center_nz_index = ptr_nz_index_radial_Asc(i_r, i_theta); \
        left_nz_index = ptr_nz_index_radial_Asc(i_r-1, i_theta); \
        right_nz_index = ptr_nz_index_radial_Asc(i_r+1, i_theta); \
        bottom_nz_index = ptr_nz_index_radial_Asc(i_r, grid.wrap_theta_index(i_theta-1)); \
        top_nz_index = ptr_nz_index_radial_Asc(i_r, grid.wrap_theta_index(i_theta+1)); \
        \
        center_index = i_r - numberSmootherCircles; \
        left_index = i_r - numberSmootherCircles - 1; \
        right_index = i_r - numberSmootherCircles + 1; \
        bottom_index = i_r - numberSmootherCircles; \
        top_index = i_r - numberSmootherCircles; \
        \
        /* Fill matrix row of (i,j) */ \
        const Stencil& CenterStencil = get_stencil(i_r); \
        \
        assert(CenterStencil[StencilType::Center] != -1); \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        center_Asc_matrix.row_index(nz_index) = center_index + 1; \
        center_Asc_matrix.col_index(nz_index) = center_index + 1; \
        center_Asc_matrix.value(nz_index) += 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */ \
        \
        assert(CenterStencil[StencilType::Left] != -1); \
        nz_index = center_nz_index + CenterStencil[StencilType::Left]; \
        center_Asc_matrix.row_index(nz_index) = center_index + 1; \
        center_Asc_matrix.col_index(nz_index) = left_index + 1; \
        center_Asc_matrix.value(nz_index) += - coeff1 * arr; /* Left */ \
        \
        assert(CenterStencil[StencilType::Right] != -1); \
        nz_index = center_nz_index + CenterStencil[StencilType::Right]; \
        center_Asc_matrix.row_index(nz_index) = center_index + 1; \
        center_Asc_matrix.col_index(nz_index) = right_index + 1; \
        center_Asc_matrix.value(nz_index) += - coeff2 * arr; /* Right */ \
        \
        assert(CenterStencil[StencilType::Center] != -1); \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        center_Asc_matrix.row_index(nz_index) = center_index + 1; \
        center_Asc_matrix.col_index(nz_index) = center_index + 1; \
        /* Center: (Left, Right, Bottom, Top) */ \
        center_Asc_matrix.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; \
        \
        /* Fill matrix row of (i-1,j) */ \
        const Stencil& LeftStencil = get_stencil(i_r-1); \
        \
        assert(LeftStencil[StencilType::Right] != -1); \
        nz_index = left_nz_index + LeftStencil[StencilType::Right]; \
        center_Asc_matrix.row_index(nz_index) = left_index + 1; \
        center_Asc_matrix.col_index(nz_index) = center_index + 1; \
        center_Asc_matrix.value(nz_index) += - coeff1 * arr; /* Right */ \
        \
        assert(LeftStencil[StencilType::Center] != -1); \
        nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
        center_Asc_matrix.row_index(nz_index) = left_index + 1; \
        center_Asc_matrix.col_index(nz_index) = left_index + 1; \
        center_Asc_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) */ \
        \
        /* Fill matrix row of (i+1,j) */ \
        const Stencil& RightStencil = get_stencil(i_r+1); \
        \
        assert(RightStencil[StencilType::Left] != -1); \
        nz_index = right_nz_index + RightStencil[StencilType::Left]; \
        center_Asc_matrix.row_index(nz_index) = right_index + 1; \
        center_Asc_matrix.col_index(nz_index) = center_index + 1; \
        center_Asc_matrix.value(nz_index) += - coeff2 * arr; /* Left */ \
        \
        assert(RightStencil[StencilType::Center] != -1); \
        nz_index = right_nz_index + RightStencil[StencilType::Center]; \
        center_Asc_matrix.row_index(nz_index) = right_index + 1; \
        center_Asc_matrix.col_index(nz_index) = right_index + 1; \
        center_Asc_matrix.value(nz_index) += coeff2 * arr; /* Center: (Left) */ \
        \
        /* Fill matrix row of (i,j-1) */ \
        const Stencil& BottomStencil = CenterStencil; \
        \
        assert(BottomStencil[StencilType::Center] != -1); \
        nz_index = bottom_nz_index + BottomStencil[StencilType::Center]; \
        bottom_Asc_matrix.row_index(nz_index) = bottom_index + 1; \
        bottom_Asc_matrix.col_index(nz_index) = bottom_index + 1; \
        bottom_Asc_matrix.value(nz_index) += coeff3 * att; /* Center: (Top) */ \
        \
        /* Fill matrix row of (i,j+1) */ \
        const Stencil& TopStencil = CenterStencil; \
        \
        assert(TopStencil[StencilType::Center] != -1); \
        nz_index = top_nz_index + TopStencil[StencilType::Center]; \
        top_Asc_matrix.row_index(nz_index) = top_index + 1; \
        top_Asc_matrix.col_index(nz_index) = top_index + 1; \
        top_Asc_matrix.value(nz_index) += coeff4 * att; /* Center: (Bottom) */ \
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
            double h2 = grid.r_dist(i_r); \
            double k1 = grid.theta_dist(i_theta-1); \
            double k2 = grid.theta_dist(i_theta); \
            double coeff2 = 0.5*(k1+k2)/h2; \
            \
            auto& center_Asc_matrix = circle_Asc_matrix[i_r]; \
            auto& right_Asc_matrix = circle_Asc_matrix[i_r+1]; \
            \
            center_nz_index = ptr_nz_index_circle_Asc(i_r, i_theta); \
            right_nz_index = ptr_nz_index_circle_Asc(i_r+1, i_theta); \
            \
            center_index = i_theta; \
            right_index = i_theta; \
            bottom_index = grid.wrap_theta_index(i_theta-1); \
            top_index = grid.wrap_theta_index(i_theta+1); \
            \
            /* Fill matrix row of (i,j) */ \
            const Stencil& CenterStencil = get_stencil(i_r); \
            \
            assert(CenterStencil[StencilType::Center] != -1); \
            nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
            center_Asc_matrix.row_index(nz_index) = center_index + 1; \
            center_Asc_matrix.col_index(nz_index) = center_index + 1; \
            center_Asc_matrix.value(nz_index) += 1.0; \
            \
            /* Give values to the interior nodes! */ \
            /* Fill matrix row of (i+1,j) */ \
            const Stencil& RightStencil = get_stencil(i_r+1); \
            \
            assert(RightStencil[StencilType::Center] != -1); \
            nz_index = right_nz_index + RightStencil[StencilType::Center]; \
            right_Asc_matrix.row_index(nz_index) = right_index + 1; \
            right_Asc_matrix.col_index(nz_index) = right_index + 1; \
            right_Asc_matrix.value(nz_index) += coeff2 * arr; /* Center: (Left) */ \
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
            auto& center_Asc_matrix = circle_Asc_matrix[i_r]; \
            auto& right_Asc_matrix = circle_Asc_matrix[i_r+1]; \
            auto& left_Asc_matrix = center_Asc_matrix; \
            \
            center_nz_index = ptr_nz_index_circle_Asc(i_r, i_theta); \
            left_nz_index = ptr_nz_index_circle_Asc(i_r, grid.wrap_theta_index(i_theta + (grid.ntheta()>>1))); \
            right_nz_index = ptr_nz_index_circle_Asc(i_r+1, i_theta); \
            bottom_nz_index = ptr_nz_index_circle_Asc(i_r, grid.wrap_theta_index(i_theta-1)); \
            top_nz_index = ptr_nz_index_circle_Asc(i_r, grid.wrap_theta_index(i_theta+1)); \
            \
            center_index = i_theta; \
            left_index = grid.wrap_theta_index(i_theta + (grid.ntheta()>>1)); \
            right_index = i_theta; \
            bottom_index = grid.wrap_theta_index(i_theta-1); \
            top_index = grid.wrap_theta_index(i_theta+1); \
            \
            /* Fill matrix row of (i,j) */ \
            const Stencil& CenterStencil = get_stencil(i_r); \
            \
            assert(CenterStencil[StencilType::Center] != -1); \
            nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
            center_Asc_matrix.row_index(nz_index) = center_index + 1; \
            center_Asc_matrix.col_index(nz_index) = center_index + 1; \
            center_Asc_matrix.value(nz_index) += 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */ \
            \
            assert(CenterStencil[StencilType::Left] != -1); \
            nz_index = center_nz_index + CenterStencil[StencilType::Left]; \
            center_Asc_matrix.row_index(nz_index) = center_index + 1; \
            center_Asc_matrix.col_index(nz_index) = left_index + 1; \
            center_Asc_matrix.value(nz_index) += - coeff1 * arr; /* Left */ \
            \
            assert(CenterStencil[StencilType::Bottom] != -1); \
            nz_index = center_nz_index + CenterStencil[StencilType::Bottom]; \
            center_Asc_matrix.row_index(nz_index) = center_index + 1; \
            center_Asc_matrix.col_index(nz_index) = bottom_index + 1; \
            center_Asc_matrix.value(nz_index) += - coeff3 * att; /* Bottom */ \
            \
            assert(CenterStencil[StencilType::Top] != -1); \
            nz_index = center_nz_index + CenterStencil[StencilType::Top]; \
            center_Asc_matrix.row_index(nz_index) = center_index + 1; \
            center_Asc_matrix.col_index(nz_index) = top_index + 1; \
            center_Asc_matrix.value(nz_index) += - coeff4 * att; /* Top */ \
            \
            assert(CenterStencil[StencilType::Center] != -1); \
            nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
            center_Asc_matrix.row_index(nz_index) = center_index + 1; \
            center_Asc_matrix.col_index(nz_index) = center_index + 1; \
            /* Center: (Left, Right, Bottom, Top) */ \
            center_Asc_matrix.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; \
            \
            /* Fill matrix row of (i-1,j) */ \
            /* From view the view of the across origin node, the directions are roatated by 180 degrees in the stencil! */ \
            const Stencil& LeftStencil = CenterStencil; \
            \
            assert(LeftStencil[StencilType::Left] != -1); \
            nz_index = left_nz_index + LeftStencil[StencilType::Left]; \
            left_Asc_matrix.row_index(nz_index) = left_index + 1; \
            left_Asc_matrix.col_index(nz_index) = center_index + 1; \
            left_Asc_matrix.value(nz_index) += - coeff1 * arr; /* Right -> Left*/ \
            \
            assert(LeftStencil[StencilType::Center] != -1); \
            nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
            left_Asc_matrix.row_index(nz_index) = left_index + 1; \
            left_Asc_matrix.col_index(nz_index) = left_index + 1; \
            left_Asc_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) -> Center: (Left) */ \
            \
            assert(LeftStencil[StencilType::BottomLeft] == -1); \
            /* Top Right -> Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = left_nz_index + LeftStencil[StencilType::BottomLeft]; */  \
            /* left_Asc_matrix.row_index(nz_index) = left_index; */  \
            /* left_Asc_matrix.col_index(nz_index) = top_index; */  \
            /* left_Asc_matrix.value(nz_index) += - 0.25 * art; // Top Right -> Bottom Left*/ \
            \
            assert(LeftStencil[StencilType::TopLeft] == -1); \
            /* Bottom Right -> Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = left_nz_index + LeftStencil[StencilType::TopLeft]; */ \
            /* left_Asc_matrix.row_index(nz_index) = left_index; */  \
            /* left_Asc_matrix.col_index(nz_index) = bottom_index; */  \
            /* left_Asc_matrix.value(nz_index) += 0.25 * art; // Bottom Right -> Top Left */ \
            \
            /* Fill matrix row of (i+1,j) */ \
            const Stencil& RightStencil = get_stencil(i_r+1); \
            \
            assert(RightStencil[StencilType::Center] != -1); \
            nz_index = right_nz_index + RightStencil[StencilType::Center]; \
            right_Asc_matrix.row_index(nz_index) = right_index + 1; \
            right_Asc_matrix.col_index(nz_index) = right_index + 1; \
            right_Asc_matrix.value(nz_index) += coeff2 * arr; /* Center: (Left) */ \
            \
            /* Fill matrix row of (i,j-1) */ \
            const Stencil& BottomStencil = CenterStencil; \
            \
            assert(BottomStencil[StencilType::Top] != -1); \
            nz_index = bottom_nz_index + BottomStencil[StencilType::Top]; \
            center_Asc_matrix.row_index(nz_index) = bottom_index + 1; \
            center_Asc_matrix.col_index(nz_index) = center_index + 1; \
            center_Asc_matrix.value(nz_index) += - coeff3 * att; /* Top */ \
            \
            assert(BottomStencil[StencilType::Center] != -1); \
            nz_index = bottom_nz_index + BottomStencil[StencilType::Center]; \
            center_Asc_matrix.row_index(nz_index) = bottom_index + 1; \
            center_Asc_matrix.col_index(nz_index) = bottom_index + 1; \
            center_Asc_matrix.value(nz_index) += coeff3 * att; /* Center: (Top) */ \
            \
            assert(BottomStencil[StencilType::TopLeft] == -1); \
            /* REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = bottom_nz_index + BottomStencil[StencilType::TopLeft]; */ \
            /* center_Asc_matrix.row_index(nz_index) = bottom_index + 1; */ \
            /* center_Asc_matrix.col_index(nz_index) = left_index + 1; */ \
            /* center_Asc_matrix.value(nz_index) += 0.25 * art; // Top Left */ \
            \
            /* Fill matrix row of (i,j+1) */ \
            const Stencil& TopStencil = CenterStencil; \
            \
            assert(TopStencil[StencilType::Bottom] != -1); \
            nz_index = top_nz_index + TopStencil[StencilType::Bottom]; \
            center_Asc_matrix.row_index(nz_index) = top_index + 1; \
            center_Asc_matrix.col_index(nz_index) = center_index + 1; \
            center_Asc_matrix.value(nz_index) += - coeff4 * att; /* Bottom */ \
            \
            assert(TopStencil[StencilType::Center] != -1); \
            nz_index = top_nz_index + TopStencil[StencilType::Center]; \
            center_Asc_matrix.row_index(nz_index) = top_index + 1; \
            center_Asc_matrix.col_index(nz_index) = top_index + 1; \
            center_Asc_matrix.value(nz_index) += coeff4 * att; /* Center: (Bottom) */ \
            \
            assert(TopStencil[StencilType::BottomLeft] == -1); \
            /* REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = top_nz_index + TopStencil[StencilType::BottomLeft]; */ \
            /* center_Asc_matrix.row_index(nz_index) = top_index + 1; */ \
            /* center_Asc_matrix.col_index(nz_index) = left_index + 1; */ \
            /* center_Asc_matrix.value(nz_index) += - 0.25 * art; // Bottom Left */ \
            \
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
        auto& bottom_Asc_matrix = radial_Asc_matrix[grid.wrap_theta_index(i_theta-1)]; \
        auto& center_Asc_matrix = radial_Asc_matrix[i_theta]; \
        auto& top_Asc_matrix = radial_Asc_matrix[grid.wrap_theta_index(i_theta+1)]; \
        auto& left_Asc_matrix = circle_Asc_matrix[i_r-1]; \
        \
        center_nz_index = ptr_nz_index_radial_Asc(i_r, i_theta); \
        left_nz_index = ptr_nz_index_circle_Asc(i_r-1, i_theta); \
        right_nz_index = ptr_nz_index_radial_Asc(i_r+1, i_theta); \
        bottom_nz_index = ptr_nz_index_radial_Asc(i_r, grid.wrap_theta_index(i_theta-1)); \
        top_nz_index = ptr_nz_index_radial_Asc(i_r, grid.wrap_theta_index(i_theta+1)); \
        \
        center_index = i_r - numberSmootherCircles; \
        left_index = i_theta; \
        right_index = i_r - numberSmootherCircles + 1; \
        bottom_index = i_r - numberSmootherCircles; \
        top_index = i_r - numberSmootherCircles; \
        \
        /* Fill matrix row of (i,j) */ \
        const Stencil& CenterStencil = get_stencil(i_r); \
        \
        assert(CenterStencil[StencilType::Center] != -1); \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        center_Asc_matrix.row_index(nz_index) = center_index + 1; \
        center_Asc_matrix.col_index(nz_index) = center_index + 1; \
        center_Asc_matrix.value(nz_index) += 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */ \
        \
        assert(CenterStencil[StencilType::Right] != -1); \
        nz_index = center_nz_index + CenterStencil[StencilType::Right]; \
        center_Asc_matrix.row_index(nz_index) = center_index + 1; \
        center_Asc_matrix.col_index(nz_index) = right_index + 1; \
        center_Asc_matrix.value(nz_index) += - coeff2 * arr; /* Right */ \
        \
        assert(CenterStencil[StencilType::Center] != -1); \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        center_Asc_matrix.row_index(nz_index) = center_index + 1; \
        center_Asc_matrix.col_index(nz_index) = center_index + 1; \
        /* Center: (Left, Right, Bottom, Top) */ \
        center_Asc_matrix.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; \
        \
        /* Fill matrix row of (i-1,j) */ \
        const Stencil& LeftStencil = get_stencil(i_r-1); \
        \
        assert(LeftStencil[StencilType::Center] != -1); \
        nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
        left_Asc_matrix.row_index(nz_index) = left_index + 1; \
        left_Asc_matrix.col_index(nz_index) = left_index + 1; \
        left_Asc_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) */ \
        \
        /* Fill matrix row of (i+1,j) */ \
        const Stencil& RightStencil = get_stencil(i_r+1); \
        \
        assert(RightStencil[StencilType::Left] != -1); \
        nz_index = right_nz_index + RightStencil[StencilType::Left]; \
        center_Asc_matrix.row_index(nz_index) = right_index + 1; \
        center_Asc_matrix.col_index(nz_index) = center_index + 1; \
        center_Asc_matrix.value(nz_index) += - coeff2 * arr; /* Left */ \
        \
        assert(RightStencil[StencilType::Center] != -1); \
        nz_index = right_nz_index + RightStencil[StencilType::Center]; \
        center_Asc_matrix.row_index(nz_index) = right_index + 1; \
        center_Asc_matrix.col_index(nz_index) = right_index + 1; \
        center_Asc_matrix.value(nz_index) += coeff2 * arr; /* Center: (Left) */ \
        \
        /* Fill matrix row of (i,j-1) */ \
        const Stencil& BottomStencil = CenterStencil; \
        \
        assert(BottomStencil[StencilType::Center] != -1); \
        nz_index = bottom_nz_index + BottomStencil[StencilType::Center]; \
        bottom_Asc_matrix.row_index(nz_index) = bottom_index + 1; \
        bottom_Asc_matrix.col_index(nz_index) = bottom_index + 1; \
        bottom_Asc_matrix.value(nz_index) += coeff3 * att; /* Center: (Top) */ \
        \
        /* Fill matrix row of (i,j+1) */ \
        const Stencil& TopStencil = CenterStencil; \
        \
        assert(TopStencil[StencilType::Center] != -1); \
        nz_index = top_nz_index + TopStencil[StencilType::Center]; \
        top_Asc_matrix.row_index(nz_index) = top_index + 1; \
        top_Asc_matrix.col_index(nz_index) = top_index + 1; \
        top_Asc_matrix.value(nz_index) += coeff4 * att; /* Center: (Bottom) */ \
    } \
    /* ------------------------------------------- */ \
    /* Radial Section: Node next to outer boundary */ \
    /* ------------------------------------------- */ \
    else if(i_r == grid.nr() - 2){ \
        double h1 = grid.r_dist(i_r-1); \
        double h2 = grid.r_dist(i_r); \
        double k1 = grid.theta_dist(i_theta-1); \
        double k2 = grid.theta_dist(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        \
        auto& bottom_Asc_matrix = radial_Asc_matrix[grid.wrap_theta_index(i_theta-1)]; \
        auto& center_Asc_matrix = radial_Asc_matrix[i_theta]; \
        auto& top_Asc_matrix = radial_Asc_matrix[grid.wrap_theta_index(i_theta+1)]; \
        \
        center_nz_index = ptr_nz_index_radial_Asc(i_r, i_theta); \
        left_nz_index = ptr_nz_index_radial_Asc(i_r-1, i_theta); \
        right_nz_index = ptr_nz_index_radial_Asc(i_r+1, i_theta); \
        bottom_nz_index = ptr_nz_index_radial_Asc(i_r, grid.wrap_theta_index(i_theta-1)); \
        top_nz_index = ptr_nz_index_radial_Asc(i_r, grid.wrap_theta_index(i_theta+1)); \
        \
        center_index = i_r - numberSmootherCircles; \
        left_index = i_r - numberSmootherCircles - 1; \
        right_index = i_r - numberSmootherCircles + 1; \
        bottom_index = i_r - numberSmootherCircles; \
        top_index = i_r - numberSmootherCircles; \
        \
        /* Fill matrix row of (i,j) */ \
        const Stencil& CenterStencil = get_stencil(i_r); \
        \
        assert(CenterStencil[StencilType::Center] != -1); \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        center_Asc_matrix.row_index(nz_index) = center_index + 1; \
        center_Asc_matrix.col_index(nz_index) = center_index + 1; \
        center_Asc_matrix.value(nz_index) += 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */ \
        \
        assert(CenterStencil[StencilType::Left] != -1); \
        nz_index = center_nz_index + CenterStencil[StencilType::Left]; \
        center_Asc_matrix.row_index(nz_index) = center_index + 1; \
        center_Asc_matrix.col_index(nz_index) = left_index + 1; \
        center_Asc_matrix.value(nz_index) += - coeff1 * arr; /* Left */ \
        \
        assert(CenterStencil[StencilType::Right] == -1); \
        /* Removed to make matrix symmetric */ \
        /* nz_index = center_nz_index + CenterStencil[StencilType::Right]; */ \
        /* center_Asc_matrix.row_index(nz_index) = center_index + 1; */ \
        /* center_Asc_matrix.col_index(nz_index) = right_index + 1; */ \
        /* center_Asc_matrix.value(nz_index) += - coeff2 * arr; // Right */ \
        \
        assert(CenterStencil[StencilType::Center] != -1); \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        center_Asc_matrix.row_index(nz_index) = center_index + 1; \
        center_Asc_matrix.col_index(nz_index) = center_index + 1; \
        /* Center: (Left, Right, Bottom, Top) */ \
        center_Asc_matrix.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; \
        \
        /* Fill matrix row of (i-1,j) */ \
        const Stencil& LeftStencil = get_stencil(i_r-1); \
        \
        assert(LeftStencil[StencilType::Right] != -1); \
        nz_index = left_nz_index + LeftStencil[StencilType::Right]; \
        center_Asc_matrix.row_index(nz_index) = left_index + 1; \
        center_Asc_matrix.col_index(nz_index) = center_index + 1; \
        center_Asc_matrix.value(nz_index) += - coeff1 * arr; /* Right */ \
        \
        assert(LeftStencil[StencilType::Center] != -1); \
        nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
        center_Asc_matrix.row_index(nz_index) = left_index + 1; \
        center_Asc_matrix.col_index(nz_index) = left_index + 1; \
        center_Asc_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) */ \
        \
        /* Fill matrix row of (i,j-1) */ \
        const Stencil& BottomStencil = CenterStencil; \
        \
        assert(BottomStencil[StencilType::Center] != -1); \
        nz_index = bottom_nz_index + BottomStencil[StencilType::Center]; \
        bottom_Asc_matrix.row_index(nz_index) = bottom_index + 1; \
        bottom_Asc_matrix.col_index(nz_index) = bottom_index + 1; \
        bottom_Asc_matrix.value(nz_index) += coeff3 * att; /* Center: (Top) */ \
        \
        /* Fill matrix row of (i,j+1) */ \
        const Stencil& TopStencil = CenterStencil; \
        \
        assert(TopStencil[StencilType::Center] != -1); \
        nz_index = top_nz_index + TopStencil[StencilType::Center]; \
        bottom_Asc_matrix.row_index(nz_index) = top_index + 1; \
        bottom_Asc_matrix.col_index(nz_index) = top_index + 1; \
        bottom_Asc_matrix.value(nz_index) += coeff4 * att; /* Center: (Bottom) */ \
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
        auto& center_Asc_matrix = radial_Asc_matrix[i_theta]; \
        \
        center_nz_index = ptr_nz_index_radial_Asc(i_r, i_theta); \
        left_nz_index = ptr_nz_index_radial_Asc(i_r-1, i_theta); \
        \
        center_index = i_r - numberSmootherCircles; \
        left_index = i_r - numberSmootherCircles - 1; \
        \
        /* Fill matrix row of (i,j) */ \
        const Stencil& CenterStencil = get_stencil(i_r); \
        \
        assert(CenterStencil[StencilType::Center] != -1); \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        center_Asc_matrix.row_index(nz_index) = center_index + 1; \
        center_Asc_matrix.col_index(nz_index) = center_index + 1; \
        center_Asc_matrix.value(nz_index) += 1.0; \
        \
        /* Give value to the interior nodes! */ \
        /* Fill matrix row of (i-1,j) */ \
        const Stencil& LeftStencil = get_stencil(i_r-1); \
        \
        assert(LeftStencil[StencilType::Center] != -1); \
        nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
        center_Asc_matrix.row_index(nz_index) = left_index + 1; \
        center_Asc_matrix.col_index(nz_index) = left_index + 1; \
        center_Asc_matrix.value(nz_index) += coeff1 * arr; /* Center: (Right) */ \
    } \
} while(0) \



#define CIRCLE_SECTION_BUILD_SMOOTHER_GIVE(i_r) \
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
        NODE_BUILD_SMOOTHER_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, \
            grid_, DirBC_Interior_, \
            circle_Asc_matrix, radial_Asc_matrix, nz_index, \
            center_nz_index, left_nz_index, right_nz_index, bottom_nz_index, top_nz_index, \
            center_index, left_index, right_index, bottom_index, top_index, \
            arr, att, art, coeff_beta, detDF); \
    } \
} while(0)



#define RADIAL_SECTION_BUILD_SMOOTHER_GIVE(i_theta) \
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
        NODE_BUILD_SMOOTHER_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, \
            grid_, DirBC_Interior_, \
            circle_Asc_matrix, radial_Asc_matrix, nz_index, \
            center_nz_index, left_nz_index, right_nz_index, bottom_nz_index, top_nz_index, \
            center_index, left_index, right_index, bottom_index, top_index, \
            arr, att, art, coeff_beta, detDF); \
    } \
} while(0)



void Smoother::build_Asc_matrices(
    std::vector<SparseMatrix<double>>& symmetric_circle_Asc_matrix, 
    std::vector<SparseMatrix<double>>& symmetric_radial_Asc_matrix)
{
    omp_set_num_threads(maxOpenMPThreads_);

    std::vector<SparseMatrix<double>> circle_Asc_matrix(grid_.numberSmootherCircles());
    const int circle_n = grid_.ntheta();
    #pragma omp parallel for
    for (int circle_Asc_index = 0; circle_Asc_index < grid_.numberSmootherCircles(); circle_Asc_index++){
        const int circle_Asc_nnz = nnz_circle_Asc(circle_Asc_index);
        circle_Asc_matrix[circle_Asc_index] = SparseMatrix<double>(circle_n, circle_n, circle_Asc_nnz);
        for (int i = 0; i < circle_Asc_nnz; i++){
            circle_Asc_matrix[circle_Asc_index].value(i) = 0.0;
        }
    }
    
    std::vector<SparseMatrix<double>> radial_Asc_matrix(grid_.ntheta());
    const int radial_n = grid_.lengthSmootherRadial();
    #pragma omp parallel for
    for (int radial_Asc_index = 0; radial_Asc_index < grid_.ntheta(); radial_Asc_index++){
        const int radial_Asc_nnz = nnz_radial_Asc(radial_Asc_index);
        radial_Asc_matrix[radial_Asc_index] = SparseMatrix<double>(radial_n, radial_n, radial_Asc_nnz);
        for (int i = 0; i < radial_Asc_nnz; i++){
            radial_Asc_matrix[radial_Asc_index].value(i) = 0.0;
        }
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

    symmetric_circle_Asc_matrix = std::move(circle_Asc_matrix);
    symmetric_radial_Asc_matrix = std::move(radial_Asc_matrix);

    // symmetric_circle_Asc_matrix.resize(grid_.numberSmootherCircles());
    // #pragma omp parallel for
    // for (int circle_Asc_index = 0; circle_Asc_index < grid_.numberSmootherCircles(); circle_Asc_index++){
    //     SparseMatrix<double>& Asc_matrix = circle_Asc_matrix[circle_Asc_index];
    //     const int circle_Asc_matrix_nnz = Asc_matrix.non_zero_size();
    //     const int symmetric_circle_Asc_matrix_nnz = circle_Asc_matrix_nnz - (circle_Asc_matrix_nnz - circle_n) / 2;
    //     symmetric_circle_Asc_matrix[circle_Asc_index] = SparseMatrix<double> (Asc_matrix.rows(), Asc_matrix.columns(), symmetric_circle_Asc_matrix_nnz);
    //     SparseMatrix<double>& symmetric_Asc_matrix = symmetric_circle_Asc_matrix[circle_Asc_index];
    //     symmetric_Asc_matrix.is_symmetric(true);
    //     int current_nz = 0;
    //     for (int nz_index = 0; nz_index < Asc_matrix.non_zero_size(); nz_index++) {
    //         int current_row = Asc_matrix.row_index(nz_index);
    //         int current_col = Asc_matrix.col_index(nz_index);
    //         if (current_row <= current_col) {
    //             symmetric_Asc_matrix.row_index(current_nz) = current_row;
    //             symmetric_Asc_matrix.col_index(current_nz) = current_col;
    //             symmetric_Asc_matrix.value(current_nz) = std::move(Asc_matrix.value(nz_index));
    //             current_nz++;
    //         }
    //     }
    // }

    // symmetric_radial_Asc_matrix.resize(grid_.ntheta());
    // #pragma omp parallel for
    // for (int radial_Asc_index = 0; radial_Asc_index < grid_.ntheta(); radial_Asc_index++){
    //     SparseMatrix<double>& Asc_matrix = radial_Asc_matrix[radial_Asc_index];
    //     const int radial_Asc_matrix_nnz = Asc_matrix.non_zero_size();
    //     const int symmetric_radial_Asc_matrix_nnz = radial_Asc_matrix_nnz - (radial_Asc_matrix_nnz - radial_n) / 2;
    //     symmetric_radial_Asc_matrix[radial_Asc_index] = SparseMatrix<double> (radial_n, radial_n, symmetric_radial_Asc_matrix_nnz);
    //     SparseMatrix<double>& symmetric_Asc_matrix = symmetric_radial_Asc_matrix[radial_Asc_index];
    //     symmetric_Asc_matrix.is_symmetric(true);
    //     int current_nz = 0;
    //     for (int nz_index = 0; nz_index < Asc_matrix.non_zero_size(); nz_index++) {
    //         int current_row = Asc_matrix.row_index(nz_index);
    //         int current_col = Asc_matrix.col_index(nz_index);
    //         if (current_row <= current_col) {
    //             symmetric_Asc_matrix.row_index(current_nz) = current_row;
    //             symmetric_Asc_matrix.col_index(current_nz) = current_col;
    //             symmetric_Asc_matrix.value(current_nz) = std::move(Asc_matrix.value(nz_index));
    //             current_nz++;
    //         }
    //     }
    // }    
}