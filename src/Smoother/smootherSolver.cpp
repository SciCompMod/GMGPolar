#include "../../include/Smoother/smoother.h"

enum class SmootherColor {
    Black = 0,
    White = 1,
};

const char* toString(SmootherColor color) {
    switch (color) {
        case SmootherColor::Black: return "Black";
        case SmootherColor::White: return "White";
        default: return "Unknown";
    }
}

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



#define NODE_APPLY_ASC_ORTHO_CIRCLE_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, \
    system_parameters, grid, DirBC_Interior, color, \
    x, result, factor, \
    arr, att, art, coeff_beta, detDF) \
do { \
    assert(i_r >= 0 && i_r <= grid_.numberSmootherCircles()); \
    bool isOddNumberSmootherCircles = (grid.numberSmootherCircles() & 1); \
    bool isOddRadialIndex = (i_r & 1); \
    SmootherColor node_color = (isOddNumberSmootherCircles == isOddRadialIndex) ? SmootherColor::White : SmootherColor::Black; \
    /* -------------------- */ \
    /* Node in the interior */ \
    /* -------------------- */ \
    if(i_r > 0 && i_r < grid.numberSmootherCircles()){ \
        double h1 = grid.r_dist(i_r-1); \
        double h2 = grid.r_dist(i_r); \
        double k1 = grid.theta_dist(i_theta-1); \
        double k2 = grid.theta_dist(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        /* -------------------- */ \
        /* Inside Section Parts */ \
        if(node_color == color){ \
            /* Fill result of (i,j) */ \
            result[grid.index(i_r,i_theta)] += \
                0.25 * (h1+h2)*(k1+k2) * system_parameters.rhs_f(r, theta, sin_theta, cos_theta) * fabs(detDF); \
            /* Fill result(i,j) */ \
            assert(get_stencil(i_r)[StencilType::Left] == -1); \
            assert(get_stencil(i_r)[StencilType::Right] == -1); \
            result[grid.index(i_r,i_theta)] += factor * ( \
                - coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Left */ \
                - coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Right */ \
            ); \
            /* Fill result(i,j-1) */ \
            assert(get_stencil(i_r)[StencilType::TopRight] == -1); \
            assert(get_stencil(i_r)[StencilType::TopLeft] == -1); \
            result[grid.index(i_r,i_theta-1)] += factor * ( \
                - 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Top Right */ \
                + 0.25 * art * x[grid.index(i_r-1,i_theta)] ); /* Top Left */ \
            /* Fill result(i,j+1) */ \
            assert(get_stencil(i_r)[StencilType::BottomRight] == -1); \
            assert(get_stencil(i_r)[StencilType::BottomLeft] == -1); \
            result[grid.index(i_r,i_theta+1)] += factor * ( \
                + 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Bottom Right */ \
                - 0.25 * art * x[grid.index(i_r-1,i_theta)] ); /* Bottom Left */ \
        } \
        /* --------------------- */ \
        /* Outside Section Parts */ \
        else if(node_color != color){ \
            /* Fill result(i-1,j) */ \
            if(!DirBC_Interior || i_r > 1) { \
                assert(get_stencil(i_r-1)[StencilType::Right] == -1); \
                assert(get_stencil(i_r-1)[StencilType::TopRight] == -1); \
                assert(get_stencil(i_r-1)[StencilType::BottomRight] == -1); \
                result[grid.index(i_r-1,i_theta)] += factor * ( \
                    - coeff1 * arr * x[grid.index(i_r,i_theta)] /* Right */ \
                    - 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Right */ \
                    + 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Right */ \
            } \
            /* Fill result(i+1,j) */ \
            if(i_r < grid.numberSmootherCircles() - 1) { \
                assert(get_stencil(i_r+1)[StencilType::Left] == -1); \
                assert(get_stencil(i_r+1)[StencilType::TopLeft] == -1); \
                assert(get_stencil(i_r+1)[StencilType::BottomLeft] == -1); \
                result[grid.index(i_r+1,i_theta)] += factor * ( \
                    - coeff2 * arr * x[grid.index(i_r,i_theta)] /* Left */ \
                    + 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Left */ \
                    - 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Left */ \
            } \
        } \
    } \
    /* -------------------- */ \
    /* Node on the boundary */ \
    /* -------------------- */ \
    else if(i_r == 0){ \
        /* ------------------------------------------------ */ \
        /* Case 1: Dirichlet boundary on the inner boundary */ \
        /* ------------------------------------------------ */ \
        if(DirBC_Interior){ \
            double h2 = grid.r_dist(i_r); \
            double k1 = grid.theta_dist(i_theta-1); \
            double k2 = grid.theta_dist(i_theta); \
            double coeff2 = 0.5*(k1+k2)/h2; \
            /* -------------------- */ \
            /* Inside Section Parts */ \
            if(node_color == color){ \
                /* Fill result of (i,j) */ \
                result[grid.index(i_r,i_theta)] += \
                    system_parameters.u_D_Interior(r, theta, sin_theta, cos_theta); \
            } \
            /* --------------------- */ \
            /* Outside Section Parts */ \
            else if(node_color != color){ \
                /* Fill result(i+1,j) */ \
                assert(get_stencil(i_r+1)[StencilType::Left] == -1); \
                assert(get_stencil(i_r+1)[StencilType::TopLeft] == -1); \
                assert(get_stencil(i_r+1)[StencilType::BottomLeft] == -1); \
                result[grid.index(i_r+1,i_theta)] += factor * ( \
                    - coeff2 * arr * x[grid.index(i_r,i_theta)] /* Left */ \
                    + 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Left */ \
                    - 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Left */ \
            } \
        } \
        else{ \
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
            /* -------------------- */ \
            /* Inside Section Parts */ \
            if(node_color == color){ \
                /* Fill result of (i,j) */ \
                result[grid.index(i_r,i_theta)] += \
                    0.25 * (h1+h2)*(k1+k2) * system_parameters.rhs_f(r, theta, sin_theta, cos_theta) * fabs(detDF); \
                /* Fill result(i,j) */ \
                assert(get_stencil(i_r)[StencilType::Right] == -1); \
                result[grid.index(i_r,i_theta)] += factor * ( \
                    /* - coeff1 * arr * x[grid.index(i_r, i_theta + (grid.ntheta()>>1))] // Left: Not in Asc_ortho */ \
                    - coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Right */ \
                ); \
                /* Fill result(i,j-1) */ \
                assert(get_stencil(i_r)[StencilType::TopRight] == -1); \
                assert(get_stencil(i_r)[StencilType::TopLeft] == -1); \
                result[grid.index(i_r,i_theta-1)] += factor * ( \
                    - 0.25 * art * x[grid.index(i_r+1,i_theta)] ); /* Top Right */ \
                /*  + 0.25 * art * x[grid.index(i_r-1,i_theta)]; // Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
                /* Fill result(i,j+1) */ \
                assert(get_stencil(i_r)[StencilType::BottomRight] == -1); \
                assert(get_stencil(i_r)[StencilType::BottomLeft] == -1); \
                result[grid.index(i_r,i_theta+1)] += factor * ( \
                    + 0.25 * art * x[grid.index(i_r+1,i_theta)] ); /* Bottom Right */ \
                /*  - 0.25 * art * x[grid.index(i_r-1,i_theta)]; // Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            } \
            /* --------------------- */ \
            /* Outside Section Parts */ \
            else if(node_color != color){ \
                /* Fill result(i+1,j) */ \
                assert(get_stencil(i_r+1)[StencilType::Left] == -1); \
                assert(get_stencil(i_r+1)[StencilType::TopLeft] == -1); \
                assert(get_stencil(i_r+1)[StencilType::BottomLeft] == -1); \
                result[grid.index(i_r+1,i_theta)] += factor * ( \
                    - coeff2 * arr * x[grid.index(i_r,i_theta)] /* Left */ \
                    + 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Left */ \
                    - 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Left */ \
            } \
        } \
    } \
    /* ----------------------------- */ \
    /* Node next to circular section */ \
    /* ----------------------------- */ \
    else if(i_r == grid.numberSmootherCircles()) { \
        assert(node_color == SmootherColor::White); \
        if(color == SmootherColor::Black){ \
            double h1 = grid.r_dist(i_r-1); \
            double h2 = grid.r_dist(i_r); \
            double k1 = grid.theta_dist(i_theta-1); \
            double k2 = grid.theta_dist(i_theta); \
            double coeff1 = 0.5*(k1+k2)/h1; \
            double coeff2 = 0.5*(k1+k2)/h2; \
            double coeff3 = 0.5*(h1+h2)/k1; \
            double coeff4 = 0.5*(h1+h2)/k2; \
            /* --------------------- */ \
            /* Outside Section Parts */ \
            /* Fill result(i-1,j) */ \
            assert(get_stencil(i_r-1)[StencilType::Right] == -1); \
            assert(get_stencil(i_r-1)[StencilType::TopRight] == -1); \
            assert(get_stencil(i_r-1)[StencilType::BottomRight] == -1); \
            result[grid.index(i_r-1,i_theta)] += factor * ( \
                - coeff1 * arr * x[grid.index(i_r,i_theta)] /* Right */ \
                - 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Right */ \
                + 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Right */ \
        } \
    } \
} while(0) \






#define NODE_APPLY_ASC_ORTHO_RADIAL_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, \
    system_parameters, grid, DirBC_Interior, color, \
    x, result, factor, \
    arr, att, art, coeff_beta, detDF) \
do { \
    assert(i_r >= grid.numberSmootherCircles()-1 && i_r < grid.nr()); \
    SmootherColor node_color = (i_theta & 1) ? SmootherColor::White : SmootherColor::Black; \
    /* -------------------- */ \
    /* Node in the interior */ \
    /* -------------------- */ \
    if(i_r > grid.numberSmootherCircles() && i_r < grid.nr()-2){ \
        double h1 = grid.r_dist(i_r-1); \
        double h2 = grid.r_dist(i_r); \
        double k1 = grid.theta_dist(i_theta-1); \
        double k2 = grid.theta_dist(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        /* -------------------- */ \
        /* Inside Section Parts */ \
        if(node_color == color){ \
            /* Fill result of (i,j) */ \
            result[grid.index(i_r,i_theta)] += \
                0.25 * (h1+h2)*(k1+k2) * system_parameters.rhs_f(r, theta, sin_theta, cos_theta) * fabs(detDF); \
            /* Fill result(i,j) */ \
            result[grid.index(i_r,i_theta)] += factor * ( \
                - coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Bottom */ \
                - coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Top */ \
            ); \
            /* Fill result(i-1,j) */ \
            result[grid.index(i_r-1,i_theta)] += factor * ( \
                - 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Right */ \
                + 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Right */ \
            /* Fill result(i+1,j) */ \
            result[grid.index(i_r+1,i_theta)] += factor * ( \
                + 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Left */ \
                - 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Left */ \
        } \
        /* --------------------- */ \
        /* Outside Section Parts */ \
        else if(node_color != color){ \
            /* Fill result(i,j-1) */ \
            result[grid.index(i_r,i_theta-1)] += factor * ( \
                - coeff3 * att * x[grid.index(i_r,i_theta)] /* Top */ \
                - 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Top Right */ \
                + 0.25 * art * x[grid.index(i_r-1,i_theta)] ); /* Top Left */ \
            /* Fill result(i,j+1) */ \
            result[grid.index(i_r,i_theta+1)] += factor * ( \
                - coeff4 * att * x[grid.index(i_r,i_theta)] /* Bottom */ \
                + 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Bottom Right */ \
                - 0.25 * art * x[grid.index(i_r-1,i_theta)] ); /* Bottom Left */ \
        } \
    } \
    else if(i_r == grid.numberSmootherCircles()-1){ \
        double h1 = grid.r_dist(i_r-1); \
        double h2 = grid.r_dist(i_r); \
        double k1 = grid.theta_dist(i_theta-1); \
        double k2 = grid.theta_dist(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        /* -------------------- */ \
        /* Inside Section Parts */ \
        if(node_color == color){ \
            /* Fill result(i+1,j) */ \
            result[grid.index(i_r+1,i_theta)] += factor * ( \
                - coeff2 * arr * x[grid.index(i_r,i_theta)] /* Left */ \
                + 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Left */ \
                - 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Left */ \
        } \
        /* --------------------- */ \
        /* Outside Section Parts */ \
        else if(node_color != color){ \
            /* Nothing to be done here */ \
        } \
    } \
    else if(i_r == grid.numberSmootherCircles()){ \
        double h1 = grid.r_dist(i_r-1); \
        double h2 = grid.r_dist(i_r); \
        double k1 = grid.theta_dist(i_theta-1); \
        double k2 = grid.theta_dist(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        /* -------------------- */ \
        /* Inside Section Parts */ \
        if(node_color == color){ \
            /* Fill result of (i,j) */ \
            result[grid.index(i_r,i_theta)] += \
                0.25 * (h1+h2)*(k1+k2) * system_parameters.rhs_f(r, theta, sin_theta, cos_theta) * fabs(detDF); \
            /* Fill result(i,j) */ \
            result[grid.index(i_r,i_theta)] += factor * ( \
                - coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Left */ \
                - coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Bottom */ \
                - coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Top */ \
            ); \
            /* Fill result(i+1,j) */ \
            result[grid.index(i_r+1,i_theta)] += factor * ( \
                + 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Left */ \
                - 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Left */ \
        } \
        /* --------------------- */ \
        /* Outside Section Parts */ \
        else if(node_color != color){ \
            /* Fill result(i,j-1) */ \
            result[grid.index(i_r,i_theta-1)] += factor * ( \
                - coeff3 * att * x[grid.index(i_r,i_theta)] /* Top */ \
                - 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Top Right */ \
                + 0.25 * art * x[grid.index(i_r-1,i_theta)] ); /* Top Left */ \
            /* Fill result(i,j+1) */ \
            result[grid.index(i_r,i_theta+1)] += factor * ( \
                - coeff4 * att * x[grid.index(i_r,i_theta)] /* Bottom */ \
                + 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Bottom Right */ \
                - 0.25 * art * x[grid.index(i_r-1,i_theta)] ); /* Bottom Left */ \
        } \
    } \
    else if(i_r == grid.nr()-2){ \
        double h1 = grid.r_dist(i_r-1); \
        double h2 = grid.r_dist(i_r); \
        double k1 = grid.theta_dist(i_theta-1); \
        double k2 = grid.theta_dist(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        /* -------------------- */ \
        /* Inside Section Parts */ \
        if(node_color == color){ \
            /* Fill result of (i,j) */ \
            result[grid.index(i_r,i_theta)] += \
                0.25 * (h1+h2)*(k1+k2) * system_parameters.rhs_f(r, theta, sin_theta, cos_theta) * fabs(detDF); \
            /* Fill result(i,j) */ \
            result[grid.index(i_r,i_theta)] += factor * ( \
            - coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Bottom */ \
            - coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Top */ \
            ); \
            /* Fill result(i-1,j) */ \
            result[grid.index(i_r-1,i_theta)] += factor * ( \
                - 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Right */ \
                + 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Right */ \
            \
            /* "Right" is part of the radial Asc^ortho matrices, */ \
            /* but is shifted over to the rhs to make the radial Asc^ortho matrices symmetric. */ \
            /* Note that the circle Asc^ortho matrices are symmetric by default. */ \
            result[grid.index(i_r,i_theta)] -= /* Right: Symmetry shift! */ \
                -coeff2 * arr * system_parameters.u_D(r+h2, theta, sin_theta, cos_theta); \
        } \
        /* --------------------- */ \
        /* Outside Section Parts */ \
        else if(node_color != color){ \
            /* Fill result(i,j-1) */ \
            result[grid.index(i_r,i_theta-1)] += factor * ( \
                - coeff3 * att * x[grid.index(i_r,i_theta)] /* Top */ \
                - 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Top Right */ \
                + 0.25 * art * x[grid.index(i_r-1,i_theta)] ); /* Top Left */ \
            /* Fill result(i,j+1) */ \
            result[grid.index(i_r,i_theta+1)] += factor * ( \
                - coeff4 * att * x[grid.index(i_r,i_theta)] /* Bottom */ \
                + 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Bottom Right */ \
                - 0.25 * art * x[grid.index(i_r-1,i_theta)] ); /* Bottom Left */ \
        } \
    } \
    else if(i_r == grid.nr()-1){ \
        /* -------------------- */ \
        /* Inside Section Parts */ \
        if(node_color == color){ \
            /* Fill result of (i,j) */ \
            result[grid.index(i_r,i_theta)] += system_parameters.u_D(r, theta, sin_theta, cos_theta); \
            /* Fill result(i-1,j) */ \
            result[grid.index(i_r-1,i_theta)] += factor * ( \
                - 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Right */ \
                + 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Right */ \
        } \
        /* --------------------- */ \
        /* Outside Section Parts */ \
        else if(node_color != color){ \
            /* Nothing to be done here */ \
        } \
    } \
} while(0) \



#define CIRCLE_SECTION_APPLY_ASC_ORTHO_GIVE(i_r, color) \
do { \
    /* std::cout << "CIRCLE ASC i_r = " << i_r << ", color = " << toString(color) << std::endl; */ \
    assert(i_r >= 0 && i_r <= grid_.numberSmootherCircles()); \
    r = grid_.radius(i_r); \
    coeff_alpha = system_parameters_.alpha(r); \
    coeff_beta = system_parameters_.beta(r); \
    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++){ \
        theta = grid_.theta(i_theta); \
        sin_theta = sin_theta_[i_theta]; \
        cos_theta = cos_theta_[i_theta]; \
        \
        ARR_ATT_ART(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, \
            arr, att, art, detDF); \
        \
        NODE_APPLY_ASC_ORTHO_CIRCLE_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, \
            system_parameters_, grid_, DirBC_Interior_, color, \
            x, temp_rhs, factor, \
            arr, att, art, coeff_beta, detDF); \
    } \
} while(0)

#define RADIAL_SECTION_APPLY_ASC_ORTHO_GIVE(i_theta, color) \
do { \
    /* std::cout << "RADIAL ASC i_theta = " << i_theta << ", color = " << toString(color) << std::endl; */ \
    theta = grid_.theta(i_theta); \
    sin_theta = sin_theta_[i_theta]; \
    cos_theta = cos_theta_[i_theta]; \
    for (int i_r = grid_.numberSmootherCircles()-1; i_r < grid_.nr(); i_r++){ \
        r = grid_.radius(i_r); \
        coeff_alpha = system_parameters_.alpha(r); \
        coeff_beta = system_parameters_.beta(r); \
        /* Get arr, att, art, detDF value at the current node */ \
        ARR_ATT_ART(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, \
            arr, att, art, detDF); \
        \
        NODE_APPLY_ASC_ORTHO_RADIAL_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, \
            system_parameters_, grid_, DirBC_Interior_, color, \
            x, temp_rhs, factor, \
            arr, att, art, coeff_beta, detDF); \
    } \
} while(0)



#define CIRCLE_SECTION_SOLVE_SMOOTER(i_r) \
do { \
    /* std::cout << "CIRCLE SOLVE i_r = " << i_r << std::endl; */ \
    const int start = grid_.index(i_r, 0); \
    const int end = start + grid_.ntheta(); \
    if(i_r == 0){ \
        inner_boundary_circle_Asc_mumps_.job = JOB_COMPUTE_SOLUTION; \
        inner_boundary_circle_Asc_mumps_.nrhs = 1; \
        inner_boundary_circle_Asc_mumps_.nz_rhs = end - start; \
        inner_boundary_circle_Asc_mumps_.rhs = temp_rhs.begin() + start; \
        inner_boundary_circle_Asc_mumps_.lrhs = end - start; \
        dmumps_c(&inner_boundary_circle_Asc_mumps_); \
        if (inner_boundary_circle_Asc_mumps_.info[0] != 0) { \
            std::cerr << "Error solving the system: " << inner_boundary_circle_Asc_mumps_.info[0] << std::endl; \
        } \
        std::move(temp_rhs.begin() + start, \
            temp_rhs.begin() + end, \
            x.begin() + start \
        );  \
    } \
    else { \
        circle_symmetric_cyclic_tridiagonal_solver_[i_r].solveInPlace(temp_rhs.begin() + start, circle_solver_storage_1.begin(), circle_solver_storage_2.begin()); \
        std::move(temp_rhs.begin() + start, \
            temp_rhs.begin() + end, \
            x.begin() + start \
        );  \
    } \
} while(0)

#define RADIAL_SECTION_SOLVE_SMOOTER(i_theta) \
do { \
    /* std::cout << "RADIAL SOLVE i_theta = " << i_theta << std::endl; */ \
    const int start = grid_.index(grid_.numberSmootherCircles(), i_theta); \
    const int end = start + grid_.lengthSmootherRadial(); \
    radial_symmetric_tridiagonal_solver_[i_theta].solveInPlace(temp_rhs.begin() + start, radial_solver_storage.begin()); \
    std::move(temp_rhs.begin() + start, \
        temp_rhs.begin() + end, \
        x.begin() + start \
    ); \
} while(0)



void Smoother::smoothing(Vector<double>& x, Vector<double>& temp_rhs){
    assert(x.size() == temp_rhs.size());

    omp_set_num_threads(maxOpenMPThreads_);

    assign(temp_rhs, 0.0);

    const double factor = -1.0;

    const int numCircleTasks = grid_.numberSmootherCircles();
    const int numRadialTasks = grid_.ntheta();

    const int additionalRadialTasks = numRadialTasks % 3;

    assert(numCircleTasks >= 2);
    assert(numRadialTasks >= 3 && numRadialTasks % 2 == 0);

    const int shift = 2; /* To stay in bounds of the dependency arrays */

    int* asc_ortho_circle_dep = new int[numCircleTasks+1 + 2*shift];
    int* asc_ortho_radial_dep = new int[numRadialTasks + 2*shift];

    int* smoother_circle_dep = new int[numCircleTasks + 2*shift];
    int* smoother_radial_dep = new int[numRadialTasks + 2*shift];

    omp_set_num_threads(openMPTaskThreads_);
    #pragma omp parallel num_threads(openMPTaskThreads_) /* Outside variable are shared by default */
    {
        /* Define thread-local variables */
        double r, theta;
        double sin_theta, cos_theta;
        double arr, att, art;
        double coeff_alpha, coeff_beta;
        double detDF;

        Vector<double> circle_solver_storage_1(grid_.ntheta());
        Vector<double> circle_solver_storage_2(grid_.ntheta());
        Vector<double> radial_solver_storage(grid_.lengthSmootherRadial());

        #pragma omp single
        {
            /* ---------------------------- */
            /* ------ CIRCLE SECTION ------ */

            /* ---------------------------- */
            /* Asc ortho Black Circle Tasks */
            /* ---------------------------- */

            /* Mod 0 Black Circles */
            for(int circle_task = -1; circle_task < numCircleTasks; circle_task += 3) {
                #pragma omp task \
                    depend(out: asc_ortho_circle_dep[circle_task + shift+1])
                {
                    int i_r = numCircleTasks - circle_task - 1;    
                    CIRCLE_SECTION_APPLY_ASC_ORTHO_GIVE(i_r, SmootherColor::Black);
                }
            }
            /* Mod 1 Black Circles */
            for(int circle_task = 0; circle_task < numCircleTasks; circle_task += 3) {
                #pragma omp task \
                    depend(out: \
                        asc_ortho_circle_dep[circle_task + shift+1]) \
                    depend(in: \
                        asc_ortho_circle_dep[circle_task-1 + shift+1], \
                        asc_ortho_circle_dep[circle_task+2 + shift+1])   
                {
                    int i_r = numCircleTasks - circle_task - 1;    
                    CIRCLE_SECTION_APPLY_ASC_ORTHO_GIVE(i_r, SmootherColor::Black);
                }
                
            }
            /* Mod 2 Black Circles */
            for(int circle_task = 1; circle_task < numCircleTasks; circle_task += 3) {
                #pragma omp task \
                    depend(out: \
                        asc_ortho_circle_dep[circle_task + shift+1]) \
                    depend(in: \
                        asc_ortho_circle_dep[circle_task-1 + shift+1], \
                        asc_ortho_circle_dep[circle_task+2 + shift+1])  
                {
                    int i_r = numCircleTasks - circle_task - 1;    
                    CIRCLE_SECTION_APPLY_ASC_ORTHO_GIVE(i_r, SmootherColor::Black);
                }
            }

            /* Black Circle Smoother */
            for(int circle_task = 0; circle_task < numCircleTasks; circle_task += 2) {
                #pragma omp task \
                    depend(out: smoother_circle_dep[circle_task + shift]) \
                    depend(in: \
                        asc_ortho_circle_dep[circle_task-1 + shift+1], \
                        asc_ortho_circle_dep[circle_task+0 + shift+1], \
                        asc_ortho_circle_dep[circle_task+1 + shift+1])   
                {
                    int i_r = numCircleTasks - circle_task - 1;    
                    CIRCLE_SECTION_SOLVE_SMOOTER(i_r);
                }
            }

            /* ---------------------------- */
            /* Asc ortho White Circle Tasks */
            /* ---------------------------- */

            /* Mod 0 White Circles */
            for(int circle_task = 0; circle_task < numCircleTasks; circle_task += 3) {
                #pragma omp task \
                    depend(out: asc_ortho_circle_dep[circle_task + shift+1]) \
                    depend(in: \
                        smoother_circle_dep[circle_task-1 + shift], \
                        smoother_circle_dep[circle_task+0 + shift], \
                        smoother_circle_dep[circle_task+1 + shift])
                {
                    int i_r = numCircleTasks - circle_task - 1;    
                    CIRCLE_SECTION_APPLY_ASC_ORTHO_GIVE(i_r, SmootherColor::White);
                }
            }
            /* Mod 1 White Circles */
            for(int circle_task = 1; circle_task < numCircleTasks; circle_task += 3) {
                #pragma omp task \
                    depend(out: asc_ortho_circle_dep[circle_task + shift+1]) \
                    depend(in: \
                        smoother_circle_dep[circle_task-1 + shift], \
                        smoother_circle_dep[circle_task+0 + shift], \
                        smoother_circle_dep[circle_task+1 + shift], \
                        asc_ortho_circle_dep[circle_task-1 + shift+1], \
                        asc_ortho_circle_dep[circle_task+2 + shift+1])  
                {
                    int i_r = numCircleTasks - circle_task - 1;    
                    CIRCLE_SECTION_APPLY_ASC_ORTHO_GIVE(i_r, SmootherColor::White);
                }
                
            }
            /* Mod 2 White Circles */
            for(int circle_task = 2; circle_task < numCircleTasks; circle_task += 3) {
                #pragma omp task \
                    depend(out: asc_ortho_circle_dep[circle_task + shift+1]) \
                    depend(in: \
                        smoother_circle_dep[circle_task-1 + shift], \
                        smoother_circle_dep[circle_task+0 + shift], \
                        smoother_circle_dep[circle_task+1 + shift], \
                        asc_ortho_circle_dep[circle_task-1 + shift+1], \
                        asc_ortho_circle_dep[circle_task+2 + shift+1])  
                {
                    int i_r = numCircleTasks - circle_task - 1;    
                    CIRCLE_SECTION_APPLY_ASC_ORTHO_GIVE(i_r, SmootherColor::White);
                }
            }

            /* White Circle Smoother */
            for(int circle_task = 1; circle_task < numCircleTasks; circle_task += 2) {
                #pragma omp task \
                    depend(out: smoother_circle_dep[circle_task + shift]) \
                    depend(in: \
                        asc_ortho_circle_dep[circle_task-1 + shift+1], \
                        asc_ortho_circle_dep[circle_task+0 + shift+1], \
                        asc_ortho_circle_dep[circle_task+1 + shift+1])   
                {
                    int i_r = numCircleTasks - circle_task - 1;    
                    CIRCLE_SECTION_SOLVE_SMOOTER(i_r);
                }
            }

            /* ---------------------------- */
            /* ------ RADIAL SECTION ------ */

            /* ---------------------------- */
            /* Asc ortho Black Radial Tasks */
            /* ---------------------------- */

            /* Mod 0 Black Radials */
            for(int radial_task = 0; radial_task < numRadialTasks - additionalRadialTasks; radial_task += 3) {
                #pragma omp task \
                    depend(out: asc_ortho_radial_dep[radial_task + shift]) \
                    depend(in : smoother_circle_dep[0 + shift])
                {
                    int i_theta = radial_task;    
                    RADIAL_SECTION_APPLY_ASC_ORTHO_GIVE(i_theta, SmootherColor::Black);
                }
            }
            /* Mod 1 Black Radials */
            for(int radial_task = 1; radial_task < numRadialTasks - additionalRadialTasks; radial_task += 3) {
                #pragma omp task \
                    depend(out: asc_ortho_radial_dep[radial_task + shift]) \
                    depend(in: \
                        asc_ortho_radial_dep[radial_task-1 + shift], \
                        asc_ortho_radial_dep[radial_task+2 + shift])   
                {
                    int i_theta = radial_task;    
                    RADIAL_SECTION_APPLY_ASC_ORTHO_GIVE(i_theta, SmootherColor::Black);
                }
                
            }
            /* Mod 2 Black Radials */
            for(int radial_task = 2; radial_task < numRadialTasks - additionalRadialTasks; radial_task += 3) {
                #pragma omp task \
                    depend(out: asc_ortho_radial_dep[radial_task + shift]) \
                    depend(in: \
                        asc_ortho_radial_dep[radial_task-1 + shift], \
                        asc_ortho_radial_dep[radial_task+2 + shift])    
                {
                    int i_theta = radial_task;    
                    RADIAL_SECTION_APPLY_ASC_ORTHO_GIVE(i_theta, SmootherColor::Black);
                }
            }

            /* First additional Radial */
            if(additionalRadialTasks >= 1){
                int radial_task = numRadialTasks - additionalRadialTasks;
                #pragma omp task \
                    depend(out: asc_ortho_radial_dep[radial_task + shift]) \
                    depend(in: \
                        asc_ortho_radial_dep[radial_task-1 + shift])    
                {
                    int i_theta = radial_task;    
                    RADIAL_SECTION_APPLY_ASC_ORTHO_GIVE(i_theta, SmootherColor::Black);
                }
            }

            /* Second additional Radial */
            if(additionalRadialTasks >= 2){
                int radial_task = numRadialTasks - additionalRadialTasks + 1;
                #pragma omp task \
                    depend(out: asc_ortho_radial_dep[radial_task + shift]) \
                    depend(in: \
                        asc_ortho_radial_dep[radial_task-1 + shift])
                {
                    int i_theta = radial_task;    
                    RADIAL_SECTION_APPLY_ASC_ORTHO_GIVE(i_theta, SmootherColor::Black);
                }
            }

            /* Black Radial Smoother */
            for(int radial_task = 0; radial_task < numRadialTasks; radial_task += 2) {
                #pragma omp task \
                    depend(out: smoother_radial_dep[radial_task + shift]) \
                    depend(in: \
                        asc_ortho_radial_dep[(radial_task-1 + numRadialTasks) % numRadialTasks + shift], \
                        asc_ortho_radial_dep[(radial_task+0 + numRadialTasks) % numRadialTasks + shift], \
                        asc_ortho_radial_dep[(radial_task+1 + numRadialTasks) % numRadialTasks + shift])     
                {
                    int i_theta = radial_task;    
                    RADIAL_SECTION_SOLVE_SMOOTER(i_theta);
                }
            }

            /* ---------------------------- */
            /* Asc ortho White Circle Tasks */
            /* ---------------------------- */

            /* Mod 0 White Radials */
            for(int radial_task = 0; radial_task < numRadialTasks - additionalRadialTasks; radial_task += 3) {
                #pragma omp task \
                    depend(out: asc_ortho_radial_dep[radial_task + shift]) \
                    depend(in : \
                        smoother_radial_dep[(radial_task-1 + numRadialTasks) % numRadialTasks + shift], \
                        smoother_radial_dep[(radial_task+0 + numRadialTasks) % numRadialTasks + shift], \
                        smoother_radial_dep[(radial_task+1 + numRadialTasks) % numRadialTasks + shift])
                {
                    int i_theta = radial_task;    
                    RADIAL_SECTION_APPLY_ASC_ORTHO_GIVE(i_theta, SmootherColor::White);
                }
            }
            /* Mod 1 White Radials */
            for(int radial_task = 1; radial_task < numRadialTasks - additionalRadialTasks; radial_task += 3) {
                #pragma omp task \
                    depend(out: asc_ortho_radial_dep[radial_task + shift]) \
                    depend(in: \
                        asc_ortho_radial_dep[radial_task-1 + shift], \
                        asc_ortho_radial_dep[radial_task+2 + shift])   
                {
                    int i_theta = radial_task;    
                    RADIAL_SECTION_APPLY_ASC_ORTHO_GIVE(i_theta, SmootherColor::White);
                }
                
            }
            /* Mod 2 White Radials */
            for(int radial_task = 2; radial_task < numRadialTasks - additionalRadialTasks; radial_task += 3) {
                #pragma omp task \
                    depend(out: asc_ortho_radial_dep[radial_task + shift]) \
                    depend(in: \
                        asc_ortho_radial_dep[radial_task-1 + shift], \
                        asc_ortho_radial_dep[radial_task+2 + shift])    
                {
                    int i_theta = radial_task;    
                    RADIAL_SECTION_APPLY_ASC_ORTHO_GIVE(i_theta, SmootherColor::White);
                }
            }

            /* First additional Radial */
            if(additionalRadialTasks >= 1){
                int radial_task = numRadialTasks - additionalRadialTasks;
                #pragma omp task \
                    depend(out: asc_ortho_radial_dep[radial_task + shift]) \
                    depend(in: \
                        asc_ortho_radial_dep[radial_task-1 + shift])    
                {
                    int i_theta = radial_task;    
                    RADIAL_SECTION_APPLY_ASC_ORTHO_GIVE(i_theta, SmootherColor::White);
                }
            }

            /* Second additional Radial */
            if(additionalRadialTasks >= 2){
                int radial_task = numRadialTasks - additionalRadialTasks + 1;
                #pragma omp task \
                    depend(out: asc_ortho_radial_dep[radial_task + shift]) \
                    depend(in: \
                        asc_ortho_radial_dep[radial_task-1 + shift])
                {
                    int i_theta = radial_task;    
                    RADIAL_SECTION_APPLY_ASC_ORTHO_GIVE(i_theta, SmootherColor::White);
                }
            }

            /* White Radial Smoother */
            for(int radial_task = 1; radial_task < numRadialTasks; radial_task += 2) {
                #pragma omp task \
                    depend(out: smoother_radial_dep[radial_task + shift]) \
                    depend(in: \
                        asc_ortho_radial_dep[(radial_task-1 + numRadialTasks) % numRadialTasks + shift], \
                        asc_ortho_radial_dep[(radial_task+0 + numRadialTasks) % numRadialTasks + shift], \
                        asc_ortho_radial_dep[(radial_task+1 + numRadialTasks) % numRadialTasks + shift])     
                {
                    int i_theta = radial_task;    
                    RADIAL_SECTION_SOLVE_SMOOTER(i_theta);
                }
            }
        }
    }

    omp_set_num_threads(maxOpenMPThreads_);

    delete[] asc_ortho_circle_dep;
    delete[] asc_ortho_radial_dep;

    delete[] smoother_circle_dep;
    delete[] smoother_radial_dep;
}
