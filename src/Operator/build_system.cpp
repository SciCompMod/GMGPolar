#include "../../include/Operator/operator.h"

#include "../../include/common/constants.h"

#define NODE_BUILD_A_GIVE(i_r, i_theta, grid, DirBC_Interior, \
    matrixA, nz_index, \
    center_nz_index, left_nz_index, right_nz_index, bottom_nz_index, top_nz_index, \
    center_index, left_index, right_index, bottom_index, top_index, \
    arr, att, art, coeff_beta, detDF) \
do { \
    /* -------------------- */ \
    /* Node in the interior */ \
    /* -------------------- */ \
    if (i_r > 1 && i_r < grid.nr() - 2) { \
        scalar_t h1 = grid.r_dist(i_r-1); \
        scalar_t h2 = grid.r_dist(i_r); \
        scalar_t k1 = grid.theta_dist(i_theta-1); \
        scalar_t k2 = grid.theta_dist(i_theta); \
        scalar_t coeff1 = 0.5*(k1+k2)/h1; \
        scalar_t coeff2 = 0.5*(k1+k2)/h2; \
        scalar_t coeff3 = 0.5*(h1+h2)/k1; \
        scalar_t coeff4 = 0.5*(h1+h2)/k2; \
        \
        center_nz_index = ptr_nz_index_matrixA(grid, i_r, i_theta, DirBC_Interior); \
        left_nz_index = ptr_nz_index_matrixA(grid, i_r-1, i_theta, DirBC_Interior); \
        right_nz_index = ptr_nz_index_matrixA(grid, i_r+1, i_theta, DirBC_Interior); \
        bottom_nz_index = ptr_nz_index_matrixA(grid, i_r, ((i_theta-1) + grid.ntheta()) % grid.ntheta(), DirBC_Interior); \
        top_nz_index = ptr_nz_index_matrixA(grid, i_r, ((i_theta+1) + grid.ntheta()) % grid.ntheta(), DirBC_Interior); \
        \
        center_index = grid.index(i_r,i_theta); \
        left_index = grid.index(i_r-1,i_theta); \
        right_index = grid.index(i_r+1,i_theta); \
        bottom_index = grid.index(i_r,i_theta-1); \
        top_index = grid.index(i_r,i_theta+1); \
        \
        /* Fill matrix row of (i,j) */ \
        const Stencil& CenterStencil = get_stencil(grid, i_r, DirBC_Interior); \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = center_index; \
        matrixA.col_index(nz_index) = center_index; \
        matrixA.value(nz_index) += 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Left]; \
        matrixA.row_index(nz_index) = center_index; \
        matrixA.col_index(nz_index) = left_index; \
        matrixA.value(nz_index) += - coeff1 * arr; /* Left */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Right]; \
        matrixA.row_index(nz_index) = center_index; \
        matrixA.col_index(nz_index) = right_index; \
        matrixA.value(nz_index) += - coeff2 * arr; /* Right */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Bottom]; \
        matrixA.row_index(nz_index) = center_index; \
        matrixA.col_index(nz_index) = bottom_index; \
        matrixA.value(nz_index) += - coeff3 * att; /* Bottom */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Top]; \
        matrixA.row_index(nz_index) = center_index; \
        matrixA.col_index(nz_index) = top_index; \
        matrixA.value(nz_index) += - coeff4 * att; /* Top */ \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = center_index; \
        matrixA.col_index(nz_index) = center_index; \
        /* Center: (Left, Right, Bottom, Top) */ \
        matrixA.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; \
        \
        /* Fill matrix row of (i-1,j) */ \
        const Stencil& LeftStencil = get_stencil(grid, i_r-1, DirBC_Interior); \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::Right]; \
        matrixA.row_index(nz_index) = left_index; \
        matrixA.col_index(nz_index) = center_index; \
        matrixA.value(nz_index) += - coeff1 * arr; /* Right */ \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = left_index; \
        matrixA.col_index(nz_index) = left_index; \
        matrixA.value(nz_index) += coeff1 * arr; /* Center: (Right) */ \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::TopRight]; \
        matrixA.row_index(nz_index) = left_index; \
        matrixA.col_index(nz_index) = top_index; \
        matrixA.value(nz_index) += - 0.25 * art; /* Top Right */ \
        \
        nz_index = left_nz_index + LeftStencil[StencilType::BottomRight]; \
        matrixA.row_index(nz_index) = left_index; \
        matrixA.col_index(nz_index) = bottom_index; \
        matrixA.value(nz_index) += 0.25 * art; /* Bottom Right */ \
        \
        /* Fill matrix row of (i+1,j) */ \
        const Stencil& RightStencil = get_stencil(grid, i_r+1, DirBC_Interior); \
        \
        nz_index = right_nz_index + RightStencil[StencilType::Left]; \
        matrixA.row_index(nz_index) = right_index; \
        matrixA.col_index(nz_index) = center_index; \
        matrixA.value(nz_index) += - coeff2 * arr; /* Left */ \
        \
        nz_index = right_nz_index + RightStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = right_index; \
        matrixA.col_index(nz_index) = right_index; \
        matrixA.value(nz_index) += coeff2 * arr; /* Center: (Left) */ \
        \
        nz_index = right_nz_index + RightStencil[StencilType::TopLeft]; \
        matrixA.row_index(nz_index) = right_index; \
        matrixA.col_index(nz_index) = top_index; \
        matrixA.value(nz_index) += 0.25 * art; /* Top Left */ \
        \
        nz_index = right_nz_index + RightStencil[StencilType::BottomLeft]; \
        matrixA.row_index(nz_index) = right_index; \
        matrixA.col_index(nz_index) = bottom_index; \
        matrixA.value(nz_index) += - 0.25 * art; /* Bottom Left */ \
        \
        /* Fill matrix row of (i,j-1) */ \
        const Stencil& BottomStencil = CenterStencil; \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::Top]; \
        matrixA.row_index(nz_index) = bottom_index; \
        matrixA.col_index(nz_index) = center_index; \
        matrixA.value(nz_index) += - coeff3 * att; /* Top */ \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = bottom_index; \
        matrixA.col_index(nz_index) = bottom_index; \
        matrixA.value(nz_index) += coeff3 * att; /* Center: (Top) */ \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::TopRight]; \
        matrixA.row_index(nz_index) = bottom_index; \
        matrixA.col_index(nz_index) = right_index; \
        matrixA.value(nz_index) += - 0.25 * art; /* Top Right */ \
        \
        nz_index = bottom_nz_index + BottomStencil[StencilType::TopLeft]; \
        matrixA.row_index(nz_index) = bottom_index; \
        matrixA.col_index(nz_index) = left_index; \
        matrixA.value(nz_index) += 0.25 * art; /* Top Left */ \
        \
        /* Fill matrix row of (i,j+1) */ \
        const Stencil& TopStencil = CenterStencil; \
        \
        nz_index = top_nz_index + TopStencil[StencilType::Bottom]; \
        matrixA.row_index(nz_index) = top_index; \
        matrixA.col_index(nz_index) = center_index; \
        matrixA.value(nz_index) += - coeff4 * att; /* Bottom */ \
        \
        nz_index = top_nz_index + TopStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = top_index; \
        matrixA.col_index(nz_index) = top_index; \
        matrixA.value(nz_index) += coeff4 * att; /* Center: (Bottom) */ \
        \
        nz_index = top_nz_index + TopStencil[StencilType::BottomRight]; \
        matrixA.row_index(nz_index) = top_index; \
        matrixA.col_index(nz_index) = right_index; \
        matrixA.value(nz_index) += 0.25 * art; /* Bottom Right */ \
        \
        nz_index = top_nz_index + TopStencil[StencilType::BottomLeft]; \
        matrixA.row_index(nz_index) = top_index; \
        matrixA.col_index(nz_index) = left_index; \
        matrixA.value(nz_index) += - 0.25 * art; /* Bottom Left */ \
        \
    /* ------------------------ */ \
    /* Node in the inner circle */ \
    /* ------------------------ */ \
    } else if (i_r == 0) { \
        /* Case 1: Dirichlet boundary on the interior ring */ \
        if(DirBC_Interior){ \
        /* Fill result(i,j) */ \
        scalar_t h2 = grid.r_dist(i_r); \
        scalar_t k1 = grid.theta_dist(i_theta-1); \
        scalar_t k2 = grid.theta_dist(i_theta); \
        scalar_t coeff2 = 0.5*(k1+k2)/h2; \
        \
        center_nz_index = ptr_nz_index_matrixA(grid, i_r, i_theta, DirBC_Interior); \
        right_nz_index = ptr_nz_index_matrixA(grid, i_r+1, i_theta, DirBC_Interior); \
        \
        center_index = grid.index(i_r,i_theta); \
        right_index = grid.index(i_r+1,i_theta); \
        bottom_index = grid.index(i_r,i_theta-1); \
        top_index = grid.index(i_r,i_theta+1); \
        \
        /* Fill matrix row of (i,j) */ \
        const Stencil& CenterStencil = get_stencil(grid, i_r, DirBC_Interior); \
        \
        nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = center_index; \
        matrixA.col_index(nz_index) = center_index; \
        matrixA.value(nz_index) = 1; \
        \
        /* Give value to the interior nodes! */ \
        /* Fill matrix row of (i+1,j) */ \
        const Stencil& RightStencil = get_stencil(grid, i_r+1, DirBC_Interior); \
        \
        nz_index = right_nz_index + RightStencil[StencilType::Left]; \
        matrixA.row_index(nz_index) = right_index; \
        matrixA.col_index(nz_index) = center_index; \
        matrixA.value(nz_index) += - coeff2 * arr; /* Left */ \
        \
        nz_index = right_nz_index + RightStencil[StencilType::Center]; \
        matrixA.row_index(nz_index) = right_index; \
        matrixA.col_index(nz_index) = right_index; \
        matrixA.value(nz_index) += coeff2 * arr; /* Center: (Left) */ \
        \
        nz_index = right_nz_index + RightStencil[StencilType::TopLeft]; \
        matrixA.row_index(nz_index) = right_index; \
        matrixA.col_index(nz_index) = top_index; \
        matrixA.value(nz_index) += 0.25 * art; /* Top Left */ \
        \
        nz_index = right_nz_index + RightStencil[StencilType::BottomLeft]; \
        matrixA.row_index(nz_index) = right_index; \
        matrixA.col_index(nz_index) = bottom_index; \
        matrixA.value(nz_index) += - 0.25 * art; /* Bottom Left */ \
        \
        /* Case 2: Across origin discretization */ \
        } else{ \
            /* h1 gets replaced with 2 * R0. */ \
            /* (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()>>1)). */ \
            /* Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil. */ \
            scalar_t h1 = 2 * grid.radius(0); \
            scalar_t h2 = grid.r_dist(i_r); \
            scalar_t k1 = grid.theta_dist(i_theta-1); \
            scalar_t k2 = grid.theta_dist(i_theta); \
            scalar_t coeff1 = 0.5*(k1+k2)/h1; \
            scalar_t coeff2 = 0.5*(k1+k2)/h2; \
            scalar_t coeff3 = 0.5*(h1+h2)/k1; \
            scalar_t coeff4 = 0.5*(h1+h2)/k2; \
            \
            center_nz_index = ptr_nz_index_matrixA(grid, i_r, i_theta, DirBC_Interior); \
            left_nz_index = ptr_nz_index_matrixA(grid, i_r, (i_theta + (grid.ntheta()>>1)) % grid.ntheta(), DirBC_Interior); \
            right_nz_index = ptr_nz_index_matrixA(grid, i_r+1, i_theta, DirBC_Interior); \
            bottom_nz_index = ptr_nz_index_matrixA(grid, i_r, ((i_theta-1) + grid.ntheta()) % grid.ntheta(), DirBC_Interior); \
            top_nz_index = ptr_nz_index_matrixA(grid, i_r, ((i_theta+1) + grid.ntheta()) % grid.ntheta(), DirBC_Interior); \
            \
            center_index = grid.index(i_r,i_theta); \
            left_index = grid.index(i_r, i_theta + (grid.ntheta()>>1)); \
            right_index = grid.index(i_r+1,i_theta); \
            bottom_index = grid.index(i_r,i_theta-1); \
            top_index = grid.index(i_r,i_theta+1); \
            \
            /* Fill matrix row of (i,j) */ \
            const Stencil& CenterStencil = get_stencil(grid, i_r, DirBC_Interior); \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
            matrixA.row_index(nz_index) = center_index; \
            matrixA.col_index(nz_index) = center_index; \
            matrixA.value(nz_index) += 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Left]; \
            matrixA.row_index(nz_index) = center_index; \
            matrixA.col_index(nz_index) = left_index; \
            matrixA.value(nz_index) += - coeff1 * arr; /* Left */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Right]; \
            matrixA.row_index(nz_index) = center_index; \
            matrixA.col_index(nz_index) = right_index; \
            matrixA.value(nz_index) += - coeff2 * arr; /* Right */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Bottom]; \
            matrixA.row_index(nz_index) = center_index; \
            matrixA.col_index(nz_index) = bottom_index; \
            matrixA.value(nz_index) += - coeff3 * att; /* Bottom */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Top]; \
            matrixA.row_index(nz_index) = center_index; \
            matrixA.col_index(nz_index) = top_index; \
            matrixA.value(nz_index) += - coeff4 * att; /* Top */ \
            \
            nz_index = center_nz_index + CenterStencil[StencilType::Center]; \
            matrixA.row_index(nz_index) = center_index; \
            matrixA.col_index(nz_index) = center_index; \
            /* Center: (Left, Right, Bottom, Top) */ \
            matrixA.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; \
            \
            /* Fill matrix row of (i-1,j) */ \
            const Stencil& LeftStencil = get_stencil(grid, i_r-1, DirBC_Interior); \
            \
            nz_index = left_nz_index + LeftStencil[StencilType::Right]; \
            matrixA.row_index(nz_index) = left_index; \
            matrixA.col_index(nz_index) = center_index; \
            matrixA.value(nz_index) += - coeff1 * arr; /* Right */ \
            \
            nz_index = left_nz_index + LeftStencil[StencilType::Center]; \
            matrixA.row_index(nz_index) = left_index; \
            matrixA.col_index(nz_index) = left_index; \
            matrixA.value(nz_index) += coeff1 * arr; /* Center: (Right) */ \
            \
            /* Top Right: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = left_nz_index + LeftStencil[StencilType::TopRight]; */  \
            /* matrixA.row_index(nz_index) = left_index; */  \
            /* matrixA.col_index(nz_index) = top_index; */  \
            /* matrixA.value(nz_index) += - 0.25 * art; // Top Right */ \
            \
            /* Bottom Right: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* nz_index = left_nz_index + LeftStencil[StencilType::BottomRight]; */ \
            /* matrixA.row_index(nz_index) = left_index; */  \
            /* matrixA.col_index(nz_index) = bottom_index; */  \
            /* matrixA.value(nz_index) += 0.25 * art; // Bottom Right */ \
            \
            /* Fill matrix row of (i+1,j) */ \
            const Stencil& RightStencil = get_stencil(grid, i_r+1, DirBC_Interior); \
            \
            nz_index = right_nz_index + RightStencil[StencilType::Left]; \
            matrixA.row_index(nz_index) = right_index; \
            matrixA.col_index(nz_index) = center_index; \
            matrixA.value(nz_index) += - coeff2 * arr; /* Left */ \
            \
            nz_index = right_nz_index + RightStencil[StencilType::Center]; \
            matrixA.row_index(nz_index) = right_index; \
            matrixA.col_index(nz_index) = right_index; \
            matrixA.value(nz_index) += coeff2 * arr; /* Center: (Left) */ \
            \
            nz_index = right_nz_index + RightStencil[StencilType::TopLeft]; \
            matrixA.row_index(nz_index) = right_index; \
            matrixA.col_index(nz_index) = top_index; \
            matrixA.value(nz_index) += 0.25 * art; /* Top Left */ \
            \
            nz_index = right_nz_index + RightStencil[StencilType::BottomLeft]; \
            matrixA.row_index(nz_index) = right_index; \
            matrixA.col_index(nz_index) = bottom_index; \
            matrixA.value(nz_index) += - 0.25 * art; /* Bottom Left */ \
            \
            /* Fill matrix row of (i,j-1) */ \
            const Stencil& BottomStencil = CenterStencil; \
            \
            nz_index = bottom_nz_index + BottomStencil[StencilType::Top]; \
            matrixA.row_index(nz_index) = bottom_index; \
            matrixA.col_index(nz_index) = center_index; \
            matrixA.value(nz_index) += - coeff3 * att; /* Top */ \
            \
            nz_index = bottom_nz_index + BottomStencil[StencilType::Center]; \
            matrixA.row_index(nz_index) = bottom_index; \
            matrixA.col_index(nz_index) = bottom_index; \
            matrixA.value(nz_index) += coeff3 * att; /* Center: (Top) */ \
            \
            nz_index = bottom_nz_index + BottomStencil[StencilType::TopRight]; \
            matrixA.row_index(nz_index) = bottom_index; \
            matrixA.col_index(nz_index) = right_index; \
            matrixA.value(nz_index) += - 0.25 * art; /* Top Right */ \
            \
            nz_index = bottom_nz_index + BottomStencil[StencilType::TopLeft]; \
            matrixA.row_index(nz_index) = bottom_index; \
            matrixA.col_index(nz_index) = left_index; \
            matrixA.value(nz_index) += 0.25 * art; /* Top Left */ \
            \
            /* Fill matrix row of (i,j+1) */ \
            const Stencil& TopStencil = CenterStencil; \
            \
            nz_index = top_nz_index + TopStencil[StencilType::Bottom]; \
            matrixA.row_index(nz_index) = top_index; \
            matrixA.col_index(nz_index) = center_index; \
            matrixA.value(nz_index) += - coeff4 * att; /* Bottom */ \
            \
            nz_index = top_nz_index + TopStencil[StencilType::Center]; \
            matrixA.row_index(nz_index) = top_index; \
            matrixA.col_index(nz_index) = top_index; \
            matrixA.value(nz_index) += coeff4 * att; /* Center: (Bottom) */ \
            \
            nz_index = top_nz_index + TopStencil[StencilType::BottomRight]; \
            matrixA.row_index(nz_index) = top_index; \
            matrixA.col_index(nz_index) = right_index; \
            matrixA.value(nz_index) += 0.25 * art; /* Bottom Right */ \
            \
            nz_index = top_nz_index + TopStencil[StencilType::BottomLeft]; \
            matrixA.row_index(nz_index) = top_index; \
            matrixA.col_index(nz_index) = left_index; \
            matrixA.value(nz_index) += - 0.25 * art; /* Bottom Left */ \
            \





            /* Fill result(i,j) */ \
            result[grid.index(i_r,i_theta)] += factor * ( \
                0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF) * x[grid.index(i_r,i_theta)] /* beta_{i,j} */ \
                - coeff1 * arr * x[grid.index(i_r, i_theta + (grid.ntheta()>>1))] /* Left */ \
                - coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Right */ \
                - coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Bottom */ \
                - coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Top */ \
                /* Center: (Left, Right, Bottom, Top) */ \
                + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) * x[grid.index(i_r,i_theta)] ); \
            \
            /* Fill result(i-1,j) */ \
            if(GiveToLeft){ \
            result[grid.index(i_r, i_theta + (grid.ntheta()>>1))] += factor * ( \
                - coeff1 * arr * x[grid.index(i_r,i_theta)] /* Right */ \
                + coeff1 * arr * x[grid.index(i_r, i_theta + (grid.ntheta()>>1))] ); /* Center: (Right) */ \
            /*  + 0.25 * art * x[grid.index(i_r,i_theta+1)]; // Top Right: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /*  - 0.25 * art * x[grid.index(i_r,i_theta-1)]; // Bottom Right: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            } \
            /* Fill result(i+1,j) */ \
            if(GiveToRight) { \
            result[grid.index(i_r+1,i_theta)] += factor * ( \
                - coeff2 * arr * x[grid.index(i_r,i_theta)] /* Left */ \
                + coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Center: (Left) */ \
                + 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Left */ \
                - 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Left */ \
            } \
            /* Fill result(i,j-1) */ \
            if(GiveToBottom) { \
            result[grid.index(i_r,i_theta-1)] += factor * ( \
                - coeff3 * att * x[grid.index(i_r,i_theta)] /* Top */ \
                + coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Center: (Top) */ \
                - 0.25 * art * x[grid.index(i_r+1,i_theta)] ); /* Top Right */ \
            /*  + 0.25 * art * x[grid.index(i_r-1,i_theta)]; // Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            } \
            /* Fill result(i,j+1) */ \
            if(GiveToTop) { \
            result[grid.index(i_r,i_theta+1)] += factor * ( \
                - coeff4 * att * x[grid.index(i_r,i_theta)] /* Bottom */ \
                + coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Center: (Bottom) */ \
                + 0.25 * art * x[grid.index(i_r+1,i_theta)] ); /* Bottom Right */ \
            /*  - 0.25 * art * x[grid.index(i_r-1,i_theta)]; // Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            } \
        } \
    /* ---------------------------- */ \
    /* Node in the 2nd inner circle */ \
    /* ---------------------------- */ \
    } else if (i_r == 1) { \
        scalar_t h1 = grid.r_dist(i_r-1); \
        scalar_t h2 = grid.r_dist(i_r); \
        scalar_t k1 = grid.theta_dist(i_theta-1); \
        scalar_t k2 = grid.theta_dist(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        /* Fill result(i,j) */ \
        if(GiveToCenter) { \
        result[grid.index(i_r,i_theta)] += factor * ( \
            0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF) * x[grid.index(i_r,i_theta)] /* beta_{i,j} */ \
            - coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Left */ \
            - coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Right */ \
            - coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Bottom */ \
            - coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Top */ \
            /* Center: (Left, Right, Bottom, Top) */ \
            + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) * x[grid.index(i_r,i_theta)] ); \
        } \
        /* Fill result(i-1,j) */ \
        if(!DirBC_Interior){ /* Don't give to the inner dirichlet boundary! */ \
            if(GiveToLeft){ \
            result[grid.index(i_r-1,i_theta)] += factor * ( \
                - coeff1 * arr * x[grid.index(i_r,i_theta)] /* Right */ \
                + coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Center: (Right) */ \
                - 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Right */ \
                + 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Right */ \
            } \
        } \
        /* Fill result(i+1,j) */ \
        if(GiveToRight) { \
        result[grid.index(i_r+1,i_theta)] += factor * ( \
            - coeff2 * arr * x[grid.index(i_r,i_theta)] /* Left */ \
            + coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Center: (Left) */ \
            + 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Left */ \
            - 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Left */ \
        } \
        /* Fill result(i,j-1) */ \
        if(GiveToBottom) { \
        result[grid.index(i_r,i_theta-1)] += factor * ( \
            - coeff3 * att * x[grid.index(i_r,i_theta)] /* Top */ \
            + coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Center: (Top) */ \
            - 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Top Right */ \
            + 0.25 * art * x[grid.index(i_r-1,i_theta)] ); /* Top Left */ \
        } \
        /* Fill result(i,j+1) */ \
        if(GiveToTop) { \
        result[grid.index(i_r,i_theta+1)] += factor * ( \
            - coeff4 * att * x[grid.index(i_r,i_theta)] /* Bottom */ \
            + coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Center: (Bottom) */ \
            + 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Bottom Right */ \
            - 0.25 * art * x[grid.index(i_r-1,i_theta)] ); /* Bottom Left */ \
        } \
    /* ---------------------------- */ \
    /* Node in the 2nd outer circle */ \
    /* ---------------------------- */ \
    } else if (i_r == grid.nr() - 2) { \
        scalar_t h1 = grid.r_dist(i_r-1); \
        scalar_t h2 = grid.r_dist(i_r); \
        scalar_t k1 = grid.theta_dist(i_theta-1); \
        scalar_t k2 = grid.theta_dist(i_theta); \
        scalar_t coeff1 = 0.5*(k1+k2)/h1; \
        scalar_t coeff2 = 0.5*(k1+k2)/h2; \
        scalar_t coeff3 = 0.5*(h1+h2)/k1; \
        scalar_t coeff4 = 0.5*(h1+h2)/k2; \
        /* Fill result(i,j) */ \
        if(GiveToCenter) { \
        result[grid.index(i_r,i_theta)] += factor * ( \
            (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF) / 4 * x[grid.index(i_r,i_theta)] /* beta_{i,j} */ \
            - coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Left */ \
            - coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Right */ \
            - coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Bottom */ \
            - coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Top */ \
            /* Center: (Left, Right, Bottom, Top) */ \
            + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) * x[grid.index(i_r,i_theta)] ); \
        } \
        /* Fill result(i-1,j) */ \
        if(GiveToLeft){ \
        result[grid.index(i_r-1,i_theta)] += factor * ( \
            - coeff1 * arr * x[grid.index(i_r,i_theta)] /* Right */ \
            + coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Center: (Right) */ \
            - 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Right */ \
            + 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Right */ \
        } \
        /* Don't give to the outer dirichlet boundary! */ \
        /* Fill result(i+1,j) */ \
        /* if(GiveToRight) { */ \
        /* result[grid.index(i_r+1,i_theta)] += factor * ( */ \
        /*     - coeff2 * arr * x[grid.index(i_r,i_theta)] // Left */ \
        /*     + coeff2 * arr * x[grid.index(i_r+1,i_theta)] // Center: (Left) */ \
        /*     + 0.25 * art * x[grid.index(i_r,i_theta+1)] // Top Left */ \
        /*     - 0.25 * art * x[grid.index(i_r,i_theta-1)] ); // Bottom Left */ \
        /* } */ \
        /* Fill result(i,j-1) */ \
        if(GiveToBottom) { \
        result[grid.index(i_r,i_theta-1)] += factor * ( \
            - coeff3 * att * x[grid.index(i_r,i_theta)] /* Top */ \
            + coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Center: (Top) */ \
            - 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Top Right */ \
            + 0.25 * art * x[grid.index(i_r-1,i_theta)] ); /* Top Left */ \
        } \
        /* Fill result(i,j+1) */ \
        if(GiveToTop) { \
        result[grid.index(i_r,i_theta+1)] += factor * ( \
            - coeff4 * att * x[grid.index(i_r,i_theta)] /* Bottom */ \
            + coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Center: (Bottom) */ \
            + 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Bottom Right */ \
            - 0.25 * art * x[grid.index(i_r-1,i_theta)] ); /* Bottom Left */ \
        } \
    /* ------------------------ */ \
    /* Node in the outer circle */ \
    /* ------------------------ */ \
    } else if (i_r == grid.nr() - 1) { \
        /* Dirichlet boundary */ \
        if(GiveToCenter) { \
        result[grid.index(i_r,i_theta)] += factor * x[grid.index(i_r,i_theta)]; \
        } \
        /* Give value to the interior nodes! */ \
        scalar_t h1 = grid.r_dist(i_r-1); \
        scalar_t k1 = grid.theta_dist(i_theta-1); \
        scalar_t k2 = grid.theta_dist(i_theta); \
        scalar_t coeff1 = 0.5*(k1+k2)/h1; \
        /* Fill result(i-1,j) */ \
        if(GiveToLeft){ \
        result[grid.index(i_r-1,i_theta)] += factor * ( \
            - coeff1 * arr * x[grid.index(i_r,i_theta)] /* Right */ \
            + coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Center: (Right) */ \
            - 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Right */ \
            + 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Right */ \
        } \
    } \
} while(0)


int ptr_nz_index_matrixA(const PolarGrid& grid, const int i_r, const int i_theta, const int DirBC_Interior){
    const int inner_circle_stencil = DirBC_Interior ? 1 : 7;
    const int interior_stencil = 9;
    const int outer_circle_stencil = 1;
    if(0 < i_r && i_r < grid.numberSmootherCircles()){
        // Interior. Circle index section
        const int number_prior_inner_boundary_nodes = grid.ntheta();
        const int number_prior_interior_circle_nodes = (i_r-1) * grid.ntheta() + i_theta;
        return inner_circle_stencil * number_prior_inner_boundary_nodes + 
            interior_stencil * number_prior_interior_circle_nodes;
    } else if(grid.numberSmootherCircles() <= i_r && i_r < grid.nr()-1){
        // Interior. Radial index section
        const int number_prior_inner_boundary_nodes = grid.ntheta();
        const int number_prior_interior_circle_nodes = grid.ntheta() * (grid.numberSmootherCircles()-1);
        const int number_prior_interior_radial_nodes = i_theta * (grid.lengthSmootherRadial()-1) + i_r - grid.numberSmootherCircles();
        const int number_prior_outer_boundary_nodes = i_theta;
        return inner_circle_stencil * number_prior_inner_boundary_nodes + 
            interior_stencil * number_prior_interior_circle_nodes + 
            interior_stencil * number_prior_interior_radial_nodes + 
            outer_circle_stencil * number_prior_outer_boundary_nodes;
    } else if(i_r == 0){
        // Boundary on the interior, inner circle.
        const int number_prior_inner_boundary_nodes = i_theta;
        return inner_circle_stencil * number_prior_inner_boundary_nodes;
    } else{
        // Boundary on the outside, outer circle.
        const int number_prior_inner_boundary_nodes = grid.ntheta();
        const int number_prior_interior_circle_nodes = grid.ntheta() * (grid.numberSmootherCircles()-1);
        const int number_prior_interior_radial_nodes = (i_theta+1) * (grid.lengthSmootherRadial()-1);
        const int number_prior_outer_boundary_nodes = i_theta;
        return inner_circle_stencil * number_prior_inner_boundary_nodes + 
            interior_stencil * number_prior_interior_circle_nodes + 
            interior_stencil * number_prior_interior_radial_nodes + 
            outer_circle_stencil * number_prior_outer_boundary_nodes;
    }
}


#include <array>
#include <stdexcept>

#include <array>
#include <stdexcept>
#include <initializer_list>
#include <algorithm>

enum class StencilType
{
    TopLeft,
    Top,
    TopRight,
    Left,
    Center,
    Right,
    BottomLeft,
    Bottom,
    BottomRight,
};

struct Stencil {
    std::array<int, 9> values;

    Stencil(std::initializer_list<int> init) : values{} {
        std::copy(init.begin(), init.end(), values.begin());
    }

    int operator[](StencilType type) const {
        return values[static_cast<size_t>(type)];
    }
};

const Stencil& get_stencil(const PolarGrid& grid, int i_r, int DirBC_Interior) {
    static const Stencil stencil_DB = 
        {-1, -1, -1,
         -1,  0, -1,
         -1, -1, -1};
         
    static const Stencil stencil_across = 
        {-1, 4, 6,
          1, 0, 2,
         -1, 3, 5};
         
    static const Stencil stencil_interior = 
        {7, 4, 8,
         1, 0, 2,
         5, 3, 6};

    if (i_r > 0 && i_r < grid.nr() - 1) {
        return stencil_interior;
    } else if (i_r == grid.nr() - 1 || (i_r == 0 && DirBC_Interior)) {
        return stencil_DB;
    } else if (i_r == 0 && !DirBC_Interior) {
        return stencil_across;
    }
    
    throw std::out_of_range("Invalid index for stencil");
}

int nnz_matrixA(const PolarGrid& grid, const int DirBC_Interior){
    const int inner_circle_stencil = DirBC_Interior ? 1 : 7;
    const int interior_stencil = 9;
    const int outer_circle_stencil = 1;
    return inner_circle_stencil * grid.ntheta() + 
        interior_stencil * (grid.nr()-2) * grid.ntheta() +
        outer_circle_stencil * grid.ntheta();
}

std::pair<SparseMatrix<scalar_t>, Vector<scalar_t>> Operator::build_system(const Level& onLevel) const{
    const PolarGrid& grid = onLevel.grid();

    const int n = grid.number_of_nodes();
    const int matrixA_nnz = nnz_matrixA(grid, DirBC_Interior_);

    SparseMatrix<scalar_t> matrixA(n, n, matrixA_nnz);
    Vector<scalar_t> rhs(n);

    #pragma omp parallel for
    for (int i = 0; i < matrixA_nnz; i++){
        matrixA.value(i) = 0.0;
    }


    // ---------------------- //
    // OpenMP Parallelization //
    // ---------------------- //
    const int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);
    const int minimalChunkSize = 4;
    const int zone = 2;

    // Distribute Tasks to each thread
    TaskDistribution CircleSmootherTasks(grid.numberSmootherCircles(), minimalChunkSize, numThreads);
    TaskDistribution RadialSmootherTasks(grid.ntheta(), minimalChunkSize, numThreads);

    // The acroos origin neighbor nodes need to be in the same TaskDistribution to prevent race condition.
    assert(grid.numberSmootherCircles() >= 1);

    #pragma omp parallel num_threads(numThreads)
    {   
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

        const int threadID = omp_get_thread_num();
        // ---------------------------------------------------------- //
        // Take care of the separation strips of the circular smoother //
        // ---------------------------------------------------------- //
        const int i_r_start = CircleSmootherTasks.getStart(threadID);
        const int i_r_end = CircleSmootherTasks.getEnd(threadID);
        const int i_r_separate = std::min(i_r_end - i_r_start, zone);

        // For loop matches circular access pattern
        for (int i_r = i_r_end - i_r_separate; i_r < i_r_end; i_r++){
            r = grid.radius(i_r);
            coeff_alpha = (*alpha_)(r);
            coeff_beta = (*beta_)(r);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
                theta = grid.theta(i_theta);
                sin_theta = sin_theta_[i_theta];
                cos_theta = cos_theta_[i_theta];
                arr_att_art(r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDF);
                // center_nz_index = ptr_nz_index_matrixA(grid, i_r, i_theta, DirBC_Interior_);







            }
        }
        
        #pragma omp barrier

        // -------------------------------------------------------- //
        // Take care of the separation strips of the radial smoother //
        // -------------------------------------------------------- //
        const int i_theta_start = RadialSmootherTasks.getStart(threadID);
        const int i_theta_end = RadialSmootherTasks.getEnd(threadID);
        const int i_theta_seperate = std::min(i_theta_end-i_theta_start, zone);

        // For loop matches radial access pattern
        for (int i_theta = i_theta_start; i_theta < i_theta_start + i_theta_seperate; i_theta++){
            theta = grid.theta(i_theta);
            sin_theta = sin_theta_[i_theta];
            cos_theta = cos_theta_[i_theta];
            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                r = grid.radius(i_r);
                coeff_alpha = (*alpha_)(r);
                coeff_beta = (*beta_)(r);
                arr_att_art(r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDF);
                // center_nz_index = ptr_nz_index_matrixA(grid, i_r, i_theta, DirBC_Interior_);


                // NODE_APPLY_A_GIVE(i_r, i_theta, grid, DirBC_Interior_, result, x, scaleAx, 
                //     arr, att, art, coeff_beta, detDF, 
                //     true, true, true, true, true);
            }
        }

        #pragma omp barrier

        // ------------------------------------------ //
        // Take care of the circular smoother section //
        // ------------------------------------------ //
        // For loop matches circular access pattern
        for (int i_r = i_r_start; i_r < i_r_end - i_r_separate; i_r++){
            r = grid.radius(i_r);
            coeff_alpha = (*alpha_)(r);
            coeff_beta = (*beta_)(r);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
                theta = grid.theta(i_theta);
                sin_theta = sin_theta_[i_theta];
                cos_theta = cos_theta_[i_theta];
                arr_att_art(r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDF);

                /* -------------------- */ 
                /* Node in the interior */ 
                /* -------------------- */ 
                if (i_r > 1 && i_r < grid.nr() - 2) { 
                    scalar_t h1 = grid.r_dist(i_r-1); 
                    scalar_t h2 = grid.r_dist(i_r); 
                    scalar_t k1 = grid.theta_dist(i_theta-1); 
                    scalar_t k2 = grid.theta_dist(i_theta); 
                    scalar_t coeff1 = 0.5*(k1+k2)/h1; 
                    scalar_t coeff2 = 0.5*(k1+k2)/h2; 
                    scalar_t coeff3 = 0.5*(h1+h2)/k1; 
                    scalar_t coeff4 = 0.5*(h1+h2)/k2;

                    center_nz_index = ptr_nz_index_matrixA(grid, i_r, i_theta, DirBC_Interior_);
                    left_nz_index = ptr_nz_index_matrixA(grid, i_r-1, i_theta, DirBC_Interior_);
                    right_nz_index = ptr_nz_index_matrixA(grid, i_r+1, i_theta, DirBC_Interior_);
                    bottom_nz_index = ptr_nz_index_matrixA(grid, i_r, ((i_theta-1) + grid.ntheta()) % grid.ntheta(), DirBC_Interior_);
                    top_nz_index = ptr_nz_index_matrixA(grid, i_r, ((i_theta+1) + grid.ntheta()) % grid.ntheta(), DirBC_Interior_);

                    center_index = grid.index(i_r,i_theta);
                    left_index = grid.index(i_r-1,i_theta);
                    right_index = grid.index(i_r+1,i_theta);
                    bottom_index = grid.index(i_r,i_theta-1);
                    top_index = grid.index(i_r,i_theta+1);

                    /* Fill result(i,j) */
                    const Stencil& CenterStencil = get_stencil(grid, i_r, DirBC_Interior_);

                    nz_index = center_nz_index + CenterStencil[StencilType::Center];
                    matrixA.row_index(nz_index) = center_index;
                    matrixA.col_index(nz_index) = center_index;
                    matrixA.value(nz_index) += 0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF); /* beta_{i,j} */ 

                    nz_index = center_nz_index + CenterStencil[StencilType::Left];
                    matrixA.row_index(nz_index) = center_index;
                    matrixA.col_index(nz_index) = left_index;
                    matrixA.value(nz_index) += - coeff1 * arr; /* Left */

                    nz_index = center_nz_index + CenterStencil[StencilType::Right];
                    matrixA.row_index(nz_index) = center_index;
                    matrixA.col_index(nz_index) = right_index;
                    matrixA.value(nz_index) += - coeff2 * arr; /* Right */

                    nz_index = center_nz_index + CenterStencil[StencilType::Bottom];
                    matrixA.row_index(nz_index) = center_index;
                    matrixA.col_index(nz_index) = bottom_index;
                    matrixA.value(nz_index) += - coeff3 * att; /* Bottom */

                    nz_index = center_nz_index + CenterStencil[StencilType::Top];
                    matrixA.row_index(nz_index) = center_index;
                    matrixA.col_index(nz_index) = top_index;
                    matrixA.value(nz_index) += - coeff4 * att; /* Top */

                    nz_index = center_nz_index + CenterStencil[StencilType::Center];
                    matrixA.row_index(nz_index) = center_index;
                    matrixA.col_index(nz_index) = center_index;
                    /* Center: (Left, Right, Bottom, Top) */
                    matrixA.value(nz_index) += (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att;

                    /* Fill result(i-1,j) */ 
                    const Stencil& LeftStencil = get_stencil(grid, i_r-1, DirBC_Interior_);

                    nz_index = left_nz_index + LeftStencil[StencilType::Right];
                    matrixA.row_index(nz_index) = left_index;
                    matrixA.col_index(nz_index) = center_index;
                    matrixA.value(nz_index) += - coeff1 * arr; /* Right */

                    nz_index = left_nz_index + LeftStencil[StencilType::Center];
                    matrixA.row_index(nz_index) = left_index;
                    matrixA.col_index(nz_index) = left_index;
                    matrixA.value(nz_index) += coeff1 * arr; /* Center: (Right) */

                    nz_index = left_nz_index + LeftStencil[StencilType::TopRight];
                    matrixA.row_index(nz_index) = left_index;
                    matrixA.col_index(nz_index) = top_index;
                    matrixA.value(nz_index) += - 0.25 * art; /* Top Right */

                    nz_index = left_nz_index + LeftStencil[StencilType::BottomRight];
                    matrixA.row_index(nz_index) = left_index;
                    matrixA.col_index(nz_index) = bottom_index;
                    matrixA.value(nz_index) += 0.25 * art; /* Bottom Right */ 

                    /* Fill result(i+1,j) */ 
                    const Stencil& RightStencil = get_stencil(grid, i_r+1, DirBC_Interior_);

                    nz_index = right_nz_index + RightStencil[StencilType::Left];
                    matrixA.row_index(nz_index) = right_index;
                    matrixA.col_index(nz_index) = center_index;
                    matrixA.value(nz_index) += - coeff2 * arr; /* Left */

                    nz_index = right_nz_index + RightStencil[StencilType::Center];
                    matrixA.row_index(nz_index) = right_index;
                    matrixA.col_index(nz_index) = right_index;
                    matrixA.value(nz_index) += coeff2 * arr; /* Center: (Left) */

                    nz_index = right_nz_index + RightStencil[StencilType::TopLeft];
                    matrixA.row_index(nz_index) = right_index;
                    matrixA.col_index(nz_index) = top_index;
                    matrixA.value(nz_index) += 0.25 * art; /* Top Left */

                    nz_index = right_nz_index + RightStencil[StencilType::BottomLeft];
                    matrixA.row_index(nz_index) = right_index;
                    matrixA.col_index(nz_index) = bottom_index;
                    matrixA.value(nz_index) += - 0.25 * art; /* Bottom Left */

                    /* Fill result(i,j-1) */
                    const Stencil& BottomStencil = CenterStencil;

                    nz_index = bottom_nz_index + BottomStencil[StencilType::Top];
                    matrixA.row_index(nz_index) = bottom_index;
                    matrixA.col_index(nz_index) = center_index;
                    matrixA.value(nz_index) += - coeff3 * att; /* Top */

                    nz_index = bottom_nz_index + BottomStencil[StencilType::Center];
                    matrixA.row_index(nz_index) = bottom_index;
                    matrixA.col_index(nz_index) = bottom_index;
                    matrixA.value(nz_index) += coeff3 * att; /* Center: (Top) */

                    nz_index = bottom_nz_index + BottomStencil[StencilType::TopRight];
                    matrixA.row_index(nz_index) = bottom_index;
                    matrixA.col_index(nz_index) = right_index;
                    matrixA.value(nz_index) += - 0.25 * art; /* Top Right */

                    nz_index = bottom_nz_index + BottomStencil[StencilType::TopLeft];
                    matrixA.row_index(nz_index) = bottom_index;
                    matrixA.col_index(nz_index) = left_index;
                    matrixA.value(nz_index) += 0.25 * art; /* Top Left */

                    /* Fill result(i,j+1) */
                    const Stencil& TopStencil = CenterStencil;

                    nz_index = top_nz_index + TopStencil[StencilType::Bottom];
                    matrixA.row_index(nz_index) = top_index;
                    matrixA.col_index(nz_index) = center_index;
                    matrixA.value(nz_index) += - coeff4 * att; /* Bottom */

                    nz_index = top_nz_index + TopStencil[StencilType::Bottom];
                    matrixA.row_index(nz_index) = top_index;
                    matrixA.col_index(nz_index) = center_index;
                    matrixA.value(nz_index) += coeff4 * att; /* Center: (Bottom) */

                    nz_index = top_nz_index + TopStencil[StencilType::BottomRight];
                    matrixA.row_index(nz_index) = top_index;
                    matrixA.col_index(nz_index) = right_index;
                    matrixA.value(nz_index) += 0.25 * art; /* Bottom Right */

                    nz_index = top_nz_index + TopStencil[StencilType::BottomLeft];
                    matrixA.row_index(nz_index) = top_index;
                    matrixA.col_index(nz_index) = left_index;
                    matrixA.value(nz_index) += - 0.25 * art; /* Bottom Left */


                /* ------------------------ */ \
                /* Node in the inner circle */ \
                /* ------------------------ */ \
                } 

            }
        }

        // ---------------------------------------- //
        // Take care of the radial smoother section //
        // ---------------------------------------- //
        // For loop matches radial access pattern
        for (int i_theta = i_theta_start + i_theta_seperate; i_theta < i_theta_end; i_theta++){
            theta = grid.theta(i_theta);
            sin_theta = sin_theta_[i_theta];
            cos_theta = cos_theta_[i_theta];
            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                r = grid.radius(i_r);
                coeff_alpha = (*alpha_)(r);
                coeff_beta = (*beta_)(r);
                arr_att_art(r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDF);
                // center_nz_index = ptr_nz_index_matrixA(grid, i_r, i_theta, DirBC_Interior_);
                // NODE_APPLY_A_GIVE(i_r, i_theta, grid, DirBC_Interior_, result, x, scaleAx, 
                //     arr, att, art, coeff_beta, detDF, 
                //     true, true, true, true, true);
            }
        }
    }



    std::cout<<matrixA<<std::endl;

    return std::make_pair(std::move(matrixA), std::move(rhs));
}