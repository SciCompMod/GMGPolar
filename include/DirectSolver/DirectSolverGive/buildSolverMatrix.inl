#pragma once

namespace direct_solver_give
{

#ifdef GMGPOLAR_USE_MUMPS
// When using the MUMPS solver, the matrix is assembled in COO format.
static KOKKOS_INLINE_FUNCTION void updateMatrixElement(const SparseMatrixCOO<double, Kokkos::HostSpace>& matrix,
                                                       int ptr, const int offset, const int row, const int column,
                                                       const double value)
{
    matrix.set_row_index(ptr + offset, row);
    matrix.set_col_index(ptr + offset, column);
    matrix.increase_value(ptr + offset, value);
}
#else
// When using the in-house solver, the matrix is stored in CSR format.
static KOKKOS_INLINE_FUNCTION void updateMatrixElement(const SparseMatrixCSR<double, Kokkos::HostSpace>& matrix,
                                                       int ptr, const int offset, const int row, const int column,
                                                       const double value)
{
    matrix.set_row_nz_index(row, offset, column);
    matrix.increase_row_nz_entry(row, offset, value);
}
#endif

template <typename LevelCacheType, typename SystemMatrix>
static KOKKOS_INLINE_FUNCTION void
nodeBuildSolverMatrixGive(const int i_r, const int i_theta, const PolarGrid& grid, const LevelCacheType& level_cache,
                          const bool DirBC_Interior, const SystemMatrix& solver_matrix)
{
    using direct_solver_give::getSolverMatrixIndex;
    using direct_solver_give::getStencil;
    using direct_solver_give::updateMatrixElement;

    /* ---------------------------------------- */
    /* Compute or retrieve stencil coefficients */
    /* ---------------------------------------- */
    const int center    = grid.index(i_r, i_theta);
    const double radius = grid.radius(i_r);
    const double theta  = grid.theta(i_theta);

    double coeff_beta, arr, att, art, detDF;
    level_cache.obtainValues(i_r, i_theta, center, radius, theta, coeff_beta, arr, att, art, detDF);

    int ptr, offset;
    int row, column;
    double value;
    /* -------------------- */
    /* Node in the interior */
    /* -------------------- */
    if (i_r > 1 && i_r < grid.nr() - 2) {
        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const double h1 = grid.radialSpacing(i_r - 1);
        const double h2 = grid.radialSpacing(i_r);
        const double k1 = grid.angularSpacing(i_theta_M1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff1 = 0.5 * (k1 + k2) / h1;
        const double coeff2 = 0.5 * (k1 + k2) / h2;
        const double coeff3 = 0.5 * (h1 + h2) / k1;
        const double coeff4 = 0.5 * (h1 + h2) / k2;
        const double coeff5 = 0.25 * (h1 + h2) * (k1 + k2);

        const int center_nz_index = getSolverMatrixIndex(i_r, i_theta, grid, DirBC_Interior);
        const int left_nz_index   = getSolverMatrixIndex(i_r - 1, i_theta, grid, DirBC_Interior);
        const int right_nz_index  = getSolverMatrixIndex(i_r + 1, i_theta, grid, DirBC_Interior);
        const int bottom_nz_index = getSolverMatrixIndex(i_r, i_theta_M1, grid, DirBC_Interior);
        const int top_nz_index    = getSolverMatrixIndex(i_r, i_theta_P1, grid, DirBC_Interior);

        const int center_index = grid.index(i_r, i_theta);
        const int left_index   = grid.index(i_r - 1, i_theta);
        const int right_index  = grid.index(i_r + 1, i_theta);
        const int bottom_index = grid.index(i_r, i_theta_M1);
        const int top_index    = grid.index(i_r, i_theta_P1);

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r, grid, DirBC_Interior);

        offset = CenterStencil[StencilPosition::Center];
        column = center_index;
        value  = coeff5 * coeff_beta * Kokkos::fabs(detDF); /* beta_{i,j} */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Left];
        column = left_index;
        value  = -coeff1 * arr; /* Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Right];
        column = right_index;
        value  = -coeff2 * arr; /* Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Bottom];
        column = bottom_index;
        value  = -coeff3 * att; /* Bottom */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Top];
        column = top_index;
        value  = -coeff4 * att; /* Top */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Center];
        column = center_index;
        value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* Fill matrix row of (i-1,j) */
        row = left_index;
        ptr = left_nz_index;

        const Stencil& LeftStencil = getStencil(i_r - 1, grid, DirBC_Interior);

        offset = LeftStencil[StencilPosition::Right];
        column = center_index;
        value  = -coeff1 * arr; /* Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = LeftStencil[StencilPosition::Center];
        column = left_index;
        value  = +coeff1 * arr; /* Center: (Right) */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = LeftStencil[StencilPosition::TopRight];
        column = top_index;
        value  = -0.25 * art; /* Top Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = LeftStencil[StencilPosition::BottomRight];
        column = bottom_index;
        value  = +0.25 * art; /* Bottom Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* Fill matrix row of (i+1,j) */
        row = right_index;
        ptr = right_nz_index;

        const Stencil& RightStencil = getStencil(i_r + 1, grid, DirBC_Interior);

        offset = RightStencil[StencilPosition::Left];
        column = center_index;
        value  = -coeff2 * arr; /* Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = RightStencil[StencilPosition::Center];
        column = right_index;
        value  = +coeff2 * arr; /* Center: (Left) */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = RightStencil[StencilPosition::TopLeft];
        column = top_index;
        value  = +0.25 * art; /* Top Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = RightStencil[StencilPosition::BottomLeft];
        column = bottom_index;
        value  = -0.25 * art; /* Bottom Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* Fill matrix row of (i,j-1) */
        row = bottom_index;
        ptr = bottom_nz_index;

        const Stencil& BottomStencil = CenterStencil;

        offset = BottomStencil[StencilPosition::Top];
        column = center_index;
        value  = -coeff3 * att; /* Top */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = BottomStencil[StencilPosition::Center];
        column = bottom_index;
        value  = +coeff3 * att; /* Center: (Top) */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = BottomStencil[StencilPosition::TopRight];
        column = right_index;
        value  = -0.25 * art; /* Top Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = BottomStencil[StencilPosition::TopLeft];
        column = left_index;
        value  = +0.25 * art; /* Top Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* Fill matrix row of (i,j+1) */
        row = top_index;
        ptr = top_nz_index;

        const Stencil& TopStencil = CenterStencil;

        offset = TopStencil[StencilPosition::Bottom];
        column = center_index;
        value  = -coeff4 * att; /* Bottom */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = TopStencil[StencilPosition::Center];
        column = top_index;
        value  = +coeff4 * att; /* Center: (Bottom) */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = TopStencil[StencilPosition::BottomRight];
        column = right_index;
        value  = +0.25 * art; /* Bottom Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = TopStencil[StencilPosition::BottomLeft];
        column = left_index;
        value  = -0.25 * art; /* Bottom Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);
    }
    /* -------------------------- */
    /* Node on the inner boundary */
    /* -------------------------- */
    else if (i_r == 0) {
        /* ------------------------------------------------ */
        /* Case 1: Dirichlet boundary on the inner boundary */
        /* ------------------------------------------------ */
        if (DirBC_Interior) {
            const double h2     = grid.radialSpacing(i_r);
            const double k1     = grid.angularSpacing(i_theta - 1);
            const double k2     = grid.angularSpacing(i_theta);
            const double coeff2 = 0.5 * (k1 + k2) / h2;

            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

            const int center_nz_index = getSolverMatrixIndex(i_r, i_theta, grid, DirBC_Interior);
            const int right_nz_index  = getSolverMatrixIndex(i_r + 1, i_theta, grid, DirBC_Interior);

            const int center_index = grid.index(i_r, i_theta);
            const int right_index  = grid.index(i_r + 1, i_theta);
            const int bottom_index = grid.index(i_r, i_theta_M1);
            const int top_index    = grid.index(i_r, i_theta_P1);

            /* Fill matrix row of (i,j) */
            row = center_index;
            ptr = center_nz_index;

            const Stencil& CenterStencil = getStencil(i_r, grid, DirBC_Interior);

            offset = CenterStencil[StencilPosition::Center];
            column = center_index;
            value  = 1.0;
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            /* Fill matrix row of (i+1,j) */
            row = right_index;
            ptr = right_nz_index;

            const Stencil& RightStencil = getStencil(i_r + 1, grid, DirBC_Interior);

            /* Left REMOVED: Moved to the right hand side to make the matrix symmetric */

            offset = RightStencil[StencilPosition::Center];
            column = right_index;
            value  = +coeff2 * arr; /* Center: (Left) */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            /* TopLeft REMOVED: Moved to the right hand side to make the matrix symmetric */

            /* BottomLeft REMOVED: Moved to the right hand side to make the matrix symmetric */
        }
        else {
            /* ------------------------------------------------------------- */
            /* Case 2: Across origin discretization on the interior boundary */
            /* ------------------------------------------------------------- */
            /* h1 gets replaced with 2 * R0. */
            /* (i_r-1,i_theta) gets replaced with (i_r, i_theta + grid.ntheta()/2). */
            /* Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil. */
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

            assert(grid.ntheta() % 2 == 0);
            const int i_theta_AcrossOrigin = grid.wrapThetaIndex(i_theta + grid.ntheta() / 2);

            double h1 = 2.0 * grid.radius(0);
            double h2 = grid.radialSpacing(i_r);
            double k1 = grid.angularSpacing(i_theta_M1);
            double k2 = grid.angularSpacing(i_theta);

            double coeff1       = 0.5 * (k1 + k2) / h1;
            double coeff2       = 0.5 * (k1 + k2) / h2;
            double coeff3       = 0.5 * (h1 + h2) / k1;
            double coeff4       = 0.5 * (h1 + h2) / k2;
            const double coeff5 = 0.25 * (h1 + h2) * (k1 + k2);

            const int center_nz_index = getSolverMatrixIndex(i_r, i_theta, grid, DirBC_Interior);
            const int left_nz_index   = getSolverMatrixIndex(i_r, i_theta_AcrossOrigin, grid, DirBC_Interior);
            const int right_nz_index  = getSolverMatrixIndex(i_r + 1, i_theta, grid, DirBC_Interior);
            const int bottom_nz_index = getSolverMatrixIndex(i_r, i_theta_M1, grid, DirBC_Interior);
            const int top_nz_index    = getSolverMatrixIndex(i_r, i_theta_P1, grid, DirBC_Interior);

            const int center_index = grid.index(i_r, i_theta);
            const int left_index   = grid.index(i_r, i_theta_AcrossOrigin);
            const int right_index  = grid.index(i_r + 1, i_theta);
            const int bottom_index = grid.index(i_r, i_theta_M1);
            const int top_index    = grid.index(i_r, i_theta_P1);

            /* Fill matrix row of (i,j) */
            row                          = center_index;
            ptr                          = center_nz_index;
            const Stencil& CenterStencil = getStencil(i_r, grid, DirBC_Interior);

            offset = CenterStencil[StencilPosition::Center];
            column = center_index;
            value  = coeff5 * coeff_beta * Kokkos::fabs(detDF); /* beta_{i,j} */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = CenterStencil[StencilPosition::Left];
            column = left_index;
            value  = -coeff1 * arr; /* Left */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = CenterStencil[StencilPosition::Right];
            column = right_index;
            value  = -coeff2 * arr; /* Right */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = CenterStencil[StencilPosition::Bottom];
            column = bottom_index;
            value  = -coeff3 * att; /* Bottom */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = CenterStencil[StencilPosition::Top];
            column = top_index;
            value  = -coeff4 * att; /* Top */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = CenterStencil[StencilPosition::Center];
            column = center_index;
            value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            /* Fill matrix row of (i-1,j) */
            /* From view the view of the across origin node, */
            /* the directions are roatated by 180 degrees in the stencil! */
            row = left_index;
            ptr = left_nz_index;

            const Stencil& LeftStencil = CenterStencil;

            offset = LeftStencil[StencilPosition::Left];
            column = center_index;
            value  = -coeff1 * arr; /* Right -> Left*/
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = LeftStencil[StencilPosition::Center];
            column = left_index;
            value  = +coeff1 * arr; /* Center: (Right) -> Center: (Left) */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            /* Top Right -> Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */

            /* Bottom Right -> Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */

            /* Fill matrix row of (i+1,j) */
            row                         = right_index;
            ptr                         = right_nz_index;
            const Stencil& RightStencil = getStencil(i_r + 1, grid, DirBC_Interior);

            offset = RightStencil[StencilPosition::Left];
            column = center_index;
            value  = -coeff2 * arr; /* Left */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = RightStencil[StencilPosition::Center];
            column = right_index;
            value  = +coeff2 * arr; /* Center: (Left) */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = RightStencil[StencilPosition::TopLeft];
            column = top_index;
            value  = +0.25 * art; /* Top Left */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = RightStencil[StencilPosition::BottomLeft];
            column = bottom_index;
            value  = -0.25 * art; /* Bottom Left */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            /* Fill matrix row of (i,j-1) */
            row = bottom_index;
            ptr = bottom_nz_index;

            const Stencil& BottomStencil = CenterStencil;

            offset = BottomStencil[StencilPosition::Top];
            column = center_index;
            value  = -coeff3 * att; /* Top */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = BottomStencil[StencilPosition::Center];
            column = bottom_index;
            value  = +coeff3 * att; /* Center: (Top) */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = BottomStencil[StencilPosition::TopRight];
            column = right_index;
            value  = -0.25 * art; /* Top Right */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            /* TopLeft REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */

            /* Fill matrix row of (i,j+1) */
            row = top_index;
            ptr = top_nz_index;

            const Stencil& TopStencil = CenterStencil;

            offset = TopStencil[StencilPosition::Bottom];
            column = center_index;
            value  = -coeff4 * att; /* Bottom */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = TopStencil[StencilPosition::Center];
            column = top_index;
            value  = +coeff4 * att; /* Center: (Bottom) */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = TopStencil[StencilPosition::BottomRight];
            column = right_index;
            value  = +0.25 * art; /* Bottom Right */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            /* BottomLeft REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
        }
    }
    /* ------------------------------- */
    /* Node next to the inner boundary */
    /* ------------------------------- */
    else if (i_r == 1) {
        const double h1 = grid.radialSpacing(i_r - 1);
        const double h2 = grid.radialSpacing(i_r);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff1 = 0.5 * (k1 + k2) / h1;
        const double coeff2 = 0.5 * (k1 + k2) / h2;
        const double coeff3 = 0.5 * (h1 + h2) / k1;
        const double coeff4 = 0.5 * (h1 + h2) / k2;
        const double coeff5 = 0.25 * (h1 + h2) * (k1 + k2);

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int center_nz_index = getSolverMatrixIndex(i_r, i_theta, grid, DirBC_Interior);
        const int left_nz_index   = getSolverMatrixIndex(i_r - 1, i_theta, grid, DirBC_Interior);
        const int right_nz_index  = getSolverMatrixIndex(i_r + 1, i_theta, grid, DirBC_Interior);
        const int bottom_nz_index = getSolverMatrixIndex(i_r, i_theta_M1, grid, DirBC_Interior);
        const int top_nz_index    = getSolverMatrixIndex(i_r, i_theta_P1, grid, DirBC_Interior);

        const int center_index = grid.index(i_r, i_theta);
        const int left_index   = grid.index(i_r - 1, i_theta);
        const int right_index  = grid.index(i_r + 1, i_theta);
        const int bottom_index = grid.index(i_r, i_theta_M1);
        const int top_index    = grid.index(i_r, i_theta_P1);

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r, grid, DirBC_Interior);

        offset = CenterStencil[StencilPosition::Center];
        column = center_index;
        value  = coeff5 * coeff_beta * Kokkos::fabs(detDF); /* beta_{i,j} */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* REMOVED: Moved to the right hand side to make the matrix symmetric */
        if (!DirBC_Interior) {
            offset = CenterStencil[StencilPosition::Left];
            column = left_index;
            value  = -coeff1 * arr; /* Left */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);
        }

        offset = CenterStencil[StencilPosition::Right];
        column = right_index;
        value  = -coeff2 * arr; /* Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Bottom];
        column = bottom_index;
        value  = -coeff3 * att; /* Bottom */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Top];
        column = top_index;
        value  = -coeff4 * att; /* Top */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Center];
        column = center_index;
        value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        if (!DirBC_Interior) { /* Don't give to the inner Dirichlet boundary! */
            /* Fill matrix row of (i-1,j) */
            row = left_index;
            ptr = left_nz_index;

            const Stencil& LeftStencil = getStencil(i_r - 1, grid, DirBC_interior);

            offset = LeftStencil[StencilPosition::Right];
            column = center_index;
            value  = -coeff1 * arr; /* Right */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = LeftStencil[StencilPosition::Center];
            column = left_index;
            value  = +coeff1 * arr; /* Center: (Right) */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = LeftStencil[StencilPosition::TopRight];
            column = top_index;
            value  = -0.25 * art; /* Top Right */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = LeftStencil[StencilPosition::BottomRight];
            column = bottom_index;
            value  = +0.25 * art; /* Bottom Right */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);
        }

        /* Fill matrix row of (i+1,j) */
        row = right_index;
        ptr = right_nz_index;

        const Stencil& RightStencil = getStencil(i_r + 1, grid, DirBC_interior);

        offset = RightStencil[StencilPosition::Left];
        column = center_index;
        value  = -coeff2 * arr; /* Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = RightStencil[StencilPosition::Center];
        column = right_index;
        value  = +coeff2 * arr; /* Center: (Left) */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = RightStencil[StencilPosition::TopLeft];
        column = top_index;
        value  = +0.25 * art; /* Top Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = RightStencil[StencilPosition::BottomLeft];
        column = bottom_index;
        value  = -0.25 * art; /* Bottom Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* Fill matrix row of (i,j-1) */
        row = bottom_index;
        ptr = bottom_nz_index;

        const Stencil& BottomStencil = CenterStencil;

        offset = BottomStencil[StencilPosition::Top];
        column = center_index;
        value  = -coeff3 * att; /* Top */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = BottomStencil[StencilPosition::Center];
        column = bottom_index;
        value  = +coeff3 * att; /* Center: (Top) */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = BottomStencil[StencilPosition::TopRight];
        column = right_index;
        value  = -0.25 * art; /* Top Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* REMOVED: Moved to the right hand side to make the matrix symmetric */
        if (!DirBC_Interior) {
            offset = BottomStencil[StencilPosition::TopLeft];
            column = left_index;
            value  = +0.25 * art; /* Top Left */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);
        }

        /* Fill matrix row of (i,j+1) */
        row = top_index;
        ptr = top_nz_index;

        const Stencil& TopStencil = CenterStencil;

        offset = TopStencil[StencilPosition::Bottom];
        column = center_index;
        value  = -coeff4 * att; /* Bottom */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = TopStencil[StencilPosition::Center];
        column = top_index;
        value  = +coeff4 * att; /* Center: (Bottom) */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = TopStencil[StencilPosition::BottomRight];
        column = right_index;
        value  = +0.25 * art; /* Bottom Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* REMOVED: Moved to the right hand side to make the matrix symmetric */
        if (!DirBC_Interior) {
            offset = TopStencil[StencilPosition::BottomLeft];
            column = left_index;
            value  = -0.25 * art; /* Bottom Left */
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);
        }
    }
    /* ------------------------------- */
    /* Node next to the outer boundary */
    /* ------------------------------- */
    else if (i_r == grid.nr() - 2) {
        const double h1 = grid.radialSpacing(i_r - 1);
        const double h2 = grid.radialSpacing(i_r);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff1 = 0.5 * (k1 + k2) / h1;
        const double coeff2 = 0.5 * (k1 + k2) / h2;
        const double coeff3 = 0.5 * (h1 + h2) / k1;
        const double coeff4 = 0.5 * (h1 + h2) / k2;
        const double coeff5 = 0.25 * (h1 + h2) * (k1 + k2);

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int center_nz_index = getSolverMatrixIndex(i_r, i_theta, grid, DirBC_Interior);
        const int left_nz_index   = getSolverMatrixIndex(i_r - 1, i_theta, grid, DirBC_Interior);
        const int right_nz_index  = getSolverMatrixIndex(i_r + 1, i_theta, grid, DirBC_Interior);
        const int bottom_nz_index = getSolverMatrixIndex(i_r, i_theta_M1, grid, DirBC_Interior);
        const int top_nz_index    = getSolverMatrixIndex(i_r, i_theta_P1, grid, DirBC_Interior);

        const int center_index = grid.index(i_r, i_theta);
        const int left_index   = grid.index(i_r - 1, i_theta);
        const int right_index  = grid.index(i_r + 1, i_theta);
        const int bottom_index = grid.index(i_r, i_theta_M1);
        const int top_index    = grid.index(i_r, i_theta_P1);

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r, grid, DirBC_Interior);

        offset = CenterStencil[StencilPosition::Center];
        column = center_index;
        value  = coeff5 * coeff_beta * Kokkos::fabs(detDF); /* beta_{i,j} */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Left];
        column = left_index;
        value  = -coeff1 * arr; /* Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* Right REMOVED: Moved to the right hand side to make the matrix symmetric */

        offset = CenterStencil[StencilPosition::Bottom];
        column = bottom_index;
        value  = -coeff3 * att; /* Bottom */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Top];
        column = top_index;
        value  = -coeff4 * att; /* Top */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Center];
        column = center_index;
        value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* Fill matrix row of (i-1,j) */
        row = left_index;
        ptr = left_nz_index;

        const Stencil& LeftStencil = getStencil(i_r - 1, grid, DirBC_Interior);

        offset = LeftStencil[StencilPosition::Right];
        column = center_index;
        value  = -coeff1 * arr; /* Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = LeftStencil[StencilPosition::Center];
        column = left_index;
        value  = coeff1 * arr; /* Center: (Right) */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = LeftStencil[StencilPosition::TopRight];
        column = top_index;
        value  = -0.25 * art; /* Top Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = LeftStencil[StencilPosition::BottomRight];
        column = bottom_index;
        value  = 0.25 * art; /* Bottom Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* Fill matrix row of (i+1,j) */
        /* Don't give to the outer dirichlet boundary! */

        /* Fill matrix row of (i,j-1) */
        row = bottom_index;
        ptr = bottom_nz_index;

        const Stencil& BottomStencil = CenterStencil;

        offset = BottomStencil[StencilPosition::Top];
        column = center_index;
        value  = -coeff3 * att; /* Top */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = BottomStencil[StencilPosition::Center];
        column = bottom_index;
        value  = coeff3 * att; /* Center: (Top) */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* TopRight REMOVED: Moved to the right hand side to make the matrix symmetric */

        offset = BottomStencil[StencilPosition::TopLeft];
        column = left_index;
        value  = 0.25 * art; /* Top Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* Fill matrix row of (i,j+1) */
        row = top_index;
        ptr = top_nz_index;

        const Stencil& TopStencil = CenterStencil;

        offset = TopStencil[StencilPosition::Bottom];
        column = center_index;
        value  = -coeff4 * att; /* Bottom */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = TopStencil[StencilPosition::Center];
        column = top_index;
        value  = coeff4 * att; /* Center: (Bottom) */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* BottomRight REMOVED: Moved to the right hand side to make the matrix symmetric */

        offset = TopStencil[StencilPosition::BottomLeft];
        column = left_index;
        value  = -0.25 * art; /* Bottom Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);
    }
    /* ------------------------------------ */
    /* Node on the outer dirichlet boundary */
    /* ------------------------------------ */
    else if (i_r == grid.nr() - 1) {
        const double h1 = grid.radialSpacing(i_r - 1);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff1 = 0.5 * (k1 + k2) / h1;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int center_nz_index = getSolverMatrixIndex(i_r, i_theta, grid, DirBC_Interior);
        const int left_nz_index   = getSolverMatrixIndex(i_r - 1, i_theta, grid, DirBC_Interior);

        const int center_index = grid.index(i_r, i_theta);
        const int left_index   = grid.index(i_r - 1, i_theta);
        const int bottom_index = grid.index(i_r, i_theta_M1);
        const int top_index    = grid.index(i_r, i_theta_P1);

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r, grid, DirBC_interior);

        offset = CenterStencil[StencilPosition::Center];
        column = center_index;
        value  = 1.0;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* Give value to the interior nodes! */
        /* Fill matrix row of (i-1,j) */
        row                        = left_index;
        ptr                        = left_nz_index;
        const Stencil& LeftStencil = getStencil(i_r - 1, grid, DirBC_Interior);

        /* Right REMOVED: Moved to the right hand side to make the matrix symmetric */

        offset = LeftStencil[StencilPosition::Center];
        column = left_index;
        value  = coeff1 * arr; /* Center: (Right) */
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* TopRight REMOVED: Moved to the right hand side to make the matrix symmetric */

        /* BottomRight REMOVED: Moved to the right hand side to make the matrix symmetric */
    }
}

} // namespace direct_solver_give

template <class LevelCacheType>
typename DirectSolverGive<LevelCacheType>::SystemMatrix DirectSolverGive<LevelCacheType>::buildSolverMatrix()
{
    using direct_solver_give::getNonZeroCountSolverMatrix;
    using direct_solver_give::getStencilSize;
    using direct_solver_give::nodeBuildSolverMatrixTake;
    using direct_solver_give::validateSolverMatrixIndexing;

    const PolarGrid& grid             = DirectSolver<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = DirectSolver<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = DirectSolver<LevelCacheType>::DirBC_Interior_;

    assert(validateSolverMatrixIndexing(grid, DirBC_Interior) && "Solver matrix indexing is inconsistent");

    const int n = grid.numberOfNodes(grid, DirBC_Interior);

#ifdef GMGPOLAR_USE_MUMPS
    const int nnz = getNonZeroCountSolverMatrix(grid, DirBC_Interior);
    SparseMatrixCOO<double, Kokkos::HostSpace> solver_matrix(n, n, nnz);
    solver_matrix.is_symmetric(true);
#else
    std::function<int(int)> nnz_per_row = [&](int global_index) {
        return getStencilSize(global_index, grid, DirBC_Interior);
    };

    SparseMatrixCSR<double, Kokkos::HostSpace> solver_matrix(n, n, nnz_per_row);
#endif

    /* ---------------- */
    /* Circular section */
    /* ---------------- */
    // We parallelize over i_r (step 3) to avoid data race conditions between adjacent circles.
    // The i_theta loop is sequential inside the kernel.
    const int num_circle_tasks = grid.numberSmootherCircles();

    for (int start_circle = 0; start_circle < 3; ++start_circle) {
        const int num_circular_tasks = (num_circle_tasks - start_circle + 2) / 3;
        Kokkos::parallel_for(
            "DirectSolverGive: BuildSolverMatrix (Circular)",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, num_circular_tasks),
            KOKKOS_LAMBDA(const int circle_task) {
                const int i_r = start_circle + circle_task * 3;
                for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                    nodeBuildSolverMatrixGive(i_r, i_theta, grid, level_cache, DirBC_Interior, solver_matrix);
                }
            });
        Kokkos::fence();
    }

    /* -------------- */
    /* Radial section */
    /* -------------- */
    // We parallelize over i_theta (step 3) to avoid data race conditions between adjacent radial lines.
    // The i_r loop is sequential inside the kernel.
    // Due to periodicity in the angular direction, handle up to 2 additional
    // radial lines (i_theta = 0 and 1) before the parallel passes.
    const int additional_radial_tasks = grid.ntheta() % 3;
    const int num_radial_tasks        = grid.ntheta() - additional_radial_tasks;

    for (int i_theta = 0; i_theta < additional_radial_tasks; i_theta++) {
        Kokkos::parallel_for(
            "DirectSolverGive: BuildSolverMatrix (Radial, additional)",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, 1), KOKKOS_LAMBDA(const int) {
                for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                    nodeBuildSolverMatrixGive(i_r, i_theta, grid, level_cache, DirBC_Interior, solver_matrix);
                }
            });
        Kokkos::fence();
    }

    for (int start_radial = 0; start_radial < 3; ++start_radial) {
        const int num_radial_batches = (num_radial_tasks - start_radial + 2) / 3;
        Kokkos::parallel_for(
            "DirectSolverGive: BuildSolverMatrix (Radial)",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, num_radial_batches),
            KOKKOS_LAMBDA(const int radial_task) {
                const int i_theta = additional_radial_tasks + start_radial + radial_task * 3;
                for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                    nodeBuildSolverMatrixGive(i_r, i_theta, grid, level_cache, DirBC_Interior, solver_matrix);
                }
            });
        Kokkos::fence();
    }

    return solver_matrix;
}
