#include "../../../include/DirectSolver/DirectSolver-COO-MUMPS-Give/directSolverGive.h"

#include "../../../include/Definitions/geometry_helper.h"

#ifdef GMGPOLAR_USE_MUMPS

static inline void updateMatrixElement(SparseMatrixCOO<double>& matrix, int ptr, int offset, int row, int col,
                                       double val)
{
    matrix.row_index(ptr + offset) = row;
    matrix.col_index(ptr + offset) = col;
    matrix.value(ptr + offset) += val;
}

void DirectSolver_COO_MUMPS_Give::nodeBuildSolverMatrixGive(int i_r, int i_theta, const PolarGrid& grid,
                                                            bool DirBC_Interior, SparseMatrixCOO<double>& solver_matrix,
                                                            double arr, double att, double art, double detDF,
                                                            double coeff_beta)
{
    int ptr, offset;
    int row, col;
    double val;
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

        const int center_nz_index = getSolverMatrixIndex(i_r, i_theta);
        const int left_nz_index   = getSolverMatrixIndex(i_r - 1, i_theta);
        const int right_nz_index  = getSolverMatrixIndex(i_r + 1, i_theta);
        const int bottom_nz_index = getSolverMatrixIndex(i_r, i_theta_M1);
        const int top_nz_index    = getSolverMatrixIndex(i_r, i_theta_P1);

        const int center_index = grid.index(i_r, i_theta);
        const int left_index   = grid.index(i_r - 1, i_theta);
        const int right_index  = grid.index(i_r + 1, i_theta);
        const int bottom_index = grid.index(i_r, i_theta_M1);
        const int top_index    = grid.index(i_r, i_theta_P1);

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r);

        offset = CenterStencil[StencilPosition::Center];
        col    = center_index;
        val    = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* beta_{i,j} */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Left];
        col    = left_index;
        val    = -coeff1 * arr; /* Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Right];
        col    = right_index;
        val    = -coeff2 * arr; /* Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Bottom];
        col    = bottom_index;
        val    = -coeff3 * att; /* Bottom */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Top];
        col    = top_index;
        val    = -coeff4 * att; /* Top */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Center];
        col    = center_index;
        val    = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* Fill matrix row of (i-1,j) */
        row = left_index;
        ptr = left_nz_index;

        const Stencil& LeftStencil = getStencil(i_r - 1);

        offset = LeftStencil[StencilPosition::Right];
        col    = center_index;
        val    = -coeff1 * arr; /* Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = LeftStencil[StencilPosition::Center];
        col    = left_index;
        val    = +coeff1 * arr; /* Center: (Right) */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = LeftStencil[StencilPosition::TopRight];
        col    = top_index;
        val    = -0.25 * art; /* Top Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = LeftStencil[StencilPosition::BottomRight];
        col    = bottom_index;
        val    = +0.25 * art; /* Bottom Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* Fill matrix row of (i+1,j) */
        row = right_index;
        ptr = right_nz_index;

        const Stencil& RightStencil = getStencil(i_r + 1);

        offset = RightStencil[StencilPosition::Left];
        col    = center_index;
        val    = -coeff2 * arr; /* Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = RightStencil[StencilPosition::Center];
        col    = right_index;
        val    = +coeff2 * arr; /* Center: (Left) */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = RightStencil[StencilPosition::TopLeft];
        col    = top_index;
        val    = +0.25 * art; /* Top Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = RightStencil[StencilPosition::BottomLeft];
        col    = bottom_index;
        val    = -0.25 * art; /* Bottom Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* Fill matrix row of (i,j-1) */
        row = bottom_index;
        ptr = bottom_nz_index;

        const Stencil& BottomStencil = CenterStencil;

        offset = BottomStencil[StencilPosition::Top];
        col    = center_index;
        val    = -coeff3 * att; /* Top */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = BottomStencil[StencilPosition::Center];
        col    = bottom_index;
        val    = +coeff3 * att; /* Center: (Top) */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = BottomStencil[StencilPosition::TopRight];
        col    = right_index;
        val    = -0.25 * art; /* Top Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = BottomStencil[StencilPosition::TopLeft];
        col    = left_index;
        val    = +0.25 * art; /* Top Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* Fill matrix row of (i,j+1) */
        row = top_index;
        ptr = top_nz_index;

        const Stencil& TopStencil = CenterStencil;

        offset = TopStencil[StencilPosition::Bottom];
        col    = center_index;
        val    = -coeff4 * att; /* Bottom */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = TopStencil[StencilPosition::Center];
        col    = top_index;
        val    = +coeff4 * att; /* Center: (Bottom) */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = TopStencil[StencilPosition::BottomRight];
        col    = right_index;
        val    = +0.25 * art; /* Bottom Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = TopStencil[StencilPosition::BottomLeft];
        col    = left_index;
        val    = -0.25 * art; /* Bottom Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);
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

            const int center_nz_index = getSolverMatrixIndex(i_r, i_theta);
            const int right_nz_index  = getSolverMatrixIndex(i_r + 1, i_theta);

            const int center_index = grid.index(i_r, i_theta);
            const int right_index  = grid.index(i_r + 1, i_theta);
            const int bottom_index = grid.index(i_r, i_theta_M1);
            const int top_index    = grid.index(i_r, i_theta_P1);

            /* Fill matrix row of (i,j) */
            row = center_index;
            ptr = center_nz_index;

            const Stencil& CenterStencil = getStencil(i_r);

            offset = CenterStencil[StencilPosition::Center];
            col    = center_index;
            val    = 1.0;
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            /* Fill matrix row of (i+1,j) */
            row = right_index;
            ptr = right_nz_index;

            const Stencil& RightStencil = getStencil(i_r + 1);

            /* Left REMOVED: Moved to the right hand side to make the matrix symmetric */

            offset = RightStencil[StencilPosition::Center];
            col    = right_index;
            val    = +coeff2 * arr; /* Center: (Left) */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

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

            assert(grid_.ntheta() % 2 == 0);
            const int i_theta_AcrossOrigin = grid.wrapThetaIndex(i_theta + grid.ntheta() / 2);

            double h1 = 2.0 * grid.radius(0);
            double h2 = grid.radialSpacing(i_r);
            double k1 = grid.angularSpacing(i_theta_M1);
            double k2 = grid.angularSpacing(i_theta);

            double coeff1 = 0.5 * (k1 + k2) / h1;
            double coeff2 = 0.5 * (k1 + k2) / h2;
            double coeff3 = 0.5 * (h1 + h2) / k1;
            double coeff4 = 0.5 * (h1 + h2) / k2;

            const int center_nz_index = getSolverMatrixIndex(i_r, i_theta);
            const int left_nz_index   = getSolverMatrixIndex(i_r, i_theta_AcrossOrigin);
            const int right_nz_index  = getSolverMatrixIndex(i_r + 1, i_theta);
            const int bottom_nz_index = getSolverMatrixIndex(i_r, i_theta_M1);
            const int top_nz_index    = getSolverMatrixIndex(i_r, i_theta_P1);

            const int center_index = grid.index(i_r, i_theta);
            const int left_index   = grid.index(i_r, i_theta_AcrossOrigin);
            const int right_index  = grid.index(i_r + 1, i_theta);
            const int bottom_index = grid.index(i_r, i_theta_M1);
            const int top_index    = grid.index(i_r, i_theta_P1);

            /* Fill matrix row of (i,j) */
            row                          = center_index;
            ptr                          = center_nz_index;
            const Stencil& CenterStencil = getStencil(i_r);

            offset = CenterStencil[StencilPosition::Center];
            col    = center_index;
            val    = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* beta_{i,j} */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = CenterStencil[StencilPosition::Left];
            col    = left_index;
            val    = -coeff1 * arr; /* Left */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = CenterStencil[StencilPosition::Right];
            col    = right_index;
            val    = -coeff2 * arr; /* Right */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = CenterStencil[StencilPosition::Bottom];
            col    = bottom_index;
            val    = -coeff3 * att; /* Bottom */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = CenterStencil[StencilPosition::Top];
            col    = top_index;
            val    = -coeff4 * att; /* Top */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = CenterStencil[StencilPosition::Center];
            col    = center_index;
            val    = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            /* Fill matrix row of (i-1,j) */
            /* From view the view of the across origin node, */
            /* the directions are roatated by 180 degrees in the stencil! */
            row = left_index;
            ptr = left_nz_index;

            const Stencil& LeftStencil = CenterStencil;

            offset = LeftStencil[StencilPosition::Left];
            col    = center_index;
            val    = -coeff1 * arr; /* Right -> Left*/
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = LeftStencil[StencilPosition::Center];
            col    = left_index;
            val    = +coeff1 * arr; /* Center: (Right) -> Center: (Left) */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            /* Top Right -> Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */

            /* Bottom Right -> Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */

            /* Fill matrix row of (i+1,j) */
            row                         = right_index;
            ptr                         = right_nz_index;
            const Stencil& RightStencil = getStencil(i_r + 1);

            offset = RightStencil[StencilPosition::Left];
            col    = center_index;
            val    = -coeff2 * arr; /* Left */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = RightStencil[StencilPosition::Center];
            col    = right_index;
            val    = +coeff2 * arr; /* Center: (Left) */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = RightStencil[StencilPosition::TopLeft];
            col    = top_index;
            val    = +0.25 * art; /* Top Left */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = RightStencil[StencilPosition::BottomLeft];
            col    = bottom_index;
            val    = -0.25 * art; /* Bottom Left */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            /* Fill matrix row of (i,j-1) */
            row = bottom_index;
            ptr = bottom_nz_index;

            const Stencil& BottomStencil = CenterStencil;

            offset = BottomStencil[StencilPosition::Top];
            col    = center_index;
            val    = -coeff3 * att; /* Top */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = BottomStencil[StencilPosition::Center];
            col    = bottom_index;
            val    = +coeff3 * att; /* Center: (Top) */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = BottomStencil[StencilPosition::TopRight];
            col    = right_index;
            val    = -0.25 * art; /* Top Right */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            /* TopLeft REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */

            /* Fill matrix row of (i,j+1) */
            row = top_index;
            ptr = top_nz_index;

            const Stencil& TopStencil = CenterStencil;

            offset = TopStencil[StencilPosition::Bottom];
            col    = center_index;
            val    = -coeff4 * att; /* Bottom */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = TopStencil[StencilPosition::Center];
            col    = top_index;
            val    = +coeff4 * att; /* Center: (Bottom) */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = TopStencil[StencilPosition::BottomRight];
            col    = right_index;
            val    = +0.25 * art; /* Bottom Right */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

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

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int center_nz_index = getSolverMatrixIndex(i_r, i_theta);
        const int left_nz_index   = getSolverMatrixIndex(i_r - 1, i_theta);
        const int right_nz_index  = getSolverMatrixIndex(i_r + 1, i_theta);
        const int bottom_nz_index = getSolverMatrixIndex(i_r, i_theta_M1);
        const int top_nz_index    = getSolverMatrixIndex(i_r, i_theta_P1);

        const int center_index = grid.index(i_r, i_theta);
        const int left_index   = grid.index(i_r - 1, i_theta);
        const int right_index  = grid.index(i_r + 1, i_theta);
        const int bottom_index = grid.index(i_r, i_theta_M1);
        const int top_index    = grid.index(i_r, i_theta_P1);

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r);

        offset = CenterStencil[StencilPosition::Center];
        col    = center_index;
        val    = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* beta_{i,j} */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* REMOVED: Moved to the right hand side to make the matrix symmetric */
        if (!DirBC_Interior) {
            offset = CenterStencil[StencilPosition::Left];
            col    = left_index;
            val    = -coeff1 * arr; /* Left */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);
        }

        offset = CenterStencil[StencilPosition::Right];
        col    = right_index;
        val    = -coeff2 * arr; /* Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Bottom];
        col    = bottom_index;
        val    = -coeff3 * att; /* Bottom */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Top];
        col    = top_index;
        val    = -coeff4 * att; /* Top */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Center];
        col    = center_index;
        val    = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        if (!DirBC_Interior) { /* Don't give to the inner Dirichlet boundary! */
            /* Fill matrix row of (i-1,j) */
            row = left_index;
            ptr = left_nz_index;

            const Stencil& LeftStencil = getStencil(i_r - 1);

            offset = LeftStencil[StencilPosition::Right];
            col    = center_index;
            val    = -coeff1 * arr; /* Right */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = LeftStencil[StencilPosition::Center];
            col    = left_index;
            val    = +coeff1 * arr; /* Center: (Right) */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = LeftStencil[StencilPosition::TopRight];
            col    = top_index;
            val    = -0.25 * art; /* Top Right */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = LeftStencil[StencilPosition::BottomRight];
            col    = bottom_index;
            val    = +0.25 * art; /* Bottom Right */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);
        }

        /* Fill matrix row of (i+1,j) */
        row = right_index;
        ptr = right_nz_index;

        const Stencil& RightStencil = getStencil(i_r + 1);

        offset = RightStencil[StencilPosition::Left];
        col    = center_index;
        val    = -coeff2 * arr; /* Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = RightStencil[StencilPosition::Center];
        col    = right_index;
        val    = +coeff2 * arr; /* Center: (Left) */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = RightStencil[StencilPosition::TopLeft];
        col    = top_index;
        val    = +0.25 * art; /* Top Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = RightStencil[StencilPosition::BottomLeft];
        col    = bottom_index;
        val    = -0.25 * art; /* Bottom Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* Fill matrix row of (i,j-1) */
        row = bottom_index;
        ptr = bottom_nz_index;

        const Stencil& BottomStencil = CenterStencil;

        offset = BottomStencil[StencilPosition::Top];
        col    = center_index;
        val    = -coeff3 * att; /* Top */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = BottomStencil[StencilPosition::Center];
        col    = bottom_index;
        val    = +coeff3 * att; /* Center: (Top) */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = BottomStencil[StencilPosition::TopRight];
        col    = right_index;
        val    = -0.25 * art; /* Top Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* REMOVED: Moved to the right hand side to make the matrix symmetric */
        if (!DirBC_Interior) {
            offset = BottomStencil[StencilPosition::TopLeft];
            col    = left_index;
            val    = +0.25 * art; /* Top Left */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);
        }

        /* Fill matrix row of (i,j+1) */
        row = top_index;
        ptr = top_nz_index;

        const Stencil& TopStencil = CenterStencil;

        offset = TopStencil[StencilPosition::Bottom];
        col    = center_index;
        val    = -coeff4 * att; /* Bottom */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = TopStencil[StencilPosition::Center];
        col    = top_index;
        val    = +coeff4 * att; /* Center: (Bottom) */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = TopStencil[StencilPosition::BottomRight];
        col    = right_index;
        val    = +0.25 * art; /* Bottom Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* REMOVED: Moved to the right hand side to make the matrix symmetric */
        if (!DirBC_Interior) {
            offset = TopStencil[StencilPosition::BottomLeft];
            col    = left_index;
            val    = -0.25 * art; /* Bottom Left */
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);
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

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int center_nz_index = getSolverMatrixIndex(i_r, i_theta);
        const int left_nz_index   = getSolverMatrixIndex(i_r - 1, i_theta);
        const int right_nz_index  = getSolverMatrixIndex(i_r + 1, i_theta);
        const int bottom_nz_index = getSolverMatrixIndex(i_r, i_theta_M1);
        const int top_nz_index    = getSolverMatrixIndex(i_r, i_theta_P1);

        const int center_index = grid.index(i_r, i_theta);
        const int left_index   = grid.index(i_r - 1, i_theta);
        const int right_index  = grid.index(i_r + 1, i_theta);
        const int bottom_index = grid.index(i_r, i_theta_M1);
        const int top_index    = grid.index(i_r, i_theta_P1);

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r);

        offset = CenterStencil[StencilPosition::Center];
        col    = center_index;
        val    = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* beta_{i,j} */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Left];
        col    = left_index;
        val    = -coeff1 * arr; /* Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* Right REMOVED: Moved to the right hand side to make the matrix symmetric */

        offset = CenterStencil[StencilPosition::Bottom];
        col    = bottom_index;
        val    = -coeff3 * att; /* Bottom */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Top];
        col    = top_index;
        val    = -coeff4 * att; /* Top */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Center];
        col    = center_index;
        val    = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* Fill matrix row of (i-1,j) */
        row = left_index;
        ptr = left_nz_index;

        const Stencil& LeftStencil = getStencil(i_r - 1);

        offset = LeftStencil[StencilPosition::Right];
        col    = center_index;
        val    = -coeff1 * arr; /* Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = LeftStencil[StencilPosition::Center];
        col    = left_index;
        val    = coeff1 * arr; /* Center: (Right) */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = LeftStencil[StencilPosition::TopRight];
        col    = top_index;
        val    = -0.25 * art; /* Top Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = LeftStencil[StencilPosition::BottomRight];
        col    = bottom_index;
        val    = 0.25 * art; /* Bottom Right */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* Fill matrix row of (i+1,j) */
        /* Don't give to the outer dirichlet boundary! */

        /* Fill matrix row of (i,j-1) */
        row = bottom_index;
        ptr = bottom_nz_index;

        const Stencil& BottomStencil = CenterStencil;

        offset = BottomStencil[StencilPosition::Top];
        col    = center_index;
        val    = -coeff3 * att; /* Top */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = BottomStencil[StencilPosition::Center];
        col    = bottom_index;
        val    = coeff3 * att; /* Center: (Top) */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* TopRight REMOVED: Moved to the right hand side to make the matrix symmetric */

        offset = BottomStencil[StencilPosition::TopLeft];
        col    = left_index;
        val    = 0.25 * art; /* Top Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* Fill matrix row of (i,j+1) */
        row = top_index;
        ptr = top_nz_index;

        const Stencil& TopStencil = CenterStencil;

        offset = TopStencil[StencilPosition::Bottom];
        col    = center_index;
        val    = -coeff4 * att; /* Bottom */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = TopStencil[StencilPosition::Center];
        col    = top_index;
        val    = coeff4 * att; /* Center: (Bottom) */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* BottomRight REMOVED: Moved to the right hand side to make the matrix symmetric */

        offset = TopStencil[StencilPosition::BottomLeft];
        col    = left_index;
        val    = -0.25 * art; /* Bottom Left */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);
    }
    /* ------------------------------------ */
    /* Node on the outer dirichlet boundary */
    /* ------------------------------------ */
    else if (i_r == grid.nr() - 1) {
        double h1 = grid.radialSpacing(i_r - 1);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int center_nz_index = getSolverMatrixIndex(i_r, i_theta);
        const int left_nz_index   = getSolverMatrixIndex(i_r - 1, i_theta);

        const int center_index = grid.index(i_r, i_theta);
        const int left_index   = grid.index(i_r - 1, i_theta);
        const int bottom_index = grid.index(i_r, i_theta_M1);
        const int top_index    = grid.index(i_r, i_theta_P1);

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r);

        offset = CenterStencil[StencilPosition::Center];
        col    = center_index;
        val    = 1.0;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* Give value to the interior nodes! */
        /* Fill matrix row of (i-1,j) */
        row                        = left_index;
        ptr                        = left_nz_index;
        const Stencil& LeftStencil = getStencil(i_r - 1);

        /* Right REMOVED: Moved to the right hand side to make the matrix symmetric */

        offset = LeftStencil[StencilPosition::Center];
        col    = left_index;
        val    = coeff1 * arr; /* Center: (Right) */
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* TopRight REMOVED: Moved to the right hand side to make the matrix symmetric */

        /* BottomRight REMOVED: Moved to the right hand side to make the matrix symmetric */
    }
}

void DirectSolver_COO_MUMPS_Give::buildSolverMatrixCircleSection(const int i_r, SparseMatrixCOO<double>& solver_matrix)
{
    const double r = grid_.radius(i_r);
    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        const int global_index = grid_.index(i_r, i_theta);
        const double theta     = grid_.theta(i_theta);

        double coeff_beta, arr, att, art, detDF;
        level_cache_.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art, detDF);

        // Build solver matrix at the current node
        nodeBuildSolverMatrixGive(i_r, i_theta, grid_, DirBC_Interior_, solver_matrix, arr, att, art, detDF,
                                  coeff_beta);
    }
}

void DirectSolver_COO_MUMPS_Give::buildSolverMatrixRadialSection(const int i_theta,
                                                                 SparseMatrixCOO<double>& solver_matrix)
{
    const double theta = grid_.theta(i_theta);
    for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
        const int global_index = grid_.index(i_r, i_theta);
        const double r         = grid_.radius(i_r);

        double coeff_beta, arr, att, art, detDF;
        level_cache_.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art, detDF);

        // Build solver matrix at the current node
        nodeBuildSolverMatrixGive(i_r, i_theta, grid_, DirBC_Interior_, solver_matrix, arr, att, art, detDF,
                                  coeff_beta);
    }
}

// clang-format off

/* ------------------------------------------------------------------------ */
/* If the indexing is not smoother-based, please adjust the access patterns */
SparseMatrixCOO<double> DirectSolver_COO_MUMPS_Give::buildSolverMatrix()
{
    const int n   = grid_.numberOfNodes();
    const int nnz = getNonZeroCountSolverMatrix();

    // Although the matrix is symmetric, we need to store all its entries, so we disable the symmetry.
    SparseMatrixCOO<double> solver_matrix(n, n, nnz);
    solver_matrix.is_symmetric(false);

    #pragma omp parallel for num_threads(num_omp_threads_) if (nnz > 10'000)
    for (int i = 0; i < nnz; i++) {
        solver_matrix.value(i) = 0.0;
    }

    if (num_omp_threads_ == 1) {
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

        #pragma omp parallel num_threads(num_omp_threads_)
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

#endif
