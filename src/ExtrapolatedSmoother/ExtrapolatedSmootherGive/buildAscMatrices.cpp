#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherGive/extrapolatedSmootherGive.h"

#include "../../../include/Definitions/geometry_helper.h"

/* Tridiagonal matrices */
static inline void updateTridiagonalElement(SymmetricTridiagonalSolver<double>& matrix, int row, int column,
                                            double value)
{
    if (row == column)
        matrix.main_diagonal(row) += value;
    else if (row == column - 1)
        matrix.sub_diagonal(row) += value;
    else if (row == 0 && column == matrix.columns() - 1)
        matrix.cyclic_corner_element() += value;
}

/* Diagonal matrices */
static inline void updateDiagonalElement(DiagonalSolver<double>& matrix, int row, int column, double value)
{
    matrix.diagonal(row) += value;
}

/* Inner Boundary COO/CSR matrix */
#ifdef GMGPOLAR_USE_MUMPS
static inline void updateCOOCSRMatrixElement(SparseMatrixCOO<double>& matrix, int ptr, int offset, int row, int col,
                                             double val)
{
    matrix.row_index(ptr + offset) = row;
    matrix.col_index(ptr + offset) = col;
    matrix.value(ptr + offset) += val;
}
#else
static inline void updateCOOCSRMatrixElement(SparseMatrixCSR<double>& matrix, int ptr, int offset, int row, int col,
                                             double val)
{
    matrix.row_nz_index(row, offset) = col;
    matrix.row_nz_entry(row, offset) += val;
}
#endif

void ExtrapolatedSmootherGive::nodeBuildSmootherGive(
    int i_r, int i_theta, const PolarGrid& grid, bool DirBC_Interior, MatrixType& inner_boundary_circle_matrix,
    std::vector<DiagonalSolver<double>>& circle_diagonal_solver,
    std::vector<DiagonalSolver<double>>& radial_diagonal_solver,
    std::vector<SymmetricTridiagonalSolver<double>>& circle_tridiagonal_solver,
    std::vector<SymmetricTridiagonalSolver<double>>& radial_tridiagonal_solver, double arr, double att, double art,
    double detDF, double coeff_beta)
{
    assert(i_r >= 0 && i_r < grid.nr());
    assert(i_theta >= 0 && i_theta < grid.ntheta());

    const int numberSmootherCircles = grid.numberSmootherCircles();
    const int lengthSmootherRadial  = grid.lengthSmootherRadial();

    assert(numberSmootherCircles >= 3);
    assert(lengthSmootherRadial >= 3);

    int ptr, offset;
    int row, column, col;
    double value, val;
    /* ------------------------------------------ */
    /* Node in the interior of the Circle Section */
    /* ------------------------------------------ */
    if (i_r > 0 && i_r < numberSmootherCircles - 1) {
        double h1 = grid.radialSpacing(i_r - 1);
        double h2 = grid.radialSpacing(i_r);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;

        int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        int center_index = i_theta;
        int left_index   = i_theta;
        int right_index  = i_theta;
        int bottom_index = i_theta_M1;
        int top_index    = i_theta_P1;
        /* -------------------------- */
        /* Cyclic Tridiagonal Section */
        /* i_r % 2 == 1               */
        if (i_r & 1) {
            /* i_theta % 2 == 1 */
            /* | x | o | x | */
            /* |   |   |   | */
            /* | o | O | o | */
            /* |   |   |   | */
            /* | x | o | x | */
            /* or */
            /* i_theta % 2 == 0 */
            /* | o | o | o | */
            /* |   |   |   | */
            /* | x | O | x | */
            /* |   |   |   | */
            /* | o | o | o | */

            auto& center_matrix = circle_tridiagonal_solver[i_r / 2];
            auto& left_matrix   = circle_diagonal_solver[(i_r - 1) / 2];
            auto& right_matrix  = circle_diagonal_solver[(i_r + 1) / 2];

            /* Fill matrix row of (i,j) */
            row    = center_index;
            column = center_index;
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */
            updateTridiagonalElement(center_matrix, row, column, value);

            row    = center_index;
            column = bottom_index;
            value  = -coeff3 * att; /* Bottom */
            updateTridiagonalElement(center_matrix, row, column, value);

            row    = center_index;
            column = top_index;
            value  = -coeff4 * att; /* Top */
            updateTridiagonalElement(center_matrix, row, column, value);

            row    = center_index;
            column = center_index;
            value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
            updateTridiagonalElement(center_matrix, row, column, value);

            /* Fill matrix row of (i,j-1) */
            row    = bottom_index;
            column = center_index;
            value  = -coeff3 * att; /* Top */
            updateTridiagonalElement(center_matrix, row, column, value);

            row    = bottom_index;
            column = bottom_index;
            value  = coeff3 * att; /* Center: (Top) */
            updateTridiagonalElement(center_matrix, row, column, value);

            /* Fill matrix row of (i,j+1) */
            row    = top_index;
            column = center_index;
            value  = -coeff4 * att; /* Bottom */
            updateTridiagonalElement(center_matrix, row, column, value);

            row    = top_index;
            column = top_index;
            value  = coeff4 * att; /* Center: (Bottom) */
            updateTridiagonalElement(center_matrix, row, column, value);

            if (i_theta & 1) {
                /* i_theta % 2 == 1 */
                /* | x | o | x | */
                /* |   |   |   | */
                /* | o | O | o | */
                /* |   |   |   | */
                /* | x | o | x | */

                /* Fill matrix row of (i-1,j) */
                if (i_r == 1) {
                    /* Only in the case of AcrossOrigin */
                    if (!DirBC_Interior) {
                        row = left_index;
                        ptr = getCircleAscIndex(i_r - 1, i_theta);

                        const Stencil& LeftStencil = getStencil(i_r - 1, i_theta);

                        offset = LeftStencil[StencilPosition::Center];
                        col    = left_index;
                        val    = +coeff1 * arr; /* Center: (Right) */
                        updateCOOCSRMatrixElement(inner_boundary_circle_matrix, ptr, offset, row, col, val);
                    }
                }
                else {
                    row    = left_index;
                    column = left_index;
                    value  = coeff1 * arr; /* Center: (Right) */
                    updateDiagonalElement(left_matrix, row, column, value);
                }

                /* Fill matrix row of (i+1,j) */
                row    = right_index;
                column = right_index;
                value  = coeff2 * arr; /* Center: (Left) */
                updateDiagonalElement(right_matrix, row, column, value);
            }
        }
        /* ---------------- */
        /* Diagonal Section */
        /* i_r % 2 == 0     */
        else {
            /* i_theta % 2 == 1 */
            /* | o | x | o | */
            /* |   |   |   | */
            /* | o | O | o | */
            /* |   |   |   | */
            /* | o | x | o | */
            /* or */
            /* i_theta % 2 == 0 */
            /* | o | o | o | */
            /* |   |   |   | */
            /* | o | X | o | */
            /* |   |   |   | */
            /* | o | o | o | */

            auto& center_matrix = circle_diagonal_solver[i_r / 2];
            auto& left_matrix   = circle_tridiagonal_solver[(i_r - 1) / 2];
            auto& right_matrix  = circle_tridiagonal_solver[(i_r + 1) / 2];

            if (i_theta & 1) { /* i_theta % 2 == 1 */
                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
                updateDiagonalElement(center_matrix, row, column, value);

                row    = center_index;
                column = center_index;
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
                updateDiagonalElement(center_matrix, row, column, value);
            }
            else { /* i_theta % 2 == 0 */
                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 1.0;
                updateDiagonalElement(center_matrix, row, column, value);

                /* Fill matrix row of (i,j-1) */
                row    = bottom_index;
                column = bottom_index;
                value  = coeff3 * att; /* Center: (Top) */
                updateDiagonalElement(center_matrix, row, column, value);

                /* Fill matrix row of (i,j+1) */
                row    = top_index;
                column = top_index;
                value  = coeff4 * att; /* Center: (Bottom) */
                updateDiagonalElement(center_matrix, row, column, value);
            }
            /* Fill matrix row of (i-1,j) */
            row    = left_index;
            column = left_index;
            value  = coeff1 * arr; /* Center: (Right) */
            updateTridiagonalElement(left_matrix, row, column, value);

            /* Fill matrix row of (i+1,j) */
            row    = right_index;
            column = right_index;
            value  = coeff2 * arr; /* Center: (Left) */
            updateTridiagonalElement(right_matrix, row, column, value);
        }
    }
    /* ------------------------------------------ */
    /* Node in the interior of the Radial Section */
    /* ------------------------------------------ */
    else if (i_r > numberSmootherCircles && i_r < grid.nr() - 2) {
        double h1 = grid.radialSpacing(i_r - 1);
        double h2 = grid.radialSpacing(i_r);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;

        int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        int center_index = i_r - numberSmootherCircles;
        int left_index   = i_r - numberSmootherCircles - 1;
        int right_index  = i_r - numberSmootherCircles + 1;
        int bottom_index = i_r - numberSmootherCircles;
        int top_index    = i_r - numberSmootherCircles;
        /* ------------------- */
        /* Tridiagonal Section */
        /* i_theta % 2 == 1    */
        if (i_theta & 1) {
            /* i_r % 2 == 1 */
            /* ---------- */
            /* x   o   x  */
            /* ---------- */
            /* o   O   o  */
            /* ---------- */
            /* x   o   x  */
            /* ---------- */
            /* or */
            /* i_r % 2 == 0 */
            /* ---------- */
            /* o   x   o  */
            /* ---------- */
            /* o   O   o  */
            /* ---------- */
            /* o   x   o  */
            /* ---------- */

            auto& center_matrix = radial_tridiagonal_solver[i_theta / 2];
            auto& bottom_matrix = radial_diagonal_solver[i_theta_M1 / 2];
            auto& top_matrix    = radial_diagonal_solver[i_theta_P1 / 2];

            /* Fill matrix row of (i,j) */
            row    = center_index;
            column = center_index;
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * fabs(detDF); /* Center: beta_{i,j} */
            updateTridiagonalElement(center_matrix, row, column, value);

            row    = center_index;
            column = left_index;
            value  = -coeff1 * arr; /* Left */
            updateTridiagonalElement(center_matrix, row, column, value);

            row    = center_index;
            column = right_index;
            value  = -coeff2 * arr; /* Right */
            updateTridiagonalElement(center_matrix, row, column, value);

            row    = center_index;
            column = center_index;
            value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
            updateTridiagonalElement(center_matrix, row, column, value);

            /* Fill matrix row of (i-1,j) */
            row    = left_index;
            column = center_index;
            value  = -coeff1 * arr; /* Right */
            updateTridiagonalElement(center_matrix, row, column, value);

            row    = left_index;
            column = left_index;
            value  = coeff1 * arr; /* Center: (Right) */
            updateTridiagonalElement(center_matrix, row, column, value);

            /* Fill matrix row of (i+1,j) */
            row    = right_index;
            column = center_index;
            value  = -coeff2 * arr; /* Left */
            updateTridiagonalElement(center_matrix, row, column, value);

            row    = right_index;
            column = right_index;
            value  = coeff2 * arr; /* Center: (Left) */
            updateTridiagonalElement(center_matrix, row, column, value);

            if (i_r & 1) { /* i_r % 2 == 1 */
                /* ---------- */
                /* x   o   x  */
                /* ---------- */
                /* o   O   o  */
                /* ---------- */
                /* x   o   x  */
                /* ---------- */
                /* Fill matrix row of (i,j-1) */
                row    = bottom_index;
                column = bottom_index;
                value  = coeff3 * att; /* Center: (Top) */
                updateDiagonalElement(bottom_matrix, row, column, value);

                /* Fill matrix row of (i,j+1) */
                row    = top_index;
                column = top_index;
                value  = coeff4 * att; /* Center: (Bottom) */
                updateDiagonalElement(top_matrix, row, column, value);
            }
        }
        /* ---------------- */
        /* Diagonal Section */
        /* i_theta % 2 == 0 */
        else {
            /* i_r % 2 == 1 */
            /* ---------- */
            /* o   o   o  */
            /* ---------- */
            /* x   O   x  */
            /* ---------- */
            /* o   o   o  */
            /* ---------- */
            /* or */
            /* i_r % 2 == 0 */
            /* ---------- */
            /* o   o   o  */
            /* ---------- */
            /* o   X   o  */
            /* ---------- */
            /* o   o   o  */
            /* ---------- */

            auto& center_matrix = radial_diagonal_solver[i_theta / 2];
            auto& bottom_matrix = radial_tridiagonal_solver[i_theta_M1 / 2];
            auto& top_matrix    = radial_tridiagonal_solver[i_theta_P1 / 2];
            if (i_r & 1) { /* i_r % 2 == 1 */
                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
                updateDiagonalElement(center_matrix, row, column, value);

                row    = center_index;
                column = center_index;
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
                updateDiagonalElement(center_matrix, row, column, value);
            }
            else { /* i_r % 2 == 0 */
                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 1.0;
                updateDiagonalElement(center_matrix, row, column, value);

                /* Fill matrix row of (i-1,j) */
                row    = left_index;
                column = left_index;
                value  = coeff1 * arr; /* Center: (Right) */
                updateDiagonalElement(center_matrix, row, column, value);

                /* Fill matrix row of (i+1,j) */
                row    = right_index;
                column = right_index;
                value  = coeff2 * arr; /* Center: (Left) */
                updateDiagonalElement(center_matrix, row, column, value);
            }
            /* Fill matrix row of (i,j-1) */
            row    = bottom_index;
            column = bottom_index;
            value  = coeff3 * att; /* Center: (Top) */
            updateTridiagonalElement(bottom_matrix, row, column, value);

            /* Fill matrix row of (i,j+1) */
            row    = top_index;
            column = top_index;
            value  = coeff4 * att; /* Center: (Bottom) */
            updateTridiagonalElement(top_matrix, row, column, value);
        }
    }
    /* ------------------------------------------ */
    /* Circle Section: Node in the inner boundary */
    /* ------------------------------------------ */
    else if (i_r == 0) {
        /* ------------------------------------------------ */
        /* Case 1: Dirichlet boundary on the inner boundary */
        /* ------------------------------------------------ */
        if (DirBC_Interior) {
            /* Fill result(i,j) */
            double h2 = grid.radialSpacing(i_r);
            double k1 = grid.angularSpacing(i_theta - 1);
            double k2 = grid.angularSpacing(i_theta);

            double coeff2 = 0.5 * (k1 + k2) / h2;

            int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
            int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

            auto& center_matrix = inner_boundary_circle_matrix;
            auto& right_matrix  = circle_tridiagonal_solver[(i_r + 1) / 2];

            int center_index = i_theta;
            int right_index  = i_theta;
            int bottom_index = i_theta_M1;
            int top_index    = i_theta_P1;

            /* Fill matrix row of (i,j) */
            row = center_index;
            ptr = getCircleAscIndex(i_r, i_theta);

            const Stencil& CenterStencil = getStencil(i_r, i_theta);

            offset = CenterStencil[StencilPosition::Center];
            col    = center_index;
            val    = 1.0;
            updateCOOCSRMatrixElement(center_matrix, ptr, offset, row, col, val);

            /* Fill matrix row of (i+1,j) */
            row    = right_index;
            column = right_index;
            value  = coeff2 * arr; /* Center: (Left) */
            updateTridiagonalElement(right_matrix, row, column, value);
        }
        else {
            /* ------------------------------------------------------------- */
            /* Case 2: Across origin discretization on the interior boundary */
            /* ------------------------------------------------------------- */
            // h1 gets replaced with 2 * R0.
            // (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()/2)).
            // Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil.
            double h1 = 2.0 * grid.radius(0);
            double h2 = grid.radialSpacing(i_r);
            double k1 = grid.angularSpacing(i_theta - 1);
            double k2 = grid.angularSpacing(i_theta);

            double coeff1 = 0.5 * (k1 + k2) / h1;
            double coeff2 = 0.5 * (k1 + k2) / h2;
            double coeff3 = 0.5 * (h1 + h2) / k1;
            double coeff4 = 0.5 * (h1 + h2) / k2;

            const int i_theta_M1           = grid.wrapThetaIndex(i_theta - 1);
            const int i_theta_P1           = grid.wrapThetaIndex(i_theta + 1);
            const int i_theta_AcrossOrigin = grid.wrapThetaIndex(i_theta + (grid.ntheta() / 2));

            const int center_index = i_theta;
            const int left_index   = i_theta_AcrossOrigin;
            const int right_index  = i_theta;
            const int bottom_index = i_theta_M1;
            const int top_index    = i_theta_P1;

            const int center_nz_index = getCircleAscIndex(i_r, i_theta);
            const int bottom_nz_index = getCircleAscIndex(i_r, i_theta_M1);
            const int top_nz_index    = getCircleAscIndex(i_r, i_theta_P1);
            const int left_nz_index   = getCircleAscIndex(i_r, i_theta_AcrossOrigin);

            int nz_index;
            const Stencil& CenterStencil = getStencil(i_r, i_theta);

            if (i_theta & 1) {
                /* i_theta % 2 == 1 */
                /* -| x | o | x | */
                /* -|   |   |   | */
                /* -| O | o | o | */
                /* -|   |   |   | */
                /* -| x | o | x | */

                auto& center_matrix = inner_boundary_circle_matrix;
                auto& right_matrix  = circle_tridiagonal_solver[(i_r + 1) / 2];
                auto& left_matrix   = inner_boundary_circle_matrix;
                /* Fill matrix row of (i,j) */
                row = center_index;
                ptr = center_nz_index;

                offset = CenterStencil[StencilPosition::Center];
                col    = center_index;
                val    = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* beta_{i,j} */
                updateCOOCSRMatrixElement(center_matrix, ptr, offset, row, col, val);

                offset = CenterStencil[StencilPosition::Left];
                col    = left_index;
                val    = -coeff1 * arr; /* Left */
                updateCOOCSRMatrixElement(center_matrix, ptr, offset, row, col, val);

                offset = CenterStencil[StencilPosition::Center];
                col    = center_index;
                val    = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
                updateCOOCSRMatrixElement(center_matrix, ptr, offset, row, col, val);

                /* Fill matrix row of (i-1,j) */
                /* From view the view of the across origin node, */
                /* the directions are roatated by 180 degrees in the stencil! */
                row = left_index;
                ptr = left_nz_index;

                const Stencil& LeftStencil = CenterStencil;

                offset = LeftStencil[StencilPosition::Left];
                col    = center_index;
                val    = -coeff1 * arr; /* Right -> Left*/
                updateCOOCSRMatrixElement(left_matrix, ptr, offset, row, col, val);

                offset = LeftStencil[StencilPosition::Center];
                col    = left_index;
                val    = +coeff1 * arr; /* Center: (Right) -> Center: (Left) */
                updateCOOCSRMatrixElement(left_matrix, ptr, offset, row, col, val);

                /* Fill matrix row of (i+1,j) */
                row    = right_index;
                column = right_index;
                value  = coeff2 * arr; /* Center: (Left) */
                updateTridiagonalElement(right_matrix, row, column, value);
            }
            else {
                /* i_theta % 2 == 0 */
                /* -| o | o | o | */
                /* -|   |   |   | */
                /* -| X | o | x | */
                /* -|   |   |   | */
                /* -| o | o | o | */

                auto& center_matrix = inner_boundary_circle_matrix;
                auto& right_matrix  = circle_tridiagonal_solver[(i_r + 1) / 2];
                auto& left_matrix   = inner_boundary_circle_matrix;
                /* Fill matrix row of (i,j) */
                row = center_index;
                ptr = center_nz_index;

                offset = CenterStencil[StencilPosition::Center];
                col    = center_index;
                val    = 1.0;
                updateCOOCSRMatrixElement(center_matrix, ptr, offset, row, col, val);

                /* Fill matrix row of (i,j-1) */
                row = bottom_index;
                ptr = bottom_nz_index;

                const Stencil& BottomStencil = CenterStencil;

                offset = BottomStencil[StencilPosition::Center];
                col    = bottom_index;
                val    = +coeff3 * att; /* Center: (Top) */
                updateCOOCSRMatrixElement(center_matrix, ptr, offset, row, col, val);

                /* Fill matrix row of (i,j+1) */
                row = top_index;
                ptr = top_nz_index;

                const Stencil& TopStencil = CenterStencil;

                offset = TopStencil[StencilPosition::Center];
                col    = top_index;
                val    = +coeff4 * att; /* Center: (Bottom) */
                updateCOOCSRMatrixElement(center_matrix, ptr, offset, row, col, val);

                /* Fill matrix row of (i+1,j) */
                row    = right_index;
                column = right_index;
                value  = coeff2 * arr; /* Center: (Left) */
                updateTridiagonalElement(right_matrix, row, column, value);
            }
        }
    }
    /* ------------------------------------------- */
    /* Circle Section: Node next to radial section */
    /* ------------------------------------------- */
    else if (i_r == numberSmootherCircles - 1) {
        assert(i_r > 1);

        double h1 = grid.radialSpacing(i_r - 1);
        double h2 = grid.radialSpacing(i_r);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;

        int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        int center_index = i_theta;
        int left_index   = i_theta;
        int right_index  = 0;
        int bottom_index = i_theta_M1;
        int top_index    = i_theta_P1;

        if (i_r & 1) {
            if (i_theta & 1) {
                /* i_r % 2 == 1 and i_theta % 2 == 1 */
                /* | o | x | o || x   o   x   o  */
                /* |   |   |   || -------------- */
                /* | o | o | O || o   o   o   o  */
                /* |   |   |   || -------------- */
                /* | o | x | o || x   o   x   o  */

                auto& center_matrix = circle_tridiagonal_solver[i_r / 2];
                auto& left_matrix   = circle_diagonal_solver[(i_r - 1) / 2];
                auto& right_matrix  = radial_tridiagonal_solver[i_theta / 2];

                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
                updateTridiagonalElement(center_matrix, row, column, value);

                row    = center_index;
                column = bottom_index;
                value  = -coeff3 * att; /* Bottom */
                updateTridiagonalElement(center_matrix, row, column, value);

                row    = center_index;
                column = top_index;
                value  = -coeff4 * att; /* Top */
                updateTridiagonalElement(center_matrix, row, column, value);

                row    = center_index;
                column = center_index;
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
                updateTridiagonalElement(center_matrix, row, column, value);

                /* Fill matrix row of (i,j-1) */
                row    = bottom_index;
                column = center_index;
                value  = -coeff3 * att; /* Top */
                updateTridiagonalElement(center_matrix, row, column, value);

                row    = bottom_index;
                column = bottom_index;
                value  = coeff3 * att; /* Center: (Top) */
                updateTridiagonalElement(center_matrix, row, column, value);

                /* Fill matrix row of (i,j+1) */
                row    = top_index;
                column = center_index;
                value  = -coeff4 * att; /* Bottom */
                updateTridiagonalElement(center_matrix, row, column, value);

                row    = top_index;
                column = top_index;
                value  = coeff4 * att; /* Center: (Bottom) */
                updateTridiagonalElement(center_matrix, row, column, value);

                /* Fill matrix row of (i-1,j) */
                row    = left_index;
                column = left_index;
                value  = coeff1 * arr; /* Center: (Right) */
                updateDiagonalElement(left_matrix, row, column, value);

                /* Fill matrix row of (i+1,j) */
                row    = right_index;
                column = right_index;
                value  = coeff2 * arr; /* Center: (Left) */
                updateTridiagonalElement(right_matrix, row, column, value);
            }
            else {
                /* i_r % 2 == 1 and i_theta % 2 == 0 */
                /* | o | o | o || o   o   o   o  */
                /* |   |   |   || -------------- */
                /* | o | x | O || x   o   x   o  */
                /* |   |   |   || -------------- */
                /* | o | o | o || o   o   o   o  */

                auto& center_matrix = circle_tridiagonal_solver[i_r / 2];
                auto& left_matrix   = circle_diagonal_solver[(i_r - 1) / 2];
                auto& right_matrix  = radial_diagonal_solver[i_theta / 2];

                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
                updateTridiagonalElement(center_matrix, row, column, value);

                row    = center_index;
                column = bottom_index;
                value  = -coeff3 * att; /* Bottom */
                updateTridiagonalElement(center_matrix, row, column, value);

                row    = center_index;
                column = top_index;
                value  = -coeff4 * att; /* Top */
                updateTridiagonalElement(center_matrix, row, column, value);

                row    = center_index;
                column = center_index;
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
                updateTridiagonalElement(center_matrix, row, column, value);

                /* Fill matrix row of (i,j-1) */
                row    = bottom_index;
                column = center_index;
                value  = -coeff3 * att; /* Top */
                updateTridiagonalElement(center_matrix, row, column, value);

                row    = bottom_index;
                column = bottom_index;
                value  = coeff3 * att; /* Center: (Top) */
                updateTridiagonalElement(center_matrix, row, column, value);

                /* Fill matrix row of (i,j+1) */
                row    = top_index;
                column = center_index;
                value  = -coeff4 * att; /* Bottom */
                updateTridiagonalElement(center_matrix, row, column, value);

                row    = top_index;
                column = top_index;
                value  = coeff4 * att; /* Center: (Bottom) */
                updateTridiagonalElement(center_matrix, row, column, value);
            }
        }
        else {
            if (i_theta & 1) {
                /* i_r % 2 == 0 and i_theta % 2 == 1 */
                /* | x | o | x || o   x   o   x  */
                /* |   |   |   || -------------- */
                /* | o | o | O || o   o   o   o  */
                /* |   |   |   || -------------- */
                /* | x | o | x || o   x   o   x  */

                auto& center_matrix = circle_diagonal_solver[i_r / 2];
                auto& left_matrix   = circle_tridiagonal_solver[(i_r - 1) / 2];
                auto& right_matrix  = radial_tridiagonal_solver[i_theta / 2];

                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
                updateDiagonalElement(center_matrix, row, column, value);

                row    = center_index;
                column = center_index;
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
                updateDiagonalElement(center_matrix, row, column, value);

                /* Fill matrix row of (i-1,j) */
                row    = left_index;
                column = left_index;
                value  = coeff1 * arr; /* Center: (Right) */
                updateTridiagonalElement(left_matrix, row, column, value);

                /* Fill matrix row of (i+1,j) */
                row    = right_index;
                column = right_index;
                value  = coeff2 * arr; /* Center: (Left) */
                updateTridiagonalElement(right_matrix, row, column, value);
            }
            else {
                /* i_r % 2 == 0 and i_theta % 2 == 0 */
                /* | o | o | o || o   o   o   o  */
                /* |   |   |   || -------------- */
                /* | x | o | X || o   x   o   x  */
                /* |   |   |   || -------------- */
                /* | o | o | o || o   o   o   o  */

                auto& center_matrix = circle_diagonal_solver[i_r / 2];
                auto& left_matrix   = circle_tridiagonal_solver[(i_r - 1) / 2];
                auto& right_matrix  = radial_diagonal_solver[i_theta / 2];

                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 1.0;
                updateDiagonalElement(center_matrix, row, column, value);

                /* Fill matrix row of (i,j-1) */
                row    = bottom_index;
                column = bottom_index;
                value  = coeff3 * att; /* Center: (Top) */
                updateDiagonalElement(center_matrix, row, column, value);

                /* Fill matrix row of (i,j+1) */
                row    = top_index;
                column = top_index;
                value  = coeff4 * att; /* Center: (Bottom) */
                updateDiagonalElement(center_matrix, row, column, value);

                /* Fill matrix row of (i-1,j) */
                row    = left_index;
                column = left_index;
                value  = coeff1 * arr; /* Center: (Right) */
                updateTridiagonalElement(left_matrix, row, column, value);

                /* Fill matrix row of (i+1,j) */
                row    = right_index;
                column = right_index;
                value  = coeff2 * arr; /* Center: (Left) */
                updateDiagonalElement(right_matrix, row, column, value);
            }
        }
    }
    /* --------------------------------------------- */
    /* Radial Section: Node next to circular section */
    /* --------------------------------------------- */
    else if (i_r == numberSmootherCircles) {
        double h1 = grid.radialSpacing(i_r - 1);
        double h2 = grid.radialSpacing(i_r);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int center_index = i_r - numberSmootherCircles;
        const int left_index   = i_theta;
        const int right_index  = i_r - numberSmootherCircles + 1;
        const int bottom_index = i_r - numberSmootherCircles;
        const int top_index    = i_r - numberSmootherCircles;

        if (i_theta & 1) {
            if (i_r & 1) {
                /* i_theta % 2 == 1 and i_r % 2 == 1 */
                /* | x | o | x || o   x   o   x  */
                /* |   |   |   || -------------- */
                /* | o | o | o || O   o   o   o  */
                /* |   |   |   || -------------- */
                /* | x | o | x || o   x   o   x  */

                auto& center_matrix = radial_tridiagonal_solver[i_theta / 2];
                auto& bottom_matrix = radial_diagonal_solver[i_theta_M1 / 2];
                auto& top_matrix    = radial_diagonal_solver[i_theta_P1 / 2];
                auto& left_matrix   = circle_diagonal_solver[(i_r - 1) / 2];

                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
                updateTridiagonalElement(center_matrix, row, column, value);

                row    = center_index;
                column = right_index;
                value  = -coeff2 * arr; /* Right */
                updateTridiagonalElement(center_matrix, row, column, value);

                row    = center_index;
                column = center_index;
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
                updateTridiagonalElement(center_matrix, row, column, value);

                /* Fill matrix row of (i-1,j) */
                row    = left_index;
                column = left_index;
                value  = coeff1 * arr; /* Center: (Right) */
                updateDiagonalElement(left_matrix, row, column, value);

                /* Fill matrix row of (i+1,j) */
                row    = right_index;
                column = center_index;
                value  = -coeff2 * arr; /* Left */
                updateTridiagonalElement(center_matrix, row, column, value);

                row    = right_index;
                column = right_index;
                value  = coeff2 * arr; /* Center: (Left) */
                updateTridiagonalElement(center_matrix, row, column, value);

                /* Fill matrix row of (i,j-1) */
                row    = bottom_index;
                column = bottom_index;
                value  = coeff3 * att; /* Center: (Top) */
                updateDiagonalElement(bottom_matrix, row, column, value);

                /* Fill matrix row of (i,j+1) */
                row    = top_index;
                column = top_index;
                value  = coeff4 * att; /* Center: (Bottom) */
                updateDiagonalElement(top_matrix, row, column, value);
            }
            else {
                /* i_theta % 2 == 1 and i_r % 2 == 0 */
                /* | o | x | o || x   o   x   o  */
                /* |   |   |   || -------------- */
                /* | o | o | o || O   o   o   o  */
                /* |   |   |   || -------------- */
                /* | o | x | o || x   o   x   o  */

                auto& center_matrix = radial_tridiagonal_solver[i_theta / 2];
                auto& bottom_matrix = radial_diagonal_solver[i_theta_M1 / 2];
                auto& top_matrix    = radial_diagonal_solver[i_theta_P1 / 2];
                auto& left_matrix   = circle_tridiagonal_solver[(i_r - 1) / 2];

                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
                updateTridiagonalElement(center_matrix, row, column, value);

                row    = center_index;
                column = right_index;
                value  = -coeff2 * arr; /* Right */
                updateTridiagonalElement(center_matrix, row, column, value);

                row    = center_index;
                column = center_index;
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
                updateTridiagonalElement(center_matrix, row, column, value);

                /* Fill matrix row of (i-1,j) */
                row    = left_index;
                column = left_index;
                value  = coeff1 * arr; /* Center: (Right) */
                updateTridiagonalElement(left_matrix, row, column, value);

                /* Fill matrix row of (i+1,j) */
                row    = right_index;
                column = center_index;
                value  = -coeff2 * arr; /* Left */
                updateTridiagonalElement(center_matrix, row, column, value);

                row    = right_index;
                column = right_index;
                value  = coeff2 * arr; /* Center: (Left) */
                updateTridiagonalElement(center_matrix, row, column, value);
            }
        }
        else {
            if (i_r & 1) {
                /* i_theta % 2 == 0 and i_r % 2 == 1 */
                /* | o | o | o || o   o   o   o  */
                /* |   |   |   || -------------- */
                /* | x | o | x || O   x   o   x  */
                /* |   |   |   || -------------- */
                /* | o | o | o || o   o   o   o  */

                auto& center_matrix = radial_diagonal_solver[i_theta / 2];
                auto& bottom_matrix = radial_tridiagonal_solver[i_theta_M1 / 2];
                auto& top_matrix    = radial_tridiagonal_solver[i_theta_P1 / 2];
                auto& left_matrix   = circle_diagonal_solver[(i_r - 1) / 2];

                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
                updateDiagonalElement(center_matrix, row, column, value);

                row    = center_index;
                column = center_index;
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
                updateDiagonalElement(center_matrix, row, column, value);

                /* Fill matrix row of (i,j-1) */
                row    = bottom_index;
                column = bottom_index;
                value  = coeff3 * att; /* Center: (Top) */
                updateTridiagonalElement(bottom_matrix, row, column, value);

                /* Fill matrix row of (i,j+1) */
                row    = top_index;
                column = top_index;
                value  = coeff4 * att; /* Center: (Bottom) */
                updateTridiagonalElement(top_matrix, row, column, value);
            }
            else {
                /* i_theta % 2 == 0 and i_r % 2 == 0 */
                /* | o | o | o || o   o   o   o  */
                /* |   |   |   || -------------- */
                /* | o | x | o || X   o   x   o  */
                /* |   |   |   || -------------- */
                /* | o | o | o || o   o   o   o  */

                auto& center_matrix = radial_diagonal_solver[i_theta / 2];
                auto& bottom_matrix = radial_tridiagonal_solver[i_theta_M1 / 2];
                auto& top_matrix    = radial_tridiagonal_solver[i_theta_P1 / 2];
                auto& left_matrix   = circle_tridiagonal_solver[(i_r - 1) / 2];

                /* Fill matrix row of (i,j) */

                row    = center_index;
                column = center_index;
                value  = 1.0;
                updateDiagonalElement(center_matrix, row, column, value);

                /* Fill matrix row of (i-1,j) */
                row    = left_index;
                column = left_index;
                value  = coeff1 * arr; /* Center: (Right) */
                updateTridiagonalElement(left_matrix, row, column, value);

                /* Fill matrix row of (i+1,j) */
                row    = right_index;
                column = right_index;
                value  = coeff2 * arr; /* Center: (Left) */
                updateDiagonalElement(center_matrix, row, column, value);

                /* Fill matrix row of (i,j-1) */
                row    = bottom_index;
                column = bottom_index;
                value  = coeff3 * att; /* Center: (Top) */
                updateTridiagonalElement(bottom_matrix, row, column, value);

                /* Fill matrix row of (i,j+1) */
                row    = top_index;
                column = top_index;
                value  = coeff4 * att; /* Center: (Bottom) */
                updateTridiagonalElement(top_matrix, row, column, value);
            }
        }
    }
    /* ------------------------------------------- */
    /* Radial Section: Node next to outer boundary */
    /* ------------------------------------------- */
    else if (i_r == grid.nr() - 2) {
        assert(i_r % 2 == 1);

        double h1 = grid.radialSpacing(i_r - 1);
        double h2 = grid.radialSpacing(i_r);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;

        int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        int center_index = i_r - numberSmootherCircles;
        int left_index   = i_r - numberSmootherCircles - 1;
        int right_index  = i_r - numberSmootherCircles + 1;
        int bottom_index = i_r - numberSmootherCircles;
        int top_index    = i_r - numberSmootherCircles;

        if (i_theta & 1) {
            /* i_theta % 2 == 1 */
            /* ---------------|| */
            /* o   x   o   x  || */
            /* ---------------|| */
            /* o   o   O   o  || */
            /* ---------------|| */
            /* o   x   o   x  || */
            /* ---------------|| */
            auto& center_matrix = radial_tridiagonal_solver[i_theta / 2];
            auto& bottom_matrix = radial_diagonal_solver[i_theta_M1 / 2];
            auto& top_matrix    = radial_diagonal_solver[i_theta_P1 / 2];
            /* Fill matrix row of (i,j) */
            row    = center_index;
            column = center_index;
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
            updateTridiagonalElement(center_matrix, row, column, value);

            row    = center_index;
            column = left_index;
            value  = -coeff1 * arr; /* Left */
            updateTridiagonalElement(center_matrix, row, column, value);

            /* Remark: Right is not included here due to the symmetry shift */

            row    = center_index;
            column = center_index;
            value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
            updateTridiagonalElement(center_matrix, row, column, value);

            /* Fill matrix row of (i-1,j) */
            row    = left_index;
            column = center_index;
            value  = -coeff1 * arr; /* Right */
            updateTridiagonalElement(center_matrix, row, column, value);

            row    = left_index;
            column = left_index;
            value  = coeff1 * arr; /* Center: (Right) */
            updateTridiagonalElement(center_matrix, row, column, value);

            /* Fill matrix row of (i+1,j) */
            /* Nothing to be done */

            /* Fill matrix row of (i,j-1) */
            row    = bottom_index;
            column = bottom_index;
            value  = coeff3 * att; /* Center: (Top) */
            updateDiagonalElement(bottom_matrix, row, column, value);

            /* Fill matrix row of (i,j+1) */
            row    = top_index;
            column = top_index;
            value  = coeff4 * att; /* Center: (Bottom) */
            updateDiagonalElement(top_matrix, row, column, value);
        }
        else {
            /* i_theta % 2 == 0 */
            /* ---------------|| */
            /* o   o   o   o  || */
            /* ---------------|| */
            /* o   x   O   x  || */
            /* ---------------|| */
            /* o   o   o   o  || */
            /* ---------------|| */
            auto& center_matrix = radial_diagonal_solver[i_theta / 2];
            auto& bottom_matrix = radial_tridiagonal_solver[i_theta_M1 / 2];
            auto& top_matrix    = radial_tridiagonal_solver[i_theta_P1 / 2];
            /* Fill matrix row of (i,j) */
            row    = center_index;
            column = center_index;
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
            updateDiagonalElement(center_matrix, row, column, value);

            row    = center_index;
            column = center_index;
            value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
            updateDiagonalElement(center_matrix, row, column, value);

            /* Fill matrix row of (i,j-1) */
            row    = bottom_index;
            column = bottom_index;
            value  = coeff3 * att; /* Center: (Top) */
            updateTridiagonalElement(bottom_matrix, row, column, value);

            /* Fill matrix row of (i,j+1) */
            row    = top_index;
            column = top_index;
            value  = coeff4 * att; /* Center: (Bottom) */
            updateTridiagonalElement(top_matrix, row, column, value);
        }
    }
    /* ------------------------------------------ */
    /* Radial Section: Node on the outer boundary */
    /* ------------------------------------------ */
    else if (i_r == grid.nr() - 1) {
        assert(!i_r % 2 == 0);

        double h1 = grid.radialSpacing(i_r - 1);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;

        int center_index = i_r - numberSmootherCircles;
        int left_index   = i_r - numberSmootherCircles - 1;

        if (i_theta & 1) {
            /* i_theta % 2 == 1 */
            /* -----------|| */
            /* x   o   x  || */
            /* -----------|| */
            /* o   o   O  || */
            /* -----------|| */
            /* x   o   x  || */
            /* -----------|| */
            auto& center_matrix = radial_tridiagonal_solver[i_theta / 2];

            /* Fill matrix row of (i,j) */
            row    = center_index;
            column = center_index;
            value  = 1.0;
            updateTridiagonalElement(center_matrix, row, column, value);

            /* Fill matrix row of (i-1,j) */
            row    = left_index;
            column = left_index;
            value  = coeff1 * arr; /* Center: (Right) */
            updateTridiagonalElement(center_matrix, row, column, value);
        }
        else {
            /* i_theta % 2 == 0 */
            /* -----------|| */
            /* o   o   o  || */
            /* -----------|| */
            /* x   o   X  || */
            /* -----------|| */
            /* o   o   o  || */
            /* -----------|| */
            auto& center_matrix = radial_diagonal_solver[i_theta / 2];

            /* Fill matrix row of (i,j) */
            row    = center_index;
            column = center_index;
            value  = 1.0;
            updateDiagonalElement(center_matrix, row, column, value);

            /* Fill matrix row of (i-1,j) */
            row    = left_index;
            column = left_index;
            value  = coeff1 * arr; /* Center: (Right) */
            updateDiagonalElement(center_matrix, row, column, value);
        }
    }
}

void ExtrapolatedSmootherGive::buildAscCircleSection(const int i_r)
{
    const double r = grid_.radius(i_r);
    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        const int global_index = grid_.index(i_r, i_theta);
        const double theta     = grid_.theta(i_theta);

        double coeff_beta, arr, att, art, detDF;
        level_cache_.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art, detDF);

        // Build Asc at the current node
        nodeBuildSmootherGive(i_r, i_theta, grid_, DirBC_Interior_, inner_boundary_circle_matrix_,
                              circle_diagonal_solver_, radial_diagonal_solver_, circle_tridiagonal_solver_,
                              radial_tridiagonal_solver_, arr, att, art, detDF, coeff_beta);
    }
}

void ExtrapolatedSmootherGive::buildAscRadialSection(const int i_theta)
{
    const double theta = grid_.theta(i_theta);
    for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
        const int global_index = grid_.index(i_r, i_theta);
        const double r         = grid_.radius(i_r);

        double coeff_beta, arr, att, art, detDF;
        level_cache_.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art, detDF);

        // Build Asc at the current node
        nodeBuildSmootherGive(i_r, i_theta, grid_, DirBC_Interior_, inner_boundary_circle_matrix_,
                              circle_diagonal_solver_, radial_diagonal_solver_, circle_tridiagonal_solver_,
                              radial_tridiagonal_solver_, arr, att, art, detDF, coeff_beta);
    }
}

// clang-format off
void ExtrapolatedSmootherGive::buildAscMatrices()
{
    /* -------------------------------------- */
    /* Part 1: Allocate Asc Smoother matrices */
    /* -------------------------------------- */

    const int number_smoother_circles = grid_.numberSmootherCircles();
    const int length_smoother_radial  = grid_.lengthSmootherRadial();

    const int num_circle_nodes = grid_.ntheta();
    circle_tridiagonal_solver_.resize(number_smoother_circles / 2);
    circle_diagonal_solver_.resize(number_smoother_circles - number_smoother_circles / 2);

    assert((grid_.ntheta() / 2) % 2 == 0);
    const int num_radial_nodes = length_smoother_radial;
    radial_tridiagonal_solver_.resize(grid_.ntheta() / 2);
    radial_diagonal_solver_.resize(grid_.ntheta() / 2);

    // Remark: circle_diagonal_solver_[0] is undefnied.
    // Use inner_boundary_circle_matrix_ instead.
    #pragma omp parallel num_threads(num_omp_threads_) if (grid_.numberOfNodes() > 10'000)
    {
        // ---------------- //
        // Circular Section //
        #pragma omp for nowait
        for (int circle_Asc_index = 0; circle_Asc_index < number_smoother_circles; circle_Asc_index++) {

            /* Inner boundary circle */
            if (circle_Asc_index == 0) {
                #ifdef GMGPOLAR_USE_MUMPS
                // Although the matrix is symmetric, we need to store all its entries, so we disable the symmetry.
                const int nnz                 = getNonZeroCountCircleAsc(circle_Asc_index);
                inner_boundary_circle_matrix_ = SparseMatrixCOO<double>(num_circle_nodes, num_circle_nodes, nnz);
                inner_boundary_circle_matrix_.is_symmetric(false);
                for (int i = 0; i < nnz; i++) {
                    inner_boundary_circle_matrix_.value(i) = 0.0;
                }
                #else
                std::function<int(int)> nnz_per_row = [&](int i_theta) {
                    if(DirBC_Interior_) return 1;
                    else return i_theta % 2 == 0 ? 1 : 2;
                };
                inner_boundary_circle_matrix_ = SparseMatrixCSR<double>(num_circle_nodes, num_circle_nodes, nnz_per_row);
                for (int i = 0; i < inner_boundary_circle_matrix_.non_zero_size(); i++) {
                    inner_boundary_circle_matrix_.values_data()[i] = 0.0;
                }
                #endif
            }

            /* Interior Circle Section */
            else {
                if (circle_Asc_index & 1) {
                    const int circle_tridiagonal_solver_index = circle_Asc_index / 2;
                    auto& solver_matrix = circle_tridiagonal_solver_[circle_tridiagonal_solver_index];

                    solver_matrix = SymmetricTridiagonalSolver<double>(num_circle_nodes);
                    solver_matrix.is_cyclic(true);
                    solver_matrix.cyclic_corner_element() = 0.0;

                    for (int i = 0; i < num_circle_nodes; i++) {
                        solver_matrix.main_diagonal(i) = 0.0;
                        if (i < num_circle_nodes - 1) {
                            solver_matrix.sub_diagonal(i) = 0.0;
                        }
                    }
                }
                else {
                    const int circle_diagonal_solver_index = circle_Asc_index / 2;
                    auto& solver_matrix                    = circle_diagonal_solver_[circle_diagonal_solver_index];

                    solver_matrix = DiagonalSolver<double>(num_circle_nodes);

                    for (int i = 0; i < num_circle_nodes; i++) {
                        solver_matrix.diagonal(i) = 0.0;
                    }
                }
            }
        }

        // -------------- //
        // Radial Section //
        #pragma omp for nowait
        for (int radial_Asc_index = 0; radial_Asc_index < grid_.ntheta(); radial_Asc_index++) {

            if (radial_Asc_index & 1) {
                const int radial_tridiagonal_solver_index = radial_Asc_index / 2;
                auto& solver_matrix                       = radial_tridiagonal_solver_[radial_tridiagonal_solver_index];

                solver_matrix = SymmetricTridiagonalSolver<double>(num_radial_nodes);
                solver_matrix.is_cyclic(false);

                for (int i = 0; i < num_radial_nodes; i++) {
                    solver_matrix.main_diagonal(i) = 0.0;
                    if (i < num_radial_nodes - 1) {
                        solver_matrix.sub_diagonal(i) = 0.0;
                    }
                }
            }
            else {
                const int radial_diagonal_solver_index = radial_Asc_index / 2;
                auto& solver_matrix                    = radial_diagonal_solver_[radial_diagonal_solver_index];

                solver_matrix = DiagonalSolver<double>(num_radial_nodes);

                for (int i = 0; i < num_radial_nodes; i++) {
                    solver_matrix.diagonal(i) = 0.0;
                }
            }
        }
    }

    /* ---------------------------------- */
    /* Part 2: Fill Asc Smoother matrices */
    /* ---------------------------------- */

    if (num_omp_threads_ == 1) {
        /* Single-threaded execution */
        for (int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
            buildAscCircleSection(i_r);
        }
        for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
            buildAscRadialSection(i_theta);
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
                buildAscCircleSection(i_r);
            }
            #pragma omp for
            for (int circle_task = 1; circle_task < num_circle_tasks; circle_task += 3) {
                int i_r = grid_.numberSmootherCircles() - circle_task - 1;
                buildAscCircleSection(i_r);
            }
            #pragma omp for nowait
            for (int circle_task = 2; circle_task < num_circle_tasks; circle_task += 3) {
                int i_r = grid_.numberSmootherCircles() - circle_task - 1;
                buildAscCircleSection(i_r);
            }

            #pragma omp for
            for (int radial_task = 0; radial_task < num_radial_tasks; radial_task += 3) {
                if (radial_task > 0) {
                    int i_theta = radial_task + additional_radial_tasks;
                    buildAscRadialSection(i_theta);
                }
                else {
                    if (additional_radial_tasks == 0) {
                        buildAscRadialSection(0);
                    }
                    else if (additional_radial_tasks >= 1) {
                        buildAscRadialSection(0);
                        buildAscRadialSection(1);
                    }
                }
            }
            #pragma omp for
            for (int radial_task = 1; radial_task < num_radial_tasks; radial_task += 3) {
                if (radial_task > 1) {
                    int i_theta = radial_task + additional_radial_tasks;
                    buildAscRadialSection(i_theta);
                }
                else {
                    if (additional_radial_tasks == 0) {
                        buildAscRadialSection(1);
                    }
                    else if (additional_radial_tasks == 1) {
                        buildAscRadialSection(2);
                    }
                    else if (additional_radial_tasks == 2) {
                        buildAscRadialSection(2);
                        buildAscRadialSection(3);
                    }
                }
            }
            #pragma omp for
            for (int radial_task = 2; radial_task < num_radial_tasks; radial_task += 3) {
                int i_theta = radial_task + additional_radial_tasks;
                buildAscRadialSection(i_theta);
            }
        }
    }

    #ifdef GMGPOLAR_USE_MUMPS
    /* ------------------------------------------------------------------- */
    /* Part 3: Convert inner_boundary_circle_matrix_ to a symmetric matrix */
    /* ------------------------------------------------------------------- */

    SparseMatrixCOO<double> full_matrix = std::move(inner_boundary_circle_matrix_);

    const int nnz           = full_matrix.non_zero_size();
    const int numRows       = full_matrix.rows();
    const int numColumns    = full_matrix.columns();
    const int symmetric_nnz = nnz - (nnz - numRows) / 2;

    inner_boundary_circle_matrix_ = SparseMatrixCOO<double>(numRows, numColumns, symmetric_nnz);
    inner_boundary_circle_matrix_.is_symmetric(true);

    int current_nz = 0;
    for (int nz_index = 0; nz_index < full_matrix.non_zero_size(); nz_index++) {
        int current_row = full_matrix.row_index(nz_index);
        int current_col = full_matrix.col_index(nz_index);
        if (current_row <= current_col) {
            inner_boundary_circle_matrix_.row_index(current_nz) = current_row;
            inner_boundary_circle_matrix_.col_index(current_nz) = current_col;
            inner_boundary_circle_matrix_.value(current_nz)     = std::move(full_matrix.value(nz_index));
            current_nz++;
        }
    }
    #endif
}
// clang-format on
