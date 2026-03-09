#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherGive/extrapolatedSmootherGive.h"

#include "../../../include/Definitions/geometry_helper.h"

#ifdef GMGPOLAR_USE_MUMPS
// When using the MUMPS solver, the matrix is assembled in COO format.
static inline void updateMatrixElement(SparseMatrixCOO<double>& matrix, int ptr, int offset, int row, int column,
                                       double value)
{
    matrix.row_index(ptr + offset) = row;
    matrix.col_index(ptr + offset) = column;
    matrix.value(ptr + offset) += value;
}
#else
// When using the in-house solver, the matrix is stored in CSR format.
static inline void updateMatrixElement(SparseMatrixCSR<double>& matrix, int ptr, int offset, int row, int column,
                                       double value)
{
    matrix.row_nz_index(row, offset) = column;
    matrix.row_nz_entry(row, offset) += value;
}
#endif

void ExtrapolatedSmootherGive::nodeBuildInteriorBoundarySolverMatrix_i_r_0(int i_theta, const PolarGrid& grid,
                                                                           bool DirBC_Interior, MatrixType& matrix,
                                                                           double arr, double att, double art,
                                                                           double detDF, double coeff_beta)
{
    assert(i_theta >= 0 && i_theta < grid.ntheta());

    int ptr, offset;
    int row, column;
    double value;

    const int i_r = 0;

    /* ------------------------------------------------ */
    /* Case 1: Dirichlet boundary on the inner boundary */
    /* ------------------------------------------------ */
    if (DirBC_Interior) {
        /* Fill result(i,j) */
        int center_index = i_theta;

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = getCircleAscIndex(i_r, i_theta);

        const Stencil& CenterStencil = getStencil(i_r, i_theta);

        offset = CenterStencil[StencilPosition::Center];
        column = center_index;
        value  = 1.0;
        updateMatrixElement(matrix, ptr, offset, row, column, value);
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

            /* Fill matrix row of (i,j) */
            row = center_index;
            ptr = center_nz_index;

            offset = CenterStencil[StencilPosition::Center];
            column = center_index;
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* beta_{i,j} */
            updateMatrixElement(matrix, ptr, offset, row, column, value);

            offset = CenterStencil[StencilPosition::Left];
            column = left_index;
            value  = -coeff1 * arr; /* Left */
            updateMatrixElement(matrix, ptr, offset, row, column, value);

            offset = CenterStencil[StencilPosition::Center];
            column = center_index;
            value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
            updateMatrixElement(matrix, ptr, offset, row, column, value);

            /* Fill matrix row of (i-1,j) */
            /* From view the view of the across origin node, */
            /* the directions are roatated by 180 degrees in the stencil! */
            row = left_index;
            ptr = left_nz_index;

            const Stencil& LeftStencil = CenterStencil;

            offset = LeftStencil[StencilPosition::Left];
            column = center_index;
            value  = -coeff1 * arr; /* Right -> Left*/
            updateMatrixElement(matrix, ptr, offset, row, column, value);

            offset = LeftStencil[StencilPosition::Center];
            column = left_index;
            value  = +coeff1 * arr; /* Center: (Right) -> Center: (Left) */
            updateMatrixElement(matrix, ptr, offset, row, column, value);
        }
        else {
            /* i_theta % 2 == 0 */
            /* -| o | o | o | */
            /* -|   |   |   | */
            /* -| X | o | x | */
            /* -|   |   |   | */
            /* -| o | o | o | */

            /* Fill matrix row of (i,j) */
            row = center_index;
            ptr = center_nz_index;

            offset = CenterStencil[StencilPosition::Center];
            column = center_index;
            value  = 1.0;
            updateMatrixElement(matrix, ptr, offset, row, column, value);

            /* Fill matrix row of (i,j-1) */
            row = bottom_index;
            ptr = bottom_nz_index;

            const Stencil& BottomStencil = CenterStencil;

            offset = BottomStencil[StencilPosition::Center];
            column = bottom_index;
            value  = +coeff3 * att; /* Center: (Top) */
            updateMatrixElement(matrix, ptr, offset, row, column, value);

            /* Fill matrix row of (i,j+1) */
            row = top_index;
            ptr = top_nz_index;

            const Stencil& TopStencil = CenterStencil;

            offset = TopStencil[StencilPosition::Center];
            column = top_index;
            value  = +coeff4 * att; /* Center: (Bottom) */
            updateMatrixElement(matrix, ptr, offset, row, column, value);
        }
    }
}

void ExtrapolatedSmootherGive::nodeBuildInteriorBoundarySolverMatrix_i_r_1(int i_theta, const PolarGrid& grid,
                                                                           bool DirBC_Interior, MatrixType& matrix,
                                                                           double arr, double att, double art,
                                                                           double detDF, double coeff_beta)
{
    assert(i_theta >= 0 && i_theta < grid.ntheta());

    int ptr, offset;
    int row, column;
    double value;

    const int i_r = 1;

    const double h1 = grid.radialSpacing(i_r - 1);
    const double h2 = grid.radialSpacing(i_r);
    const double k1 = grid.angularSpacing(i_theta - 1);
    const double k2 = grid.angularSpacing(i_theta);

    const double coeff1 = 0.5 * (k1 + k2) / h1;

    const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
    const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

    const int left_index = i_theta;

    /* -------------------------- */
    /* Cyclic Tridiagonal Section */
    /* i_r % 2 == 1               */
    if (i_r & 1) {
        if (i_theta & 1) {
            /* i_theta % 2 == 1 */
            /* | x | o | x | */
            /* |   |   |   | */
            /* | o | O | o | */
            /* |   |   |   | */
            /* | x | o | x | */

            /* Fill matrix row of (i-1,j) */

            /* Only in the case of AcrossOrigin */
            if (!DirBC_Interior) {
                row = left_index;
                ptr = getCircleAscIndex(i_r - 1, i_theta);

                const Stencil& LeftStencil = getStencil(i_r - 1, i_theta);

                offset = LeftStencil[StencilPosition::Center];
                column = left_index;
                value  = +coeff1 * arr; /* Center: (Right) */
                updateMatrixElement(matrix, ptr, offset, row, column, value);
            }
        }
    }
}

ExtrapolatedSmootherGive::MatrixType ExtrapolatedSmootherGive::buildInteriorBoundarySolverMatrix()
{
    const int ntheta = grid_.ntheta();

#ifdef GMGPOLAR_USE_MUMPS
    // Although the matrix is symmetric, we need to store all its entries, so we disable the symmetry.
    const int nnz = getNonZeroCountCircleAsc(0);
    SparseMatrixCOO<double> inner_boundary_solver_matrix(ntheta, ntheta, nnz);
    inner_boundary_solver_matrix.is_symmetric(false);
#else
    std::function<int(int)> nnz_per_row = [&](int i_theta) {
        if (DirBC_Interior_)
            return 1;
        else
            return i_theta % 2 == 0 ? 1 : 2;
    };
    SparseMatrixCSR<double> inner_boundary_solver_matrix(ntheta, ntheta, nnz_per_row);
#endif

    {
        const int i_r  = 0;
        const double r = grid_.radius(i_r);
        for (int i_theta = 0; i_theta < ntheta; i_theta++) {
            {
                const int global_index = grid_.index(i_r, i_theta);
                const double theta     = grid_.theta(i_theta);

                double coeff_beta, arr, att, art, detDF;
                level_cache_.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art, detDF);

                nodeBuildInteriorBoundarySolverMatrix_i_r_0(
                    i_theta, grid_, DirBC_Interior_, inner_boundary_solver_matrix, arr, att, art, detDF, coeff_beta);
            }
        }
    }

    {
        const int i_r  = 1;
        const double r = grid_.radius(i_r);
        for (int i_theta = 0; i_theta < ntheta; i_theta++) {
            {
                const int global_index = grid_.index(i_r, i_theta);
                const double theta     = grid_.theta(i_theta);

                double coeff_beta, arr, att, art, detDF;
                level_cache_.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art, detDF);

                nodeBuildInteriorBoundarySolverMatrix_i_r_1(
                    i_theta, grid_, DirBC_Interior_, inner_boundary_solver_matrix, arr, att, art, detDF, coeff_beta);
            }
        }
    }

#ifdef GMGPOLAR_USE_MUMPS
    /* Mumps: In the case of symmetric matrices, only half of the matrix should be provided. */
    const bool construct_symmetric = true;
    if (!construct_symmetric) {
        return inner_boundary_solver_matrix;
    }

    const int full_nnz      = inner_boundary_solver_matrix.non_zero_size();
    const int numRows       = inner_boundary_solver_matrix.rows();
    const int numColumns    = inner_boundary_solver_matrix.columns();
    const int symmetric_nnz = full_nnz - (full_nnz - numRows) / 2;

    SparseMatrixCOO<double> inner_boundary_solver_matrix_symmetric(numRows, numColumns, symmetric_nnz);
    inner_boundary_solver_matrix_symmetric.is_symmetric(true);

    int current_nz = 0; // Current non-zero index in the symmetric matrix
    for (int nz_index = 0; nz_index < full_nnz; nz_index++) {
        const int current_row    = inner_boundary_solver_matrix.row_index(nz_index);
        const int current_column = inner_boundary_solver_matrix.col_index(nz_index);
        if (current_row <= current_column) {
            inner_boundary_solver_matrix_symmetric.row_index(current_nz) = current_row;
            inner_boundary_solver_matrix_symmetric.col_index(current_nz) = current_column;
            inner_boundary_solver_matrix_symmetric.value(current_nz)     = inner_boundary_solver_matrix.value(nz_index);
            current_nz++;
        }
    }
    return inner_boundary_solver_matrix_symmetric;
#else
    return inner_boundary_solver_matrix;
#endif
}
