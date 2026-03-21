#pragma once

namespace smoother_give
{

#ifdef GMGPOLAR_USE_MUMPS
// When using the MUMPS solver, the matrix is assembled in COO format.
static inline void update_CSR_COO_MatrixElement(SparseMatrixCOO<double>& matrix, int ptr, int offset, int row,
                                                int column, double value)
{
    matrix.row_index(ptr + offset) = row;
    matrix.col_index(ptr + offset) = column;
    matrix.value(ptr + offset) += value;
}
#else
// When using the in-house solver, the matrix is stored in CSR format.
static inline void update_CSR_COO_MatrixElement(SparseMatrixCSR<double>& matrix, int ptr, int offset, int row,
                                                int column, double value)
{
    matrix.row_nz_index(row, offset) = column;
    matrix.row_nz_entry(row, offset) += value;
}
#endif

} // namespace smoother_give

template <concepts::DomainGeometry DomainGeometry>
void SmootherGive<DomainGeometry>::nodeBuildInteriorBoundarySolverMatrix_i_r_0(int i_theta, const PolarGrid& grid,
                                                                               bool DirBC_Interior,
                                                                               InnerBoundaryMatrix& matrix, double arr,
                                                                               double att, double art, double detDF,
                                                                               double coeff_beta)
{
    using smoother_give::update_CSR_COO_MatrixElement;

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
        const int center_index = i_theta;

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = getCircleAscIndex(i_r, i_theta);

        const Stencil& CenterStencil = getStencil(i_r);

        offset = CenterStencil[StencilPosition::Center];
        column = center_index;
        value  = 1.0;
        update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);
    }
    else {
        /* ------------------------------------------------------------- */
        /* Case 2: Across origin discretization on the interior boundary */
        /* ------------------------------------------------------------- */
        // h1 gets replaced with 2 * R0.
        // (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()>>1)).
        // Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil.
        const double h1 = 2 * grid.radius(0);
        const double h2 = grid.radialSpacing(i_r);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff1 = 0.5 * (k1 + k2) / h1;
        const double coeff2 = 0.5 * (k1 + k2) / h2;
        const double coeff3 = 0.5 * (h1 + h2) / k1;
        const double coeff4 = 0.5 * (h1 + h2) / k2;

        /* -| x | o | x | */
        /* -|   |   |   | */
        /* -| O | o | o | */
        /* -|   |   |   | */
        /* -| x | o | x | */
        const int i_theta_M1           = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1           = grid.wrapThetaIndex(i_theta + 1);
        const int i_theta_AcrossOrigin = grid.wrapThetaIndex(i_theta + grid.ntheta() / 2);

        const int center_index = i_theta;
        const int left_index   = i_theta_AcrossOrigin;
        const int bottom_index = i_theta_M1;
        const int top_index    = i_theta_P1;

        const int center_nz_index = getCircleAscIndex(i_r, i_theta);
        const int bottom_nz_index = getCircleAscIndex(i_r, i_theta_M1);
        const int top_nz_index    = getCircleAscIndex(i_r, i_theta_P1);
        const int left_nz_index   = getCircleAscIndex(i_r, i_theta_AcrossOrigin);

        int nz_index; /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r);

        offset = CenterStencil[StencilPosition::Center];
        column = center_index;
        value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* beta_{i,j} */
        update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Left];
        column = left_index;
        value  = -coeff1 * arr; /* Left */
        update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Bottom];
        column = bottom_index;
        value  = -coeff3 * att; /* Bottom */
        update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Top];
        column = top_index;
        value  = -coeff4 * att; /* Top */
        update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Center];
        column = center_index;
        value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
        update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

        /* Fill matrix row of (i-1,j) */
        /* From view the view of the across origin node, */ /* the directions are roatated by 180 degrees in the stencil! */
        row = left_index;
        ptr = left_nz_index;

        const Stencil& LeftStencil = CenterStencil;

        offset = LeftStencil[StencilPosition::Left];
        column = center_index;
        value  = -coeff1 * arr; /* Right -> Left*/
        update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

        offset = LeftStencil[StencilPosition::Center];
        column = left_index;
        value  = +coeff1 * arr; /* Center: (Right) -> Center: (Left) */
        update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

        /* Top Right -> Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */

        /* Bottom Right -> Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */

        /* Fill matrix row of (i,j-1) */
        row = bottom_index;
        ptr = bottom_nz_index;

        const Stencil& BottomStencil = CenterStencil;

        offset = BottomStencil[StencilPosition::Top];
        column = center_index;
        value  = -coeff3 * att; /* Top */
        update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

        offset = BottomStencil[StencilPosition::Center];
        column = bottom_index;
        value  = +coeff3 * att; /* Center: (Top) */
        update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

        /* TopLeft: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */

        /* Fill matrix row of (i,j+1) */
        row = top_index;
        ptr = top_nz_index;

        const Stencil& TopStencil = CenterStencil;

        offset = TopStencil[StencilPosition::Bottom];
        column = center_index;
        value  = -coeff4 * att; /* Bottom */
        update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

        offset = TopStencil[StencilPosition::Center];
        column = top_index;
        value  = +coeff4 * att; /* Center: (Bottom) */
        update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

        /* BottomLeft: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */
    }
}

template <concepts::DomainGeometry DomainGeometry>
void SmootherGive<DomainGeometry>::nodeBuildInteriorBoundarySolverMatrix_i_r_1(int i_theta, const PolarGrid& grid,
                                                                               bool DirBC_Interior,
                                                                               InnerBoundaryMatrix& matrix, double arr,
                                                                               double att, double art, double detDF,
                                                                               double coeff_beta)
{
    using smoother_give::update_CSR_COO_MatrixElement;

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

    /* Fill matrix row of (i-1,j) */
    if (!DirBC_Interior) {
        row = left_index;
        ptr = getCircleAscIndex(i_r - 1, i_theta);

        const Stencil& LeftStencil = getStencil(i_r - 1);

        offset = LeftStencil[StencilPosition::Center];
        column = left_index;
        value  = coeff1 * arr; /* Center: (Right) */
        update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);
    }
}

template <concepts::DomainGeometry DomainGeometry>
typename SmootherGive<DomainGeometry>::InnerBoundaryMatrix
SmootherGive<DomainGeometry>::buildInteriorBoundarySolverMatrix()
{
    const PolarGrid& grid                         = Smoother<DomainGeometry>::grid_;
    const LevelCache<DomainGeometry>& level_cache = Smoother<DomainGeometry>::level_cache_;
    const bool DirBC_Interior                     = Smoother<DomainGeometry>::DirBC_Interior_;

    const int i_r    = 0;
    const int ntheta = grid.ntheta();

#ifdef GMGPOLAR_USE_MUMPS
    // Although the matrix is symmetric, we need to store all its entries, so we disable the symmetry.
    const int nnz = getNonZeroCountCircleAsc(i_r);
    SparseMatrixCOO<double> inner_boundary_solver_matrix(ntheta, ntheta, nnz);
    inner_boundary_solver_matrix.is_symmetric(true);
#else
    // The stencils size for the inner boundary matrix is either 1 (Dirichlet BC) or 4 (across-origin discretization).
    std::function<int(int)> nnz_per_row = [&](int i_theta) {
        return DirBC_Interior ? 1 : 4;
    };
    SparseMatrixCSR<double> inner_boundary_solver_matrix(ntheta, ntheta, nnz_per_row);
#endif

    {
        const int i_r  = 0;
        const double r = grid.radius(i_r);
        for (int i_theta = 0; i_theta < ntheta; i_theta++) {
            {
                const int global_index = grid.index(i_r, i_theta);
                const double theta     = grid.theta(i_theta);

                double coeff_beta, arr, att, art, detDF;
                level_cache.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art, detDF);

                nodeBuildInteriorBoundarySolverMatrix_i_r_0(i_theta, grid, DirBC_Interior, inner_boundary_solver_matrix,
                                                            arr, att, art, detDF, coeff_beta);
            }
        }
    }

    {
        const int i_r  = 1;
        const double r = grid.radius(i_r);
        for (int i_theta = 0; i_theta < ntheta; i_theta++) {
            {
                const int global_index = grid.index(i_r, i_theta);
                const double theta     = grid.theta(i_theta);

                double coeff_beta, arr, att, art, detDF;
                level_cache.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art, detDF);

                nodeBuildInteriorBoundarySolverMatrix_i_r_1(i_theta, grid, DirBC_Interior, inner_boundary_solver_matrix,
                                                            arr, att, art, detDF, coeff_beta);
            }
        }
    }

    return inner_boundary_solver_matrix;
}
