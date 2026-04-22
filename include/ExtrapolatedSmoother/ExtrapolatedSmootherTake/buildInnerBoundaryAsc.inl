#pragma once

namespace extrapolated_smoother_take
{

#ifdef GMGPOLAR_USE_MUMPS
// When using the MUMPS solver, the matrix is assembled in COO format.
static inline void update_CSR_COO_MatrixElement(SparseMatrixCOO<double>& matrix, int ptr, int offset, int row,
                                                int column, double value)
{
    matrix.row_index(ptr + offset) = row;
    matrix.col_index(ptr + offset) = column;
    matrix.value(ptr + offset)     = value;
}
#else
// When using the in-house solver, the matrix is stored in CSR format.
static inline void update_CSR_COO_MatrixElement(SparseMatrixCSR<double>& matrix, int ptr, int offset, int row,
                                                int column, double value)
{
    matrix.row_nz_index(row, offset) = column;
    matrix.row_nz_entry(row, offset) = value;
}
#endif

} // namespace extrapolated_smoother_take

template <class LevelCacheType>
void ExtrapolatedSmootherTake<LevelCacheType>::nodeBuildInteriorBoundarySolverMatrix(
    int i_theta, const PolarGrid& grid, bool DirBC_Interior, InnerBoundaryMatrix& matrix, ConstVector<double>& arr,
    ConstVector<double>& att, ConstVector<double>& art, ConstVector<double>& detDF, ConstVector<double>& coeff_beta)
{
    using extrapolated_smoother_take::update_CSR_COO_MatrixElement;

    assert(i_theta >= 0 && i_theta < grid.ntheta());

    /* ------------------------------------------ */
    /* Circle Section: Node in the inner boundary */
    /* ------------------------------------------ */
    const int i_r = 0;

    int ptr, offset;
    int row, column;
    double value;

    /* ------------------------------------------------ */
    /* Case 1: Dirichlet boundary on the inner boundary */
    /* ------------------------------------------------ */
    if (DirBC_Interior) {
        const int center_index    = i_theta;
        const int center_nz_index = getCircleAscIndex(i_r, i_theta);

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r, i_theta);

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
        double coeff5 = 0.25 * (h1 + h2) * (k1 + k2);

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

            const int left   = grid.index(i_r, i_theta_AcrossOrigin);
            const int bottom = grid.index(i_r, i_theta_M1);
            const int center = grid.index(i_r, i_theta);
            const int top    = grid.index(i_r, i_theta_P1);
            const int right  = grid.index(i_r + 1, i_theta);

            const double center_value = coeff5 * coeff_beta[center] * fabs(detDF[center]) +
                                        coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +
                                        coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);
            const double left_value   = -coeff1 * (arr[center] + arr[left]);

            /* Fill matrix row of (i,j) */
            row = center_index;
            ptr = center_nz_index;

            const Stencil& CenterStencil = getStencil(i_r, i_theta);

            offset = CenterStencil[StencilPosition::Center];
            column = center_index;
            value  = center_value;
            update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

            offset = CenterStencil[StencilPosition::Left];
            column = left_index;
            value  = left_value;
            update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);
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

            const Stencil& CenterStencil = getStencil(i_r, i_theta);

            offset = CenterStencil[StencilPosition::Center];
            column = center_index;
            value  = 1.0;
            update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);
        }
    }
}

template <class LevelCacheType>
typename ExtrapolatedSmootherTake<LevelCacheType>::InnerBoundaryMatrix
ExtrapolatedSmootherTake<LevelCacheType>::buildInteriorBoundarySolverMatrix()
{
    const PolarGrid& grid             = ExtrapolatedSmootherTake<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = ExtrapolatedSmootherTake<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = ExtrapolatedSmootherTake<LevelCacheType>::DirBC_Interior_;
    const int num_omp_threads         = ExtrapolatedSmootherTake<LevelCacheType>::num_omp_threads_;

    const int i_r    = 0;
    const int ntheta = grid.ntheta();

#ifdef GMGPOLAR_USE_MUMPS
    const int nnz = getNonZeroCountCircleAsc(i_r);
    SparseMatrixCOO<double> inner_boundary_solver_matrix(ntheta, ntheta, nnz);
    inner_boundary_solver_matrix.is_symmetric(true);
#else
    std::function<int(int)> nnz_per_row = [&](int i_theta) {
        if (DirBC_Interior)
            return 1;
        else
            return i_theta % 2 == 0 ? 1 : 2;
    };
    SparseMatrixCSR<double> inner_boundary_solver_matrix(ntheta, ntheta, nnz_per_row);
#endif

    assert(level_cache.cacheDensityProfileCoefficients());
    assert(level_cache.cacheDomainGeometry());

    ConstVector<double> arr        = level_cache.arr();
    ConstVector<double> att        = level_cache.att();
    ConstVector<double> art        = level_cache.art();
    ConstVector<double> detDF      = level_cache.detDF();
    ConstVector<double> coeff_beta = level_cache.coeff_beta();

#pragma omp parallel for num_threads(num_omp_threads)
    for (int i_theta = 0; i_theta < ntheta; i_theta++) {
        nodeBuildInteriorBoundarySolverMatrix(i_theta, grid, DirBC_Interior, inner_boundary_solver_matrix, arr, att,
                                              art, detDF, coeff_beta);
    }

    return inner_boundary_solver_matrix;
}
