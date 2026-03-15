#pragma once

namespace smoother_take
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

} // namespace smoother_take

template <concepts::DomainGeometry DomainGeometry>
void SmootherTake<DomainGeometry>::nodeBuildInteriorBoundarySolverMatrix(
    int i_theta, const PolarGrid& grid, bool DirBC_Interior, MatrixType& matrix, ConstVector<double>& arr,
    ConstVector<double>& att, ConstVector<double>& art, ConstVector<double>& detDF, ConstVector<double>& coeff_beta)
{
    using smoother_take::update_CSR_COO_MatrixElement;

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
        // (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()/2)).
        // Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil.
        const double h1 = 2.0 * grid.radius(0);
        const double h2 = grid.radialSpacing(i_r);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff1 = 0.5 * (k1 + k2) / h1;
        const double coeff2 = 0.5 * (k1 + k2) / h2;
        const double coeff3 = 0.5 * (h1 + h2) / k1;
        const double coeff4 = 0.5 * (h1 + h2) / k2;

        const int i_theta_M1           = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1           = grid.wrapThetaIndex(i_theta + 1);
        const int i_theta_AcrossOrigin = grid.wrapThetaIndex(i_theta + grid.ntheta() / 2);

        const int left   = grid.index(i_r, i_theta_AcrossOrigin);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int center = grid.index(i_r, i_theta);
        const int top    = grid.index(i_r, i_theta_P1);
        const int right  = grid.index(i_r + 1, i_theta);

        const int center_index = i_theta;
        const int left_index   = i_theta_AcrossOrigin;
        const int bottom_index = i_theta_M1;
        const int top_index    = i_theta_P1;

        const int center_nz_index = getCircleAscIndex(i_r, i_theta);

        const double left_value   = -coeff1 * (arr[center] + arr[left]);
        const double right_value  = -coeff2 * (arr[center] + arr[right]);
        const double bottom_value = -coeff3 * (att[center] + att[bottom]);
        const double top_value    = -coeff4 * (att[center] + att[top]);

        const double center_value = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center] * std::fabs(detDF[center]) -
                                    (left_value + right_value + bottom_value + top_value);

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r);

        offset = CenterStencil[StencilPosition::Center];
        column = center_index;
        value  = center_value;
        update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Left];
        column = left_index;
        value  = left_value;
        update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Bottom];
        column = bottom_index;
        value  = bottom_value;
        update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Top];
        column = top_index;
        value  = top_value;
        update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);
    }
}

template <concepts::DomainGeometry DomainGeometry>
typename SmootherTake<DomainGeometry>::MatrixType SmootherTake<DomainGeometry>::buildInteriorBoundarySolverMatrix()
{
    const PolarGrid& grid                         = Smoother<DomainGeometry>::grid_;
    const LevelCache<DomainGeometry>& level_cache = Smoother<DomainGeometry>::level_cache_;
    const bool DirBC_Interior                     = Smoother<DomainGeometry>::DirBC_Interior_;
    const int num_omp_threads                     = Smoother<DomainGeometry>::num_omp_threads_;

    const int i_r    = 0;
    const int ntheta = grid.ntheta();

#ifdef GMGPOLAR_USE_MUMPS
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
