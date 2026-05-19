#pragma once

#include "smootherStencil.inl"

namespace extrapolated_smoother_take
{

#ifdef GMGPOLAR_USE_MUMPS
// When using the MUMPS solver, the matrix is assembled in COO format.
static KOKKOS_INLINE_FUNCTION void
update_CSR_COO_MatrixElement(const SparseMatrixCOO<double>& matrix, const int ptr, const int offset,
                             const int row, const int column, const double value)
{
    matrix.set_row_index(ptr + offset, row);
    matrix.set_col_index(ptr + offset, column);
    matrix.set_value(ptr + offset, value);
}
#else
// When using the in-house solver, the matrix is stored in CSR format.
static KOKKOS_INLINE_FUNCTION void
update_CSR_COO_MatrixElement(const SparseMatrixCSR<double>& matrix, const int ptr, const int offset,
                             const int row, const int column, const double value)
{
    matrix.set_row_nz_index(row, offset, column);
    matrix.set_row_nz_entry(row, offset, value);
}
#endif

template <class InnerBoundaryMatrix>
static KOKKOS_INLINE_FUNCTION void
nodeBuildInteriorBoundarySolverMatrix(const int i_theta, const PolarGrid<DefaultMemorySpace>& grid, const bool DirBC_Interior,
                                      const InnerBoundaryMatrix& matrix, ConstVector<double>& arr,
                                      ConstVector<double>& att, ConstVector<double>& art,
                                      ConstVector<double>& detDF, ConstVector<double>& coeff_beta)
{
    using extrapolated_smoother_take::getCircleAscIndex;
    using extrapolated_smoother_take::getStencil;
    using extrapolated_smoother_take::update_CSR_COO_MatrixElement;

    KOKKOS_ASSERT(i_theta >= 0 && i_theta < grid.ntheta());

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
        const int center_nz_index = getCircleAscIndex(i_r, i_theta, grid, DirBC_Interior);

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r, i_theta, grid, DirBC_Interior);

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
        const double coeff5 = 0.25 * (h1 + h2) * (k1 + k2);

        const int i_theta_M1           = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1           = grid.wrapThetaIndex(i_theta + 1);
        const int i_theta_AcrossOrigin = grid.wrapThetaIndex(i_theta + (grid.ntheta() / 2));

        const int center_index = i_theta;
        const int left_index   = i_theta_AcrossOrigin;
        const int right_index  = i_theta;
        const int bottom_index = i_theta_M1;
        const int top_index    = i_theta_P1;

        const int center_nz_index = getCircleAscIndex(i_r, i_theta, grid, DirBC_Interior);
        const int bottom_nz_index = getCircleAscIndex(i_r, i_theta_M1, grid, DirBC_Interior);
        const int top_nz_index    = getCircleAscIndex(i_r, i_theta_P1, grid, DirBC_Interior);
        const int left_nz_index   = getCircleAscIndex(i_r, i_theta_AcrossOrigin, grid, DirBC_Interior);

        int nz_index;
        const Stencil& CenterStencil = getStencil(i_r, i_theta, grid, DirBC_Interior);

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

            const double center_value = coeff5 * coeff_beta[center] * Kokkos::fabs(detDF[center]) +
                                        coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +
                                        coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);
            const double left_value = -coeff1 * (arr[center] + arr[left]);

            /* Fill matrix row of (i,j) */
            row = center_index;
            ptr = center_nz_index;

            const Stencil& CenterStencil = getStencil(i_r, i_theta, grid, DirBC_Interior);

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

            const Stencil& CenterStencil = getStencil(i_r, i_theta, grid, DirBC_Interior);

            offset = CenterStencil[StencilPosition::Center];
            column = center_index;
            value  = 1.0;
            update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);
        }
    }
}

} // namespace extrapolated_smoother_take

template <class LevelCacheType>
typename ExtrapolatedSmootherTake<LevelCacheType>::InnerBoundaryMatrix
ExtrapolatedSmootherTake<LevelCacheType>::buildInteriorBoundarySolverMatrix()
{
    using extrapolated_smoother_take::getNonZeroCountCircleAsc;
    using extrapolated_smoother_take::nodeBuildInteriorBoundarySolverMatrix;

    const PolarGrid<DefaultMemorySpace>& grid             = ExtrapolatedSmoother<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = ExtrapolatedSmoother<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = ExtrapolatedSmoother<LevelCacheType>::DirBC_Interior_;

    const int i_r    = 0;
    const int ntheta = grid.ntheta();

    // The interior boundary matrix is symmetric due to the periodicity in the theta direction
    // and the assumption that ntheta is even, which is required for the across-origin discretization.
    // We store all non-zero entries of the matrix, both in COO format (for MUMPS)
    // and in CSR format (for the in-house solver). If the COO matrix is marked as symmetric,
    // the COO_Mumps_Solver optimizes the factorization by only using the upper triangular part of the matrix,
    // which is extracted by the COO_Mumps_Solver internally.
#ifdef GMGPOLAR_USE_MUMPS
    const int nnz = getNonZeroCountCircleAsc(i_r, grid, DirBC_Interior);
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

    Kokkos::parallel_for(
        "ExtrapolatedSmootherTake: BuildInnerBoundaryMatrix",
        Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, ntheta), KOKKOS_LAMBDA(const int i_theta) {
            nodeBuildInteriorBoundarySolverMatrix(i_theta, grid, DirBC_Interior, inner_boundary_solver_matrix, arr, att,
                                                  art, detDF, coeff_beta);
        });

    Kokkos::fence();

    return inner_boundary_solver_matrix;
}
