#include <ExtrapolatedSmoother/ExtrapolatedSmootherGive/extrapolatedSmootherGive.h>

#include <Definitions/geometry_helper.h>

#include "smootherStencil.inl"

namespace extrapolated_smoother_give
{

#ifdef GMGPOLAR_USE_MUMPS
// When using the MUMPS solver, the matrix is assembled in COO format.
static KOKKOS_INLINE_FUNCTION void
update_CSR_COO_MatrixElement(const SparseMatrixCOO<double, Kokkos::HostSpace>& matrix, const int ptr, const int offset,
                             const int row, const int column, const double value)
{
    matrix.set_row_index(ptr + offset, row);
    matrix.set_col_index(ptr + offset, column);
    matrix.increase_value(ptr + offset, value);
}
#else
// When using the in-house solver, the matrix is stored in CSR format.
static KOKKOS_INLINE_FUNCTION void
update_CSR_COO_MatrixElement(const SparseMatrixCSR<double, Kokkos::HostSpace>& matrix, const int ptr, const int offset,
                             const int row, const int column, const double value)
{
    matrix.set_row_nz_index(row, offset, column);
    matrix.increase_row_nz_entry(row, offset, value);
}
#endif

template <typename LevelCacheType, typename InnerBoundaryMatrix>
static KOKKOS_INLINE_FUNCTION void
nodeBuildInteriorBoundarySolverMatrix_i_r_0(const int i_theta, const PolarGrid& grid, const LevelCacheType& level_cache,
                                            const bool DirBC_Interior, const InnerBoundaryMatrix& matrix)
{
    using extrapolated_smoother_give::update_CSR_COO_MatrixElement;

    KOKKOS_ASSERT(i_theta >= 0 && i_theta < grid.ntheta());

    int ptr, offset;
    int row, column;
    double value;

    const int i_r = 0;

    /* ---------------------------------------- */
    /* Compute or retrieve stencil coefficients */
    /* ---------------------------------------- */
    const int center    = grid.index(i_r, i_theta);
    const double radius = grid.radius(i_r);
    const double theta  = grid.theta(i_theta);

    double coeff_beta, arr, att, art, detDF;
    level_cache.obtainValues(i_r, i_theta, center, radius, theta, coeff_beta, arr, att, art, detDF);

    /* ------------------------------------------------ */
    /* Case 1: Dirichlet boundary on the inner boundary */
    /* ------------------------------------------------ */
    if (DirBC_Interior) {
        /* Fill result(i,j) */
        int center_index = i_theta;

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = getCircleAscIndex(i_r, i_theta, grid, DirBC_Interior);

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
        const int bottom_index = i_theta_M1;
        const int top_index    = i_theta_P1;

        const int center_nz_index = getCircleAscIndex(i_r, i_theta, grid, DirBC_Interior);
        const int bottom_nz_index = getCircleAscIndex(i_r, i_theta_M1, grid, DirBC_Interior);
        const int top_nz_index    = getCircleAscIndex(i_r, i_theta_P1, grid, DirBC_Interior);
        const int left_nz_index   = getCircleAscIndex(i_r, i_theta_AcrossOrigin, grid, DirBC_Interior);

        const Stencil& CenterStencil = getStencil(i_r, i_theta, grid, DirBC_Interior);

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
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * Kokkos::fabs(detDF); /* beta_{i,j} */
            update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

            offset = CenterStencil[StencilPosition::Left];
            column = left_index;
            value  = -coeff1 * arr; /* Left */
            update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

            offset = CenterStencil[StencilPosition::Center];
            column = center_index;
            value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
            update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

            /* Fill matrix row of (i-1,j) */
            /* From view the view of the across origin node, */
            /* the directions are roatated by 180 degrees in the stencil! */
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
            update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

            /* Fill matrix row of (i,j-1) */
            row = bottom_index;
            ptr = bottom_nz_index;

            const Stencil& BottomStencil = CenterStencil;

            offset = BottomStencil[StencilPosition::Center];
            column = bottom_index;
            value  = +coeff3 * att; /* Center: (Top) */
            update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);

            /* Fill matrix row of (i,j+1) */
            row = top_index;
            ptr = top_nz_index;

            const Stencil& TopStencil = CenterStencil;

            offset = TopStencil[StencilPosition::Center];
            column = top_index;
            value  = +coeff4 * att; /* Center: (Bottom) */
            update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);
        }
    }
}

template <typename LevelCacheType, typename InnerBoundaryMatrix>
static KOKKOS_INLINE_FUNCTION void
nodeBuildInteriorBoundarySolverMatrix_i_r_1(const int i_theta, const PolarGrid& grid, const LevelCacheType& level_cache,
                                            const bool DirBC_Interior, const InnerBoundaryMatrix& matrix)
{
    using extrapolated_smoother_give::update_CSR_COO_MatrixElement;

    assert(i_theta >= 0 && i_theta < grid.ntheta());

    int ptr, offset;
    int row, column;
    double value;

    const int i_r = 1;

    /* ---------------------------------------- */
    /* Compute or retrieve stencil coefficients */
    /* ---------------------------------------- */
    const int center    = grid.index(i_r, i_theta);
    const double radius = grid.radius(i_r);
    const double theta  = grid.theta(i_theta);

    double coeff_beta, arr, att, art, detDF;
    level_cache.obtainValues(i_r, i_theta, center, radius, theta, coeff_beta, arr, att, art, detDF);

    const double h1 = grid.radialSpacing(i_r - 1);
    const double k1 = grid.angularSpacing(i_theta - 1);
    const double k2 = grid.angularSpacing(i_theta);

    const double coeff1 = 0.5 * (k1 + k2) / h1;

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
                ptr = getCircleAscIndex(i_r - 1, i_theta, grid, DirBC_Interior);

                const Stencil& LeftStencil = getStencil(i_r - 1, i_theta, grid, DirBC_Interior);

                offset = LeftStencil[StencilPosition::Center];
                column = left_index;
                value  = +coeff1 * arr; /* Center: (Right) */
                update_CSR_COO_MatrixElement(matrix, ptr, offset, row, column, value);
            }
        }
    }
}

} // namespace extrapolated_smoother_give

template <class LevelCacheType>
typename ExtrapolatedSmootherGive<LevelCacheType>::InnerBoundaryMatrix
ExtrapolatedSmootherGive<LevelCacheType>::buildInteriorBoundarySolverMatrix()
{
    using extrapolated_smoother_give::getNonZeroCountCircleAsc;
    using extrapolated_smoother_give::nodeBuildInteriorBoundarySolverMatrix_i_r_0;
    using extrapolated_smoother_give::nodeBuildInteriorBoundarySolverMatrix_i_r_1;

    const PolarGrid& grid             = ExtrapolatedSmoother<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = ExtrapolatedSmoother<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = ExtrapolatedSmoother<LevelCacheType>::DirBC_Interior_;

    const int ntheta = grid.ntheta();

    // The interior boundary matrix is symmetric due to the periodicity in the theta direction
    // and the assumption that ntheta is even, which is required for the across-origin discretization.
    // We store all non-zero entries of the matrix, both in COO format (for MUMPS)
    // and in CSR format (for the in-house solver). If the COO matrix is marked as symmetric,
    // the COO_Mumps_Solver optimizes the factorization by only using the upper triangular part of the matrix,
    // which is extracted by the COO_Mumps_Solver internally.
#ifdef GMGPOLAR_USE_MUMPS
    const int nnz = getNonZeroCountCircleAsc(0, grid, DirBC_Interior);
    SparseMatrixCOO<double, Kokkos::HostSpace> inner_boundary_solver_matrix(ntheta, ntheta, nnz);
    inner_boundary_solver_matrix.is_symmetric(true);
#else
    std::function<int(int)> nnz_per_row = [&](int i_theta) {
        if (DirBC_Interior)
            return 1;
        else
            return i_theta % 2 == 0 ? 1 : 2;
    };
    SparseMatrixCSR<double, Kokkos::HostSpace> inner_boundary_solver_matrix(ntheta, ntheta, nnz_per_row);
#endif

    {
        Kokkos::parallel_for(
            "ExtrapolatedSmootherGive: BuildInnerBoundaryMatrix",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, 1), KOKKOS_LAMBDA(const int) {
                for (int i_theta = 0; i_theta < ntheta; i_theta++) {
                    nodeBuildInteriorBoundarySolverMatrix_i_r_0(i_theta, grid, level_cache, DirBC_Interior,
                                                                inner_boundary_solver_matrix);
                }
            });
        Kokkos::fence();
    }
    {
        Kokkos::parallel_for(
            "ExtrapolatedSmootherGive: BuildInnerBoundaryMatrix",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, 1), KOKKOS_LAMBDA(const int) {
                for (int i_theta = 0; i_theta < ntheta; i_theta++) {
                    nodeBuildInteriorBoundarySolverMatrix_i_r_1(i_theta, grid, level_cache, DirBC_Interior,
                                                                inner_boundary_solver_matrix);
                }
            });
        Kokkos::fence();
    }

    return inner_boundary_solver_matrix;
}
