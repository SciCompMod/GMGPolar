#pragma once

#include "matrixStencil.inl"

namespace direct_solver_take
{

#ifdef GMGPOLAR_USE_MUMPS
// When using the MUMPS solver, the matrix is assembled in COO format.
static KOKKOS_INLINE_FUNCTION void updateMatrixElement(const SparseMatrixCOO<double, Kokkos::HostSpace>& matrix,
                                                       const int ptr, const int offset, const int row, const int column,
                                                       const double value)
{
    matrix.set_row_index(ptr + offset, row);
    matrix.set_col_index(ptr + offset, column);
    matrix.set_value(ptr + offset, value);
}
#else
// When using the in-house solver, the matrix is stored in CSR format.
static KOKKOS_INLINE_FUNCTION void updateMatrixElement(const SparseMatrixCSR<double, Kokkos::HostSpace>& matrix,
                                                       const int ptr, const int offset, const int row, const int column,
                                                       const double value)
{
    matrix.set_row_nz_index(row, offset, column);
    matrix.set_row_nz_entry(row, offset, value);
}
#endif

template <typename SystemMatrix>
static KOKKOS_INLINE_FUNCTION void
nodeBuildSolverMatrixTake(const int i_r, const int i_theta, const PolarGrid& grid, const bool DirBC_Interior,
                          const SystemMatrix& solver_matrix, HostConstVector<double>& arr, HostConstVector<double>& att,
                          HostConstVector<double>& art, HostConstVector<double>& detDF, const double coeff_beta)
{
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

        const int center_index       = grid.index(i_r, i_theta);
        const int left_index         = grid.index(i_r - 1, i_theta);
        const int right_index        = grid.index(i_r + 1, i_theta);
        const int bottom_index       = grid.index(i_r, i_theta_M1);
        const int top_index          = grid.index(i_r, i_theta_P1);
        const int bottom_left_index  = grid.index(i_r - 1, i_theta_M1);
        const int bottom_right_index = grid.index(i_r + 1, i_theta_M1);
        const int top_left_index     = grid.index(i_r - 1, i_theta_P1);
        const int top_right_index    = grid.index(i_r + 1, i_theta_P1);

        const double left_value   = -coeff1 * (arr(center_index) + arr(left_index)); /* Left */
        const double right_value  = -coeff2 * (arr(center_index) + arr(right_index)); /* Right */
        const double bottom_value = -coeff3 * (att(center_index) + att(bottom_index)); /* Bottom */
        const double top_value    = -coeff4 * (att(center_index) + att(top_index)); /* Top */

        const double center_value = (coeff5 * coeff_beta * Kokkos::fabs(detDF(center_index)) /* beta_{i,j} */
                                     - left_value /* Center: (Left) */
                                     - right_value /* Center: (Right) */
                                     - bottom_value /* Center: (Bottom) */
                                     - top_value /* Center: (Top) */
        );

        const double bottom_left_value  = -0.25 * (art(left_index) + art(bottom_index)); /* Bottom Left */
        const double bottom_right_value = +0.25 * (art(right_index) + art(bottom_index)); /* Bottom Right */
        const double top_left_value     = +0.25 * (art(left_index) + art(top_index)); /* Top Left */
        const double top_right_value    = -0.25 * (art(right_index) + art(top_index)); /* Top Right */

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r, grid, DirBC_Interior);

        offset = CenterStencil[StencilPosition::Center];
        column = center_index;
        value  = center_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Left];
        column = left_index;
        value  = left_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Right];
        column = right_index;
        value  = right_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Bottom];
        column = bottom_index;
        value  = bottom_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Top];
        column = top_index;
        value  = top_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::BottomLeft];
        column = bottom_left_index;
        value  = bottom_left_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::BottomRight];
        column = bottom_right_index;
        value  = bottom_right_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::TopLeft];
        column = top_left_index;
        value  = top_left_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::TopRight];
        column = top_right_index;
        value  = top_right_value;
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
            const int center_nz_index = getSolverMatrixIndex(i_r, i_theta, grid, DirBC_Interior);

            const int center_index = grid.index(i_r, i_theta);

            /* Fill matrix row of (i,j) */
            row = center_index;
            ptr = center_nz_index;

            const Stencil& CenterStencil = getStencil(i_r, grid, DirBC_Interior);

            offset = CenterStencil[StencilPosition::Center];
            column = center_index;
            value  = 1.0;
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);
        }
        else {
            /* ------------------------------------------------------------- */
            /* Case 2: Across origin discretization on the interior boundary */
            /* ------------------------------------------------------------- */
            // h1 gets replaced with 2 * R0.
            // (i_r-1,i_theta) gets replaced with (i_r, i_theta + grid.ntheta()/2).
            // Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil.
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

            KOKKOS_ASSERT(grid.ntheta() % 2 == 0);
            const int i_theta_AcrossOrigin = grid.wrapThetaIndex(i_theta + grid.ntheta() / 2);

            const double h1 = 2.0 * grid.radius(0);
            const double h2 = grid.radialSpacing(i_r);
            const double k1 = grid.angularSpacing(i_theta_M1);
            const double k2 = grid.angularSpacing(i_theta);

            const double coeff1 = 0.5 * (k1 + k2) / h1;
            const double coeff2 = 0.5 * (k1 + k2) / h2;
            const double coeff3 = 0.5 * (h1 + h2) / k1;
            const double coeff4 = 0.5 * (h1 + h2) / k2;
            const double coeff5 = 0.25 * (h1 + h2) * (k1 + k2);

            const int center_nz_index = getSolverMatrixIndex(i_r, i_theta, grid, DirBC_Interior);

            const int center_index       = grid.index(i_r, i_theta);
            const int left_index         = grid.index(i_r, i_theta_AcrossOrigin);
            const int right_index        = grid.index(i_r + 1, i_theta);
            const int bottom_index       = grid.index(i_r, i_theta_M1);
            const int top_index          = grid.index(i_r, i_theta_P1);
            const int bottom_right_index = grid.index(i_r + 1, i_theta_M1);
            const int top_right_index    = grid.index(i_r + 1, i_theta_P1);

            const double left_value   = -coeff1 * (arr(center_index) + arr(left_index)); /* Left */
            const double right_value  = -coeff2 * (arr(center_index) + arr(right_index)); /* Right */
            const double bottom_value = -coeff3 * (att(center_index) + att(bottom_index)); /* Bottom */
            const double top_value    = -coeff4 * (att(center_index) + att(top_index)); /* Top */

            const double center_value = (coeff5 * coeff_beta * Kokkos::fabs(detDF(center_index)) /* beta_{i,j} */
                                         - left_value /* Center: (Left) */
                                         - right_value /* Center: (Right) */
                                         - bottom_value /* Center: (Bottom) */
                                         - top_value /* Center: (Top) */
            );

            const double bottom_right_value = +0.25 * (art(right_index) + art(bottom_index)); /* Bottom Right */
            const double top_right_value    = -0.25 * (art(right_index) + art(top_index)); /* Top Right */

            /* Fill matrix row of (i,j) */
            row = center_index;
            ptr = center_nz_index;

            const Stencil& CenterStencil = getStencil(i_r, grid, DirBC_Interior);

            offset = CenterStencil[StencilPosition::Center];
            column = center_index;
            value  = center_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = CenterStencil[StencilPosition::Left];
            column = left_index;
            value  = left_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = CenterStencil[StencilPosition::Right];
            column = right_index;
            value  = right_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = CenterStencil[StencilPosition::Bottom];
            column = bottom_index;
            value  = bottom_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            offset = CenterStencil[StencilPosition::Top];
            column = top_index;
            value  = top_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            /* BottomLeft: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */

            offset = CenterStencil[StencilPosition::BottomRight];
            column = bottom_right_index;
            value  = bottom_right_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

            /* TopLeft: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */

            offset = CenterStencil[StencilPosition::TopRight];
            column = top_right_index;
            value  = top_right_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);
        }
    }
    /* ------------------------------- */
    /* Node next to the inner boundary */
    /* ------------------------------- */
    else if (i_r == 1) {
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

        const int center_index       = grid.index(i_r, i_theta);
        const int left_index         = grid.index(i_r - 1, i_theta);
        const int right_index        = grid.index(i_r + 1, i_theta);
        const int bottom_index       = grid.index(i_r, i_theta_M1);
        const int top_index          = grid.index(i_r, i_theta_P1);
        const int bottom_left_index  = grid.index(i_r - 1, i_theta_M1);
        const int bottom_right_index = grid.index(i_r + 1, i_theta_M1);
        const int top_left_index     = grid.index(i_r - 1, i_theta_P1);
        const int top_right_index    = grid.index(i_r + 1, i_theta_P1);

        const double left_value   = -coeff1 * (arr(center_index) + arr(left_index)); /* Left */
        const double right_value  = -coeff2 * (arr(center_index) + arr(right_index)); /* Right */
        const double bottom_value = -coeff3 * (att(center_index) + att(bottom_index)); /* Bottom */
        const double top_value    = -coeff4 * (att(center_index) + att(top_index)); /* Top */

        const double center_value = (coeff5 * coeff_beta * Kokkos::fabs(detDF(center_index)) /* beta_{i,j} */
                                     - left_value /* Center: (Left) */
                                     - right_value /* Center: (Right) */
                                     - bottom_value /* Center: (Bottom) */
                                     - top_value /* Center: (Top) */
        );

        const double bottom_left_value  = -0.25 * (art(left_index) + art(bottom_index)); /* Bottom Left */
        const double bottom_right_value = +0.25 * (art(right_index) + art(bottom_index)); /* Bottom Right */
        const double top_left_value     = +0.25 * (art(left_index) + art(top_index)); /* Top Left */
        const double top_right_value    = -0.25 * (art(right_index) + art(top_index)); /* Top Right */

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r, grid, DirBC_Interior);

        offset = CenterStencil[StencilPosition::Center];
        column = center_index;
        value  = center_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* REMOVED: Moved to the right hand side to make the matrix symmetric */
        if (!DirBC_Interior) {
            offset = CenterStencil[StencilPosition::Left];
            column = left_index;
            value  = left_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);
        }

        offset = CenterStencil[StencilPosition::Right];
        column = right_index;
        value  = right_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Bottom];
        column = bottom_index;
        value  = bottom_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Top];
        column = top_index;
        value  = top_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* REMOVED: Moved to the right hand side to make the matrix symmetric */
        if (!DirBC_Interior) {
            offset = CenterStencil[StencilPosition::BottomLeft];
            column = bottom_left_index;
            value  = bottom_left_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);
        }

        offset = CenterStencil[StencilPosition::BottomRight];
        column = bottom_right_index;
        value  = bottom_right_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* REMOVED: Moved to the right hand side to make the matrix symmetric */
        if (!DirBC_Interior) {
            offset = CenterStencil[StencilPosition::TopLeft];
            column = top_left_index;
            value  = top_left_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, column, value);
        }

        offset = CenterStencil[StencilPosition::TopRight];
        column = top_right_index;
        value  = top_right_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);
    }
    /* ------------------------------- */
    /* Node next to the outer boundary */
    /* ------------------------------- */
    else if (i_r == grid.nr() - 2) {
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

        const int center_index       = grid.index(i_r, i_theta);
        const int left_index         = grid.index(i_r - 1, i_theta);
        const int right_index        = grid.index(i_r + 1, i_theta);
        const int bottom_index       = grid.index(i_r, i_theta_M1);
        const int top_index          = grid.index(i_r, i_theta_P1);
        const int bottom_left_index  = grid.index(i_r - 1, i_theta_M1);
        const int bottom_right_index = grid.index(i_r + 1, i_theta_M1);
        const int top_left_index     = grid.index(i_r - 1, i_theta_P1);
        const int top_right_index    = grid.index(i_r + 1, i_theta_P1);

        const double left_value   = -coeff1 * (arr(center_index) + arr(left_index)); /* Left */
        const double right_value  = -coeff2 * (arr(center_index) + arr(right_index)); /* Right */
        const double bottom_value = -coeff3 * (att(center_index) + att(bottom_index)); /* Bottom */
        const double top_value    = -coeff4 * (att(center_index) + att(top_index)); /* Top */

        const double center_value = (coeff5 * coeff_beta * Kokkos::fabs(detDF(center_index)) /* beta_{i,j} */
                                     - left_value /* Center: (Left) */
                                     - right_value /* Center: (Right) */
                                     - bottom_value /* Center: (Bottom) */
                                     - top_value /* Center: (Top) */
        );

        const double bottom_left_value  = -0.25 * (art(left_index) + art(bottom_index)); /* Bottom Left */
        const double bottom_right_value = +0.25 * (art(right_index) + art(bottom_index)); /* Bottom Right */
        const double top_left_value     = +0.25 * (art(left_index) + art(top_index)); /* Top Left */
        const double top_right_value    = -0.25 * (art(right_index) + art(top_index)); /* Top Right */

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r, grid, DirBC_Interior);

        offset = CenterStencil[StencilPosition::Center];
        column = center_index;
        value  = center_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Left];
        column = left_index;
        value  = left_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* Right REMOVED: Moved to the right hand side to make the matrix symmetric */

        offset = CenterStencil[StencilPosition::Bottom];
        column = bottom_index;
        value  = bottom_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::Top];
        column = top_index;
        value  = top_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        offset = CenterStencil[StencilPosition::BottomLeft];
        column = bottom_left_index;
        value  = bottom_left_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* BottomRight REMOVED: Moved to the right hand side to make the matrix symmetric */

        offset = CenterStencil[StencilPosition::TopLeft];
        column = top_left_index;
        value  = top_left_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);

        /* TopRight REMOVED: Moved to the right hand side to make the matrix symmetric */
    }
    /* ------------------------------------ */
    /* Node on the outer dirichlet boundary */
    /* ------------------------------------ */
    else if (i_r == grid.nr() - 1) {
        const int center_nz_index = getSolverMatrixIndex(i_r, i_theta, grid, DirBC_Interior);

        const int center_index = grid.index(i_r, i_theta);

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r, grid, DirBC_Interior);

        offset = CenterStencil[StencilPosition::Center];
        column = center_index;
        value  = 1.0;
        updateMatrixElement(solver_matrix, ptr, offset, row, column, value);
    }
}

} // namespace direct_solver_take

template <class LevelCacheType>
typename DirectSolverTake<LevelCacheType>::SystemMatrix DirectSolverTake<LevelCacheType>::buildSolverMatrix()
{
    using direct_solver_take::getNonZeroCountSolverMatrix;
    using direct_solver_take::getStencilSize;
    using direct_solver_take::nodeBuildSolverMatrixTake;
    using direct_solver_take::validateSolverMatrixIndexing;

    const PolarGrid& grid             = DirectSolver<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = DirectSolver<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = DirectSolver<LevelCacheType>::DirBC_Interior_;

    assert(validateSolverMatrixIndexing(grid, DirBC_Interior) && "Solver matrix indexing is inconsistent");

    const int n = grid.numberOfNodes();

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

    assert(level_cache.cacheDomainGeometry());
    HostConstVector<double> arr   = level_cache.arr();
    HostConstVector<double> att   = level_cache.att();
    HostConstVector<double> art   = level_cache.art();
    HostConstVector<double> detDF = level_cache.detDF();

    /* We split the loops into two regions to better respect the */
    /* access patterns of the smoother and improve cache locality. */

    // The For loop matches circular access pattern */
    Kokkos::parallel_for(
        "Residual Take: Apply System Operator (Circular)",
        Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
            {0, 0}, // Starting point of the index space
            {grid.numberSmootherCircles(), grid.ntheta()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_r, const int i_theta) {
            const double coeff_beta = level_cache.obtainBeta(i_r, i_theta, grid);
            nodeBuildSolverMatrixTake(i_r, i_theta, grid, DirBC_Interior, solver_matrix, arr, att, art, detDF,
                                      coeff_beta);
        });

    /* For loop matches radial access pattern */
    Kokkos::parallel_for(
        "Residual Take: Apply System Operator (Radial)",
        Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
            {0, grid.numberSmootherCircles()}, // Starting point of the index space
            {grid.ntheta(), grid.nr()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_theta, const int i_r) {
            const double coeff_beta = level_cache.obtainBeta(i_r, i_theta, grid);
            nodeBuildSolverMatrixTake(i_r, i_theta, grid, DirBC_Interior, solver_matrix, arr, att, art, detDF,
                                      coeff_beta);
        });

    Kokkos::fence();

    return solver_matrix;
}
