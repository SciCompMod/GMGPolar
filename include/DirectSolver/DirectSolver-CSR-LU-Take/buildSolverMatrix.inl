#pragma once

namespace direct_solver_csr_lu_take
{

static inline void updateMatrixElement(SparseMatrixCSR<double>& matrix, int offset, int row, int col, double val)
{
    matrix.row_nz_index(row, offset) = col;
    matrix.row_nz_entry(row, offset) = val;
}

} // namespace direct_solver_csr_lu_take

template <concepts::DomainGeometry DomainGeometry>
void DirectSolver_CSR_LU_Take<DomainGeometry>::nodeBuildSolverMatrixTake(
    int i_r, int i_theta, const PolarGrid& grid, bool DirBC_Interior, SparseMatrixCSR<double>& solver_matrix,
    ConstVector<double>& arr, ConstVector<double>& att, ConstVector<double>& art, ConstVector<double>& detDF,
    ConstVector<double>& coeff_beta)
{
    using direct_solver_csr_lu_take::updateMatrixElement;
    int offset;
    int row, col;
    double val;
    /* -------------------- */
    /* Node in the interior */
    /* -------------------- */
    if (i_r > 1 && i_r < grid.nr() - 2) {
        int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        double h1 = grid.radialSpacing(i_r - 1);
        double h2 = grid.radialSpacing(i_r);
        double k1 = grid.angularSpacing(i_theta_M1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;

        int center_index       = grid.index(i_r, i_theta);
        int left_index         = grid.index(i_r - 1, i_theta);
        int right_index        = grid.index(i_r + 1, i_theta);
        int bottom_index       = grid.index(i_r, i_theta_M1);
        int top_index          = grid.index(i_r, i_theta_P1);
        int bottom_left_index  = grid.index(i_r - 1, i_theta_M1);
        int bottom_right_index = grid.index(i_r + 1, i_theta_M1);
        int top_left_index     = grid.index(i_r - 1, i_theta_P1);
        int top_right_index    = grid.index(i_r + 1, i_theta_P1);

        double left_value   = -coeff1 * (arr[center_index] + arr[left_index]); /* Left */
        double right_value  = -coeff2 * (arr[center_index] + arr[right_index]); /* Right */
        double bottom_value = -coeff3 * (att[center_index] + att[bottom_index]); /* Bottom */
        double top_value    = -coeff4 * (att[center_index] + att[top_index]); /* Top */

        double center_value =
            (+0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center_index] * std::fabs(detDF[center_index]) /* beta_{i,j} */
             - left_value /* Center: (Left) */
             - right_value /* Center: (Right) */
             - bottom_value /* Center: (Bottom) */
             - top_value /* Center: (Top) */
            );

        double bottom_left_value  = -0.25 * (art[left_index] + art[bottom_index]); /* Bottom Left */
        double bottom_right_value = +0.25 * (art[right_index] + art[bottom_index]); /* Bottom Right */
        double top_left_value     = +0.25 * (art[left_index] + art[top_index]); /* Top Left */
        double top_right_value    = -0.25 * (art[right_index] + art[top_index]); /* Top Right */

        /* Fill matrix row of (i,j) */
        row = center_index;

        const Stencil& CenterStencil = getStencil(i_r);

        offset = CenterStencil[StencilPosition::Center];
        col    = center_index;
        val    = center_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Left];
        col    = left_index;
        val    = left_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Right];
        col    = right_index;
        val    = right_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Bottom];
        col    = bottom_index;
        val    = bottom_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Top];
        col    = top_index;
        val    = top_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);

        offset = CenterStencil[StencilPosition::BottomLeft];
        col    = bottom_left_index;
        val    = bottom_left_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);

        offset = CenterStencil[StencilPosition::BottomRight];
        col    = bottom_right_index;
        val    = bottom_right_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);

        offset = CenterStencil[StencilPosition::TopLeft];
        col    = top_left_index;
        val    = top_left_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);

        offset = CenterStencil[StencilPosition::TopRight];
        col    = top_right_index;
        val    = top_right_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);
    }
    /* -------------------------- */
    /* Node on the inner boundary */
    /* -------------------------- */
    else if (i_r == 0) {
        /* ------------------------------------------------ */
        /* Case 1: Dirichlet boundary on the inner boundary */
        /* ------------------------------------------------ */
        if (DirBC_Interior) {
            int center_index = grid.index(i_r, i_theta);

            /* Fill matrix row of (i,j) */
            row = center_index;

            const Stencil& CenterStencil = getStencil(i_r);

            offset = CenterStencil[StencilPosition::Center];
            col    = center_index;
            val    = 1.0;
            updateMatrixElement(solver_matrix, offset, row, col, val);
        }
        else {
            /* ------------------------------------------------------------- */
            /* Case 2: Across origin discretization on the interior boundary */
            /* ------------------------------------------------------------- */
            // h1 gets replaced with 2 * R0.
            // (i_r-1,i_theta) gets replaced with (i_r, i_theta + grid.ntheta()/2).
            // Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil.
            assert(grid.ntheta() % 2 == 0);

            int i_theta_M1           = grid.wrapThetaIndex(i_theta - 1);
            int i_theta_P1           = grid.wrapThetaIndex(i_theta + 1);
            int i_theta_AcrossOrigin = grid.wrapThetaIndex(i_theta + grid.ntheta() / 2);

            double h1 = 2.0 * grid.radius(0);
            double h2 = grid.radialSpacing(i_r);
            double k1 = grid.angularSpacing(i_theta_M1);
            double k2 = grid.angularSpacing(i_theta);

            double coeff1 = 0.5 * (k1 + k2) / h1;
            double coeff2 = 0.5 * (k1 + k2) / h2;
            double coeff3 = 0.5 * (h1 + h2) / k1;
            double coeff4 = 0.5 * (h1 + h2) / k2;

            int center_index       = grid.index(i_r, i_theta);
            int left_index         = grid.index(i_r, i_theta_AcrossOrigin);
            int right_index        = grid.index(i_r + 1, i_theta);
            int bottom_index       = grid.index(i_r, i_theta_M1);
            int top_index          = grid.index(i_r, i_theta_P1);
            int bottom_right_index = grid.index(i_r + 1, i_theta_M1);
            int top_right_index    = grid.index(i_r + 1, i_theta_P1);

            double left_value   = -coeff1 * (arr[center_index] + arr[left_index]); /* Left */
            double right_value  = -coeff2 * (arr[center_index] + arr[right_index]); /* Right */
            double bottom_value = -coeff3 * (att[center_index] + att[bottom_index]); /* Bottom */
            double top_value    = -coeff4 * (att[center_index] + att[top_index]); /* Top */

            double center_value =
                (+0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center_index] * fabs(detDF[center_index]) /* beta_{i,j} */
                 - left_value /* Center: (Left) */
                 - right_value /* Center: (Right) */
                 - bottom_value /* Center: (Bottom) */
                 - top_value /* Center: (Top) */
                );

            double bottom_right_value = +0.25 * (art[right_index] + art[bottom_index]); /* Bottom Right */
            double top_right_value    = -0.25 * (art[right_index] + art[top_index]); /* Top Right */

            /* Fill matrix row of (i,j) */
            row = center_index;

            const Stencil& CenterStencil = getStencil(i_r);

            offset = CenterStencil[StencilPosition::Center];
            col    = center_index;
            val    = center_value;
            updateMatrixElement(solver_matrix, offset, row, col, val);

            offset = CenterStencil[StencilPosition::Left];
            col    = left_index;
            val    = left_value;
            updateMatrixElement(solver_matrix, offset, row, col, val);

            offset = CenterStencil[StencilPosition::Right];
            col    = right_index;
            val    = right_value;
            updateMatrixElement(solver_matrix, offset, row, col, val);

            offset = CenterStencil[StencilPosition::Bottom];
            col    = bottom_index;
            val    = bottom_value;
            updateMatrixElement(solver_matrix, offset, row, col, val);

            offset = CenterStencil[StencilPosition::Top];
            col    = top_index;
            val    = top_value;
            updateMatrixElement(solver_matrix, offset, row, col, val);

            /* BottomLeft: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */

            offset = CenterStencil[StencilPosition::BottomRight];
            col    = bottom_right_index;
            val    = bottom_right_value;
            updateMatrixElement(solver_matrix, offset, row, col, val);

            /* TopLeft: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */

            offset = CenterStencil[StencilPosition::TopRight];
            col    = top_right_index;
            val    = top_right_value;
            updateMatrixElement(solver_matrix, offset, row, col, val);
        }
    }
    /* ------------------------------- */
    /* Node next to the inner boundary */
    /* ------------------------------- */
    else if (i_r == 1) {
        int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        double h1 = grid.radialSpacing(i_r - 1);
        double h2 = grid.radialSpacing(i_r);
        double k1 = grid.angularSpacing(i_theta_M1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;

        int center_index       = grid.index(i_r, i_theta);
        int left_index         = grid.index(i_r - 1, i_theta);
        int right_index        = grid.index(i_r + 1, i_theta);
        int bottom_index       = grid.index(i_r, i_theta_M1);
        int top_index          = grid.index(i_r, i_theta_P1);
        int bottom_left_index  = grid.index(i_r - 1, i_theta_M1);
        int bottom_right_index = grid.index(i_r + 1, i_theta_M1);
        int top_left_index     = grid.index(i_r - 1, i_theta_P1);
        int top_right_index    = grid.index(i_r + 1, i_theta_P1);

        double left_value   = -coeff1 * (arr[center_index] + arr[left_index]); /* Left */
        double right_value  = -coeff2 * (arr[center_index] + arr[right_index]); /* Right */
        double bottom_value = -coeff3 * (att[center_index] + att[bottom_index]); /* Bottom */
        double top_value    = -coeff4 * (att[center_index] + att[top_index]); /* Top */

        double center_value =
            (+0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center_index] * std::fabs(detDF[center_index]) /* beta_{i,j} */
             - left_value /* Center: (Left) */
             - right_value /* Center: (Right) */
             - bottom_value /* Center: (Bottom) */
             - top_value /* Center: (Top) */
            );

        double bottom_left_value  = -0.25 * (art[left_index] + art[bottom_index]); /* Bottom Left */
        double bottom_right_value = +0.25 * (art[right_index] + art[bottom_index]); /* Bottom Right */
        double top_left_value     = +0.25 * (art[left_index] + art[top_index]); /* Top Left */
        double top_right_value    = -0.25 * (art[right_index] + art[top_index]); /* Top Right */

        /* Fill matrix row of (i,j) */
        row = center_index;

        const Stencil& CenterStencil = getStencil(i_r);

        offset = CenterStencil[StencilPosition::Center];
        col    = center_index;
        val    = center_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);

        /* REMOVED: Moved to the right hand side to make the matrix symmetric */
        if (!DirBC_Interior) {
            offset = CenterStencil[StencilPosition::Left];
            col    = left_index;
            val    = left_value;
            updateMatrixElement(solver_matrix, offset, row, col, val);
        }

        offset = CenterStencil[StencilPosition::Right];
        col    = right_index;
        val    = right_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Bottom];
        col    = bottom_index;
        val    = bottom_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Top];
        col    = top_index;
        val    = top_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);

        /* REMOVED: Moved to the right hand side to make the matrix symmetric */
        if (!DirBC_Interior) {
            offset = CenterStencil[StencilPosition::BottomLeft];
            col    = bottom_left_index;
            val    = bottom_left_value;
            updateMatrixElement(solver_matrix, offset, row, col, val);
        }

        offset = CenterStencil[StencilPosition::BottomRight];
        col    = bottom_right_index;
        val    = bottom_right_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);

        /* REMOVED: Moved to the right hand side to make the matrix symmetric */
        if (!DirBC_Interior) {
            offset = CenterStencil[StencilPosition::TopLeft];
            col    = top_left_index;
            val    = top_left_value;
            updateMatrixElement(solver_matrix, offset, row, col, val);
        }

        offset = CenterStencil[StencilPosition::TopRight];
        col    = top_right_index;
        val    = top_right_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);
    }
    /* ------------------------------- */
    /* Node next to the outer boundary */
    /* ------------------------------- */
    else if (i_r == grid.nr() - 2) {
        int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        double h1 = grid.radialSpacing(i_r - 1);
        double h2 = grid.radialSpacing(i_r);
        double k1 = grid.angularSpacing(i_theta_M1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;

        int center_index       = grid.index(i_r, i_theta);
        int left_index         = grid.index(i_r - 1, i_theta);
        int right_index        = grid.index(i_r + 1, i_theta);
        int bottom_index       = grid.index(i_r, i_theta_M1);
        int top_index          = grid.index(i_r, i_theta_P1);
        int bottom_left_index  = grid.index(i_r - 1, i_theta_M1);
        int bottom_right_index = grid.index(i_r + 1, i_theta_M1);
        int top_left_index     = grid.index(i_r - 1, i_theta_P1);
        int top_right_index    = grid.index(i_r + 1, i_theta_P1);

        double left_value   = -coeff1 * (arr[center_index] + arr[left_index]); /* Left */
        double right_value  = -coeff2 * (arr[center_index] + arr[right_index]); /* Right */
        double bottom_value = -coeff3 * (att[center_index] + att[bottom_index]); /* Bottom */
        double top_value    = -coeff4 * (att[center_index] + att[top_index]); /* Top */

        double center_value =
            (+0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center_index] * std::fabs(detDF[center_index]) /* beta_{i,j} */
             - left_value /* Center: (Left) */
             - right_value /* Center: (Right) */
             - bottom_value /* Center: (Bottom) */
             - top_value /* Center: (Top) */
            );

        double bottom_left_value  = -0.25 * (art[left_index] + art[bottom_index]); /* Bottom Left */
        double bottom_right_value = +0.25 * (art[right_index] + art[bottom_index]); /* Bottom Right */
        double top_left_value     = +0.25 * (art[left_index] + art[top_index]); /* Top Left */
        double top_right_value    = -0.25 * (art[right_index] + art[top_index]); /* Top Right */

        /* Fill matrix row of (i,j) */
        row = center_index;

        const Stencil& CenterStencil = getStencil(i_r);

        offset = CenterStencil[StencilPosition::Center];
        col    = center_index;
        val    = center_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Left];
        col    = left_index;
        val    = left_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);

        /* Right REMOVED: Moved to the right hand side to make the matrix symmetric */

        offset = CenterStencil[StencilPosition::Bottom];
        col    = bottom_index;
        val    = bottom_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Top];
        col    = top_index;
        val    = top_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);

        offset = CenterStencil[StencilPosition::BottomLeft];
        col    = bottom_left_index;
        val    = bottom_left_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);

        /* BottomRight REMOVED: Moved to the right hand side to make the matrix symmetric */

        offset = CenterStencil[StencilPosition::TopLeft];
        col    = top_left_index;
        val    = top_left_value;
        updateMatrixElement(solver_matrix, offset, row, col, val);

        /* TopRight REMOVED: Moved to the right hand side to make the matrix symmetric */
    }
    /* ------------------------------------ */
    /* Node on the outer dirichlet boundary */
    /* ------------------------------------ */
    else if (i_r == grid.nr() - 1) {
        int center_index = grid.index(i_r, i_theta);

        /* Fill matrix row of (i,j) */
        row = center_index;

        const Stencil& CenterStencil = getStencil(i_r);

        offset = CenterStencil[StencilPosition::Center];
        col    = center_index;
        val    = 1.0;
        updateMatrixElement(solver_matrix, offset, row, col, val);
    }
}

template <concepts::DomainGeometry DomainGeometry>
SparseMatrixCSR<double> DirectSolver_CSR_LU_Take<DomainGeometry>::buildSolverMatrix()
{
    const PolarGrid& grid                         = DirectSolver<DomainGeometry>::grid_;
    const LevelCache<DomainGeometry>& level_cache = DirectSolver<DomainGeometry>::level_cache_;
    const int num_omp_threads                     = DirectSolver<DomainGeometry>::num_omp_threads_;
    const bool DirBC_Interior                     = DirectSolver<DomainGeometry>::DirBC_Interior_;

    const int n = grid.numberOfNodes();

    std::function<int(int)> nnz_per_row = [&](int global_index) {
        return getStencilSize(global_index);
    };

    SparseMatrixCSR<double> solver_matrix(n, n, nnz_per_row);

    assert(level_cache.cacheDensityProfileCoefficients());
    assert(level_cache.cacheDomainGeometry());

    ConstVector<double> arr        = level_cache.arr();
    ConstVector<double> att        = level_cache.att();
    ConstVector<double> art        = level_cache.art();
    ConstVector<double> detDF      = level_cache.detDF();
    ConstVector<double> coeff_beta = level_cache.coeff_beta();

#pragma omp parallel num_threads(num_omp_threads)
    {
        /* Circle Section */
#pragma omp for nowait
        for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++) {
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                nodeBuildSolverMatrixTake(i_r, i_theta, grid, DirBC_Interior, solver_matrix, arr, att, art, detDF,
                                          coeff_beta);
            }
        }
        /* Radial Section */
#pragma omp for nowait
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                nodeBuildSolverMatrixTake(i_r, i_theta, grid, DirBC_Interior, solver_matrix, arr, att, art, detDF,
                                          coeff_beta);
            }
        }
    }

    return solver_matrix;
}
