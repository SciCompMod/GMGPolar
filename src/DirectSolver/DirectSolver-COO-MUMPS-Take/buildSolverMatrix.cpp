#include "../../../include/DirectSolver/DirectSolver-COO-MUMPS-Take/directSolverTake.h"

#ifdef GMGPOLAR_USE_MUMPS

inline void updateMatrixElement(SparseMatrixCOO<double>& matrix, int ptr, int offset, int row, int col, double val)
{
    matrix.row_index(ptr + offset) = row;
    matrix.col_index(ptr + offset) = col;
    matrix.value(ptr + offset)     = val;
}

void DirectSolver_COO_MUMPS_Take::nodeBuildSolverMatrixTake(int i_r, int i_theta, const PolarGrid& grid,
                                                            bool DirBC_Interior, SparseMatrixCOO<double>& solver_matrix,
                                                            ConstVector<double>& arr, ConstVector<double>& att,
                                                            ConstVector<double>& art, ConstVector<double>& detDF,
                                                            ConstVector<double>& coeff_beta)
{
    int ptr, offset;
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

        int center_nz_index = getSolverMatrixIndex(i_r, i_theta);

        int center_index       = grid.index(i_r, i_theta);
        int left_index         = grid.index(i_r - 1, i_theta);
        int right_index        = grid.index(i_r + 1, i_theta);
        int bottom_index       = grid.index(i_r, i_theta_M1);
        int top_index          = grid.index(i_r, i_theta_P1);
        int bottom_left_index  = grid.index(i_r - 1, i_theta_M1);
        int bottom_right_index = grid.index(i_r + 1, i_theta_M1);
        int top_left_index     = grid.index(i_r - 1, i_theta_P1);
        int top_right_index    = grid.index(i_r + 1, i_theta_P1);

        double left_value   = -coeff1 * (arr(center_index) + arr(left_index)); /* Left */
        double right_value  = -coeff2 * (arr(center_index) + arr(right_index)); /* Right */
        double bottom_value = -coeff3 * (att(center_index) + att(bottom_index)); /* Bottom */
        double top_value    = -coeff4 * (att(center_index) + att(top_index)); /* Top */

        double center_value =
            (+0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center_index] * fabs(detDF(center_index)) /* beta_{i,j} */
             - left_value /* Center: (Left) */
             - right_value /* Center: (Right) */
             - bottom_value /* Center: (Bottom) */
             - top_value /* Center: (Top) */
            );

        double bottom_left_value  = -0.25 * (art(left_index) + art(bottom_index)); /* Bottom Left */
        double bottom_right_value = +0.25 * (art(right_index) + art(bottom_index)); /* Bottom Right */
        double top_left_value     = +0.25 * (art(left_index) + art(top_index)); /* Top Left */
        double top_right_value    = -0.25 * (art(right_index) + art(top_index)); /* Top Right */

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r);

        offset = CenterStencil[StencilPosition::Center];
        col    = center_index;
        val    = center_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Left];
        col    = left_index;
        val    = left_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Right];
        col    = right_index;
        val    = right_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Bottom];
        col    = bottom_index;
        val    = bottom_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Top];
        col    = top_index;
        val    = top_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::BottomLeft];
        col    = bottom_left_index;
        val    = bottom_left_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::BottomRight];
        col    = bottom_right_index;
        val    = bottom_right_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::TopLeft];
        col    = top_left_index;
        val    = top_left_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::TopRight];
        col    = top_right_index;
        val    = top_right_value;
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
            int center_nz_index = getSolverMatrixIndex(i_r, i_theta);

            int center_index = grid.index(i_r, i_theta);

            /* Fill matrix row of (i,j) */
            row = center_index;
            ptr = center_nz_index;

            const Stencil& CenterStencil = getStencil(i_r);

            offset = CenterStencil[StencilPosition::Center];
            col    = center_index;
            val    = 1.0;
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);
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

            int center_nz_index = getSolverMatrixIndex(i_r, i_theta);

            int center_index       = grid.index(i_r, i_theta);
            int left_index         = grid.index(i_r, i_theta_AcrossOrigin);
            int right_index        = grid.index(i_r + 1, i_theta);
            int bottom_index       = grid.index(i_r, i_theta_M1);
            int top_index          = grid.index(i_r, i_theta_P1);
            int bottom_right_index = grid.index(i_r + 1, i_theta_M1);
            int top_right_index    = grid.index(i_r + 1, i_theta_P1);

            double left_value   = -coeff1 * (arr(center_index) + arr(left_index)); /* Left */
            double right_value  = -coeff2 * (arr(center_index) + arr(right_index)); /* Right */
            double bottom_value = -coeff3 * (att(center_index) + att(bottom_index)); /* Bottom */
            double top_value    = -coeff4 * (att(center_index) + att(top_index)); /* Top */

            double center_value =
                (+0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center_index] * fabs(detDF(center_index)) /* beta_{i,j} */
                 - left_value /* Center: (Left) */
                 - right_value /* Center: (Right) */
                 - bottom_value /* Center: (Bottom) */
                 - top_value /* Center: (Top) */
                );

            double bottom_right_value = +0.25 * (art(right_index) + art(bottom_index)); /* Bottom Right */
            double top_right_value    = -0.25 * (art(right_index) + art(top_index)); /* Top Right */

            /* Fill matrix row of (i,j) */
            row = center_index;
            ptr = center_nz_index;

            const Stencil& CenterStencil = getStencil(i_r);

            offset = CenterStencil[StencilPosition::Center];
            col    = center_index;
            val    = center_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = CenterStencil[StencilPosition::Left];
            col    = left_index;
            val    = left_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = CenterStencil[StencilPosition::Right];
            col    = right_index;
            val    = right_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = CenterStencil[StencilPosition::Bottom];
            col    = bottom_index;
            val    = bottom_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            offset = CenterStencil[StencilPosition::Top];
            col    = top_index;
            val    = top_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            /* BottomLeft: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */

            offset = CenterStencil[StencilPosition::BottomRight];
            col    = bottom_right_index;
            val    = bottom_right_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

            /* TopLeft: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */

            offset = CenterStencil[StencilPosition::TopRight];
            col    = top_right_index;
            val    = top_right_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);
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

        int center_nz_index = getSolverMatrixIndex(i_r, i_theta);

        int center_index             = grid.index(i_r, i_theta);
        int left_index               = grid.index(i_r - 1, i_theta);
        int right_index              = grid.index(i_r + 1, i_theta);
        int bottom_index             = grid.index(i_r, i_theta_M1);
        int top_index                = grid.index(i_r, i_theta_P1);
        int bottom_left_index        = grid.index(i_r - 1, i_theta_M1);
        const int bottom_right_index = grid.index(i_r + 1, i_theta_M1);
        int top_left_index           = grid.index(i_r - 1, i_theta_P1);
        int top_right_index          = grid.index(i_r + 1, i_theta_P1);

        double left_value   = -coeff1 * (arr(center_index) + arr(left_index)); /* Left */
        double right_value  = -coeff2 * (arr(center_index) + arr(right_index)); /* Right */
        double bottom_value = -coeff3 * (att(center_index) + att(bottom_index)); /* Bottom */
        double top_value    = -coeff4 * (att(center_index) + att(top_index)); /* Top */

        double center_value =
            (+0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center_index] * fabs(detDF(center_index)) /* beta_{i,j} */
             - left_value /* Center: (Left) */
             - right_value /* Center: (Right) */
             - bottom_value /* Center: (Bottom) */
             - top_value /* Center: (Top) */
            );

        double bottom_left_value  = -0.25 * (art(left_index) + art(bottom_index)); /* Bottom Left */
        double bottom_right_value = +0.25 * (art(right_index) + art(bottom_index)); /* Bottom Right */
        double top_left_value     = +0.25 * (art(left_index) + art(top_index)); /* Top Left */
        double top_right_value    = -0.25 * (art(right_index) + art(top_index)); /* Top Right */

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r);

        offset = CenterStencil[StencilPosition::Center];
        col    = center_index;
        val    = center_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* REMOVED: Moved to the right hand side to make the matrix symmetric */
        if (!DirBC_Interior) {
            offset = CenterStencil[StencilPosition::Left];
            col    = left_index;
            val    = left_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);
        }

        offset = CenterStencil[StencilPosition::Right];
        col    = right_index;
        val    = right_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Bottom];
        col    = bottom_index;
        val    = bottom_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Top];
        col    = top_index;
        val    = top_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* REMOVED: Moved to the right hand side to make the matrix symmetric */
        if (!DirBC_Interior) {
            offset = CenterStencil[StencilPosition::BottomLeft];
            col    = bottom_left_index;
            val    = bottom_left_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);
        }

        offset = CenterStencil[StencilPosition::BottomRight];
        col    = bottom_right_index;
        val    = bottom_right_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* REMOVED: Moved to the right hand side to make the matrix symmetric */
        if (!DirBC_Interior) {
            offset = CenterStencil[StencilPosition::TopLeft];
            col    = top_left_index;
            val    = top_left_value;
            updateMatrixElement(solver_matrix, ptr, offset, row, col, val);
        }

        offset = CenterStencil[StencilPosition::TopRight];
        col    = top_right_index;
        val    = top_right_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);
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

        int center_nz_index = getSolverMatrixIndex(i_r, i_theta);

        int center_index       = grid.index(i_r, i_theta);
        int left_index         = grid.index(i_r - 1, i_theta);
        int right_index        = grid.index(i_r + 1, i_theta);
        int bottom_index       = grid.index(i_r, i_theta_M1);
        int top_index          = grid.index(i_r, i_theta_P1);
        int bottom_left_index  = grid.index(i_r - 1, i_theta_M1);
        int bottom_right_index = grid.index(i_r + 1, i_theta_M1);
        int top_left_index     = grid.index(i_r - 1, i_theta_P1);
        int top_right_index    = grid.index(i_r + 1, i_theta_P1);

        double left_value   = -coeff1 * (arr(center_index) + arr(left_index)); /* Left */
        double right_value  = -coeff2 * (arr(center_index) + arr(right_index)); /* Right */
        double bottom_value = -coeff3 * (att(center_index) + att(bottom_index)); /* Bottom */
        double top_value    = -coeff4 * (att(center_index) + att(top_index)); /* Top */

        double center_value =
            (+0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center_index] * fabs(detDF(center_index)) /* beta_{i,j} */
             - left_value /* Center: (Left) */
             - right_value /* Center: (Right) */
             - bottom_value /* Center: (Bottom) */
             - top_value /* Center: (Top) */
            );

        double bottom_left_value  = -0.25 * (art(left_index) + art(bottom_index)); /* Bottom Left */
        double bottom_right_value = +0.25 * (art(right_index) + art(bottom_index)); /* Bottom Right */
        double top_left_value     = +0.25 * (art(left_index) + art(top_index)); /* Top Left */
        double top_right_value    = -0.25 * (art(right_index) + art(top_index)); /* Top Right */

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r);

        offset = CenterStencil[StencilPosition::Center];
        col    = center_index;
        val    = center_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Left];
        col    = left_index;
        val    = left_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* Right REMOVED: Moved to the right hand side to make the matrix symmetric */

        offset = CenterStencil[StencilPosition::Bottom];
        col    = bottom_index;
        val    = bottom_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::Top];
        col    = top_index;
        val    = top_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        offset = CenterStencil[StencilPosition::BottomLeft];
        col    = bottom_left_index;
        val    = bottom_left_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* BottomRight REMOVED: Moved to the right hand side to make the matrix symmetric */

        offset = CenterStencil[StencilPosition::TopLeft];
        col    = top_left_index;
        val    = top_left_value;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);

        /* TopRight REMOVED: Moved to the right hand side to make the matrix symmetric */
    }
    /* ------------------------------------ */
    /* Node on the outer dirichlet boundary */
    /* ------------------------------------ */
    else if (i_r == grid.nr() - 1) {
        int center_nz_index = getSolverMatrixIndex(i_r, i_theta);

        int center_index = grid.index(i_r, i_theta);

        /* Fill matrix row of (i,j) */
        row = center_index;
        ptr = center_nz_index;

        const Stencil& CenterStencil = getStencil(i_r);

        offset = CenterStencil[StencilPosition::Center];
        col    = center_index;
        val    = 1.0;
        updateMatrixElement(solver_matrix, ptr, offset, row, col, val);
    }
}

void DirectSolver_COO_MUMPS_Take::buildSolverMatrixCircleSection(const int i_r, SparseMatrixCOO<double>& solver_matrix)
{
    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    ConstVector<double> arr        = level_cache_.arr();
    ConstVector<double> att        = level_cache_.att();
    ConstVector<double> art        = level_cache_.art();
    ConstVector<double> detDF      = level_cache_.detDF();
    ConstVector<double> coeff_beta = level_cache_.coeff_beta();

    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        // Build solver matrix at the current node
        nodeBuildSolverMatrixTake(i_r, i_theta, grid_, DirBC_Interior_, solver_matrix, arr, att, art, detDF,
                                  coeff_beta);
    }
}

void DirectSolver_COO_MUMPS_Take::buildSolverMatrixRadialSection(const int i_theta,
                                                                 SparseMatrixCOO<double>& solver_matrix)
{
    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    ConstVector<double> arr        = level_cache_.arr();
    ConstVector<double> att        = level_cache_.att();
    ConstVector<double> art        = level_cache_.art();
    ConstVector<double> detDF      = level_cache_.detDF();
    ConstVector<double> coeff_beta = level_cache_.coeff_beta();

    for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
        // Build solver matrix at the current node
        nodeBuildSolverMatrixTake(i_r, i_theta, grid_, DirBC_Interior_, solver_matrix, arr, att, art, detDF,
                                  coeff_beta);
    }
}

// clang-format off

/* ------------------------------------------------------------------------ */
/* If the indexing is not smoother-based, please adjust the access patterns */
SparseMatrixCOO<double> DirectSolver_COO_MUMPS_Take::buildSolverMatrix()
{
    const int n   = grid_.numberOfNodes();
    const int nnz = getNonZeroCountSolverMatrix();

    // Although the matrix is symmetric, we need to store all its entries, so we disable the symmetry.
    SparseMatrixCOO<double> solver_matrix(n, n, nnz);
    solver_matrix.is_symmetric(false);

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
        /* Multi-threaded execution */
        #pragma omp parallel num_threads(num_omp_threads_)
        {
            /* Circle Section */
            #pragma omp for nowait
            for (int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
                buildSolverMatrixCircleSection(i_r, solver_matrix);
            }
            /* Radial Section */
            #pragma omp for nowait
            for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
                buildSolverMatrixRadialSection(i_theta, solver_matrix);
            }
        }
    }

    /* Mumps: In the case of symmetric matrices, only half of the matrix should be provided. */
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
