#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTake/extrapolatedSmootherTake.h"

/* Tridiagonal matrices */
static inline void updateMatrixElement(BatchedTridiagonalSolver<double>& solver, int batch, int row, int column,
                                       double value)
{
    if (row == column)
        solver.main_diagonal(batch, row) = value;
    else if (row == column - 1)
        solver.sub_diagonal(batch, row) = value;
    else if (row == 0 && column == solver.matrixDimension() - 1)
        solver.cyclic_corner(batch) = value;
}

/* Inner Boundary COO/CSR matrix */
#ifdef GMGPOLAR_USE_MUMPS
static inline void updateCOOCSRMatrixElement(SparseMatrixCOO<double>& matrix, int ptr, int offset, int row, int col,
                                             double val)
{
    matrix.row_index(ptr + offset) = row;
    matrix.col_index(ptr + offset) = col;
    matrix.value(ptr + offset)     = val;
}
#else
static inline void updateCOOCSRMatrixElement(SparseMatrixCSR<double>& matrix, int ptr, int offset, int row, int col,
                                             double val)
{
    matrix.row_nz_index(row, offset) = col;
    matrix.row_nz_entry(row, offset) = val;
}
#endif

void ExtrapolatedSmootherTake::nodeBuildAscTake(int i_r, int i_theta, const PolarGrid& grid, bool DirBC_Interior,
                                                MatrixType& inner_boundary_circle_matrix,
                                                BatchedTridiagonalSolver<double>& circle_tridiagonal_solver,
                                                BatchedTridiagonalSolver<double>& radial_tridiagonal_solver,
                                                ConstVector<double>& arr, ConstVector<double>& att,
                                                ConstVector<double>& art, ConstVector<double>& detDF,
                                                ConstVector<double>& coeff_beta)
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
    if (i_r > 0 && i_r < numberSmootherCircles) {
        /* i_r = numberSmootherCircles-1 is included here! */
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

        const int left   = grid.index(i_r - 1, i_theta);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int center = grid.index(i_r, i_theta);
        const int top    = grid.index(i_r, i_theta_P1);
        const int right  = grid.index(i_r + 1, i_theta);

        int center_index = i_theta;
        int bottom_index = i_theta_M1;
        int top_index    = i_theta_P1;

        auto& solver = circle_tridiagonal_solver;
        int batch    = i_r;

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
            /* or */ /* i_theta % 2 == 0 */
            /* | o | o | o | */
            /* |   |   |   | */
            /* | x | O | x | */
            /* |   |   |   | */
            /* | o | o | o | */

            /* Center: (Left, Right, Bottom, Top) */
            row    = center_index;
            column = center_index;
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center] * fabs(detDF[center]) +
                    coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +
                    coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);
            updateMatrixElement(solver, batch, row, column, value);

            /* Bottom */
            row    = center_index;
            column = bottom_index;
            value  = -coeff3 * (att[center] + att[bottom]);
            updateMatrixElement(solver, batch, row, column, value);

            /* Top */
            row    = center_index;
            column = top_index;
            value  = -coeff4 * (att[center] + att[top]);
            updateMatrixElement(solver, batch, row, column, value);
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
            /* or */ /* i_theta % 2 == 0 */ /* | o | o | o | */
            /* |   |   |   | */
            /* | o | X | o | */
            /* |   |   |   | */
            /* | o | o | o | */

            if (i_theta & 1) {
                /* i_theta % 2 == 1 */

                /* Center: (Left, Right, Bottom, Top) */
                row    = center_index;
                column = center_index;
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center] * fabs(detDF[center]) +
                        coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +
                        coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);
                updateMatrixElement(solver, batch, row, column, value);
            }
            else {
                /* i_theta % 2 == 0 */

                /* Center: Coarse */
                row    = center_index;
                column = center_index;
                value  = 1.0;
                updateMatrixElement(solver, batch, row, column, value);
            }
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
            auto& matrix              = inner_boundary_circle_matrix;
            const int center_index    = i_theta;
            const int center_nz_index = getCircleAscIndex(i_r, i_theta);

            /* Fill matrix row of (i,j) */
            row = center_index;
            ptr = center_nz_index;

            const Stencil& CenterStencil = getStencil(i_r, i_theta);

            offset = CenterStencil[StencilPosition::Center];
            col    = center_index;
            val    = 1.0;
            updateCOOCSRMatrixElement(matrix, ptr, offset, row, col, val);
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

            auto& matrix = inner_boundary_circle_matrix;

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

                const double center_value = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center] * fabs(detDF[center]) +
                                            coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +
                                            coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);
                const double left_value = -coeff1 * (arr[center] + arr[left]);

                /* Fill matrix row of (i,j) */
                row = center_index;
                ptr = center_nz_index;

                const Stencil& CenterStencil = getStencil(i_r, i_theta);

                offset = CenterStencil[StencilPosition::Center];
                col    = center_index;
                val    = center_value;
                updateCOOCSRMatrixElement(matrix, ptr, offset, row, col, val);

                offset = CenterStencil[StencilPosition::Left];
                col    = left_index;
                val    = left_value;
                updateCOOCSRMatrixElement(matrix, ptr, offset, row, col, val);
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
                col    = center_index;
                val    = 1.0;
                updateCOOCSRMatrixElement(matrix, ptr, offset, row, col, val);
            }
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

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int left   = grid.index(i_r - 1, i_theta);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int center = grid.index(i_r, i_theta);
        const int top    = grid.index(i_r, i_theta_P1);
        const int right  = grid.index(i_r + 1, i_theta);

        const int center_index = i_r - numberSmootherCircles;
        const int left_index   = i_r - numberSmootherCircles - 1;
        const int right_index  = i_r - numberSmootherCircles + 1;

        auto& solver = radial_tridiagonal_solver;
        int batch    = i_theta;

        /* ------------------- */
        /* Tridiagonal Section */
        /* i_theta % 2 == 1 */
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

            /* Center: (Left, Right, Bottom, Top) */
            row    = center_index;
            column = center_index;
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center] * fabs(detDF[center]) +
                    coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +
                    coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);
            updateMatrixElement(solver, batch, row, column, value);

            /* Left */
            row    = center_index;
            column = left_index;
            value  = -coeff1 * (arr[center] + arr[left]);
            updateMatrixElement(solver, batch, row, column, value);

            /* Right */
            row    = center_index;
            column = right_index;
            value  = -coeff2 * (arr[center] + arr[right]);
            updateMatrixElement(solver, batch, row, column, value);
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
            /* or */ /* i_r % 2 == 0 */
            /* ---------- */
            /* o   o   o  */
            /* ---------- */
            /* o   X   o  */
            /* ---------- */
            /* o   o   o  */
            /* ---------- */

            if (i_r & 1) {
                /* i_r % 2 == 1 */

                /* Center: (Left, Right, Bottom, Top) */
                row    = center_index;
                column = center_index;
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center] * fabs(detDF[center]) +
                        coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +
                        coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);
                updateMatrixElement(solver, batch, row, column, value);
            }
            else {
                /* i_r % 2 == 0 */

                /* Center: Coarse */
                row    = center_index;
                column = center_index;
                value  = 1.0;
                updateMatrixElement(solver, batch, row, column, value);
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

        const int left   = grid.index(i_r - 1, i_theta);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int center = grid.index(i_r, i_theta);
        const int top    = grid.index(i_r, i_theta_P1);
        const int right  = grid.index(i_r + 1, i_theta);

        const int center_index = i_r - numberSmootherCircles;
        const int right_index  = i_r - numberSmootherCircles + 1;

        auto& solver = radial_tridiagonal_solver;
        int batch    = i_theta;

        if (i_theta & 1) {
            /* i_theta % 2 == 1 and i_r % 2 == 1 */
            /* | x | o | x || o   x   o   x  */
            /* |   |   |   || -------------- */
            /* | o | o | o || O   o   o   o  */
            /* |   |   |   || -------------- */
            /* | x | o | x || o   x   o   x  */
            /* or */
            /* i_theta % 2 == 1 and i_r % 2 == 0 */
            /* | o | x | o || x   o   x   o  */
            /* |   |   |   || -------------- */
            /* | o | o | o || O   o   o   o  */
            /* |   |   |   || -------------- */
            /* | o | x | o || x   o   x   o  */

            /* Center: (Left, Right, Bottom, Top) */
            row    = center_index;
            column = center_index;
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center] * fabs(detDF[center]) +
                    coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +
                    coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);
            updateMatrixElement(solver, batch, row, column, value);

            /* Right */
            row    = center_index;
            column = right_index;
            value  = -coeff2 * (arr[center] + arr[right]);
            updateMatrixElement(solver, batch, row, column, value);
        }
        else {
            if (i_r & 1) {
                /* i_theta % 2 == 0 and i_r % 2 == 1 */
                /* | o | o | o || o   o   o   o  */
                /* |   |   |   || -------------- */
                /* | x | o | x || O   x   o   x  */
                /* |   |   |   || -------------- */
                /* | o | o | o || o   o   o   o  */

                /* Center: (Left, Right, Bottom, Top) */
                row    = center_index;
                column = center_index;
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center] * fabs(detDF[center]) +
                        coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +
                        coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);
                updateMatrixElement(solver, batch, row, column, value);
            }
            else {
                /* i_theta % 2 == 0 and i_r % 2 == 0 */
                /* | o | o | o || o   o   o   o  */
                /* |   |   |   || -------------- */
                /* | o | x | o || X   o   x   o  */
                /* |   |   |   || -------------- */
                /* | o | o | o || o   o   o   o  */

                /* Center: Coarse */
                row    = center_index;
                column = center_index;
                value  = 1.0;
                updateMatrixElement(solver, batch, row, column, value);
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

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int left   = grid.index(i_r - 1, i_theta);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int center = grid.index(i_r, i_theta);
        const int top    = grid.index(i_r, i_theta_P1);
        const int right  = grid.index(i_r + 1, i_theta);

        const int center_index = i_r - numberSmootherCircles;
        const int left_index   = i_r - numberSmootherCircles - 1;
        const int right_index  = i_r - numberSmootherCircles + 1;

        auto& solver = radial_tridiagonal_solver;
        int batch    = i_theta;

        if (i_theta & 1) {
            /* i_theta % 2 == 1 */
            /* ---------------|| */
            /* o   x   o   x  || */
            /* ---------------|| */
            /* o   o   O   o  || */
            /* ---------------|| */
            /* o   x   o   x  || */
            /* ---------------|| */

            /* Center: (Left, Right, Bottom, Top) */
            row    = center_index;
            column = center_index;
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center] * fabs(detDF[center]) +
                    coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +
                    coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);
            updateMatrixElement(solver, batch, row, column, value);

            /* Left */
            row    = center_index;
            column = left_index;
            value  = -coeff1 * (arr[center] + arr[left]);
            updateMatrixElement(solver, batch, row, column, value); /* Right */
            row    = center_index;
            column = right_index;
            value  = 0.0; /* Make tridiagonal matrix symmetric */
            updateMatrixElement(solver, batch, row, column, value);
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

            /* Center: (Left, Right, Bottom, Top) */
            row    = center_index;
            column = center_index;
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center] * fabs(detDF[center]) +
                    coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) +
                    coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top]);
            updateMatrixElement(solver, batch, row, column, value);
        }
    }
    /* ------------------------------------------ */
    /* Radial Section: Node on the outer boundary */
    /* ------------------------------------------ */
    else if (i_r == grid.nr() - 1) {
        assert(!i_r % 2 == 0);

        int center_index = i_r - numberSmootherCircles;
        int left_index   = i_r - numberSmootherCircles - 1;

        auto& solver = radial_tridiagonal_solver;
        int batch    = i_theta;

        if (i_theta & 1) {
            /* i_theta % 2 == 1 */
            /* -----------|| */
            /* x   o   x  || */
            /* -----------|| */
            /* o   o   O  || */
            /* -----------|| */
            /* x   o   x  || */
            /* -----------|| */

            /* Fill matrix row of (i,j) */
            row    = center_index;
            column = center_index;
            value  = 1.0;
            updateMatrixElement(solver, batch, row, column, value);

            row    = center_index;
            column = left_index;
            value  = 0.0; /* Make tridiagonal matrix symmetric */
            updateMatrixElement(solver, batch, row, column, value);
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

            /* Fill matrix row of (i,j) */
            row    = center_index;
            column = center_index;
            value  = 1.0;
            updateMatrixElement(solver, batch, row, column, value);
        }
    }
}

void ExtrapolatedSmootherTake::buildAscCircleSection(int i_r)
{
    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    ConstVector<double> arr        = level_cache_.arr();
    ConstVector<double> att        = level_cache_.att();
    ConstVector<double> art        = level_cache_.art();
    ConstVector<double> detDF      = level_cache_.detDF();
    ConstVector<double> coeff_beta = level_cache_.coeff_beta();

    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        // Build Asc at the current node
        nodeBuildAscTake(i_r, i_theta, grid_, DirBC_Interior_, inner_boundary_circle_matrix_,
                         circle_tridiagonal_solver_, radial_tridiagonal_solver_, arr, att, art, detDF, coeff_beta);
    }
}

void ExtrapolatedSmootherTake::buildAscRadialSection(int i_theta)
{
    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    ConstVector<double> arr        = level_cache_.arr();
    ConstVector<double> att        = level_cache_.att();
    ConstVector<double> art        = level_cache_.art();
    ConstVector<double> detDF      = level_cache_.detDF();
    ConstVector<double> coeff_beta = level_cache_.coeff_beta();

    for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
        // Build Asc at the current node
        nodeBuildAscTake(i_r, i_theta, grid_, DirBC_Interior_, inner_boundary_circle_matrix_,
                         circle_tridiagonal_solver_, radial_tridiagonal_solver_, arr, att, art, detDF, coeff_beta);
    }
}

void ExtrapolatedSmootherTake::buildAscMatrices()
{
    /* -------------------------------------- */
    /* Part 1: Allocate Asc Smoother matrices */
    /* -------------------------------------- */
    // BatchedTridiagonalSolvers allocations are handled in the SmootherTake constructor.
    // circle_tridiagonal_solver_[batch_index=0] is unitialized. Use inner_boundary_circle_matrix_ instead.

#ifdef GMGPOLAR_USE_MUMPS
    // Although the matrix is symmetric, we need to store all its entries, so we disable the symmetry.
    const int inner_i_r           = 0;
    const int inner_nnz           = getNonZeroCountCircleAsc(inner_i_r);
    const int num_circle_nodes    = grid_.ntheta();
    inner_boundary_circle_matrix_ = SparseMatrixCOO<double>(num_circle_nodes, num_circle_nodes, inner_nnz);
    inner_boundary_circle_matrix_.is_symmetric(false);
#else
    std::function<int(int)> nnz_per_row = [&](int i_theta) {
        if (DirBC_Interior_)
            return 1;
        else
            return i_theta % 2 == 0 ? 1 : 2;
    };
    const int num_circle_nodes    = grid_.ntheta();
    inner_boundary_circle_matrix_ = SparseMatrixCSR<double>(num_circle_nodes, num_circle_nodes, nnz_per_row);
#endif

    /* ---------------------------------- */
    /* Part 2: Fill Asc Smoother matrices */
    /* ---------------------------------- */

#pragma omp parallel num_threads(num_omp_threads_)
    {
#pragma omp for nowait
        for (int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
            buildAscCircleSection(i_r);
        }

#pragma omp for nowait
        for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
            buildAscRadialSection(i_theta);
        }
    }

    circle_tridiagonal_solver_.setup();
    radial_tridiagonal_solver_.setup();

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
