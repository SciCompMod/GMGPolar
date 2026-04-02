#pragma once

namespace extrapolated_smoother_take
{

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

} // namespace extrapolated_smoother_take

template <class LevelCacheType>
void ExtrapolatedSmootherTake<LevelCacheType>::nodeBuildTridiagonalSolverMatrices(
    int i_r, int i_theta, const PolarGrid& grid, bool DirBC_Interior,
    BatchedTridiagonalSolver<double>& circle_tridiagonal_solver,
    BatchedTridiagonalSolver<double>& radial_tridiagonal_solver, ConstVector<double>& arr, ConstVector<double>& att,
    ConstVector<double>& art, ConstVector<double>& detDF, ConstVector<double>& coeff_beta)
{
    using extrapolated_smoother_take::updateMatrixElement;

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
        // The inner boundary circle line are is handled by the inner_boundary_solver, so we fill in the identity matrix.
        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        auto& solver    = circle_tridiagonal_solver;
        const int batch = i_r;

        const int center_index = i_theta;
        const int bottom_index = i_theta_M1;
        const int top_index    = i_theta_P1;

        /* Center: (Left, Right, Bottom, Top) */
        row    = center_index;
        column = center_index;
        value  = 1.0;
        updateMatrixElement(solver, batch, row, column, value);

        /* Bottom */
        row    = center_index;
        column = bottom_index;
        value  = 0.0;
        updateMatrixElement(solver, batch, row, column, value);

        /* Top */
        row    = center_index;
        column = top_index;
        value  = 0.0;
        updateMatrixElement(solver, batch, row, column, value);
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

template <class LevelCacheType>
void ExtrapolatedSmootherTake<LevelCacheType>::buildTridiagonalSolverMatrices()
{
    const PolarGrid& grid                         = ExtrapolatedSmootherTake<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = ExtrapolatedSmootherTake<LevelCacheType>::level_cache_;
    const bool DirBC_Interior                     = ExtrapolatedSmootherTake<LevelCacheType>::DirBC_Interior_;
    const int num_omp_threads                     = ExtrapolatedSmootherTake<LevelCacheType>::num_omp_threads_;

    assert(level_cache.cacheDensityProfileCoefficients());
    assert(level_cache.cacheDomainGeometry());

    ConstVector<double> arr        = level_cache.arr();
    ConstVector<double> att        = level_cache.att();
    ConstVector<double> art        = level_cache.art();
    ConstVector<double> detDF      = level_cache.detDF();
    ConstVector<double> coeff_beta = level_cache.coeff_beta();

#pragma omp parallel num_threads(num_omp_threads)
    {
#pragma omp for nowait
        for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++) {
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                nodeBuildTridiagonalSolverMatrices(i_r, i_theta, grid, DirBC_Interior, circle_tridiagonal_solver_,
                                                   radial_tridiagonal_solver_, arr, att, art, detDF, coeff_beta);
            }
        }

#pragma omp for nowait
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                nodeBuildTridiagonalSolverMatrices(i_r, i_theta, grid, DirBC_Interior, circle_tridiagonal_solver_,
                                                   radial_tridiagonal_solver_, arr, att, art, detDF, coeff_beta);
            }
        }
    }
}
