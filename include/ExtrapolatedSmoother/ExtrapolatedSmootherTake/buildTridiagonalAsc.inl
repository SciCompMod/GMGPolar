#pragma once

namespace extrapolated_smoother_take
{

static KOKKOS_INLINE_FUNCTION void updateMatrixElement(const BatchedTridiagonalSolver<double>& solver, const int batch,
                                                       const int row, const int column, const double value)
{
    if (row == column)
        solver.set_main_diagonal(batch, row, value);
    else if (row == column - 1)
        solver.set_sub_diagonal(batch, row, value);
    else if (row == 0 && column == solver.matrixDimension() - 1)
        solver.set_cyclic_corner(batch, value);
}

static KOKKOS_INLINE_FUNCTION void nodeBuildTridiagonalSolverMatricesCircleSection(
    const int i_r, const int i_theta, const PolarGrid<Kokkos::HostSpace>& grid, const bool DirBC_Interior,
    const BatchedTridiagonalSolver<double>& circle_tridiagonal_solver,
    const BatchedTridiagonalSolver<double>& radial_tridiagonal_solver, HostConstVector<double>& arr,
    HostConstVector<double>& att, HostConstVector<double>& art, HostConstVector<double>& detDF,
    HostConstVector<double>& coeff_beta)
{
    using extrapolated_smoother_take::updateMatrixElement;

    KOKKOS_ASSERT(i_r >= 0 && i_r < grid.nr());
    KOKKOS_ASSERT(i_theta >= 0 && i_theta < grid.ntheta());

    const int numberSmootherCircles = grid.numberSmootherCircles();
    const int lengthRadialSmoother  = grid.lengthRadialSmoother();

    KOKKOS_ASSERT(numberSmootherCircles >= 3);
    KOKKOS_ASSERT(lengthRadialSmoother >= 3);

    int ptr, offset;
    int row, column, col;
    double value, val;

    /* ------------------------------------------ */
    /* Node in the interior of the Circle Section */
    /* ------------------------------------------ */
    if (i_r > 0 && i_r < numberSmootherCircles) {
        /* i_r = numberSmootherCircles-1 is included here! */
        const double h1 = grid.radialSpacing(i_r - 1);
        const double h2 = grid.radialSpacing(i_r);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff1 = 0.5 * (k1 + k2) / h1;
        const double coeff2 = 0.5 * (k1 + k2) / h2;
        const double coeff3 = 0.5 * (h1 + h2) / k1;
        const double coeff4 = 0.5 * (h1 + h2) / k2;
        const double coeff5 = 0.25 * (h1 + h2) * (k1 + k2);

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int left   = grid.index(i_r - 1, i_theta);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int center = grid.index(i_r, i_theta);
        const int top    = grid.index(i_r, i_theta_P1);
        const int right  = grid.index(i_r + 1, i_theta);

        const int center_index = i_theta;
        const int bottom_index = i_theta_M1;
        const int top_index    = i_theta_P1;

        const auto& solver = circle_tridiagonal_solver;
        const int batch    = i_r;

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
            /* or */
            /* i_theta % 2 == 0 */
            /* | o | o | o | */
            /* |   |   |   | */
            /* | x | O | x | */
            /* |   |   |   | */
            /* | o | o | o | */

            /* Center: (Left, Right, Bottom, Top) */
            row    = center_index;
            column = center_index;
            value  = coeff5 * coeff_beta[center] * Kokkos::fabs(detDF[center]) + coeff1 * (arr[center] + arr[left]) +
                    coeff2 * (arr[center] + arr[right]) + coeff3 * (att[center] + att[bottom]) +
                    coeff4 * (att[center] + att[top]);
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
            /* or */
            /* i_theta % 2 == 0 */
            /* | o | o | o | */
            /* |   |   |   | */
            /* | o | X | o | */
            /* |   |   |   | */
            /* | o | o | o | */

            if (i_theta & 1) {
                /* i_theta % 2 == 1 */

                /* Center: (Left, Right, Bottom, Top) */
                row    = center_index;
                column = center_index;
                value = coeff5 * coeff_beta[center] * Kokkos::fabs(detDF[center]) + coeff1 * (arr[center] + arr[left]) +
                        coeff2 * (arr[center] + arr[right]) + coeff3 * (att[center] + att[bottom]) +
                        coeff4 * (att[center] + att[top]);
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
        const auto& solver = circle_tridiagonal_solver;
        const int batch    = i_r;

        const int center_index = i_theta;

        /* Center: (Left, Right, Bottom, Top) */
        row    = center_index;
        column = center_index;
        value  = 1.0;
        updateMatrixElement(solver, batch, row, column, value);
    }
}

static KOKKOS_INLINE_FUNCTION void nodeBuildTridiagonalSolverMatricesRadialSection(
    const int i_r, const int i_theta, const PolarGrid<Kokkos::HostSpace>& grid, const bool DirBC_Interior,
    const BatchedTridiagonalSolver<double>& circle_tridiagonal_solver,
    const BatchedTridiagonalSolver<double>& radial_tridiagonal_solver, HostConstVector<double>& arr,
    HostConstVector<double>& att, HostConstVector<double>& art, HostConstVector<double>& detDF,
    HostConstVector<double>& coeff_beta)
{
    using extrapolated_smoother_take::updateMatrixElement;

    KOKKOS_ASSERT(i_r >= 0 && i_r < grid.nr());
    KOKKOS_ASSERT(i_theta >= 0 && i_theta < grid.ntheta());

    const int numberSmootherCircles = grid.numberSmootherCircles();
    const int lengthRadialSmoother  = grid.lengthRadialSmoother();

    KOKKOS_ASSERT(numberSmootherCircles >= 3);
    KOKKOS_ASSERT(lengthRadialSmoother >= 3);

    int ptr, offset;
    int row, column, col;
    double value, val;

    /* ------------------------------------------ */
    /* Node in the interior of the Radial Section */
    /* ------------------------------------------ */
    if (i_r > numberSmootherCircles && i_r < grid.nr() - 2) {
        const double h1 = grid.radialSpacing(i_r - 1);
        const double h2 = grid.radialSpacing(i_r);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff1 = 0.5 * (k1 + k2) / h1;
        const double coeff2 = 0.5 * (k1 + k2) / h2;
        const double coeff3 = 0.5 * (h1 + h2) / k1;
        const double coeff4 = 0.5 * (h1 + h2) / k2;
        const double coeff5 = 0.25 * (h1 + h2) * (k1 + k2);

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

        const auto& solver = radial_tridiagonal_solver;
        const int batch    = i_theta;

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
            value  = coeff5 * coeff_beta[center] * Kokkos::fabs(detDF[center]) + coeff1 * (arr[center] + arr[left]) +
                    coeff2 * (arr[center] + arr[right]) + coeff3 * (att[center] + att[bottom]) +
                    coeff4 * (att[center] + att[top]);
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
            /* or */
            /* i_r % 2 == 0 */
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
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center] * Kokkos::fabs(detDF[center]) +
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
        const double h1 = grid.radialSpacing(i_r - 1);
        const double h2 = grid.radialSpacing(i_r);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff1 = 0.5 * (k1 + k2) / h1;
        const double coeff2 = 0.5 * (k1 + k2) / h2;
        const double coeff3 = 0.5 * (h1 + h2) / k1;
        const double coeff4 = 0.5 * (h1 + h2) / k2;
        const double coeff5 = 0.25 * (h1 + h2) * (k1 + k2);

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int left   = grid.index(i_r - 1, i_theta);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int center = grid.index(i_r, i_theta);
        const int top    = grid.index(i_r, i_theta_P1);
        const int right  = grid.index(i_r + 1, i_theta);

        const int center_index = i_r - numberSmootherCircles;
        const int right_index  = i_r - numberSmootherCircles + 1;

        const auto& solver = radial_tridiagonal_solver;
        const int batch    = i_theta;

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
            value  = coeff5 * coeff_beta[center] * Kokkos::fabs(detDF[center]) + coeff1 * (arr[center] + arr[left]) +
                    coeff2 * (arr[center] + arr[right]) + coeff3 * (att[center] + att[bottom]) +
                    coeff4 * (att[center] + att[top]);
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
                value = coeff5 * coeff_beta[center] * Kokkos::fabs(detDF[center]) + coeff1 * (arr[center] + arr[left]) +
                        coeff2 * (arr[center] + arr[right]) + coeff3 * (att[center] + att[bottom]) +
                        coeff4 * (att[center] + att[top]);
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
        KOKKOS_ASSERT(i_r % 2 == 1);

        const double h1 = grid.radialSpacing(i_r - 1);
        const double h2 = grid.radialSpacing(i_r);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff1 = 0.5 * (k1 + k2) / h1;
        const double coeff2 = 0.5 * (k1 + k2) / h2;
        const double coeff3 = 0.5 * (h1 + h2) / k1;
        const double coeff4 = 0.5 * (h1 + h2) / k2;
        const double coeff5 = 0.25 * (h1 + h2) * (k1 + k2);

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

        const auto& solver = radial_tridiagonal_solver;
        const int batch    = i_theta;

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
            value  = coeff5 * coeff_beta[center] * Kokkos::fabs(detDF[center]) + coeff1 * (arr[center] + arr[left]) +
                    coeff2 * (arr[center] + arr[right]) + coeff3 * (att[center] + att[bottom]) +
                    coeff4 * (att[center] + att[top]);
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
            value  = coeff5 * coeff_beta[center] * Kokkos::fabs(detDF[center]) + coeff1 * (arr[center] + arr[left]) +
                    coeff2 * (arr[center] + arr[right]) + coeff3 * (att[center] + att[bottom]) +
                    coeff4 * (att[center] + att[top]);
            updateMatrixElement(solver, batch, row, column, value);
        }
    }
    /* ------------------------------------------ */
    /* Radial Section: Node on the outer boundary */
    /* ------------------------------------------ */
    else if (i_r == grid.nr() - 1) {
        KOKKOS_ASSERT(i_r % 2 == 0);

        const int center_index = i_r - numberSmootherCircles;
        const int left_index   = i_r - numberSmootherCircles - 1;

        const auto& solver = radial_tridiagonal_solver;
        const int batch    = i_theta;

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

} // namespace extrapolated_smoother_take

template <class LevelCacheType>
void ExtrapolatedSmootherTake<LevelCacheType>::buildTridiagonalSolverMatrices()
{
    using extrapolated_smoother_take::nodeBuildTridiagonalSolverMatricesCircleSection;
    using extrapolated_smoother_take::nodeBuildTridiagonalSolverMatricesRadialSection;

    const PolarGrid<Kokkos::HostSpace>& grid             = ExtrapolatedSmoother<LevelCacheType>::grid_;
    const bool DirBC_Interior         = ExtrapolatedSmoother<LevelCacheType>::DirBC_Interior_;
    const LevelCacheType& level_cache = ExtrapolatedSmoother<LevelCacheType>::level_cache_;

    assert(level_cache.cacheDensityProfileCoefficients());
    assert(level_cache.cacheDomainGeometry());

    HostConstVector<double> arr        = level_cache.arr();
    HostConstVector<double> att        = level_cache.att();
    HostConstVector<double> art        = level_cache.art();
    HostConstVector<double> detDF      = level_cache.detDF();
    HostConstVector<double> coeff_beta = level_cache.coeff_beta();

    const BatchedTridiagonalSolver<double>& circle_tridiagonal_solver = circle_tridiagonal_solver_;
    const BatchedTridiagonalSolver<double>& radial_tridiagonal_solver = radial_tridiagonal_solver_;

    /* We split the loops into two regions to better respect the */
    /* access patterns of the smoother and improve cache locality. */

    // The For loop matches circular access pattern */
    Kokkos::parallel_for(
        "Extrapolated Smoother Take: Build Tridiagonal Matrices (Circular)",
        Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
            {0, 0}, // Starting point of the index space
            {grid.numberSmootherCircles(), grid.ntheta()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_r, const int i_theta) {
            nodeBuildTridiagonalSolverMatricesCircleSection(i_r, i_theta, grid, DirBC_Interior,
                                                            circle_tridiagonal_solver, radial_tridiagonal_solver, arr,
                                                            att, art, detDF, coeff_beta);
        });

    /* For loop matches radial access pattern */
    Kokkos::parallel_for(
        "Extrapolated Smoother Take: Build Tridiagonal Matrices (Radial)",
        Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
            {0, grid.numberSmootherCircles()}, // Starting point of the index space
            {grid.ntheta(), grid.nr()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_theta, const int i_r) {
            nodeBuildTridiagonalSolverMatricesRadialSection(i_r, i_theta, grid, DirBC_Interior,
                                                            circle_tridiagonal_solver, radial_tridiagonal_solver, arr,
                                                            att, art, detDF, coeff_beta);
        });

    Kokkos::fence();
}
