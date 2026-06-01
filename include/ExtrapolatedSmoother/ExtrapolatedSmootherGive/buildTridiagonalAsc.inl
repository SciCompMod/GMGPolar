#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherGive/extrapolatedSmootherGive.h"

#include "../../../include/Definitions/geometry_helper.h"

namespace extrapolated_smoother_give
{

static KOKKOS_INLINE_FUNCTION void updateMatrixElement(const BatchedTridiagonalSolver<double>& solver, const int batch,
                                                       const int row, const int column, const double value)
{
    if (row == column)
        solver.increase_main_diagonal(batch, row, value);
    else if (row == column - 1)
        solver.increase_sub_diagonal(batch, row, value);
    else if (row == 0 && column == solver.matrixDimension() - 1)
        solver.increase_cyclic_corner(batch, value);
}

template <typename LevelCacheType>
static KOKKOS_INLINE_FUNCTION void nodeBuildTridiagonalSolverMatricesCircleSection(
    const int i_r, const int i_theta, const PolarGrid<DefaultMemorySpace>& grid, const LevelCacheType& level_cache,
    const bool DirBC_Interior, const BatchedTridiagonalSolver<double>& circle_tridiagonal_solver,
    const BatchedTridiagonalSolver<double>& radial_tridiagonal_solver)
{
    using extrapolated_smoother_give::updateMatrixElement;

    KOKKOS_ASSERT(i_r >= 0 && i_r < grid.nr());
    KOKKOS_ASSERT(i_theta >= 0 && i_theta < grid.ntheta());

    const int numberSmootherCircles = grid.numberSmootherCircles();
    const int lengthRadialSmoother  = grid.lengthRadialSmoother();

    KOKKOS_ASSERT(numberSmootherCircles >= 3);
    KOKKOS_ASSERT(lengthRadialSmoother >= 3);

    /* ---------------------------------------- */
    /* Compute or retrieve stencil coefficients */
    /* ---------------------------------------- */
    const int center    = grid.index(i_r, i_theta);
    const double radius = grid.radius(i_r);
    const double theta  = grid.theta(i_theta);

    double coeff_beta, arr, att, art, detDF;
    level_cache.obtainValues(i_r, i_theta, center, radius, theta, coeff_beta, arr, att, art, detDF);

    int row, column;
    double value;
    /* ------------------------------------------ */
    /* Node in the interior of the Circle Section */
    /* ------------------------------------------ */
    if (i_r > 0 && i_r < numberSmootherCircles - 1) {
        const double h1 = grid.radialSpacing(i_r - 1);
        const double h2 = grid.radialSpacing(i_r);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff1 = 0.5 * (k1 + k2) / h1;
        const double coeff2 = 0.5 * (k1 + k2) / h2;
        const double coeff3 = 0.5 * (h1 + h2) / k1;
        const double coeff4 = 0.5 * (h1 + h2) / k2;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int center_index = i_theta;
        const int left_index   = i_theta;
        const int right_index  = i_theta;
        const int bottom_index = i_theta_M1;
        const int top_index    = i_theta_P1;

        const auto& center_solver = circle_tridiagonal_solver;
        const int center_batch    = i_r;
        const auto& left_solver   = circle_tridiagonal_solver;
        const int left_batch      = i_r - 1;
        const auto& right_solver  = circle_tridiagonal_solver;
        const int right_batch     = i_r + 1;

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

            /* Fill matrix row of (i,j) */
            row    = center_index;
            column = center_index;
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * Kokkos::fabs(detDF); /* Center: beta_{i,j} */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            row    = center_index;
            column = bottom_index;
            value  = -coeff3 * att; /* Bottom */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            row    = center_index;
            column = top_index;
            value  = -coeff4 * att; /* Top */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            row    = center_index;
            column = center_index;
            value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            /* Fill matrix row of (i,j-1) */
            row    = bottom_index;
            column = center_index;
            value  = -coeff3 * att; /* Top */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            row    = bottom_index;
            column = bottom_index;
            value  = coeff3 * att; /* Center: (Top) */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            /* Fill matrix row of (i,j+1) */
            row    = top_index;
            column = center_index;
            value  = -coeff4 * att; /* Bottom */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            row    = top_index;
            column = top_index;
            value  = coeff4 * att; /* Center: (Bottom) */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            if (i_theta & 1) {
                /* i_theta % 2 == 1 */
                /* | x | o | x | */
                /* |   |   |   | */
                /* | o | O | o | */
                /* |   |   |   | */
                /* | x | o | x | */

                /* Fill matrix row of (i-1,j) */
                // The inner boundary circle line are is handled by the inner_boundary_mumps_solver, so we fill in the identity matrix.
                if (i_r > 1) {
                    row    = left_index;
                    column = left_index;
                    value  = coeff1 * arr; /* Center: (Right) */
                    updateMatrixElement(left_solver, left_batch, row, column, value);
                }

                /* Fill matrix row of (i+1,j) */
                row    = right_index;
                column = right_index;
                value  = coeff2 * arr; /* Center: (Left) */
                updateMatrixElement(right_solver, right_batch, row, column, value);
            }
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

            if (i_theta & 1) { /* i_theta % 2 == 1 */
                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * Kokkos::fabs(detDF); /* Center: beta_{i,j} */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = center_index;
                column = center_index;
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
                updateMatrixElement(center_solver, center_batch, row, column, value);
            }
            else { /* i_theta % 2 == 0 */
                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 1.0;
                updateMatrixElement(center_solver, center_batch, row, column, value);

                /* Fill matrix row of (i,j-1) */
                row    = bottom_index;
                column = bottom_index;
                value  = coeff3 * att; /* Center: (Top) */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                /* Fill matrix row of (i,j+1) */
                row    = top_index;
                column = top_index;
                value  = coeff4 * att; /* Center: (Bottom) */
                updateMatrixElement(center_solver, center_batch, row, column, value);
            }
            /* Fill matrix row of (i-1,j) */
            row    = left_index;
            column = left_index;
            value  = coeff1 * arr; /* Center: (Right) */
            updateMatrixElement(left_solver, left_batch, row, column, value);

            /* Fill matrix row of (i+1,j) */
            row    = right_index;
            column = right_index;
            value  = coeff2 * arr; /* Center: (Left) */
            updateMatrixElement(right_solver, right_batch, row, column, value);
        }
    }
    /* ------------------------------------------ */
    /* Circle Section: Node in the inner boundary */
    /* ------------------------------------------ */
    else if (i_r == 0) {
        // The inner boundary circle line are is handled by the inner_boundary_mumps_solver, so we fill in the identity matrix.
        const auto& center_solver = circle_tridiagonal_solver;
        const int center_batch    = i_r;
        const auto& right_solver  = circle_tridiagonal_solver;
        const int right_batch     = i_r + 1;

        /* Fill result(i,j) */
        const double h2 = grid.radialSpacing(i_r);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff2 = 0.5 * (k1 + k2) / h2;

        const int center_index = i_theta;
        const int right_index  = i_theta;

        /* Fill matrix row of (i,j) */
        row    = center_index;
        column = center_index;
        value  = 1.0;
        updateMatrixElement(center_solver, center_batch, row, column, value);

        /* Fill matrix row of (i+1,j) */
        row    = right_index;
        column = right_index;
        value  = coeff2 * arr; /* Center: (Left) */
        updateMatrixElement(right_solver, right_batch, row, column, value);
    }
    /* ------------------------------------------- */
    /* Circle Section: Node next to radial section */
    /* ------------------------------------------- */
    else if (i_r == numberSmootherCircles - 1) {
        KOKKOS_ASSERT(i_r > 1);

        const double h1 = grid.radialSpacing(i_r - 1);
        const double h2 = grid.radialSpacing(i_r);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff1 = 0.5 * (k1 + k2) / h1;
        const double coeff2 = 0.5 * (k1 + k2) / h2;
        const double coeff3 = 0.5 * (h1 + h2) / k1;
        const double coeff4 = 0.5 * (h1 + h2) / k2;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int center_index = i_theta;
        const int left_index   = i_theta;
        const int right_index  = 0;
        const int bottom_index = i_theta_M1;
        const int top_index    = i_theta_P1;

        const auto& center_solver = circle_tridiagonal_solver;
        const int center_batch    = i_r;
        const auto& left_solver   = circle_tridiagonal_solver;
        const int left_batch      = i_r - 1;
        const auto& right_solver  = radial_tridiagonal_solver;
        const int right_batch     = i_theta;

        if (i_r & 1) {
            if (i_theta & 1) {
                /* i_r % 2 == 1 and i_theta % 2 == 1 */
                /* | o | x | o || x   o   x   o  */
                /* |   |   |   || -------------- */
                /* | o | o | O || o   o   o   o  */
                /* |   |   |   || -------------- */
                /* | o | x | o || x   o   x   o  */

                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * Kokkos::fabs(detDF); /* Center: beta_{i,j} */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = center_index;
                column = bottom_index;
                value  = -coeff3 * att; /* Bottom */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = center_index;
                column = top_index;
                value  = -coeff4 * att; /* Top */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = center_index;
                column = center_index;
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                /* Fill matrix row of (i,j-1) */
                row    = bottom_index;
                column = center_index;
                value  = -coeff3 * att; /* Top */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = bottom_index;
                column = bottom_index;
                value  = coeff3 * att; /* Center: (Top) */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                /* Fill matrix row of (i,j+1) */
                row    = top_index;
                column = center_index;
                value  = -coeff4 * att; /* Bottom */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = top_index;
                column = top_index;
                value  = coeff4 * att; /* Center: (Bottom) */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                /* Fill matrix row of (i-1,j) */
                row    = left_index;
                column = left_index;
                value  = coeff1 * arr; /* Center: (Right) */
                updateMatrixElement(left_solver, left_batch, row, column, value);

                /* Fill matrix row of (i+1,j) */
                row    = right_index;
                column = right_index;
                value  = coeff2 * arr; /* Center: (Left) */
                updateMatrixElement(right_solver, right_batch, row, column, value);
            }
            else {
                /* i_r % 2 == 1 and i_theta % 2 == 0 */
                /* | o | o | o || o   o   o   o  */
                /* |   |   |   || -------------- */
                /* | o | x | O || x   o   x   o  */
                /* |   |   |   || -------------- */
                /* | o | o | o || o   o   o   o  */

                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * Kokkos::fabs(detDF); /* Center: beta_{i,j} */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = center_index;
                column = bottom_index;
                value  = -coeff3 * att; /* Bottom */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = center_index;
                column = top_index;
                value  = -coeff4 * att; /* Top */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = center_index;
                column = center_index;
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                /* Fill matrix row of (i,j-1) */
                row    = bottom_index;
                column = center_index;
                value  = -coeff3 * att; /* Top */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = bottom_index;
                column = bottom_index;
                value  = coeff3 * att; /* Center: (Top) */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                /* Fill matrix row of (i,j+1) */
                row    = top_index;
                column = center_index;
                value  = -coeff4 * att; /* Bottom */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = top_index;
                column = top_index;
                value  = coeff4 * att; /* Center: (Bottom) */
                updateMatrixElement(center_solver, center_batch, row, column, value);
            }
        }
        else {
            if (i_theta & 1) {
                /* i_r % 2 == 0 and i_theta % 2 == 1 */
                /* | x | o | x || o   x   o   x  */
                /* |   |   |   || -------------- */
                /* | o | o | O || o   o   o   o  */
                /* |   |   |   || -------------- */
                /* | x | o | x || o   x   o   x  */

                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * Kokkos::fabs(detDF); /* Center: beta_{i,j} */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = center_index;
                column = center_index;
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                /* Fill matrix row of (i-1,j) */
                row    = left_index;
                column = left_index;
                value  = coeff1 * arr; /* Center: (Right) */
                updateMatrixElement(left_solver, left_batch, row, column, value);

                /* Fill matrix row of (i+1,j) */
                row    = right_index;
                column = right_index;
                value  = coeff2 * arr; /* Center: (Left) */
                updateMatrixElement(right_solver, right_batch, row, column, value);
            }
            else {
                /* i_r % 2 == 0 and i_theta % 2 == 0 */
                /* | o | o | o || o   o   o   o  */
                /* |   |   |   || -------------- */
                /* | x | o | X || o   x   o   x  */
                /* |   |   |   || -------------- */
                /* | o | o | o || o   o   o   o  */

                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 1.0;
                updateMatrixElement(center_solver, center_batch, row, column, value);

                /* Fill matrix row of (i,j-1) */
                row    = bottom_index;
                column = bottom_index;
                value  = coeff3 * att; /* Center: (Top) */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                /* Fill matrix row of (i,j+1) */
                row    = top_index;
                column = top_index;
                value  = coeff4 * att; /* Center: (Bottom) */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                /* Fill matrix row of (i-1,j) */
                row    = left_index;
                column = left_index;
                value  = coeff1 * arr; /* Center: (Right) */
                updateMatrixElement(left_solver, left_batch, row, column, value);

                /* Fill matrix row of (i+1,j) */
                row    = right_index;
                column = right_index;
                value  = coeff2 * arr; /* Center: (Left) */
                updateMatrixElement(right_solver, right_batch, row, column, value);
            }
        }
    }
}

template <class LevelCacheType>
static KOKKOS_INLINE_FUNCTION void nodeBuildTridiagonalSolverMatricesRadialSection(
    const int i_r, const int i_theta, const PolarGrid<DefaultMemorySpace>& grid, const LevelCacheType& level_cache,
    const bool DirBC_Interior, const BatchedTridiagonalSolver<double>& circle_tridiagonal_solver,
    const BatchedTridiagonalSolver<double>& radial_tridiagonal_solver)
{
    using extrapolated_smoother_give::updateMatrixElement;

    KOKKOS_ASSERT(i_r >= 0 && i_r < grid.nr());
    KOKKOS_ASSERT(i_theta >= 0 && i_theta < grid.ntheta());

    const int numberSmootherCircles = grid.numberSmootherCircles();
    const int lengthRadialSmoother  = grid.lengthRadialSmoother();

    KOKKOS_ASSERT(numberSmootherCircles >= 3);
    KOKKOS_ASSERT(lengthRadialSmoother >= 3);

    /* ---------------------------------------- */
    /* Compute or retrieve stencil coefficients */
    /* ---------------------------------------- */
    const int center    = grid.index(i_r, i_theta);
    const double radius = grid.radius(i_r);
    const double theta  = grid.theta(i_theta);

    double coeff_beta, arr, att, art, detDF;
    level_cache.obtainValues(i_r, i_theta, center, radius, theta, coeff_beta, arr, att, art, detDF);

    int row, column;
    double value;

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

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int center_index = i_r - numberSmootherCircles;
        const int left_index   = i_r - numberSmootherCircles - 1;
        const int right_index  = i_r - numberSmootherCircles + 1;
        const int bottom_index = i_r - numberSmootherCircles;
        const int top_index    = i_r - numberSmootherCircles;

        const auto& center_solver = radial_tridiagonal_solver;
        const int center_batch    = i_theta;
        const auto& bottom_solver = radial_tridiagonal_solver;
        const int bottom_batch    = i_theta_M1;
        const auto& top_solver    = radial_tridiagonal_solver;
        const int top_batch       = i_theta_P1;

        /* ------------------- */
        /* Tridiagonal Section */
        /* i_theta % 2 == 1    */
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

            /* Fill matrix row of (i,j) */
            row    = center_index;
            column = center_index;
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * Kokkos::fabs(detDF); /* Center: beta_{i,j} */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            row    = center_index;
            column = left_index;
            value  = -coeff1 * arr; /* Left */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            row    = center_index;
            column = right_index;
            value  = -coeff2 * arr; /* Right */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            row    = center_index;
            column = center_index;
            value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            /* Fill matrix row of (i-1,j) */
            row    = left_index;
            column = center_index;
            value  = -coeff1 * arr; /* Right */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            row    = left_index;
            column = left_index;
            value  = coeff1 * arr; /* Center: (Right) */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            /* Fill matrix row of (i+1,j) */
            row    = right_index;
            column = center_index;
            value  = -coeff2 * arr; /* Left */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            row    = right_index;
            column = right_index;
            value  = coeff2 * arr; /* Center: (Left) */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            if (i_r & 1) { /* i_r % 2 == 1 */
                /* ---------- */
                /* x   o   x  */
                /* ---------- */
                /* o   O   o  */
                /* ---------- */
                /* x   o   x  */
                /* ---------- */
                /* Fill matrix row of (i,j-1) */
                row    = bottom_index;
                column = bottom_index;
                value  = coeff3 * att; /* Center: (Top) */
                updateMatrixElement(bottom_solver, bottom_batch, row, column, value);

                /* Fill matrix row of (i,j+1) */
                row    = top_index;
                column = top_index;
                value  = coeff4 * att; /* Center: (Bottom) */
                updateMatrixElement(top_solver, top_batch, row, column, value);
            }
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

            if (i_r & 1) { /* i_r % 2 == 1 */
                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * Kokkos::fabs(detDF); /* Center: beta_{i,j} */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = center_index;
                column = center_index;
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
                updateMatrixElement(center_solver, center_batch, row, column, value);
            }
            else { /* i_r % 2 == 0 */
                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 1.0;
                updateMatrixElement(center_solver, center_batch, row, column, value);

                /* Fill matrix row of (i-1,j) */
                row    = left_index;
                column = left_index;
                value  = coeff1 * arr; /* Center: (Right) */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                /* Fill matrix row of (i+1,j) */
                row    = right_index;
                column = right_index;
                value  = coeff2 * arr; /* Center: (Left) */
                updateMatrixElement(center_solver, center_batch, row, column, value);
            }
            /* Fill matrix row of (i,j-1) */
            row    = bottom_index;
            column = bottom_index;
            value  = coeff3 * att; /* Center: (Top) */
            updateMatrixElement(bottom_solver, bottom_batch, row, column, value);

            /* Fill matrix row of (i,j+1) */
            row    = top_index;
            column = top_index;
            value  = coeff4 * att; /* Center: (Bottom) */
            updateMatrixElement(top_solver, top_batch, row, column, value);
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

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int center_index = i_r - numberSmootherCircles;
        const int left_index   = i_theta;
        const int right_index  = i_r - numberSmootherCircles + 1;
        const int bottom_index = i_r - numberSmootherCircles;
        const int top_index    = i_r - numberSmootherCircles;

        const auto& center_solver = radial_tridiagonal_solver;
        const int center_batch    = i_theta;
        const auto& bottom_solver = radial_tridiagonal_solver;
        const int bottom_batch    = i_theta_M1;
        const auto& top_solver    = radial_tridiagonal_solver;
        const int top_batch       = i_theta_P1;
        const auto& left_solver   = circle_tridiagonal_solver;
        const int left_batch      = i_r - 1;

        if (i_theta & 1) {
            if (i_r & 1) {
                /* i_theta % 2 == 1 and i_r % 2 == 1 */
                /* | x | o | x || o   x   o   x  */
                /* |   |   |   || -------------- */
                /* | o | o | o || O   o   o   o  */
                /* |   |   |   || -------------- */
                /* | x | o | x || o   x   o   x  */

                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * Kokkos::fabs(detDF); /* Center: beta_{i,j} */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = center_index;
                column = right_index;
                value  = -coeff2 * arr; /* Right */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = center_index;
                column = center_index;
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                /* Fill matrix row of (i-1,j) */
                row    = left_index;
                column = left_index;
                value  = coeff1 * arr; /* Center: (Right) */
                updateMatrixElement(left_solver, left_batch, row, column, value);

                /* Fill matrix row of (i+1,j) */
                row    = right_index;
                column = center_index;
                value  = -coeff2 * arr; /* Left */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = right_index;
                column = right_index;
                value  = coeff2 * arr; /* Center: (Left) */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                /* Fill matrix row of (i,j-1) */
                row    = bottom_index;
                column = bottom_index;
                value  = coeff3 * att; /* Center: (Top) */
                updateMatrixElement(bottom_solver, bottom_batch, row, column, value);

                /* Fill matrix row of (i,j+1) */
                row    = top_index;
                column = top_index;
                value  = coeff4 * att; /* Center: (Bottom) */
                updateMatrixElement(top_solver, top_batch, row, column, value);
            }
            else {
                /* i_theta % 2 == 1 and i_r % 2 == 0 */
                /* | o | x | o || x   o   x   o  */
                /* |   |   |   || -------------- */
                /* | o | o | o || O   o   o   o  */
                /* |   |   |   || -------------- */
                /* | o | x | o || x   o   x   o  */

                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * Kokkos::fabs(detDF); /* Center: beta_{i,j} */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = center_index;
                column = right_index;
                value  = -coeff2 * arr; /* Right */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = center_index;
                column = center_index;
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                /* Fill matrix row of (i-1,j) */
                row    = left_index;
                column = left_index;
                value  = coeff1 * arr; /* Center: (Right) */
                updateMatrixElement(left_solver, left_batch, row, column, value);

                /* Fill matrix row of (i+1,j) */
                row    = right_index;
                column = center_index;
                value  = -coeff2 * arr; /* Left */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = right_index;
                column = right_index;
                value  = coeff2 * arr; /* Center: (Left) */
                updateMatrixElement(center_solver, center_batch, row, column, value);
            }
        }
        else {
            if (i_r & 1) {
                /* i_theta % 2 == 0 and i_r % 2 == 1 */
                /* | o | o | o || o   o   o   o  */
                /* |   |   |   || -------------- */
                /* | x | o | x || O   x   o   x  */
                /* |   |   |   || -------------- */
                /* | o | o | o || o   o   o   o  */

                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * Kokkos::fabs(detDF); /* Center: beta_{i,j} */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                row    = center_index;
                column = center_index;
                value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                /* Fill matrix row of (i,j-1) */
                row    = bottom_index;
                column = bottom_index;
                value  = coeff3 * att; /* Center: (Top) */
                updateMatrixElement(bottom_solver, bottom_batch, row, column, value);

                /* Fill matrix row of (i,j+1) */
                row    = top_index;
                column = top_index;
                value  = coeff4 * att; /* Center: (Bottom) */
                updateMatrixElement(top_solver, top_batch, row, column, value);
            }
            else {
                /* i_theta % 2 == 0 and i_r % 2 == 0 */
                /* | o | o | o || o   o   o   o  */
                /* |   |   |   || -------------- */
                /* | o | x | o || X   o   x   o  */
                /* |   |   |   || -------------- */
                /* | o | o | o || o   o   o   o  */

                /* Fill matrix row of (i,j) */
                row    = center_index;
                column = center_index;
                value  = 1.0;
                updateMatrixElement(center_solver, center_batch, row, column, value);

                /* Fill matrix row of (i-1,j) */
                row    = left_index;
                column = left_index;
                value  = coeff1 * arr; /* Center: (Right) */
                updateMatrixElement(left_solver, left_batch, row, column, value);

                /* Fill matrix row of (i+1,j) */
                row    = right_index;
                column = right_index;
                value  = coeff2 * arr; /* Center: (Left) */
                updateMatrixElement(center_solver, center_batch, row, column, value);

                /* Fill matrix row of (i,j-1) */
                row    = bottom_index;
                column = bottom_index;
                value  = coeff3 * att; /* Center: (Top) */
                updateMatrixElement(bottom_solver, bottom_batch, row, column, value);

                /* Fill matrix row of (i,j+1) */
                row    = top_index;
                column = top_index;
                value  = coeff4 * att; /* Center: (Bottom) */
                updateMatrixElement(top_solver, top_batch, row, column, value);
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

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        const int center_index = i_r - numberSmootherCircles;
        const int left_index   = i_r - numberSmootherCircles - 1;
        //const int right_index  = i_r - numberSmootherCircles + 1;
        const int bottom_index = i_r - numberSmootherCircles;
        const int top_index    = i_r - numberSmootherCircles;

        const auto& center_solver = radial_tridiagonal_solver;
        const int center_batch    = i_theta;
        const auto& bottom_solver = radial_tridiagonal_solver;
        const int bottom_batch    = i_theta_M1;
        const auto& top_solver    = radial_tridiagonal_solver;
        const int top_batch       = i_theta_P1;

        if (i_theta & 1) {
            /* i_theta % 2 == 1 */
            /* ---------------|| */
            /* o   x   o   x  || */
            /* ---------------|| */
            /* o   o   O   o  || */
            /* ---------------|| */
            /* o   x   o   x  || */
            /* ---------------|| */

            /* Fill matrix row of (i,j) */
            row    = center_index;
            column = center_index;
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * Kokkos::fabs(detDF); /* Center: beta_{i,j} */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            row    = center_index;
            column = left_index;
            value  = -coeff1 * arr; /* Left */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            /* Remark: Right is not included here due to the symmetry shift */

            row    = center_index;
            column = center_index;
            value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            /* Fill matrix row of (i-1,j) */
            row    = left_index;
            column = center_index;
            value  = -coeff1 * arr; /* Right */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            row    = left_index;
            column = left_index;
            value  = coeff1 * arr; /* Center: (Right) */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            /* Fill matrix row of (i+1,j) */
            /* Nothing to be done */

            /* Fill matrix row of (i,j-1) */
            row    = bottom_index;
            column = bottom_index;
            value  = coeff3 * att; /* Center: (Top) */
            updateMatrixElement(bottom_solver, bottom_batch, row, column, value);

            /* Fill matrix row of (i,j+1) */
            row    = top_index;
            column = top_index;
            value  = coeff4 * att; /* Center: (Bottom) */
            updateMatrixElement(top_solver, top_batch, row, column, value);
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

            /* Fill matrix row of (i,j) */
            row    = center_index;
            column = center_index;
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * Kokkos::fabs(detDF); /* Center: beta_{i,j} */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            row    = center_index;
            column = center_index;
            value  = (coeff1 + coeff2) * arr + (coeff3 + coeff4) * att; /* Center: (Left, Right, Bottom, Top) */
            updateMatrixElement(center_solver, center_batch, row, column, value);

            /* Fill matrix row of (i,j-1) */
            row    = bottom_index;
            column = bottom_index;
            value  = coeff3 * att; /* Center: (Top) */
            updateMatrixElement(bottom_solver, bottom_batch, row, column, value);

            /* Fill matrix row of (i,j+1) */
            row    = top_index;
            column = top_index;
            value  = coeff4 * att; /* Center: (Bottom) */
            updateMatrixElement(top_solver, top_batch, row, column, value);
        }
    }
    /* ------------------------------------------ */
    /* Radial Section: Node on the outer boundary */
    /* ------------------------------------------ */
    else if (i_r == grid.nr() - 1) {
        KOKKOS_ASSERT(i_r % 2 == 0);

        const double h1 = grid.radialSpacing(i_r - 1);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff1 = 0.5 * (k1 + k2) / h1;

        const int center_index = i_r - numberSmootherCircles;
        const int left_index   = i_r - numberSmootherCircles - 1;

        const auto& center_solver = radial_tridiagonal_solver;
        const int center_batch    = i_theta;

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
            updateMatrixElement(center_solver, center_batch, row, column, value);

            /* Fill matrix row of (i-1,j) */
            row    = left_index;
            column = left_index;
            value  = coeff1 * arr; /* Center: (Right) */
            updateMatrixElement(center_solver, center_batch, row, column, value);
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
            updateMatrixElement(center_solver, center_batch, row, column, value);

            /* Fill matrix row of (i-1,j) */
            row    = left_index;
            column = left_index;
            value  = coeff1 * arr; /* Center: (Right) */
            updateMatrixElement(center_solver, center_batch, row, column, value);
        }
    }
}

} // namespace extrapolated_smoother_give

template <class LevelCacheType>
void ExtrapolatedSmootherGive<LevelCacheType>::buildTridiagonalSolverMatrices()
{
    using extrapolated_smoother_give::nodeBuildTridiagonalSolverMatricesCircleSection;
    using extrapolated_smoother_give::nodeBuildTridiagonalSolverMatricesRadialSection;

    const PolarGrid<DefaultMemorySpace>& grid = ExtrapolatedSmoother<LevelCacheType>::grid_;
    const LevelCacheType& level_cache         = ExtrapolatedSmoother<LevelCacheType>::level_cache_;
    const bool DirBC_Interior                 = ExtrapolatedSmoother<LevelCacheType>::DirBC_Interior_;

    const BatchedTridiagonalSolver<double>& circle_tridiagonal_solver = circle_tridiagonal_solver_;
    const BatchedTridiagonalSolver<double>& radial_tridiagonal_solver = radial_tridiagonal_solver_;

    /* ---------------- */
    /* Circular section */
    /* ---------------- */
    // We parallelize over i_r (step 3) to avoid data race conditions between adjacent circles.
    // The i_theta loop is sequential inside the kernel.
    const int num_circle_tasks = grid.numberSmootherCircles();

    for (int start_circle = 0; start_circle < 3; ++start_circle) {
        const int num_circular_tasks = (num_circle_tasks - start_circle + 2) / 3;
        Kokkos::parallel_for(
            "SmootherGive: buildTridiagonalSolverMatrices (Circular)",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, num_circular_tasks),
            KOKKOS_LAMBDA(const int circle_task) {
                const int i_r = start_circle + circle_task * 3;
                for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                    nodeBuildTridiagonalSolverMatricesCircleSection(i_r, i_theta, grid, level_cache, DirBC_Interior,
                                                                    circle_tridiagonal_solver,
                                                                    radial_tridiagonal_solver);
                }
            });
        Kokkos::fence();
    }

    /* -------------- */
    /* Radial section */
    /* -------------- */
    // We parallelize over i_theta (step 3) to avoid data race conditions between adjacent radial lines.
    // The i_r loop is sequential inside the kernel.
    // Due to periodicity in the angular direction, handle up to 2 additional
    // radial lines (i_theta = 0 and 1) before the parallel passes.
    const int additional_radial_tasks = grid.ntheta() % 3;
    const int num_radial_tasks        = grid.ntheta() - additional_radial_tasks;

    for (int i_theta = 0; i_theta < additional_radial_tasks; i_theta++) {
        Kokkos::parallel_for(
            "SmootherGive: buildTridiagonalSolverMatrices (Radial, additional)",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, 1), KOKKOS_LAMBDA(const int) {
                for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                    nodeBuildTridiagonalSolverMatricesRadialSection(i_r, i_theta, grid, level_cache, DirBC_Interior,
                                                                    circle_tridiagonal_solver,
                                                                    radial_tridiagonal_solver);
                }
            });
        Kokkos::fence();
    }

    for (int start_radial = 0; start_radial < 3; ++start_radial) {
        const int num_radial_batches = (num_radial_tasks - start_radial + 2) / 3;
        Kokkos::parallel_for(
            "SmootherGive: buildTridiagonalSolverMatrices (Radial)",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, num_radial_batches),
            KOKKOS_LAMBDA(const int radial_task) {
                const int i_theta = additional_radial_tasks + start_radial + radial_task * 3;
                for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                    nodeBuildTridiagonalSolverMatricesRadialSection(i_r, i_theta, grid, level_cache, DirBC_Interior,
                                                                    circle_tridiagonal_solver,
                                                                    radial_tridiagonal_solver);
                }
            });
        Kokkos::fence();
    }
}
