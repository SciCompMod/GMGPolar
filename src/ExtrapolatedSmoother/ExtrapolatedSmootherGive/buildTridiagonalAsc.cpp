#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherGive/extrapolatedSmootherGive.h"

#include "../../../include/Definitions/geometry_helper.h"

/* Tridiagonal matrices */
static inline void updateMatrixElement(BatchedTridiagonalSolver<double>& solver, int batch, int row, int column,
                                       double value)
{
    if (row == column)
        solver.main_diagonal(batch, row) += value;
    else if (row == column - 1)
        solver.sub_diagonal(batch, row) += value;
    else if (row == 0 && column == solver.matrixDimension() - 1)
        solver.cyclic_corner(batch) += value;
}

void ExtrapolatedSmootherGive::nodeBuildTridiagonalSolverMatrices(
    int i_r, int i_theta, const PolarGrid& grid, bool DirBC_Interior,
    BatchedTridiagonalSolver<double>& circle_tridiagonal_solver,
    BatchedTridiagonalSolver<double>& radial_tridiagonal_solver, double arr, double att, double art, double detDF,
    double coeff_beta)
{
    assert(i_r >= 0 && i_r < grid.nr());
    assert(i_theta >= 0 && i_theta < grid.ntheta());

    const int numberSmootherCircles = grid.numberSmootherCircles();
    const int lengthSmootherRadial  = grid.lengthSmootherRadial();

    assert(numberSmootherCircles >= 3);
    assert(lengthSmootherRadial >= 3);

    int ptr, offset;
    int row, column;
    double value;
    /* ------------------------------------------ */
    /* Node in the interior of the Circle Section */
    /* ------------------------------------------ */
    if (i_r > 0 && i_r < numberSmootherCircles - 1) {
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

        int center_index = i_theta;
        int left_index   = i_theta;
        int right_index  = i_theta;
        int bottom_index = i_theta_M1;
        int top_index    = i_theta_P1;

        auto& center_solver = circle_tridiagonal_solver;
        int center_batch    = i_r;
        auto& left_solver   = circle_tridiagonal_solver;
        int left_batch      = i_r - 1;
        auto& right_solver  = circle_tridiagonal_solver;
        int right_batch     = i_r + 1;

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
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
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
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
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

        int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        int center_index = i_r - numberSmootherCircles;
        int left_index   = i_r - numberSmootherCircles - 1;
        int right_index  = i_r - numberSmootherCircles + 1;
        int bottom_index = i_r - numberSmootherCircles;
        int top_index    = i_r - numberSmootherCircles;

        auto& center_solver = radial_tridiagonal_solver;
        int center_batch    = i_theta;
        auto& bottom_solver = radial_tridiagonal_solver;
        int bottom_batch    = i_theta_M1;
        auto& top_solver    = radial_tridiagonal_solver;
        int top_batch       = i_theta_P1;

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
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
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
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
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
    /* ------------------------------------------ */
    /* Circle Section: Node in the inner boundary */
    /* ------------------------------------------ */
    else if (i_r == 0) {
        // The inner boundary circle line are is handled by the inner_boundary_mumps_solver, so we fill in the identity matrix.
        auto& center_solver = circle_tridiagonal_solver;
        int center_batch    = i_r;
        auto& right_solver  = circle_tridiagonal_solver;
        int right_batch     = i_r + 1;

        /* Fill result(i,j) */
        double h2 = grid.radialSpacing(i_r);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff2 = 0.5 * (k1 + k2) / h2;

        int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        int center_index = i_theta;
        int right_index  = i_theta;
        int bottom_index = i_theta_M1;
        int top_index    = i_theta_P1;

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
        assert(i_r > 1);

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

        int center_index = i_theta;
        int left_index   = i_theta;
        int right_index  = 0;
        int bottom_index = i_theta_M1;
        int top_index    = i_theta_P1;

        auto& center_solver = circle_tridiagonal_solver;
        int center_batch    = i_r;
        auto& left_solver   = circle_tridiagonal_solver;
        int left_batch      = i_r - 1;
        auto& right_solver  = radial_tridiagonal_solver;
        int right_batch     = i_theta;

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
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
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
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
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
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
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

        const int center_index = i_r - numberSmootherCircles;
        const int left_index   = i_theta;
        const int right_index  = i_r - numberSmootherCircles + 1;
        const int bottom_index = i_r - numberSmootherCircles;
        const int top_index    = i_r - numberSmootherCircles;

        auto& center_solver = radial_tridiagonal_solver;
        int center_batch    = i_theta;
        auto& bottom_solver = radial_tridiagonal_solver;
        int bottom_batch    = i_theta_M1;
        auto& top_solver    = radial_tridiagonal_solver;
        int top_batch       = i_theta_P1;
        auto& left_solver   = circle_tridiagonal_solver;
        int left_batch      = i_r - 1;

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
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
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
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
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
                value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
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
        assert(i_r % 2 == 1);

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

        int center_index = i_r - numberSmootherCircles;
        int left_index   = i_r - numberSmootherCircles - 1;
        int right_index  = i_r - numberSmootherCircles + 1;
        int bottom_index = i_r - numberSmootherCircles;
        int top_index    = i_r - numberSmootherCircles;

        auto& center_solver = radial_tridiagonal_solver;
        int center_batch    = i_theta;
        auto& bottom_solver = radial_tridiagonal_solver;
        int bottom_batch    = i_theta_M1;
        auto& top_solver    = radial_tridiagonal_solver;
        int top_batch       = i_theta_P1;

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
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
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
            value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
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
        assert(!i_r % 2 == 0);

        double h1 = grid.radialSpacing(i_r - 1);
        double k1 = grid.angularSpacing(i_theta - 1);
        double k2 = grid.angularSpacing(i_theta);

        double coeff1 = 0.5 * (k1 + k2) / h1;

        int center_index = i_r - numberSmootherCircles;
        int left_index   = i_r - numberSmootherCircles - 1;

        auto& center_solver = radial_tridiagonal_solver;
        int center_batch    = i_theta;

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

void ExtrapolatedSmootherGive::buildTridiagonalCircleSection(int i_r)
{
    const double r = grid_.radius(i_r);
    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        const int global_index = grid_.index(i_r, i_theta);
        const double theta     = grid_.theta(i_theta);

        double coeff_beta, arr, att, art, detDF;
        level_cache_.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art, detDF);

        // Build Asc at the current node
        nodeBuildTridiagonalSolverMatrices(i_r, i_theta, grid_, DirBC_Interior_, circle_tridiagonal_solver_,
                                           radial_tridiagonal_solver_, arr, att, art, detDF, coeff_beta);
    }
}

void ExtrapolatedSmootherGive::buildTridiagonalRadialSection(int i_theta)
{
    const double theta = grid_.theta(i_theta);
    for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
        const int global_index = grid_.index(i_r, i_theta);
        const double r         = grid_.radius(i_r);

        double coeff_beta, arr, att, art, detDF;
        level_cache_.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art, detDF);

        // Build Asc at the current node
        nodeBuildTridiagonalSolverMatrices(i_r, i_theta, grid_, DirBC_Interior_, circle_tridiagonal_solver_,
                                           radial_tridiagonal_solver_, arr, att, art, detDF, coeff_beta);
    }
}

void ExtrapolatedSmootherGive::buildTridiagonalSolverMatrices()
{
    /*  Multi-threaded execution: */
    const int num_smoother_circles    = grid_.numberSmootherCircles();
    const int additional_radial_tasks = grid_.ntheta() % 3;
    const int num_radial_tasks        = grid_.ntheta() - additional_radial_tasks;

#pragma omp parallel num_threads(num_omp_threads_)
    {
#pragma omp for
        for (int i_r = 0; i_r < num_smoother_circles; i_r += 3) {
            buildTridiagonalCircleSection(i_r);
        }
#pragma omp for
        for (int i_r = 1; i_r < num_smoother_circles; i_r += 3) {
            buildTridiagonalCircleSection(i_r);
        }
#pragma omp for
        for (int i_r = 2; i_r < num_smoother_circles; i_r += 3) {
            buildTridiagonalCircleSection(i_r);
        }

#pragma omp for
        for (int radial_task = 0; radial_task < num_radial_tasks; radial_task += 3) {
            if (radial_task > 0) {
                int i_theta = radial_task + additional_radial_tasks;
                buildTridiagonalRadialSection(i_theta);
            }
            else {
                if (additional_radial_tasks == 0) {
                    buildTridiagonalRadialSection(0);
                }
                else if (additional_radial_tasks >= 1) {
                    buildTridiagonalRadialSection(0);
                    buildTridiagonalRadialSection(1);
                }
            }
        }
#pragma omp for
        for (int radial_task = 1; radial_task < num_radial_tasks; radial_task += 3) {
            if (radial_task > 1) {
                int i_theta = radial_task + additional_radial_tasks;
                buildTridiagonalRadialSection(i_theta);
            }
            else {
                if (additional_radial_tasks == 0) {
                    buildTridiagonalRadialSection(1);
                }
                else if (additional_radial_tasks == 1) {
                    buildTridiagonalRadialSection(2);
                }
                else if (additional_radial_tasks == 2) {
                    buildTridiagonalRadialSection(2);
                    buildTridiagonalRadialSection(3);
                }
            }
        }
#pragma omp for
        for (int radial_task = 2; radial_task < num_radial_tasks; radial_task += 3) {
            int i_theta = radial_task + additional_radial_tasks;
            buildTridiagonalRadialSection(i_theta);
        }
    }
}
