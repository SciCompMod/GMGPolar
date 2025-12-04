#include "../../../include/DirectSolver/DirectSolver-COO-MUMPS-Give/directSolverGive.h"

#include "../../../include/common/geometry_helper.h"

#ifdef GMGPOLAR_USE_MUMPS

/* ----------------------- */
/* Boundary Symmetry Shift */
/* ----------------------- */

void DirectSolverGive::applySymmetryShiftInnerBoundary(Vector<double> x) const
{
    assert(DirBC_Interior_);

    int i_r;
    double r;
    int global_index;
    double h1, h2, k1, k2;
    double coeff1, coeff2;
    double coeff_beta, arr, att, art, detDF;

    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        const double theta = grid_.theta(i_theta);
        /* -------------------------- */
        /* Node on the inner boundary */
        /* -------------------------- */
        i_r          = 0;
        r            = grid_.radius(i_r);
        global_index = grid_.index(i_r, i_theta);

        level_cache_.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art,
                                  detDF);

        h2 = grid_.radialSpacing(i_r);
        k1 = grid_.angularSpacing(i_theta - 1);
        k2 = grid_.angularSpacing(i_theta);

        coeff2 = 0.5 * (k1 + k2) / h2;

        /* Fill x(i+1,j) */
        x(grid_.index(i_r + 1, i_theta)) -= -coeff2 * arr * x(grid_.index(i_r, i_theta)) /* Left */
                                            + 0.25 * art * x(grid_.index(i_r, i_theta + 1)) /* Top Left */
                                            - 0.25 * art * x(grid_.index(i_r, i_theta - 1)); /* Bottom Left */

        /* --------------------------- */
        /* Node next to inner boundary */
        /* --------------------------- */
        i_r          = 1;
        r            = grid_.radius(i_r);
        global_index = grid_.index(i_r, i_theta);

        level_cache_.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art,
                                  detDF);

        h1 = grid_.radialSpacing(i_r - 1);
        k1 = grid_.angularSpacing(i_theta - 1);
        k2 = grid_.angularSpacing(i_theta);

        coeff1 = 0.5 * (k1 + k2) / h1;

        /* Fill x(i,j) */
        x(grid_.index(i_r, i_theta)) -= -coeff1 * arr * x(grid_.index(i_r - 1, i_theta)); /* Left */
        /* Fill x(i,j-1) */
        x(grid_.index(i_r, i_theta - 1)) -= +0.25 * art * x(grid_.index(i_r - 1, i_theta)); /* Top Left */
        /* Fill x(i,j+1) */
        x(grid_.index(i_r, i_theta + 1)) -= -0.25 * art * x(grid_.index(i_r - 1, i_theta)); /* Bottom Left */
    }
}

void DirectSolverGive::applySymmetryShiftOuterBoundary(Vector<double> x) const
{
    int i_r;
    double r;
    int global_index;
    double h1, h2, k1, k2;
    double coeff1, coeff2;
    double coeff_beta, arr, att, art, detDF;

    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        const double theta = grid_.theta(i_theta);
        /* --------------------------- */
        /* Node next to outer boundary */
        /* --------------------------- */
        i_r          = grid_.nr() - 2;
        r            = grid_.radius(i_r);
        global_index = grid_.index(i_r, i_theta);

        level_cache_.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art,
                                  detDF);

        h2 = grid_.radialSpacing(i_r);
        k1 = grid_.angularSpacing(i_theta - 1);
        k2 = grid_.angularSpacing(i_theta);

        coeff2 = 0.5 * (k1 + k2) / h2;

        /* Fill result(i,j) */
        x(grid_.index(i_r, i_theta)) -= -coeff2 * arr * x(grid_.index(i_r + 1, i_theta)); /* Right */
        /* Fill result(i,j-1) */
        x(grid_.index(i_r, i_theta - 1)) -= -0.25 * art * x(grid_.index(i_r + 1, i_theta)); /* Top Right */
        /* Fill result(i,j+1) */
        x(grid_.index(i_r, i_theta + 1)) -= +0.25 * art * x(grid_.index(i_r + 1, i_theta)); /* Bottom Right */

        /* -------------------------- */
        /* Node on the outer boundary */
        /* -------------------------- */
        i_r          = grid_.nr() - 1;
        r            = grid_.radius(i_r);
        global_index = grid_.index(i_r, i_theta);

        level_cache_.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art,
                                  detDF);

        h1 = grid_.radialSpacing(i_r - 1);
        k1 = grid_.angularSpacing(i_theta - 1);
        k2 = grid_.angularSpacing(i_theta);

        coeff1 = 0.5 * (k1 + k2) / h1;

        /* Fill result(i-1,j) */
        x(grid_.index(i_r - 1, i_theta)) -= -coeff1 * arr * x(grid_.index(i_r, i_theta)) /* Right */
                                            - 0.25 * art * x(grid_.index(i_r, i_theta + 1)) /* Top Right */
                                            + 0.25 * art * x(grid_.index(i_r, i_theta - 1)); /* Bottom Right */
    }
}

// clang-format off
void DirectSolverGive::applySymmetryShift(Vector<double> x) const
{
    assert(x.size() == grid_.numberOfNodes());
    assert(grid_.nr() >= 4);

    if (num_omp_threads_ == 1) {
        /* Single-threaded execution */
        if (DirBC_Interior_) {
            applySymmetryShiftInnerBoundary(x);
        }
        applySymmetryShiftOuterBoundary(x);
    }
    else {
        #pragma omp parallel sections num_threads(num_omp_threads_)
        {
            #pragma omp section
            {
                if (DirBC_Interior_) {
                    applySymmetryShiftInnerBoundary(x);
                }
            }

            #pragma omp section
            {
                applySymmetryShiftOuterBoundary(x);
            }
        }
    }
}
// clang-format on
#endif
