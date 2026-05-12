#pragma once

/* ----------------------- */
/* Boundary Symmetry Shift */
/* ----------------------- */

template <class LevelCacheType>
void DirectSolverGive<LevelCacheType>::applySymmetryShiftInnerBoundary(Vector<double> x) const
{
    const PolarGrid& grid             = DirectSolver<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = DirectSolver<LevelCacheType>::level_cache_;

    assert(DirectSolver<LevelCacheType>::DirBC_Interior_);

    int i_r;
    double r;
    int global_index;
    double h1, h2, k1, k2;
    double coeff1, coeff2;
    double coeff_beta, arr, att, art, detDF;

    for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
        const double theta = grid.theta(i_theta);
        /* -------------------------- */
        /* Node on the inner boundary */
        /* -------------------------- */
        i_r          = 0;
        r            = grid.radius(i_r);
        global_index = grid.index(i_r, i_theta);

        level_cache.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art, detDF);

        h2 = grid.radialSpacing(i_r);
        k1 = grid.angularSpacing(i_theta - 1);
        k2 = grid.angularSpacing(i_theta);

        coeff2 = 0.5 * (k1 + k2) / h2;

        /* Fill x(i+1,j) */
        x(grid.index(i_r + 1, i_theta)) -= -coeff2 * arr * x(grid.index(i_r, i_theta)) /* Left */
                                           + 0.25 * art * x(grid.index(i_r, i_theta + 1)) /* Top Left */
                                           - 0.25 * art * x(grid.index(i_r, i_theta - 1)); /* Bottom Left */

        /* --------------------------- */
        /* Node next to inner boundary */
        /* --------------------------- */
        i_r          = 1;
        r            = grid.radius(i_r);
        global_index = grid.index(i_r, i_theta);

        level_cache.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art, detDF);

        h1 = grid.radialSpacing(i_r - 1);
        k1 = grid.angularSpacing(i_theta - 1);
        k2 = grid.angularSpacing(i_theta);

        coeff1 = 0.5 * (k1 + k2) / h1;

        /* Fill x(i,j) */
        x(grid.index(i_r, i_theta)) -= -coeff1 * arr * x(grid.index(i_r - 1, i_theta)); /* Left */
        /* Fill x(i,j-1) */
        x(grid.index(i_r, i_theta - 1)) -= +0.25 * art * x(grid.index(i_r - 1, i_theta)); /* Top Left */
        /* Fill x(i,j+1) */
        x(grid.index(i_r, i_theta + 1)) -= -0.25 * art * x(grid.index(i_r - 1, i_theta)); /* Bottom Left */
    }
}

template <class LevelCacheType>
void DirectSolverGive<LevelCacheType>::applySymmetryShiftOuterBoundary(Vector<double> x) const
{
    const PolarGrid& grid             = DirectSolver<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = DirectSolver<LevelCacheType>::level_cache_;

    int i_r;
    double r;
    int global_index;
    double h1, h2, k1, k2;
    double coeff1, coeff2;
    double coeff_beta, arr, att, art, detDF;

    for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
        const double theta = grid.theta(i_theta);
        /* --------------------------- */
        /* Node next to outer boundary */
        /* --------------------------- */
        i_r          = grid.nr() - 2;
        r            = grid.radius(i_r);
        global_index = grid.index(i_r, i_theta);

        level_cache.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art, detDF);

        h2 = grid.radialSpacing(i_r);
        k1 = grid.angularSpacing(i_theta - 1);
        k2 = grid.angularSpacing(i_theta);

        coeff2 = 0.5 * (k1 + k2) / h2;

        /* Fill result(i,j) */
        x(grid.index(i_r, i_theta)) -= -coeff2 * arr * x(grid.index(i_r + 1, i_theta)); /* Right */
        /* Fill result(i,j-1) */
        x(grid.index(i_r, i_theta - 1)) -= -0.25 * art * x(grid.index(i_r + 1, i_theta)); /* Top Right */
        /* Fill result(i,j+1) */
        x(grid.index(i_r, i_theta + 1)) -= +0.25 * art * x(grid.index(i_r + 1, i_theta)); /* Bottom Right */

        /* -------------------------- */
        /* Node on the outer boundary */
        /* -------------------------- */
        i_r          = grid.nr() - 1;
        r            = grid.radius(i_r);
        global_index = grid.index(i_r, i_theta);

        level_cache.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art, detDF);

        h1 = grid.radialSpacing(i_r - 1);
        k1 = grid.angularSpacing(i_theta - 1);
        k2 = grid.angularSpacing(i_theta);

        coeff1 = 0.5 * (k1 + k2) / h1;

        /* Fill result(i-1,j) */
        x(grid.index(i_r - 1, i_theta)) -= -coeff1 * arr * x(grid.index(i_r, i_theta)) /* Right */
                                           - 0.25 * art * x(grid.index(i_r, i_theta + 1)) /* Top Right */
                                           + 0.25 * art * x(grid.index(i_r, i_theta - 1)); /* Bottom Right */
    }
}

template <class LevelCacheType>
void DirectSolverGive<LevelCacheType>::applySymmetryShift(Vector<double> x) const
{
    const PolarGrid& grid     = DirectSolver<LevelCacheType>::grid_;
    const bool DirBC_Interior = DirectSolver<LevelCacheType>::DirBC_Interior_;

    assert(std::ssize(x) == grid.numberOfNodes());
    assert(grid.nr() >= 4);

    if (DirBC_Interior) {
        applySymmetryShiftInnerBoundary(x);
    }

    applySymmetryShiftOuterBoundary(x);
}
