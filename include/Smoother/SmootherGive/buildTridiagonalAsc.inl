#pragma once

namespace smoother_give
{

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

} // namespace smoother_give

template <class LevelCacheType>
void SmootherGive<LevelCacheType>::nodeBuildTridiagonalSolverMatrices(
    int i_r, int i_theta, const PolarGrid& grid, bool DirBC_Interior,
    BatchedTridiagonalSolver<double>& circle_tridiagonal_solver,
    BatchedTridiagonalSolver<double>& radial_tridiagonal_solver, double arr, double att, double art, double detDF,
    double coeff_beta)
{
    using smoother_give::updateMatrixElement;

    assert(i_r >= 0 && i_r < grid.nr());
    assert(i_theta >= 0 && i_theta < grid.ntheta());

    const int numberSmootherCircles = grid.numberSmootherCircles();
    const int lengthSmootherRadial  = grid.lengthSmootherRadial();

    assert(numberSmootherCircles >= 2);
    assert(lengthSmootherRadial >= 3);

    int ptr, offset;
    int row, column, col;
    double value, val;

    /* ------------------------------------------ */
    /* Node in the interior of the Circle Section */
    /* ------------------------------------------ */
    if (i_r > 0 && i_r < numberSmootherCircles) {
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
        const int right_index  = (i_r + 1 == numberSmootherCircles) ? 0 : i_theta;
        const int bottom_index = i_theta_M1;
        const int top_index    = i_theta_P1;

        /* Visualization of the sourrounding tridiagonal matrices. */
        /* left_matrix, center_matrix, right_matrix */
        /* | o | o | o | */
        /* |   |   |   | */
        /* | o | O | o | */
        /* |   |   |   | */
        /* | o | o | o | */
        /* or */
        /* left_matrix, right_matrix */
        /* | o | o | o || o   o   o   o  */
        /* |   |   |   || -------------- */
        /* | o | o | O || o   o   o   o  <- right_matrix */
        /* |   |   |   || -------------- */
        /* | o | o | o || o   o   o   o  */
        auto& left_solver   = circle_tridiagonal_solver;
        int left_batch      = i_r - 1;
        auto& center_solver = circle_tridiagonal_solver;
        int center_batch    = i_r;
        auto& right_solver = (i_r + 1 == numberSmootherCircles) ? radial_tridiagonal_solver : circle_tridiagonal_solver;
        int right_batch    = (i_r + 1 == numberSmootherCircles) ? i_theta : i_r + 1;

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
    /* ------------------------------------------ */
    /* Node in the interior of the Radial Section */
    /* ------------------------------------------ */
    else if (i_r > numberSmootherCircles && i_r < grid.nr() - 2) {
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

        /* ---------- */
        /* o   o   o  <- top_matrix */
        /* ---------- */
        /* o   O   o  <- center_matrix */
        /* ---------- */
        /* o   o   o  <- bottom_matrix */
        /* ---------- */
        auto& bottom_solver = radial_tridiagonal_solver;
        int bottom_batch    = i_theta_M1;
        auto& center_solver = radial_tridiagonal_solver;
        int center_batch    = i_theta;
        auto& top_solver    = radial_tridiagonal_solver;
        int top_batch       = i_theta_P1;

        const int center_index = i_r - numberSmootherCircles;
        const int left_index   = i_r - numberSmootherCircles - 1;
        const int right_index  = i_r - numberSmootherCircles + 1;
        const int bottom_index = i_r - numberSmootherCircles;
        const int top_index    = i_r - numberSmootherCircles;

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
    /* ------------------------------------------ */
    /* Circle Section: Node in the inner boundary */
    /* ------------------------------------------ */
    else if (i_r == 0) {
        // The inner boundary circle line are is handled by the inner_boundary_mumps_solver, so we fill in the identity matrix.

        /* Fill result(i,j) */
        const double h2 = grid.radialSpacing(i_r);
        const double k1 = grid.angularSpacing(i_theta - 1);
        const double k2 = grid.angularSpacing(i_theta);

        const double coeff2 = 0.5 * (k1 + k2) / h2;

        const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
        const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

        auto& center_solver = circle_tridiagonal_solver;
        int center_batch    = i_r;
        auto& right_solver  = circle_tridiagonal_solver;
        int right_batch     = i_r + 1;

        const int center_index = i_theta;
        const int right_index  = i_theta;
        const int bottom_index = i_theta_M1;
        const int top_index    = i_theta_P1;

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

        /* | o | o | o || o   o   o   o  <- top_matrix */
        /* |   |   |   || -------------- */
        /* | o | o | o || O   o   o   o  <- center_matrix */
        /* |   |   |   || -------------- */
        /* | o | o | o || o   o   o   o  <- bottom_matrix */
        auto& bottom_solver = radial_tridiagonal_solver;
        int bottom_batch    = i_theta_M1;
        auto& center_solver = radial_tridiagonal_solver;
        int center_batch    = i_theta;
        auto& top_solver    = radial_tridiagonal_solver;
        int top_batch       = i_theta_P1;
        auto& left_solver   = circle_tridiagonal_solver;
        int left_batch      = i_r - 1;

        const int center_index = i_r - numberSmootherCircles;
        const int left_index   = i_theta;
        const int right_index  = i_r - numberSmootherCircles + 1;
        const int bottom_index = i_r - numberSmootherCircles;
        const int top_index    = i_r - numberSmootherCircles;

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
    /* ------------------------------------------- */
    /* Radial Section: Node next to outer boundary */
    /* ------------------------------------------- */
    else if (i_r == grid.nr() - 2) {
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

        /* ---------------|| */
        /* o   o   o   o  || <- top_matrix */
        /* ---------------|| */
        /* o   o   O   o  || <- center_matrix */
        /* ---------------|| */
        /* o   o   o   o  || <- bottom_matrix */
        /* ---------------|| */
        auto& bottom_solver = radial_tridiagonal_solver;
        int bottom_batch    = i_theta_M1;
        auto& center_solver = radial_tridiagonal_solver;
        int center_batch    = i_theta;
        auto& top_solver    = radial_tridiagonal_solver;
        int top_batch       = i_theta_P1;

        const int center_index = i_r - numberSmootherCircles;
        const int left_index   = i_r - numberSmootherCircles - 1;
        const int right_index  = i_r - numberSmootherCircles + 1;
        const int bottom_index = i_r - numberSmootherCircles;
        const int top_index    = i_r - numberSmootherCircles;

        /* ---------------------------- */ /* Give values to center matrix */
        /* ---------------------------- */ /* Fill matrix row of (i,j) */
        row    = center_index;
        column = center_index;
        value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta * std::fabs(detDF); /* Center: beta_{i,j} */
        updateMatrixElement(center_solver, center_batch, row, column, value);

        row    = center_index;
        column = left_index;
        value  = -coeff1 * arr; /* Left */
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
    /* ------------------------------------------ */
    /* Radial Section: Node on the outer boundary */
    /* ------------------------------------------ */
    else if (i_r == grid.nr() - 1) {
        double h1     = grid.radialSpacing(i_r - 1);
        double k1     = grid.angularSpacing(i_theta - 1);
        double k2     = grid.angularSpacing(i_theta);
        double coeff1 = 0.5 * (k1 + k2) / h1;

        /* -----------|| */
        /* o   o   o  || */
        /* -----------|| */
        /* o   o   O  || <- center_matrix */
        /* -----------|| */
        /* o   o   o  || */
        /* -----------|| */
        auto& center_solver = radial_tridiagonal_solver;
        int center_batch    = i_theta;

        const int center_index = i_r - numberSmootherCircles;
        const int left_index   = i_r - numberSmootherCircles - 1;

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

template <class LevelCacheType>
void SmootherGive<LevelCacheType>::buildTridiagonalCircleSection(int i_r)
{
    const PolarGrid& grid             = Smoother<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = Smoother<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = Smoother<LevelCacheType>::DirBC_Interior_;

    // Access pattern is aligned with the memory layout of the grid data to maximize cache efficiency.
    const double r = grid.radius(i_r);
    for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
        const int global_index = grid.index(i_r, i_theta);
        const double theta     = grid.theta(i_theta);

        double coeff_beta, arr, att, art, detDF;
        level_cache.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art, detDF);

        // Build Asc at the current node
        nodeBuildTridiagonalSolverMatrices(i_r, i_theta, grid, DirBC_Interior, circle_tridiagonal_solver_,
                                           radial_tridiagonal_solver_, arr, att, art, detDF, coeff_beta);
    }
}

template <class LevelCacheType>
void SmootherGive<LevelCacheType>::buildTridiagonalRadialSection(int i_theta)
{
    const PolarGrid& grid             = Smoother<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = Smoother<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = Smoother<LevelCacheType>::DirBC_Interior_;

    // Access pattern is aligned with the memory layout of the grid data to maximize cache efficiency.
    const double theta = grid.theta(i_theta);
    for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
        const int global_index = grid.index(i_r, i_theta);
        const double r         = grid.radius(i_r);

        double coeff_beta, arr, att, art, detDF;
        level_cache.obtainValues(i_r, i_theta, global_index, r, theta, coeff_beta, arr, att, art, detDF);

        // Build Asc at the current node
        nodeBuildTridiagonalSolverMatrices(i_r, i_theta, grid, DirBC_Interior, circle_tridiagonal_solver_,
                                           radial_tridiagonal_solver_, arr, att, art, detDF, coeff_beta);
    }
}

template <class LevelCacheType>
void SmootherGive<LevelCacheType>::buildTridiagonalSolverMatrices()
{
    const PolarGrid& grid             = Smoother<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = Smoother<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = Smoother<LevelCacheType>::DirBC_Interior_;
    const int num_omp_threads         = Smoother<LevelCacheType>::num_omp_threads_;

    const int num_smoother_circles    = grid.numberSmootherCircles();
    const int additional_radial_tasks = grid.ntheta() % 3;
    const int num_radial_tasks        = grid.ntheta() - additional_radial_tasks;

    /* ---------------- */
    /* Circular section */
    /* ---------------- */
    // We parallelize the loop with step 3 to avoid data race conditions between adjacent circles.
#pragma omp parallel num_threads(num_omp_threads)
    {
#pragma omp for
        for (int i_r = 0; i_r < num_smoother_circles; i_r += 3) {
            buildTridiagonalCircleSection(i_r);
        } /* Implicit barrier */
#pragma omp for
        for (int i_r = 1; i_r < num_smoother_circles; i_r += 3) {
            buildTridiagonalCircleSection(i_r);
        } /* Implicit barrier */
#pragma omp for
        for (int i_r = 2; i_r < num_smoother_circles; i_r += 3) {
            buildTridiagonalCircleSection(i_r);
        } /* Implicit barrier */
    }

    /* ---------------- */
    /* Radial section */
    /* ---------------- */
    // We parallelize the loop with step 3 to avoid data race conditions between adjacent radial lines.
    // Due to the periodicity in the angular direction, we can have at most 2 additional radial tasks
    // that are handled serially before the parallel loops.
    if (additional_radial_tasks > 0) {
        const int i_theta = 0;
        buildTridiagonalRadialSection(i_theta);
    }

    if (additional_radial_tasks > 1) {
        const int i_theta = 1;
        buildTridiagonalRadialSection(i_theta);
    }

#pragma omp parallel num_threads(num_omp_threads)
    {
#pragma omp for
        for (int radial_task = 0; radial_task < num_radial_tasks; radial_task += 3) {
            const int i_theta = radial_task + additional_radial_tasks;
            buildTridiagonalRadialSection(i_theta);
        } /* Implicit barrier */
#pragma omp for
        for (int radial_task = 1; radial_task < num_radial_tasks; radial_task += 3) {
            const int i_theta = radial_task + additional_radial_tasks;
            buildTridiagonalRadialSection(i_theta);
        } /* Implicit barrier */
#pragma omp for
        for (int radial_task = 2; radial_task < num_radial_tasks; radial_task += 3) {
            const int i_theta = radial_task + additional_radial_tasks;
            buildTridiagonalRadialSection(i_theta);
        } /* Implicit barrier */
    }
}