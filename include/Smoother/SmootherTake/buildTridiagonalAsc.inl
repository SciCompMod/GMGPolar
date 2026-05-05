#include "../../../include/Smoother/SmootherTake/smootherTake.h"

namespace smoother_take
{

static inline void updateMatrixElement(const BatchedTridiagonalSolver<double>& solver, int batch, int row, int column,
                                       double value)
{
    if (row == column)
        solver.set_main_diagonal(batch, row, value);
    else if (row == column - 1)
        solver.set_sub_diagonal(batch, row, value);
    else if (row == 0 && column == solver.matrixDimension() - 1)
        solver.set_cyclic_corner(batch, value);
}

// Build the tridiagonal solver matrices for a specific node (i_r, i_theta)
static inline void nodeBuildTridiagonalSolverMatrices(int i_r, int i_theta, const PolarGrid& grid, bool DirBC_Interior,
                                                      const BatchedTridiagonalSolver<double>& circle_tridiagonal_solver,
                                                      const BatchedTridiagonalSolver<double>& radial_tridiagonal_solver,
                                                      ConstVector<double>& arr, ConstVector<double>& att,
                                                      ConstVector<double>& art, ConstVector<double>& detDF,
                                                      ConstVector<double>& coeff_beta)
{
    using smoother_take::updateMatrixElement;

    assert(i_r >= 0 && i_r < grid.nr());
    assert(i_theta >= 0 && i_theta < grid.ntheta());

    const int numberSmootherCircles = grid.numberSmootherCircles();
    const int lengthRadialSmoother  = grid.lengthRadialSmoother();

    assert(numberSmootherCircles >= 2);
    assert(lengthRadialSmoother >= 3);

    int row, column;
    double value;
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

        const int left   = grid.index(i_r - 1, i_theta);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int center = grid.index(i_r, i_theta);
        const int top    = grid.index(i_r, i_theta_P1);
        const int right  = grid.index(i_r + 1, i_theta);

        auto& solver    = circle_tridiagonal_solver;
        const int batch = i_r;

        const int center_index = i_theta;
        const int bottom_index = i_theta_M1;
        const int top_index    = i_theta_P1;

        const double left_value   = -coeff1 * (arr[center] + arr[left]);
        const double right_value  = -coeff2 * (arr[center] + arr[right]);
        const double bottom_value = -coeff3 * (att[center] + att[bottom]);
        const double top_value    = -coeff4 * (att[center] + att[top]);

        /* Center: (Left, Right, Bottom, Top) */
        row    = center_index;
        column = center_index;
        value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center] * std::fabs(detDF[center]) -
                 (left_value + right_value + bottom_value + top_value);
        updateMatrixElement(solver, batch, row, column, value);

        /* Bottom */
        row    = center_index;
        column = bottom_index;
        value  = bottom_value;
        updateMatrixElement(solver, batch, row, column, value);

        /* Top */
        row    = center_index;
        column = top_index;
        value  = top_value;
        updateMatrixElement(solver, batch, row, column, value);
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

        const int left   = grid.index(i_r - 1, i_theta);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int center = grid.index(i_r, i_theta);
        const int top    = grid.index(i_r, i_theta_P1);
        const int right  = grid.index(i_r + 1, i_theta);

        auto& solver    = radial_tridiagonal_solver;
        const int batch = i_theta;

        const int center_index = i_r - numberSmootherCircles;
        const int left_index   = i_r - numberSmootherCircles - 1;
        const int right_index  = i_r - numberSmootherCircles + 1;

        const double left_value   = -coeff1 * (arr[center] + arr[left]);
        const double right_value  = -coeff2 * (arr[center] + arr[right]);
        const double bottom_value = -coeff3 * (att[center] + att[bottom]);
        const double top_value    = -coeff4 * (att[center] + att[top]);

        /* Center: (Left, Right, Bottom, Top) */
        row    = center_index;
        column = center_index;
        value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center] * std::fabs(detDF[center]) -
                 (left_value + right_value + bottom_value + top_value);
        updateMatrixElement(solver, batch, row, column, value);

        /* Left */
        row    = center_index;
        column = left_index;
        value  = left_value;
        updateMatrixElement(solver, batch, row, column, value);

        /* Right */
        row    = center_index;
        column = right_index;
        value  = right_value;
        updateMatrixElement(solver, batch, row, column, value);
    }
    /* ------------------------------------------ */
    /* Circle Section: Node in the inner boundary */
    /* ------------------------------------------ */
    else if (i_r == 0) {
        // The inner boundary circle line are is handled by the inner_boundary_mumps_solver, so we fill in the identity matrix.
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

        const int left   = grid.index(i_r - 1, i_theta);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int center = grid.index(i_r, i_theta);
        const int top    = grid.index(i_r, i_theta_P1);
        const int right  = grid.index(i_r + 1, i_theta);

        auto& solver    = radial_tridiagonal_solver;
        const int batch = i_theta;

        const int center_index = i_r - numberSmootherCircles;
        const int right_index  = i_r - numberSmootherCircles + 1;

        const double left_value   = -coeff1 * (arr[center] + arr[left]);
        const double right_value  = -coeff2 * (arr[center] + arr[right]);
        const double bottom_value = -coeff3 * (att[center] + att[bottom]);
        const double top_value    = -coeff4 * (att[center] + att[top]);

        /* Center: (Left, Right, Bottom, Top) */
        row    = center_index;
        column = center_index;
        value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center] * std::fabs(detDF[center]) -
                 (left_value + right_value + bottom_value + top_value);
        updateMatrixElement(solver, batch, row, column, value);

        /* Right */
        row    = center_index;
        column = right_index;
        value  = right_value;
        updateMatrixElement(solver, batch, row, column, value);
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

        const int left   = grid.index(i_r - 1, i_theta);
        const int bottom = grid.index(i_r, i_theta_M1);
        const int center = grid.index(i_r, i_theta);
        const int top    = grid.index(i_r, i_theta_P1);
        const int right  = grid.index(i_r + 1, i_theta);

        auto& solver    = radial_tridiagonal_solver;
        const int batch = i_theta;

        const int center_index = i_r - numberSmootherCircles;
        const int left_index   = i_r - numberSmootherCircles - 1;
        const int right_index  = i_r - numberSmootherCircles + 1;

        const double left_value   = -coeff1 * (arr[center] + arr[left]);
        const double right_value  = -coeff2 * (arr[center] + arr[right]);
        const double bottom_value = -coeff3 * (att[center] + att[bottom]);
        const double top_value    = -coeff4 * (att[center] + att[top]);

        /* Center: (Left, Right, Bottom, Top) */
        row    = center_index;
        column = center_index;
        value  = 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[center] * std::fabs(detDF[center]) -
                 (left_value + right_value + bottom_value + top_value);
        updateMatrixElement(solver, batch, row, column, value);

        /* Left */
        row    = center_index;
        column = left_index;
        value  = left_value;
        updateMatrixElement(solver, batch, row, column, value);

        /* Right: NOT INCLUDED! */
        row    = center_index;
        column = right_index;
        value  = 0.0;
        updateMatrixElement(solver, batch, row, column, value);
    }
    /* ------------------------------------------ */
    /* Radial Section: Node on the outer boundary */
    /* ------------------------------------------ */
    else if (i_r == grid.nr() - 1) {
        auto& solver    = radial_tridiagonal_solver;
        const int batch = i_theta;

        const int center_index = i_r - numberSmootherCircles;
        const int left_index   = i_r - numberSmootherCircles - 1;

        /* Fill matrix row of (i,j) */
        row    = center_index;
        column = center_index;
        value  = 1.0;
        updateMatrixElement(solver, batch, row, column, value);

        /* Left: NOT INCLUDED */
        row    = center_index;
        column = left_index;
        value  = 0.0;
        updateMatrixElement(solver, batch, row, column, value);
    }
}

} // namespace smoother_take

template <class LevelCacheType>
void SmootherTake<LevelCacheType>::buildTridiagonalSolverMatrices()
{
    using smoother_take::nodeBuildTridiagonalSolverMatrices;

    const PolarGrid& grid             = Smoother<LevelCacheType>::grid_;
    const LevelCacheType& level_cache = Smoother<LevelCacheType>::level_cache_;
    const bool DirBC_Interior         = Smoother<LevelCacheType>::DirBC_Interior_;
    const int num_omp_threads         = Smoother<LevelCacheType>::num_omp_threads_;

    const BatchedTridiagonalSolver<double>* circle_tridiagonal_solver_ptr = &circle_tridiagonal_solver_;
    const BatchedTridiagonalSolver<double>* radial_tridiagonal_solver_ptr = &radial_tridiagonal_solver_;

    assert(level_cache.cacheDensityProfileCoefficients());
    assert(level_cache.cacheDomainGeometry());

    ConstVector<double> arr        = level_cache.arr();
    ConstVector<double> att        = level_cache.att();
    ConstVector<double> art        = level_cache.art();
    ConstVector<double> detDF      = level_cache.detDF();
    ConstVector<double> coeff_beta = level_cache.coeff_beta();

    const PolarGrid* grid_ptr = &grid;

    /* We split the loops into two regions to better respect the */
    /* access patterns of the smoother and improve cache locality. */

    // The For loop matches circular access pattern */
    Kokkos::parallel_for(
        "Smoother Take: BuildTridiagonalAsc (Circular)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>( // Rank of the index space
            {0, 0}, // Starting point of the index space
            {grid.numberSmootherCircles(), grid.ntheta()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_r, const int i_theta) {
            nodeBuildTridiagonalSolverMatrices(i_r, i_theta, *grid_ptr, DirBC_Interior, circle_tridiagonal_solver,
                                               radial_tridiagonal_solver, arr, att, art, detDF, coeff_beta);
        });

    /* For loop matches radial access pattern */
    Kokkos::parallel_for(
        "Smoother Take: BuildTridiagonalAsc (Radial)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>( // Rank of the index space
            {0, grid.numberSmootherCircles()}, // Starting point of the index space
            {grid.ntheta(), grid.nr()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_theta, const int i_r) {
            nodeBuildTridiagonalSolverMatrices(i_r, i_theta, *grid_ptr, DirBC_Interior, circle_tridiagonal_solver,
                                               radial_tridiagonal_solver, arr, att, art, detDF, coeff_beta);
        });

    Kokkos::fence();
}