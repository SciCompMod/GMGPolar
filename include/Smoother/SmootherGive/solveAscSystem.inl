#pragma once

template <class LevelCacheType>
void SmootherGive<LevelCacheType>::solveBlackCircleSection(Vector<double> x, Vector<double> temp)
{
    const PolarGrid& grid     = Smoother<LevelCacheType>::grid_;
    const int num_omp_threads = Smoother<LevelCacheType>::num_omp_threads_;

    int start                     = 0;
    int end                       = grid.numberCircularSmootherNodes();
    Vector<double> circle_section = Kokkos::subview(temp, Kokkos::make_pair(start, end));

    bool is_inner_circle_black = grid.numberSmootherCircles() % 2 != 0;

    int batch_offset = is_inner_circle_black ? 2 : 1;
    int batch_stride = 2;
    circle_tridiagonal_solver_.solve(circle_section, batch_offset, batch_stride);

    if (is_inner_circle_black) {
        Vector<double> inner_boundary = Kokkos::subview(temp, Kokkos::make_pair(0, grid.ntheta()));
        inner_boundary_solver_.solveInPlace(inner_boundary);
    }

    // Move updated values to x
    int start_black_circles = is_inner_circle_black ? 0 : 1;
#pragma omp parallel for num_threads(num_omp_threads)
    for (int i_r = start_black_circles; i_r < grid.numberSmootherCircles(); i_r += 2) {
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
            x[grid.index(i_r, i_theta)] = temp[grid.index(i_r, i_theta)];
        }
    }
}

template <class LevelCacheType>
void SmootherGive<LevelCacheType>::solveWhiteCircleSection(Vector<double> x, Vector<double> temp)
{
    const PolarGrid& grid     = Smoother<LevelCacheType>::grid_;
    const int num_omp_threads = Smoother<LevelCacheType>::num_omp_threads_;

    int start                     = 0;
    int end                       = grid.numberCircularSmootherNodes();
    Vector<double> circle_section = Kokkos::subview(temp, Kokkos::make_pair(start, end));

    bool is_inner_circle_white = grid.numberSmootherCircles() % 2 == 0;

    int batch_offset = is_inner_circle_white ? 2 : 1;
    int batch_stride = 2;
    circle_tridiagonal_solver_.solve(circle_section, batch_offset, batch_stride);

    if (is_inner_circle_white) {
        Vector<double> inner_boundary = Kokkos::subview(temp, Kokkos::make_pair(0, grid.ntheta()));
        inner_boundary_solver_.solveInPlace(inner_boundary);
    }

    // Move updated values to x
    int start_white_circles = is_inner_circle_white ? 0 : 1;
#pragma omp parallel for num_threads(num_omp_threads)
    for (int i_r = start_white_circles; i_r < grid.numberSmootherCircles(); i_r += 2) {
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
            x[grid.index(i_r, i_theta)] = temp[grid.index(i_r, i_theta)];
        }
    }
}

template <class LevelCacheType>
void SmootherGive<LevelCacheType>::solveBlackRadialSection(Vector<double> x, Vector<double> temp)
{
    const PolarGrid& grid     = Smoother<LevelCacheType>::grid_;
    const int num_omp_threads = Smoother<LevelCacheType>::num_omp_threads_;

    int start                     = grid.numberCircularSmootherNodes();
    int end                       = grid.numberOfNodes();
    Vector<double> radial_section = Kokkos::subview(temp, Kokkos::make_pair(start, end));

    int batch_offset = 0;
    int batch_stride = 2;
    radial_tridiagonal_solver_.solve(radial_section, batch_offset, batch_stride);

// Move updated values to x
#pragma omp parallel for num_threads(num_omp_threads)
    for (int i_theta = 0; i_theta < grid.ntheta(); i_theta += 2) {
        for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
            x[grid.index(i_r, i_theta)] = temp[grid.index(i_r, i_theta)];
        }
    }
}

template <class LevelCacheType>
void SmootherGive<LevelCacheType>::solveWhiteRadialSection(Vector<double> x, Vector<double> temp)
{
    const PolarGrid& grid     = Smoother<LevelCacheType>::grid_;
    const int num_omp_threads = Smoother<LevelCacheType>::num_omp_threads_;

    int start                     = grid.numberCircularSmootherNodes();
    int end                       = grid.numberOfNodes();
    Vector<double> radial_section = Kokkos::subview(temp, Kokkos::make_pair(start, end));

    int batch_offset = 1;
    int batch_stride = 2;
    radial_tridiagonal_solver_.solve(radial_section, batch_offset, batch_stride);

// Move updated values to x
#pragma omp parallel for num_threads(num_omp_threads)
    for (int i_theta = 1; i_theta < grid.ntheta(); i_theta += 2) {
        for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
            x[grid.index(i_r, i_theta)] = temp[grid.index(i_r, i_theta)];
        }
    }
}