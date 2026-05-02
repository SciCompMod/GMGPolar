#pragma once

template <class LevelCacheType>
void SmootherTake<LevelCacheType>::solveBlackCircleSection(Vector<double> x, Vector<double> temp)
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
    const int start_black_circles = is_inner_circle_black ? 0 : 1;
    const int num_black_circles   = (grid.numberSmootherCircles() - start_black_circles + 1) / 2;
    Kokkos::parallel_for(
        "SmootherTake: moveUpdatedValues (Black Circular)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num_black_circles, grid.ntheta()}),
        KOKKOS_LAMBDA(const int circle_task, const int i_theta) {
            const int i_r               = start_black_circles + circle_task * 2;
            x[grid.index(i_r, i_theta)] = temp[grid.index(i_r, i_theta)];
        });
    Kokkos::fence();
}

template <class LevelCacheType>
void SmootherTake<LevelCacheType>::solveWhiteCircleSection(Vector<double> x, Vector<double> temp)
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
    const int start_white_circles = is_inner_circle_white ? 0 : 1;
    const int num_white_circles   = (grid.numberSmootherCircles() - start_white_circles + 1) / 2;
    Kokkos::parallel_for(
        "SmootherTake: moveUpdatedValues (White Circular)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num_white_circles, grid.ntheta()}),
        KOKKOS_LAMBDA(const int circle_task, const int i_theta) {
            const int i_r               = start_white_circles + circle_task * 2;
            x[grid.index(i_r, i_theta)] = temp[grid.index(i_r, i_theta)];
        });
    Kokkos::fence();
}

template <class LevelCacheType>
void SmootherTake<LevelCacheType>::solveBlackRadialSection(Vector<double> x, Vector<double> temp)
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
    assert(grid.ntheta() % 2 == 0);
    const int start_black_radials    = 0;
    const int num_black_radial_lines = grid.ntheta() / 2;
    Kokkos::parallel_for(
        "SmootherTake: moveUpdatedValues (Black Radial)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, grid.numberSmootherCircles()}, {num_black_radial_lines, grid.nr()}),
        KOKKOS_LAMBDA(const int radial_task, const int i_r) {
            const int i_theta           = start_black_radials + radial_task * 2;
            x[grid.index(i_r, i_theta)] = temp[grid.index(i_r, i_theta)];
        });
    Kokkos::fence();
}

template <class LevelCacheType>
void SmootherTake<LevelCacheType>::solveWhiteRadialSection(Vector<double> x, Vector<double> temp)
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
    assert(grid.ntheta() % 2 == 0);
    const int start_white_radials    = 1;
    const int num_white_radial_lines = grid.ntheta() / 2;
    Kokkos::parallel_for(
        "SmootherTake: moveUpdatedValues (White Radial)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, grid.numberSmootherCircles()}, {num_white_radial_lines, grid.nr()}),
        KOKKOS_LAMBDA(const int radial_task, const int i_r) {
            const int i_theta           = start_white_radials + radial_task * 2;
            x[grid.index(i_r, i_theta)] = temp[grid.index(i_r, i_theta)];
        });
    Kokkos::fence();
}
