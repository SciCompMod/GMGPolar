#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherGive/extrapolatedSmootherGive.h"

void ExtrapolatedSmootherGive::solveBlackCircleSection(Vector<double> x, Vector<double> temp)
{
    int start                     = 0;
    int end                       = grid_.numberCircularSmootherNodes();
    Vector<double> circle_section = Kokkos::subview(temp, Kokkos::make_pair(start, end));

    bool is_inner_circle_black = grid_.numberSmootherCircles() % 2 != 0;

    if (!is_inner_circle_black) {
        int batch_offset = 1;
        int batch_stride = 2;
        circle_tridiagonal_solver_.solve(circle_section, batch_offset, batch_stride);
    }
    else {
        int batch_offset = 2;
        int batch_stride = 2;
        circle_tridiagonal_solver_.solve_diagonal(circle_section, batch_offset, batch_stride);

#ifdef GMGPOLAR_USE_MUMPS
        inner_boundary_mumps_solver_.job    = JOB_COMPUTE_SOLUTION;
        inner_boundary_mumps_solver_.nrhs   = 1; // single rhs vector
        inner_boundary_mumps_solver_.nz_rhs = grid_.ntheta(); // non-zeros in rhs
        inner_boundary_mumps_solver_.rhs    = circle_section.data();
        inner_boundary_mumps_solver_.lrhs   = grid_.ntheta(); // leading dimension of rhs
        dmumps_c(&inner_boundary_mumps_solver_);
        if (inner_boundary_mumps_solver_.info[0] != 0) {
            std::cerr << "Error solving the system: " << inner_boundary_mumps_solver_.info[0] << std::endl;
        }
#else
        inner_boundary_lu_solver_.solveInPlace(circle_section.data());
#endif
    }

    // Move updated values to x
    int start_black_circles = is_inner_circle_black ? 0 : 1;
#pragma omp parallel for num_threads(num_omp_threads_)
    for (int i_r = start_black_circles; i_r < grid_.numberSmootherCircles(); i_r += 2) {
        for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
            x[grid_.index(i_r, i_theta)] = temp[grid_.index(i_r, i_theta)];
        }
    }
}

void ExtrapolatedSmootherGive::solveWhiteCircleSection(Vector<double> x, Vector<double> temp)
{
    int start                     = 0;
    int end                       = grid_.numberCircularSmootherNodes();
    Vector<double> circle_section = Kokkos::subview(temp, Kokkos::make_pair(start, end));

    bool is_inner_circle_white = grid_.numberSmootherCircles() % 2 == 0;

    if (!is_inner_circle_white) {
        int batch_offset = 1;
        int batch_stride = 2;
        circle_tridiagonal_solver_.solve(circle_section, batch_offset, batch_stride);
    }
    else {
        int batch_offset = 2;
        int batch_stride = 2;
        circle_tridiagonal_solver_.solve_diagonal(circle_section, batch_offset, batch_stride);

#ifdef GMGPOLAR_USE_MUMPS
        inner_boundary_mumps_solver_.job    = JOB_COMPUTE_SOLUTION;
        inner_boundary_mumps_solver_.nrhs   = 1; // single rhs vector
        inner_boundary_mumps_solver_.nz_rhs = grid_.ntheta(); // non-zeros in rhs
        inner_boundary_mumps_solver_.rhs    = circle_section.data();
        inner_boundary_mumps_solver_.lrhs   = grid_.ntheta(); // leading dimension of rhs
        dmumps_c(&inner_boundary_mumps_solver_);
        if (inner_boundary_mumps_solver_.info[0] != 0) {
            std::cerr << "Error solving the system: " << inner_boundary_mumps_solver_.info[0] << std::endl;
        }
#else
        inner_boundary_lu_solver_.solveInPlace(circle_section.data());
#endif
    }

    // Move updated values to x
    int start_white_circles = is_inner_circle_white ? 0 : 1;
#pragma omp parallel for num_threads(num_omp_threads_)
    for (int i_r = start_white_circles; i_r < grid_.numberSmootherCircles(); i_r += 2) {
        for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
            x[grid_.index(i_r, i_theta)] = temp[grid_.index(i_r, i_theta)];
        }
    }
}

void ExtrapolatedSmootherGive::solveBlackRadialSection(Vector<double> x, Vector<double> temp)
{
    int start                     = grid_.numberCircularSmootherNodes();
    int end                       = grid_.numberOfNodes();
    Vector<double> radial_section = Kokkos::subview(temp, Kokkos::make_pair(start, end));

    int batch_offset = 0;
    int batch_stride = 2;
    radial_tridiagonal_solver_.solve_diagonal(radial_section, batch_offset, batch_stride);

// Move updated values to x
#pragma omp parallel for num_threads(num_omp_threads_)
    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta += 2) {
        for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
            x[grid_.index(i_r, i_theta)] = temp[grid_.index(i_r, i_theta)];
        }
    }
}

void ExtrapolatedSmootherGive::solveWhiteRadialSection(Vector<double> x, Vector<double> temp)
{
    int start                     = grid_.numberCircularSmootherNodes();
    int end                       = grid_.numberOfNodes();
    Vector<double> radial_section = Kokkos::subview(temp, Kokkos::make_pair(start, end));

    int batch_offset = 1;
    int batch_stride = 2;
    radial_tridiagonal_solver_.solve(radial_section, batch_offset, batch_stride);

// Move updated values to x
#pragma omp parallel for num_threads(num_omp_threads_)
    for (int i_theta = 1; i_theta < grid_.ntheta(); i_theta += 2) {
        for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
            x[grid_.index(i_r, i_theta)] = temp[grid_.index(i_r, i_theta)];
        }
    }
}
