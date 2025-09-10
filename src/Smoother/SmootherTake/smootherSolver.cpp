#include "../../../include/Smoother/SmootherTake/smootherTake.h"

#define NODE_APPLY_ASC_ORTHO_CIRCLE_TAKE(i_r, i_theta, grid, DirBC_Interior, smoother_color, x, rhs, temp, arr, att,                    \
                                         art, detDF, coeff_beta)                                                                        \
    do {                                                                                                                                \
        assert(i_r >= 0 && i_r < grid_.numberSmootherCircles());                                                                        \
                                                                                                                                        \
        if (i_r > 0 && i_r < grid.numberSmootherCircles()) {                                                                            \
            double h1     = grid.radialSpacing(i_r - 1);                                                                                \
            double h2     = grid.radialSpacing(i_r);                                                                                    \
            double k1     = grid.angularSpacing(i_theta - 1);                                                                           \
            double k2     = grid.angularSpacing(i_theta);                                                                               \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                                       \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                                       \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                                       \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                                       \
                                                                                                                                        \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                                    \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                                    \
                                                                                                                                        \
            const int bottom_left  = grid.index(i_r - 1, i_theta_M1);                                                                   \
            const int left         = grid.index(i_r - 1, i_theta);                                                                      \
            const int top_left     = grid.index(i_r - 1, i_theta_P1);                                                                   \
            const int bottom       = grid.index(i_r, i_theta_M1);                                                                       \
            const int center       = grid.index(i_r, i_theta);                                                                          \
            const int top          = grid.index(i_r, i_theta_P1);                                                                       \
            const int bottom_right = grid.index(i_r + 1, i_theta_M1);                                                                   \
            const int right        = grid.index(i_r + 1, i_theta);                                                                      \
            const int top_right    = grid.index(i_r + 1, i_theta_P1);                                                                   \
                                                                                                                                        \
            temp[center] = rhs[center] - (-coeff1 * (arr[center] + arr[left]) * x[left] /* Left */                                      \
                                          - coeff2 * (arr[center] + arr[right]) * x[right] /* Right */                                  \
                                                                                                                                        \
                                          - 0.25 * (art[left] + art[bottom]) * x[bottom_left] /* Bottom Left */                         \
                                          + 0.25 * (art[right] + art[bottom]) * x[bottom_right] /* Bottom Right */                      \
                                          + 0.25 * (art[left] + art[top]) * x[top_left] /* Top Left */                                  \
                                          - 0.25 * (art[right] + art[top]) * x[top_right] /* Top Right */                               \
                                         );                                                                                             \
        }                                                                                                                               \
        else if (i_r == 0) {                                                                                                            \
            if (DirBC_Interior) {                                                                                                       \
                const int center = grid.index(i_r, i_theta);                                                                            \
                temp[center]     = rhs[center];                                                                                         \
            }                                                                                                                           \
            else {                                                                                                                      \
                /* ------------------------------------------------------------- */                                                     \
                /* Case 2: Across origin discretization on the interior boundary */                                                     \
                /* ------------------------------------------------------------- */                                                     \
                /* h1 gets replaced with 2 * R0. */                                                                                     \
                /* (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()/2)). */                                            \
                /* Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil. */                     \
                double h1 = 2.0 * grid.radius(0);                                                                                       \
                double h2 = grid.radialSpacing(i_r);                                                                                    \
                double k1 = grid.angularSpacing(i_theta - 1);                                                                           \
                double k2 = grid.angularSpacing(i_theta);                                                                               \
                                                                                                                                        \
                double coeff1 = 0.5 * (k1 + k2) / h1;                                                                                   \
                double coeff2 = 0.5 * (k1 + k2) / h2;                                                                                   \
                double coeff3 = 0.5 * (h1 + h2) / k1;                                                                                   \
                double coeff4 = 0.5 * (h1 + h2) / k2;                                                                                   \
                                                                                                                                        \
                const int i_theta_M1     = grid_.wrapThetaIndex(i_theta - 1);                                                           \
                const int i_theta_P1     = grid_.wrapThetaIndex(i_theta + 1);                                                           \
                const int i_theta_Across = grid_.wrapThetaIndex(i_theta + grid_.ntheta() / 2);                                          \
                                                                                                                                        \
                const int left         = grid_.index(i_r, i_theta_Across);                                                              \
                const int bottom       = grid_.index(i_r, i_theta_M1);                                                                  \
                const int center       = grid.index(i_r, i_theta);                                                                      \
                const int top          = grid_.index(i_r, i_theta_P1);                                                                  \
                const int bottom_right = grid_.index(i_r + 1, i_theta_M1);                                                              \
                const int right        = grid_.index(i_r + 1, i_theta);                                                                 \
                const int top_right    = grid_.index(i_r + 1, i_theta_P1);                                                              \
                                                                                                                                        \
                temp[center] =                                                                                                          \
                    rhs[center] -                                                                                                       \
                    (-coeff2 * (arr[center] + arr[right]) * x[right] /* Right */                                                        \
                                                                                                                                        \
                     /* - 0.25 * (art[left] + art[bottom]) * x[bottom_left] // Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
                     + 0.25 * (art[right] + art[bottom]) * x[bottom_right] /* Bottom Right */                                           \
                                                                                                                                        \
                     /* + 0.25 * (art[left] + art[top]) * x[top_left] // Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */          \
                     - 0.25 * (art[right] + art[top]) * x[top_right] /* Top Right */                                                    \
                    );                                                                                                                  \
            }                                                                                                                           \
        }                                                                                                                               \
    } while (0)

#define NODE_APPLY_ASC_ORTHO_RADIAL_TAKE(i_r, i_theta, grid, DirBC_Interior, smoother_color, x, rhs, temp, arr, att,   \
                                         art, detDF, coeff_beta)                                                       \
    do {                                                                                                               \
        assert(i_r >= grid.numberSmootherCircles() && i_r < grid.nr());                                                \
        /* -------------------- */                                                                                     \
        /* Node in the interior */                                                                                     \
        /* -------------------- */                                                                                     \
        if (i_r > grid.numberSmootherCircles() && i_r < grid.nr() - 2) {                                               \
            double h1 = grid.radialSpacing(i_r - 1);                                                                   \
            double h2 = grid.radialSpacing(i_r);                                                                       \
            double k1 = grid.angularSpacing(i_theta - 1);                                                              \
            double k2 = grid.angularSpacing(i_theta);                                                                  \
                                                                                                                       \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                      \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                      \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                      \
                                                                                                                       \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                   \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                   \
                                                                                                                       \
            const int bottom_left  = grid.index(i_r - 1, i_theta_M1);                                                  \
            const int left         = grid.index(i_r - 1, i_theta);                                                     \
            const int top_left     = grid.index(i_r - 1, i_theta_P1);                                                  \
            const int bottom       = grid.index(i_r, i_theta_M1);                                                      \
            const int center       = grid.index(i_r, i_theta);                                                         \
            const int top          = grid.index(i_r, i_theta_P1);                                                      \
            const int bottom_right = grid.index(i_r + 1, i_theta_M1);                                                  \
            const int right        = grid.index(i_r + 1, i_theta);                                                     \
            const int top_right    = grid.index(i_r + 1, i_theta_P1);                                                  \
                                                                                                                       \
            temp[center] = rhs[center] - (-coeff3 * (att[center] + att[bottom]) * x[bottom] /* Bottom */               \
                                          - coeff4 * (att[center] + att[top]) * x[top] /* Top */                       \
                                                                                                                       \
                                          - 0.25 * (art[left] + art[bottom]) * x[bottom_left] /* Bottom Left */        \
                                          + 0.25 * (art[right] + art[bottom]) * x[bottom_right] /* Bottom Right */     \
                                          + 0.25 * (art[left] + art[top]) * x[top_left] /* Top Left */                 \
                                          - 0.25 * (art[right] + art[top]) * x[top_right] /* Top Right */              \
                                         );                                                                            \
        }                                                                                                              \
        else if (i_r == grid.numberSmootherCircles()) {                                                                \
            double h1 = grid.radialSpacing(i_r - 1);                                                                   \
            double h2 = grid.radialSpacing(i_r);                                                                       \
            double k1 = grid.angularSpacing(i_theta - 1);                                                              \
            double k2 = grid.angularSpacing(i_theta);                                                                  \
                                                                                                                       \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                      \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                      \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                      \
                                                                                                                       \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                   \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                   \
                                                                                                                       \
            const int bottom_left  = grid.index(i_r - 1, i_theta_M1);                                                  \
            const int left         = grid.index(i_r - 1, i_theta);                                                     \
            const int top_left     = grid.index(i_r - 1, i_theta_P1);                                                  \
            const int bottom       = grid.index(i_r, i_theta_M1);                                                      \
            const int center       = grid.index(i_r, i_theta);                                                         \
            const int top          = grid.index(i_r, i_theta_P1);                                                      \
            const int bottom_right = grid.index(i_r + 1, i_theta_M1);                                                  \
            const int right        = grid.index(i_r + 1, i_theta);                                                     \
            const int top_right    = grid.index(i_r + 1, i_theta_P1);                                                  \
                                                                                                                       \
            temp[center] = rhs[center] - (-coeff1 * (arr[center] + arr[left]) * x[left] /* Left */                     \
                                          - coeff3 * (att[center] + att[bottom]) * x[bottom] /* Bottom */              \
                                          - coeff4 * (att[center] + att[top]) * x[top] /* Top */                       \
                                                                                                                       \
                                          - 0.25 * (art[left] + art[bottom]) * x[bottom_left] /* Bottom Left */        \
                                          + 0.25 * (art[right] + art[bottom]) * x[bottom_right] /* Bottom Right */     \
                                          + 0.25 * (art[left] + art[top]) * x[top_left] /* Top Left */                 \
                                          - 0.25 * (art[right] + art[top]) * x[top_right] /* Top Right */              \
                                         );                                                                            \
        }                                                                                                              \
        else if (i_r == grid.nr() - 2) {                                                                               \
            double h1 = grid.radialSpacing(i_r - 1);                                                                   \
            double h2 = grid.radialSpacing(i_r);                                                                       \
            double k1 = grid.angularSpacing(i_theta - 1);                                                              \
            double k2 = grid.angularSpacing(i_theta);                                                                  \
                                                                                                                       \
            double coeff1 = 0.5 * (k1 + k2) / h1;                                                                      \
            double coeff2 = 0.5 * (k1 + k2) / h2;                                                                      \
            double coeff3 = 0.5 * (h1 + h2) / k1;                                                                      \
            double coeff4 = 0.5 * (h1 + h2) / k2;                                                                      \
                                                                                                                       \
            const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);                                                   \
            const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);                                                   \
                                                                                                                       \
            const int bottom_left  = grid.index(i_r - 1, i_theta_M1);                                                  \
            const int left         = grid.index(i_r - 1, i_theta);                                                     \
            const int top_left     = grid.index(i_r - 1, i_theta_P1);                                                  \
            const int bottom       = grid.index(i_r, i_theta_M1);                                                      \
            const int center       = grid.index(i_r, i_theta);                                                         \
            const int top          = grid.index(i_r, i_theta_P1);                                                      \
            const int bottom_right = grid.index(i_r + 1, i_theta_M1);                                                  \
            const int right        = grid.index(i_r + 1, i_theta);                                                     \
            const int top_right    = grid.index(i_r + 1, i_theta_P1);                                                  \
            /* "Right" is part of the radial Asc smoother matrices, */                                                 \
            /* but is shifted over to the rhs to make the radial Asc smoother matrices symmetric. */                   \
            /* Note that the circle Asc smoother matrices are symmetric by default. */                                 \
            /* Note that rhs[right] contains the correct boundary value of u_D. */                                     \
            temp[center] = rhs[center] - (-coeff2 * (arr[center] + arr[right]) * rhs[right] /* Right */                \
                                          - coeff3 * (att[center] + att[bottom]) * x[bottom] /* Bottom */              \
                                          - coeff4 * (att[center] + att[top]) * x[top] /* Top */                       \
                                                                                                                       \
                                          - 0.25 * (art[left] + art[bottom]) * x[bottom_left] /* Bottom Left */        \
                                          + 0.25 * (art[right] + art[bottom]) * x[bottom_right] /* Bottom Right */     \
                                          + 0.25 * (art[left] + art[top]) * x[top_left] /* Top Left */                 \
                                          - 0.25 * (art[right] + art[top]) * x[top_right] /* Top Right */              \
                                         );                                                                            \
        }                                                                                                              \
        else if (i_r == grid.nr() - 1) {                                                                               \
            const int center = grid.index(i_r, i_theta);                                                               \
            temp[center]     = rhs[center];                                                                            \
        }                                                                                                              \
    } while (0)

void SmootherTake::applyAscOrthoCircleSection(const int i_r, const SmootherColor smoother_color,
                                              const Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> x,
                                              const Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> rhs,
                                              Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> temp)
{
    assert(i_r >= 0 && i_r < grid_.numberSmootherCircles());

    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    const auto& arr        = level_cache_.arr();
    const auto& att        = level_cache_.att();
    const auto& art        = level_cache_.art();
    const auto& detDF      = level_cache_.detDF();
    const auto& coeff_beta = level_cache_.coeff_beta();

    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
        NODE_APPLY_ASC_ORTHO_CIRCLE_TAKE(i_r, i_theta, grid_, DirBC_Interior_, smoother_color, x, rhs, temp, arr, att,
                                         art, detDF, coeff_beta);
    }
}

void SmootherTake::applyAscOrthoRadialSection(const int i_theta, const SmootherColor smoother_color,
                                              const Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> x,
                                              const Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> rhs,
                                              Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> temp)
{
    assert(i_theta >= 0 && i_theta < grid_.ntheta());

    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    const auto& arr        = level_cache_.arr();
    const auto& att        = level_cache_.att();
    const auto& art        = level_cache_.art();
    const auto& detDF      = level_cache_.detDF();
    const auto& coeff_beta = level_cache_.coeff_beta();

    for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
        NODE_APPLY_ASC_ORTHO_RADIAL_TAKE(i_r, i_theta, grid_, DirBC_Interior_, smoother_color, x, rhs, temp, arr, att,
                                         art, detDF, coeff_beta);
    }
}

void SmootherTake::solveCircleSection(const int i_r, Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> x,
                                      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> temp,
                                      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> solver_storage_1,
                                      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> solver_storage_2)
{
    const int start = grid_.index(i_r, 0);
    const int end   = start + grid_.ntheta();
    if (i_r == 0) {
#ifdef GMGPOLAR_USE_MUMPS
        inner_boundary_mumps_solver_.job    = JOB_COMPUTE_SOLUTION;
        inner_boundary_mumps_solver_.nrhs   = 1; // single rhs vector
        inner_boundary_mumps_solver_.nz_rhs = grid_.ntheta(); // non-zeros in rhs
        inner_boundary_mumps_solver_.rhs    = temp.begin() + start;
        inner_boundary_mumps_solver_.lrhs   = grid_.ntheta(); // leading dimension of rhs
        dmumps_c(&inner_boundary_mumps_solver_);
        if (inner_boundary_mumps_solver_.info[0] != 0) {
            std::cerr << "Error solving the system: " << inner_boundary_mumps_solver_.info[0] << std::endl;
        }
#else
        inner_boundary_lu_solver_.solveInPlace(temp.begin() + start);
#endif
    }
    else {
        circle_tridiagonal_solver_[i_r].solveInPlace(temp.begin() + start, solver_storage_1.begin(),
                                                     solver_storage_2.begin());
    }
    // Move updated values to x
    std::move(temp.begin() + start, temp.begin() + end, x.begin() + start);
}

void SmootherTake::solveRadialSection(const int i_theta,
                                      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> x,
                                      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> temp,
                                      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> solver_storage)
{
    const int start = grid_.index(grid_.numberSmootherCircles(), i_theta);
    const int end   = start + grid_.lengthSmootherRadial();

    radial_tridiagonal_solver_[i_theta].solveInPlace(temp.begin() + start, solver_storage.begin());
    // Move updated values to x
    std::move(temp.begin() + start, temp.begin() + end, x.begin() + start);
}

// clang-format off

// In temp we store the vector 'rhs - A_sc^ortho u_sc^ortho' and then we solve the system
// Asc * u_sc = temp in place and move the updated values into 'x'.
void SmootherTake::smoothing(Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> x, const Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> rhs, Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> temp)
{
    assert(x.size() == rhs.size());
    assert(temp.size() == rhs.size());

    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    #pragma omp parallel
    {
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> circle_solver_storage_1("circle_solver_storage_1",grid_.ntheta());
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> circle_solver_storage_2("circle_solver_storage_2",grid_.ntheta());
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> radial_solver_storage("circle_solver_storage",grid_.lengthSmootherRadial());

        /* The outer most circle next to the radial section is defined to be black. */
        /* Priority: Black -> White. */
        const int start_black_circles = (grid_.numberSmootherCircles() % 2 == 0) ? 1 : 0;
        const int start_white_circles = (grid_.numberSmootherCircles() % 2 == 0) ? 0 : 1;

        /* Black Circle Section */
        #pragma omp for
        for (int i_r = start_black_circles; i_r < grid_.numberSmootherCircles(); i_r += 2) {
            applyAscOrthoCircleSection(i_r, SmootherColor::Black, x, rhs, temp);
            solveCircleSection(i_r, x, temp, circle_solver_storage_1, circle_solver_storage_2);
        } /* Implicit barrier */

        /* White Circle Section */
        #pragma omp for nowait
        for (int i_r = start_white_circles; i_r < grid_.numberSmootherCircles(); i_r += 2) {
            applyAscOrthoCircleSection(i_r, SmootherColor::White, x, rhs, temp);
            solveCircleSection(i_r, x, temp, circle_solver_storage_1, circle_solver_storage_2);
        }
        /* Black Radial Section */
        #pragma omp for
        for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta += 2) {
            applyAscOrthoRadialSection(i_theta, SmootherColor::Black, x, rhs, temp);
            solveRadialSection(i_theta, x, temp, radial_solver_storage);
        } /* Implicit barrier */

        /* White Radial Section*/
        #pragma omp for
        for (int i_theta = 1; i_theta < grid_.ntheta(); i_theta += 2) {
            applyAscOrthoRadialSection(i_theta, SmootherColor::White, x, rhs, temp);
            solveRadialSection(i_theta, x, temp, radial_solver_storage);
        } /* Implicit barrier */
    }
}
// clang-format on
