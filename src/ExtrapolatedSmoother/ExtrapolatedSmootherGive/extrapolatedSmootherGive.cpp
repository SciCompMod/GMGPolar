#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherGive/extrapolatedSmootherGive.h"

ExtrapolatedSmootherGive::ExtrapolatedSmootherGive(const PolarGrid& grid, const LevelCache& level_cache,
                                                   const DomainGeometry& domain_geometry,
                                                   const DensityProfileCoefficients& density_profile_coefficients,
                                                   const bool DirBC_Interior, const int num_omp_threads)
    : ExtrapolatedSmoother(grid, level_cache, domain_geometry, density_profile_coefficients, DirBC_Interior,
                           num_omp_threads)
    , circle_tridiagonal_solver_(grid.ntheta(), grid.numberSmootherCircles(), true)
    , radial_tridiagonal_solver_(grid.lengthSmootherRadial(), grid.ntheta(), false)
{
    buildAscMatrices();
#ifdef GMGPOLAR_USE_MUMPS
    initializeMumpsSolver(inner_boundary_mumps_solver_, inner_boundary_circle_matrix_);
#else
    inner_boundary_lu_solver_ = SparseLUSolver<double>(inner_boundary_circle_matrix_);
#endif
}

ExtrapolatedSmootherGive::~ExtrapolatedSmootherGive()
{
#ifdef GMGPOLAR_USE_MUMPS
    finalizeMumpsSolver(inner_boundary_mumps_solver_);
#endif
}

// The smoothing solves linear systems of the form:
//   A_sc * u_sc = f_sc − A_sc^ortho * u_sc^ortho
// where:
//   - s in {Circle, Radial} denotes the smoother section type,
//   - c in {Black, White} denotes the coloring (even/odd line sub-system).
//
// The update sequence is as follows:
//   1. Black-Circle update (u_bc):
//      A_bc * u_bc = f_bc − A_bc^ortho * u_bc^ortho
//   2. White-Circle update (u_wc):
//      A_wc * u_wc = f_wc − A_wc^ortho * u_wc^ortho
//   3. Black-Radial update (u_br):
//      A_br * u_br = f_br − A_br^ortho * u_br^ortho
//   4. White-Radial update (u_wr):
//      A_wr * u_wr = f_wr − A_wr^ortho * u_wr^ortho
//
// Algorithm details:
//   - 'rhs' corresponds to the f vector, 'x' stores the final solution,
//     and 'temp' is used for temporary storage during updates.
//   - First, temp is updated with f_sc − A_sc^ortho * u_sc^ortho.
//   - The system is then solved in-place in temp, and the results
//     are copied back to x.
void ExtrapolatedSmootherGive::extrapolatedSmoothing(Vector<double> x, ConstVector<double> rhs, Vector<double> temp)
{
    assert(x.size() == rhs.size());
    assert(temp.size() == rhs.size());

#pragma omp parallel num_threads(num_omp_threads_)
    {
#pragma omp for nowait
        for (int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
            for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
                const int index = grid_.index(i_r, i_theta);
                temp[index]     = (i_r & 1 || i_theta & 1) ? rhs[index] : x[index];
            }
        }
#pragma omp for
        for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
            for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
                const int index = grid_.index(i_r, i_theta);
                temp[index]     = (i_r & 1 || i_theta & 1) ? rhs[index] : x[index];
            }
        }
    }

    /* Multi-threaded execution */
    const int num_smoother_circles = grid_.numberSmootherCircles();
    const int num_radial_lines     = grid_.ntheta();

    /* ----------------------------------------------- */
    /* 1. Black-Circle update (u_bc):                  */
    /*    A_bc * u_bc = f_bc − A_bc^ortho * u_bc^ortho */
    /* ----------------------------------------------- */
#pragma omp parallel num_threads(num_omp_threads_)
    {
        /* Inside Black Section */
#pragma omp for
        for (int circle_task = 0; circle_task < num_smoother_circles; circle_task += 2) {
            int i_r = num_smoother_circles - circle_task - 1;
            applyAscOrthoCircleSection(i_r, SmootherColor::Black, x, rhs, temp);
        }
        /* Outside Black Section (Part 1)*/
#pragma omp for
        for (int circle_task = -1; circle_task < num_smoother_circles; circle_task += 4) {
            int i_r = num_smoother_circles - circle_task - 1;
            applyAscOrthoCircleSection(i_r, SmootherColor::Black, x, rhs, temp);
        }
        /* Outside Black Section (Part 2)*/
#pragma omp for
        for (int circle_task = 1; circle_task < num_smoother_circles; circle_task += 4) {
            int i_r = num_smoother_circles - circle_task - 1;
            applyAscOrthoCircleSection(i_r, SmootherColor::Black, x, rhs, temp);
        }
    }
    solveBlackCircleSection(x, temp);

    /* ----------------------------------------------- */
    /* 2. White-Circle update (u_wc):                  */
    /*    A_wc * u_wc = f_wc − A_wc^ortho * u_wc^ortho */
    /* ----------------------------------------------- */
#pragma omp parallel num_threads(num_omp_threads_)
    {
        /* Inside White Section */
#pragma omp for
        for (int circle_task = 1; circle_task < num_smoother_circles; circle_task += 2) {
            int i_r = num_smoother_circles - circle_task - 1;
            applyAscOrthoCircleSection(i_r, SmootherColor::White, x, rhs, temp);
        }
        /* Outside White Section (Part 1)*/
#pragma omp for
        for (int circle_task = 0; circle_task < num_smoother_circles; circle_task += 4) {
            int i_r = num_smoother_circles - circle_task - 1;
            applyAscOrthoCircleSection(i_r, SmootherColor::White, x, rhs, temp);
        }
        /* Outside White Section (Part 2)*/
#pragma omp for
        for (int circle_task = 2; circle_task < num_smoother_circles; circle_task += 4) {
            int i_r = num_smoother_circles - circle_task - 1;
            applyAscOrthoCircleSection(i_r, SmootherColor::White, x, rhs, temp);
        }
    }
    solveWhiteCircleSection(x, temp);

    /* ----------------------------------------------- */
    /* 3. Black-Radial update (u_br):                  */
    /*    A_br * u_br = f_br − A_br^ortho * u_br^ortho */
    /* ----------------------------------------------- */
#pragma omp parallel num_threads(num_omp_threads_)
    {
        /* Inside Black Section */
#pragma omp for
        for (int i_theta = 0; i_theta < num_radial_lines; i_theta += 2) {
            applyAscOrthoRadialSection(i_theta, SmootherColor::Black, x, rhs, temp);
        }
        /* Outside Black Section (Part 1) */
#pragma omp for
        for (int i_theta = 1; i_theta < num_radial_lines; i_theta += 4) {
            applyAscOrthoRadialSection(i_theta, SmootherColor::Black, x, rhs, temp);
        }
        /* Outside Black Section (Part 2) */
#pragma omp for
        for (int i_theta = 3; i_theta < num_radial_lines; i_theta += 4) {
            applyAscOrthoRadialSection(i_theta, SmootherColor::Black, x, rhs, temp);
        }
    }
    solveBlackRadialSection(x, temp);

    /* ----------------------------------------------- */
    /* 4. White-Radial update (u_wr):                  */
    /*    A_wr * u_wr = f_wr − A_wr^ortho * u_wr^ortho */
    /* ----------------------------------------------- */
#pragma omp parallel num_threads(num_omp_threads_)
    {
        /* Inside Black Section */
#pragma omp for
        for (int i_theta = 1; i_theta < num_radial_lines; i_theta += 2) {
            applyAscOrthoRadialSection(i_theta, SmootherColor::White, x, rhs, temp);
        }
        /* Outside Black Section (Part 1) */
#pragma omp for
        for (int i_theta = 0; i_theta < num_radial_lines; i_theta += 4) {
            applyAscOrthoRadialSection(i_theta, SmootherColor::White, x, rhs, temp);
        }
        /* Outside Black Section (Part 2) */
#pragma omp for
        for (int i_theta = 2; i_theta < num_radial_lines; i_theta += 4) {
            applyAscOrthoRadialSection(i_theta, SmootherColor::White, x, rhs, temp);
        }
    }
    solveWhiteRadialSection(x, temp);
}
