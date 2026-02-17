#include "../../../include/Smoother/SmootherTake/smootherTake.h"

SmootherTake::SmootherTake(const PolarGrid& grid, const LevelCache& level_cache, const DomainGeometry& domain_geometry,
                           const DensityProfileCoefficients& density_profile_coefficients, bool DirBC_Interior,
                           int num_omp_threads)
    : Smoother(grid, level_cache, domain_geometry, density_profile_coefficients, DirBC_Interior, num_omp_threads)
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

SmootherTake::~SmootherTake()
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
void SmootherTake::smoothing(Vector<double> x, ConstVector<double> rhs, Vector<double> temp)
{
    assert(x.size() == rhs.size());
    assert(temp.size() == rhs.size());

    assert(level_cache_.cacheDensityProfileCoefficients());
    assert(level_cache_.cacheDomainGeometry());

    /* The outer most circle next to the radial section is defined to be black. */
    /* Priority: Black -> White. */
    const int start_black_circles = (grid_.numberSmootherCircles() % 2 == 0) ? 1 : 0;
    const int start_white_circles = (grid_.numberSmootherCircles() % 2 == 0) ? 0 : 1;

    /* Black Circle Section */
#pragma omp parallel for num_threads(num_omp_threads_)
    for (int i_r = start_black_circles; i_r < grid_.numberSmootherCircles(); i_r += 2) {
        applyAscOrthoCircleSection(i_r, x, rhs, temp);
    } /* Implicit barrier */

    solveBlackCircleSection(x, temp);

    /* White Circle Section */
#pragma omp parallel for num_threads(num_omp_threads_)
    for (int i_r = start_white_circles; i_r < grid_.numberSmootherCircles(); i_r += 2) {
        applyAscOrthoCircleSection(i_r, x, rhs, temp);
    } /* Implicit barrier */

    solveWhiteCircleSection(x, temp);

    /* Black Radial Section */
#pragma omp parallel for num_threads(num_omp_threads_)
    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta += 2) {
        applyAscOrthoRadialSection(i_theta, x, rhs, temp);
    } /* Implicit barrier */

    solveBlackRadialSection(x, temp);

    /* White Radial Section*/
#pragma omp parallel for num_threads(num_omp_threads_)
    for (int i_theta = 1; i_theta < grid_.ntheta(); i_theta += 2) {
        applyAscOrthoRadialSection(i_theta, x, rhs, temp);
    } /* Implicit barrier */

    solveWhiteRadialSection(x, temp);
}
