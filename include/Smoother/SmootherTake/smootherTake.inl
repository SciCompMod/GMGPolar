#pragma once

template <class LevelCacheType>
SmootherTake<LevelCacheType>::SmootherTake(const PolarGrid& grid, const LevelCacheType& level_cache,
                                           bool DirBC_Interior, int num_omp_threads)
    : Smoother<LevelCacheType>(grid, level_cache, DirBC_Interior, num_omp_threads)
    , circle_tridiagonal_solver_(grid.ntheta(), grid.numberSmootherCircles(), true)
    , radial_tridiagonal_solver_(grid.lengthRadialSmoother(), grid.ntheta(), false)
#ifdef GMGPOLAR_USE_MUMPS
    , inner_boundary_solver_(buildInteriorBoundarySolverMatrix())
#else
    , inner_boundary_circle_matrix_(buildInteriorBoundarySolverMatrix())
    , inner_boundary_solver_(inner_boundary_circle_matrix_)
#endif
{
    buildTridiagonalSolverMatrices();
    circle_tridiagonal_solver_.setup();
    radial_tridiagonal_solver_.setup();
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
template <class LevelCacheType>
void SmootherTake<LevelCacheType>::smoothing(Vector<double> x, ConstVector<double> rhs, Vector<double> temp)
{
    assert(x.size() == rhs.size());
    assert(temp.size() == rhs.size());

    private_applyAscOrthoBlackCircleSection(x, rhs, temp);
    solveBlackCircleSection(x, temp);

    private_applyAscOrthoWhiteCircleSection(x, rhs, temp);
    solveWhiteCircleSection(x, temp);

    private_applyAscOrthoBlackRadialSection(x, rhs, temp);
    solveBlackRadialSection(x, temp);

    private_applyAscOrthoWhiteRadialSection(x, rhs, temp);
    solveWhiteRadialSection(x, temp);
}
