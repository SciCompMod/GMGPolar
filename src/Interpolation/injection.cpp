#include "../../include/Interpolation/interpolation.h"

using namespace gmgpolar;

/* Remark: This injection is not scaled. */
static inline void coarseNodeInjection(int i_r_coarse, int i_theta_coarse, const PolarGrid& fine_grid,
                                       const PolarGrid& coarse_grid, Vector<double>& coarse_result,
                                       ConstVector<double>& fine_values)
{
    const int i_r     = i_r_coarse * 2;
    const int i_theta = i_theta_coarse * 2;

    const int coarse_index = coarse_grid.index(i_r_coarse, i_theta_coarse);
    const int fine_index   = fine_grid.index(i_r, i_theta);

    coarse_result[coarse_index] = fine_values[fine_grid.index(i_r, i_theta)];
}

void Interpolation::applyInjection(const PolarGrid& fine_grid, const PolarGrid& coarse_grid,
                                   Vector<double> coarse_result, ConstVector<double> fine_values) const
{
    assert(std::ssize(fine_values) == fine_grid.numberOfNodes());
    assert(std::ssize(coarse_result) == coarse_grid.numberOfNodes());

    const int smoother_circles = coarse_grid.numberSmootherCircles();
    const int ntheta           = coarse_grid.ntheta();
    const int nr               = coarse_grid.nr();

    /* For loop matches circular access pattern */
    Kokkos::parallel_for(
        "injection_circular", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {smoother_circles, ntheta}),
        KOKKOS_LAMBDA(int i_r_coarse, int i_theta_coarse) {
            coarseNodeInjection(i_r_coarse, i_theta_coarse, fine_grid, coarse_grid, coarse_result, fine_values);
        });

    /* For loop matches radial access pattern */
    Kokkos::parallel_for(
        "injection_radial", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, smoother_circles}, {ntheta, nr}),
        KOKKOS_LAMBDA(int i_theta_coarse, int i_r_coarse) {
            coarseNodeInjection(i_r_coarse, i_theta_coarse, fine_grid, coarse_grid, coarse_result, fine_values);
        });
}
