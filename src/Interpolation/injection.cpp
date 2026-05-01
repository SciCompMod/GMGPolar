#include <Interpolation/interpolation.h>
using namespace gmgpolar;

/* Remark: This injection is not scaled. */
static inline void coarseNodeInjection(const int i_r_coarse, const int i_theta_coarse, const PolarGrid& fine_grid,
                                       const PolarGrid& coarse_grid, Vector<double>& coarse_result,
                                       ConstVector<double>& fine_values)
{
    const int i_r_fine     = i_r_coarse * 2;
    const int i_theta_fine = i_theta_coarse * 2;

    const int coarse_index = coarse_grid.index(i_r_coarse, i_theta_coarse);
    const int fine_index   = fine_grid.index(i_r_fine, i_theta_fine);

    coarse_result[coarse_index] = fine_values[fine_index];
}

void Interpolation::applyInjection(const PolarGrid& fine_grid, const PolarGrid& coarse_grid,
                                   Vector<double> coarse_result, ConstVector<double> fine_values) const
{
    assert(std::ssize(fine_values) == fine_grid.numberOfNodes());
    assert(std::ssize(coarse_result) == coarse_grid.numberOfNodes());

    // The For loop matches circular access pattern */
    Kokkos::parallel_for(
        "Interpolation: Injection (Circular)",
        // Rank of the index space, Iteration pattern between tiles, Iteration pattern within tiles
        Kokkos::MDRangePolicy<Kokkos::Rank<2, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>( // Iteration policy
            {0, 0}, // Starting point of the index space
            {coarse_grid.numberSmootherCircles(), coarse_grid.ntheta()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_r_coarse, const int i_theta_coarse) {
            coarseNodeInjection(i_r_coarse, i_theta_coarse, fine_grid, coarse_grid, coarse_result, fine_values);
        });

    /* For loop matches radial access pattern */
    Kokkos::parallel_for(
        "Interpolation: Injection (Radial)",
        // Rank of the index space, Iteration pattern between tiles, Iteration pattern within tiles
        Kokkos::MDRangePolicy<Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>( // Iteration policy
            {coarse_grid.numberSmootherCircles(), 0}, // Starting point of the index space
            {coarse_grid.nr(), coarse_grid.ntheta()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_r_coarse, const int i_theta_coarse) {
            coarseNodeInjection(i_r_coarse, i_theta_coarse, fine_grid, coarse_grid, coarse_result, fine_values);
        });

    Kokkos::fence();
}
