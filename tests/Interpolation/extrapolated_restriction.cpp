#include <gtest/gtest.h>
#include <random>

#include "../test_tools.h"

#include <GMGPolar/gmgpolar.h>
#include <Interpolation/interpolation.h>
#include <InputFunctions/DensityProfileCoefficients/poissonCoefficients.h>
using namespace gmgpolar;

// Helper that computes the mathematically expected extrapolated restriction value
static double expected_extrapolated_restriction_value(const PolarGrid<DefaultMemorySpace>& fine,
                                                      const PolarGrid<DefaultMemorySpace>& coarse,
                                                      ConstVector<double> fine_vals, int i_r_coarse, int i_theta_coarse)
{
    int i_r     = i_r_coarse * 2;
    int i_theta = i_theta_coarse * 2;

    double result = 0;
    Kokkos::parallel_reduce(
        "extrap_restriction_test", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, 1),
        KOKKOS_LAMBDA(const int idx, double& local_result) {
            // Angular indices with periodic wrapping
            int i_theta_M1 = fine.wrapThetaIndex(i_theta - 1);
            int i_theta_P1 = fine.wrapThetaIndex(i_theta + 1);

            // Center + Angular contributions (always present)
            local_result = fine_vals[fine.index(i_r, i_theta)] + 0.5 * fine_vals[fine.index(i_r, i_theta_M1)] +
                           0.5 * fine_vals[fine.index(i_r, i_theta_P1)];

            // Left contributions (if not at inner boundary)
            if (i_r_coarse > 0) {
                local_result += 0.5 * fine_vals[fine.index(i_r - 1, i_theta)] +
                                0.5 * fine_vals[fine.index(i_r - 1, i_theta_P1)]; // Top-Left diagonal
            }

            // Right contributions (if not at outer boundary)
            if (i_r_coarse < coarse.nr() - 1) {
                local_result += 0.5 * fine_vals[fine.index(i_r + 1, i_theta)] +
                                0.5 * fine_vals[fine.index(i_r + 1, i_theta_M1)]; // Bottom-Right diagonal
            }
        },
        result);

    return result;
}

TEST(ExtrapolatedRestrictionTest, ExtrapolatedRestrictionMatchesStencil)
{
    std::vector<double> fine_radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> fine_angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, 2 * M_PI};

    PolarGrid<DefaultMemorySpace> fine_grid(fine_radii, fine_angles);
    PolarGrid<DefaultMemorySpace> coarse_grid = coarseningGrid(fine_grid);

    Interpolation I(/*DirBC*/ true);

    Vector<double> fine_values = generate_random_sample_data(fine_grid, 9012, 0.0, 1.0);
    Vector<double> coarse_result("coarse_result", coarse_grid.numberOfNodes());

    I.applyExtrapolatedRestriction(fine_grid, coarse_grid, coarse_result, fine_values);

    auto h_coarse_result = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), coarse_result);

    for (int i_r_coarse = 0; i_r_coarse < coarse_grid.nr(); ++i_r_coarse) {
        for (int i_theta_coarse = 0; i_theta_coarse < coarse_grid.ntheta(); ++i_theta_coarse) {
            double expected = expected_extrapolated_restriction_value(fine_grid, coarse_grid, fine_values, i_r_coarse,
                                                                      i_theta_coarse);
            double got      = h_coarse_result[coarse_grid.index(i_r_coarse, i_theta_coarse)];
            ASSERT_NEAR(expected, got, 1e-10) << "Mismatch at (" << i_r_coarse << ", " << i_theta_coarse << ")";
        }
    }
}
