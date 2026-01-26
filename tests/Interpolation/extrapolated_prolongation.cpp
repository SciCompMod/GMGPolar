#include <gtest/gtest.h>
#include <random>

#include "../test_tools.h"

#include "../../include/GMGPolar/gmgpolar.h"
#include "../../include/Interpolation/interpolation.h"
#include "../../include/InputFunctions/DensityProfileCoefficients/poissonCoefficients.h"

// Helper that computes the mathematically expected extrapolated prolongation value
static double expected_extrapolated_value(const PolarGrid& coarse, const PolarGrid& fine,
                                          ConstVector<double> coarse_vals, int i_r, int i_theta)
{
    int i_r_coarse     = i_r / 2;
    int i_theta_coarse = i_theta / 2;

    bool r_even = (i_r % 2 == 0);
    bool t_even = (i_theta % 2 == 0);

    if (r_even && t_even) {
        // Node coincides with a coarse node
        return coarse_vals[coarse.index(i_r_coarse, i_theta_coarse)];
    }

    if (!r_even && t_even) {
        // Radial midpoint - arithmetic mean of left and right
        return 0.5 * (coarse_vals[coarse.index(i_r_coarse, i_theta_coarse)] +
                      coarse_vals[coarse.index(i_r_coarse + 1, i_theta_coarse)]);
    }

    if (r_even && !t_even) {
        // Angular midpoint - arithmetic mean of bottom and top
        return 0.5 * (coarse_vals[coarse.index(i_r_coarse, i_theta_coarse)] +
                      coarse_vals[coarse.index(i_r_coarse, i_theta_coarse + 1)]);
    }

    // Center of coarse cell - arithmetic mean of diagonal nodes (bottom-right + top-left)
    return 0.5 * (coarse_vals[coarse.index(i_r_coarse + 1, i_theta_coarse)] +
                  coarse_vals[coarse.index(i_r_coarse, i_theta_coarse + 1)]);
}

TEST(ExtrapolatedProlongationTest, ExtrapolatedProlongationMatchesStencil)
{
    std::vector<double> fine_radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> fine_angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, 2 * M_PI};

    PolarGrid fine_grid(fine_radii, fine_angles);
    PolarGrid coarse_grid = coarseningGrid(fine_grid);

    Interpolation I(/*threads*/ 16, /*DirBC*/ true);

    Vector<double> coarse_values = generate_random_sample_data(coarse_grid, 1234, 0.0, 1.0);
    Vector<double> fine_result("fine_result", fine_grid.numberOfNodes());

    I.applyExtrapolatedProlongation(coarse_grid, fine_grid, fine_result, coarse_values);

    for (int i_r = 0; i_r < fine_grid.nr(); ++i_r) {
        for (int i_theta = 0; i_theta < fine_grid.ntheta(); ++i_theta) {
            double expected = expected_extrapolated_value(coarse_grid, fine_grid, coarse_values, i_r, i_theta);
            double got      = fine_result[fine_grid.index(i_r, i_theta)];
            ASSERT_NEAR(expected, got, 1e-10) << "Mismatch at (" << i_r << ", " << i_theta << ")";
        }
    }
}
