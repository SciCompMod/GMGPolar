#include <gtest/gtest.h>
#include <random>

#include "../../include/Interpolation/interpolation.h"
#include "../../include/InputFunctions/DensityProfileCoefficients/poissonCoefficients.h"

static Vector<double> generate_random_sample_data(const PolarGrid& grid, unsigned int seed)
{
    Vector<double> vec("vec", grid.numberOfNodes());
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-20.0, 20.0);

    for (uint i = 0; i < vec.size(); ++i)
        vec[i] = dist(gen);

    return vec;
}

// Helper that computes the mathematically expected prolongation value
static double expected_value(const PolarGrid& coarse, const PolarGrid& fine, ConstVector<double> coarse_vals, int i_r,
                             int i_theta)
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
        // Radial midpoint
        double h1 = fine.radialSpacing(i_r - 1);
        double h2 = fine.radialSpacing(i_r);

        return (h1 * coarse_vals[coarse.index(i_r_coarse, i_theta_coarse)] +
                h2 * coarse_vals[coarse.index(i_r_coarse + 1, i_theta_coarse)]) /
               (h1 + h2);
    }

    if (r_even && !t_even) {
        // Angular midpoint
        double k1 = fine.angularSpacing(i_theta - 1);
        double k2 = fine.angularSpacing(i_theta);

        return (k1 * coarse_vals[coarse.index(i_r_coarse, t_even ? i_theta_coarse : i_theta_coarse)] +
                k2 * coarse_vals[coarse.index(i_r_coarse, i_theta_coarse + 1)]) /
               (k1 + k2);
    }

    // Center of coarse cell
    double h1 = fine.radialSpacing(i_r - 1);
    double h2 = fine.radialSpacing(i_r);
    double k1 = fine.angularSpacing(i_theta - 1);
    double k2 = fine.angularSpacing(i_theta);

    return (h1 * k1 * coarse_vals[coarse.index(i_r_coarse, i_theta_coarse)] +
            h2 * k1 * coarse_vals[coarse.index(i_r_coarse + 1, i_theta_coarse)] +
            h1 * k2 * coarse_vals[coarse.index(i_r_coarse, i_theta_coarse + 1)] +
            h2 * k2 * coarse_vals[coarse.index(i_r_coarse + 1, i_theta_coarse + 1)]) /
           ((h1 + h2) * (k1 + k2));
}

TEST(ProlongationTest, ProlongationMatchesStencil)
{
    std::vector<double> fine_radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> fine_angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, 2 * M_PI};

    PolarGrid fine_grid(fine_radii, fine_angles);
    PolarGrid coarse_grid = coarseningGrid(fine_grid);

    Interpolation I(/*threads*/ 16, /*DirBC*/ true);

    Vector<double> coarse_values = generate_random_sample_data(coarse_grid, 1234);
    Vector<double> fine_result("fine_result", fine_grid.numberOfNodes());

    I.applyProlongation(coarse_grid, fine_grid, fine_result, coarse_values);

    for (int i_r = 0; i_r < fine_grid.nr(); ++i_r) {
        for (int i_theta = 0; i_theta < fine_grid.ntheta(); ++i_theta) {
            double expected = expected_value(coarse_grid, fine_grid, coarse_values, i_r, i_theta);
            double got      = fine_result[fine_grid.index(i_r, i_theta)];
            ASSERT_NEAR(expected, got, 1e-10) << "Mismatch at (" << i_r << ", " << i_theta << ")";
        }
    }
}
