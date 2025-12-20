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

// Helper that computes the mathematically expected restriction value
static double expected_restriction_value(const PolarGrid& fine, const PolarGrid& coarse, ConstVector<double> fine_vals,
                                         int i_r_coarse, int i_theta_coarse)
{
    int i_r     = i_r_coarse * 2;
    int i_theta = i_theta_coarse * 2;

    // Angular indices with periodic wrapping
    int i_theta_M2 = fine.wrapThetaIndex(i_theta - 2);
    int i_theta_M1 = fine.wrapThetaIndex(i_theta - 1);
    int i_theta_P1 = fine.wrapThetaIndex(i_theta + 1);

    // Angular spacings
    double k1 = fine.angularSpacing(i_theta_M2);
    double k2 = fine.angularSpacing(i_theta_M1);
    double k3 = fine.angularSpacing(i_theta);
    double k4 = fine.angularSpacing(i_theta_P1);

    // Center + Angular contributions (always present)
    double value = fine_vals[fine.index(i_r, i_theta)] + k2 / (k1 + k2) * fine_vals[fine.index(i_r, i_theta_M1)] +
                   k3 / (k3 + k4) * fine_vals[fine.index(i_r, i_theta_P1)];

    // Left contributions (if not at inner boundary)
    if (i_r_coarse > 0) {
        double h1 = fine.radialSpacing(i_r - 2);
        double h2 = fine.radialSpacing(i_r - 1);
        value += h2 / (h1 + h2) * fine_vals[fine.index(i_r - 1, i_theta)] +
                 h2 * k2 / ((h1 + h2) * (k1 + k2)) * fine_vals[fine.index(i_r - 1, i_theta_M1)] +
                 h2 * k3 / ((h1 + h2) * (k3 + k4)) * fine_vals[fine.index(i_r - 1, i_theta_P1)];
    }

    // Right contributions (if not at outer boundary)
    if (i_r_coarse < coarse.nr() - 1) {
        double h3 = fine.radialSpacing(i_r);
        double h4 = fine.radialSpacing(i_r + 1);
        value += h3 / (h3 + h4) * fine_vals[fine.index(i_r + 1, i_theta)] +
                 h3 * k2 / ((h3 + h4) * (k1 + k2)) * fine_vals[fine.index(i_r + 1, i_theta_M1)] +
                 h3 * k3 / ((h3 + h4) * (k3 + k4)) * fine_vals[fine.index(i_r + 1, i_theta_P1)];
    }

    return value;
}

TEST(RestrictionTest, RestrictionMatchesStencil)
{
    std::vector<double> fine_radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> fine_angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, 2 * M_PI};

    PolarGrid fine_grid(fine_radii, fine_angles);
    PolarGrid coarse_grid = coarseningGrid(fine_grid);

    Interpolation I(/*threads*/ 16, /*DirBC*/ true);

    Vector<double> fine_values = generate_random_sample_data(fine_grid, 5678);
    Vector<double> coarse_result("coarse_result", coarse_grid.numberOfNodes());

    I.applyRestriction(fine_grid, coarse_grid, coarse_result, fine_values);

    for (int i_r_coarse = 0; i_r_coarse < coarse_grid.nr(); ++i_r_coarse) {
        for (int i_theta_coarse = 0; i_theta_coarse < coarse_grid.ntheta(); ++i_theta_coarse) {
            double expected =
                expected_restriction_value(fine_grid, coarse_grid, fine_values, i_r_coarse, i_theta_coarse);
            double got = coarse_result[coarse_grid.index(i_r_coarse, i_theta_coarse)];
            ASSERT_NEAR(expected, got, 1e-10) << "Mismatch at (" << i_r_coarse << ", " << i_theta_coarse << ")";
        }
    }
}
