#include <gtest/gtest.h>

#include <random>

#include "../../include/GMGPolar/gmgpolar.h"
#include "../../include/Interpolation/interpolation.h"
#include "../../include/InputFunctions/DensityProfileCoefficients/poissonCoefficients.h"

namespace ProlongationTest
{
Vector<double> generate_random_sample_data(const PolarGrid& grid, unsigned int seed)
{
    Vector<double> vector("vector", grid.numberOfNodes());
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-100.0, 100.0);
    for (uint i = 0; i < vector.size(); ++i) {
        vector[i] = dist(gen);
    }
    return vector;
}
} // namespace ProlongationTest

using namespace ProlongationTest;

TEST(ProlongationTest, ProlongationTest)
{
    std::vector<double> fine_radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> fine_angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    int maxOpenMPThreads = 16;
    bool DirBC_Interior  = true;

    PolarGrid finest_grid(fine_radii, fine_angles);
    PolarGrid coarse_grid = coarseningGrid(finest_grid);

    Interpolation interpolation_operator(maxOpenMPThreads, DirBC_Interior);

    Vector<double> x = generate_random_sample_data(coarse_grid, 42);

    // Apply prolongation to both functions
    Vector<double> result1("result1", finest_grid.numberOfNodes());
    Vector<double> result2("result2", finest_grid.numberOfNodes());

    interpolation_operator.applyProlongation0(coarse_grid, finest_grid, result1, x);
    interpolation_operator.applyProlongation(coarse_grid, finest_grid, result2, x);

    ASSERT_EQ(result1.size(), result2.size());
    for (uint i = 0; i < result1.size(); ++i) {
        ASSERT_NEAR(result1[i], result2[i], 1e-10);
    }
}
