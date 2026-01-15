#include <gtest/gtest.h>

#include <random>

#include "../test_tools.h"

#include "../../include/GMGPolar/gmgpolar.h"
#include "../../include/Interpolation/interpolation.h"
#include "../../include/InputFunctions/domainGeometry.h"
#include "../../include/InputFunctions/densityProfileCoefficients.h"

#include "../../include/InputFunctions/DomainGeometry/circularGeometry.h"
#include "../../include/InputFunctions/DensityProfileCoefficients/poissonCoefficients.h"

TEST(ExtrapolatedProlongationTest, ExtrapolatedProlongationSmoothingRadius)
{
    std::vector<double> fine_radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> fine_angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    int maxOpenMPThreads = 16;
    bool DirBC_Interior  = true;

    PolarGrid finest_grid(fine_radii, fine_angles);
    PolarGrid coarse_grid = coarseningGrid(finest_grid);

    Interpolation interpolation_operator(maxOpenMPThreads, DirBC_Interior);

    unsigned int seed = 42;
    Vector<double> x  = generate_random_sample_data(coarse_grid, seed, 0.0, 1.0);

    // Apply prolongation to both functions
    Vector<double> result1("result1", finest_grid.numberOfNodes());
    Vector<double> result2("result2", finest_grid.numberOfNodes());

    interpolation_operator.applyExtrapolatedProlongation0(coarse_grid, finest_grid, result1, x);
    interpolation_operator.applyExtrapolatedProlongation(coarse_grid, finest_grid, result2, x);

    ASSERT_EQ(result1.size(), result2.size());
    for (uint i = 0; i < result1.size(); ++i) {
        ASSERT_DOUBLE_EQ(result1[i], result2[i]);
    }
}
