#include <gtest/gtest.h>

#include <vector>
#include <random>

#include "../../include/GMGPolar/gmgpolar.h"

#include "../../include/Residual/ResidualGive/residualGive.h"
#include "../../include/Residual/ResidualTake/residualTake.h"

#include "../../include/InputFunctions/domainGeometry.h"
#include "../../include/InputFunctions/densityProfileCoefficients.h"
#include "../../include/InputFunctions/boundaryConditions.h"
#include "../../include/InputFunctions/sourceTerm.h"
/* --------- */
/* Test Case */
/* --------- */
#include "../include/InputFunctions/DomainGeometry/czarnyGeometry.h"
#include "../include/InputFunctions/BoundaryConditions/polarR6_Boundary_CzarnyGeometry.h"
#include "../include/InputFunctions/DensityProfileCoefficients/zoniShiftedCoefficients.h"
#include "../include/InputFunctions/SourceTerms/polarR6_ZoniShifted_CzarnyGeometry.h"

namespace ResidualTest
{
Vector<double> generate_random_sample_data(const PolarGrid& grid, unsigned int seed)
{
    Vector<double> x(grid.numberOfNodes());
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-100.0, 100.0);
    for (int i = 0; i < x.size(); ++i) {
        x[i] = dist(gen);
    }
    return x;
}
} // namespace ResidualTest

using namespace ResidualTest;

/* Test 1/1: */
/* Does the Take and Give Implementation match up? */

TEST(OperatorATest, applyA_DirBC_Interior)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = true;
    int maxOpenMPThreads = 16;

    // "Take" requires cached values
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, false);

    ResidualGive residualGive_operator(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualTake residualTake_operator(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);

    Vector<double> x   = generate_random_sample_data(level.grid(), 42);
    Vector<double> rhs = generate_random_sample_data(level.grid(), 69);

    Vector<double> result_Give(level.grid().numberOfNodes());
    residualGive_operator.computeResidual(result_Give, rhs, x);

    Vector<double> result_Take(level.grid().numberOfNodes());
    residualTake_operator.computeResidual(result_Take, rhs, x);

    ASSERT_EQ(result_Give.size(), result_Take.size());
    for (int index = 0; index < result_Give.size(); index++) {
        MultiIndex alpha = level.grid().multiIndex(index);
        if (alpha[0] == 0 && !DirBC_Interior)
            ASSERT_NEAR(result_Give[index], result_Take[index], 1e-8);
        else
            ASSERT_NEAR(result_Give[index], result_Take[index], 1e-11);
    }
}

TEST(OperatorATest, applyA_AcrossOrigin)
{
    std::vector<double> radii  = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

    double Rmax      = radii.back();
    double kappa_eps = 0.3;
    double delta_e   = 1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081 * Rmax;
    std::unique_ptr<DensityProfileCoefficients> coefficients =
        std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions =
        std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term =
        std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    bool DirBC_Interior  = false;
    int maxOpenMPThreads = 16;

    // "Take" requires cached values
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid>(radii, angles);
    auto levelCache = std::make_unique<LevelCache>(*grid, *coefficients, domain_geometry,
                                                   cache_density_rpofile_coefficients, cache_domain_geometry);
    Level level(0, std::move(grid), std::move(levelCache), ExtrapolationType::NONE, false);

    ResidualGive residualGive_operator(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);
    ResidualTake residualTake_operator(level.grid(), level.levelCache(), domain_geometry, *coefficients, DirBC_Interior,
                                       maxOpenMPThreads);

    Vector<double> x   = generate_random_sample_data(level.grid(), 42);
    Vector<double> rhs = generate_random_sample_data(level.grid(), 69);

    Vector<double> result_Give(level.grid().numberOfNodes());
    residualGive_operator.computeResidual(result_Give, rhs, x);

    Vector<double> result_Take(level.grid().numberOfNodes());
    residualTake_operator.computeResidual(result_Take, rhs, x);

    ASSERT_EQ(result_Give.size(), result_Take.size());
    for (int index = 0; index < result_Give.size(); index++) {
        MultiIndex alpha = level.grid().multiIndex(index);
        if (alpha[0] == 0 && !DirBC_Interior)
            ASSERT_NEAR(result_Give[index], result_Take[index], 1e-8);
        else
            ASSERT_NEAR(result_Give[index], result_Take[index], 1e-11);
    }
}

int main(int argc, char* argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}