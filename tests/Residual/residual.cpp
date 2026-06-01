#include <gtest/gtest.h>

#include <vector>
#include <random>

#include "../test_tools.h"

#include "../../include/GMGPolar/gmgpolar.h"

#include "../../include/Residual/ResidualGive/residualGive.h"
#include "../../include/Residual/ResidualTake/residualTake.h"

#include "../../include/InputFunctions/domainGeometry.h"
#include "../../include/InputFunctions/densityProfileCoefficients.h"
#include "../../include/InputFunctions/boundaryConditions.h"
/* --------- */
/* Test Case */
/* --------- */
#include "../include/InputFunctions/DomainGeometry/czarnyGeometry.h"
#include "../include/InputFunctions/DensityProfileCoefficients/zoniShiftedCoefficients.h"
using namespace gmgpolar;

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

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    ZoniShiftedCoefficients coefficients(Rmax, alpha_jump);

    bool DirBC_Interior = true;

    // "Take" requires cached values
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid<DefaultMemorySpace>>(radii, angles);
    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(grid), std::move(levelCache),
                                                                    ExtrapolationType::NONE, false);

    ResidualGive residualGive_operator(level.grid(), level.levelCache(), DirBC_Interior);
    ResidualTake residualTake_operator(level.grid(), level.levelCache(), DirBC_Interior);

    HostVector<double> x   = generate_random_sample_data(PolarGrid<Kokkos::HostSpace>(level.grid()), 42);
    HostVector<double> rhs = generate_random_sample_data(PolarGrid<Kokkos::HostSpace>(level.grid()), 69);

    HostVector<double> result_Give("result_Give", level.grid().numberOfNodes());
    residualGive_operator.computeResidual(result_Give, rhs, x);

    HostVector<double> result_Take("result_Take", level.grid().numberOfNodes());
    residualTake_operator.computeResidual(result_Take, rhs, x);

    ASSERT_EQ(result_Give.size(), result_Take.size());
    for (uint index = 0; index < result_Give.size(); index++) {
        int i_r, i_theta;
        level.grid().multiIndex(index, i_r, i_theta);
        if (i_r == 0 && !DirBC_Interior)
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

    using DomainGeometryType = CzarnyGeometry;
    DomainGeometryType domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump                    = 0.678 * Rmax;
    using DensityProfileCoefficientsType = ZoniShiftedCoefficients;
    DensityProfileCoefficientsType coefficients(Rmax, alpha_jump);

    bool DirBC_Interior = false;

    // "Take" requires cached values
    bool cache_density_rpofile_coefficients = true;
    bool cache_domain_geometry              = true;

    auto grid       = std::make_unique<PolarGrid<DefaultMemorySpace>>(radii, angles);
    auto levelCache = std::make_unique<LevelCache<DomainGeometryType, DensityProfileCoefficientsType>>(
        *grid, coefficients, domain_geometry, cache_density_rpofile_coefficients, cache_domain_geometry);
    Level<DomainGeometryType, DensityProfileCoefficientsType> level(0, std::move(grid), std::move(levelCache),
                                                                    ExtrapolationType::NONE, false);

    ResidualGive residualGive_operator(level.grid(), level.levelCache(), DirBC_Interior);
    ResidualTake residualTake_operator(level.grid(), level.levelCache(), DirBC_Interior);

    HostVector<double> x   = generate_random_sample_data(PolarGrid<Kokkos::HostSpace>(level.grid()), 42);
    HostVector<double> rhs = generate_random_sample_data(PolarGrid<Kokkos::HostSpace>(level.grid()), 69);

    HostVector<double> result_Give("result_Give", level.grid().numberOfNodes());
    residualGive_operator.computeResidual(result_Give, rhs, x);

    HostVector<double> result_Take("result_Take", level.grid().numberOfNodes());
    residualTake_operator.computeResidual(result_Take, rhs, x);

    ASSERT_EQ(result_Give.size(), result_Take.size());
    for (uint index = 0; index < result_Give.size(); index++) {
        int i_r, i_theta;
        level.grid().multiIndex(index, i_r, i_theta);
        if (i_r == 0 && !DirBC_Interior)
            ASSERT_NEAR(result_Give[index], result_Take[index], 1e-8);
        else
            ASSERT_NEAR(result_Give[index], result_Take[index], 1e-11);
    }
}
