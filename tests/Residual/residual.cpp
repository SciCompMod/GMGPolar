#include <gtest/gtest.h>

#include "../../include/GMGPolar/gmgpolar.h"

#include "../../include/InputFunctions/domainGeometry.h"
#include "../../include/InputFunctions/systemParameters.h"
/* --------- */
/* Test Case */
/* --------- */
#include "../include/InputFunctions/DomainGeometry/czarnyGeometry.h"
#include "../include/InputFunctions/BoundaryConditions/polarR6_Boundary_CzarnyGeometry.h"
#include "../include/InputFunctions/DensityProfileCoefficients/zoniShiftedCoefficients.h"
#include "../include/InputFunctions/SourceTerms/polarR6_ZoniShifted_CzarnyGeometry.h"

#include <random>

// Function to generate sample data for vector x using random values with seed
Vector<double> generate_random_sample_data_Residual(const PolarGrid& grid, unsigned int seed) {
    Vector<double> x(grid.number_of_nodes());
    std::mt19937 gen(seed);  // Standard mersenne_twister_engine seeded with seed
    std::uniform_real_distribution<double> dist(-100.0, 100.0); 
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = dist(gen);
    }
    return x;
}

TEST(OperatorATest, applyA_DirBC_Interior) {
    std::vector<double> radii = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    auto grid = std::make_unique<PolarGrid>(radii, angles);

    auto levelCache = std::make_unique<LevelCache>(*grid);

    int extrapolation = 0;
    Level level(0, std::move(grid), std::move(levelCache), extrapolation);

    double Rmax = radii.back();
    double kappa_eps=0.3;
    double delta_e=1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    SystemParameters system_parameters(std::move(coefficients), std::move(boundary_conditions), std::move(source_term));

    bool DirBC_Interior = true;
    int maxOpenMPThreads = 1;
    int openMPTaskThreads = 1;
    
    Residual residual_operator(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);

    unsigned int seed = 42;
    Vector<double> x = generate_random_sample_data_Residual(level.grid(), seed);
    Vector<double> rhs = generate_random_sample_data_Residual(level.grid(), seed);

    Vector<double> result1(level.grid().number_of_nodes());
    Vector<double> result2(level.grid().number_of_nodes());
    Vector<double> result3(level.grid().number_of_nodes());

    residual_operator.computeResidual_V1(result1, rhs, x);
    residual_operator.computeResidual_V2(result2, rhs, x);
    residual_operator.computeResidual_V3(result3, rhs, x);

    ASSERT_EQ(result1.size(), result2.size());
    for (size_t i = 0; i < result1.size(); ++i) {
        ASSERT_NEAR(result1[i], result2[i], 1e-10);
    }

    ASSERT_EQ(result2.size(), result3.size());
    for (size_t i = 0; i < result2.size(); ++i) {
        ASSERT_NEAR(result2[i], result3[i], 1e-10);
    }

    Vector<double> result4(level.grid().number_of_nodes());
    result4 = rhs;
    residual_operator.applyATake0(result4, x, -1.0);

    ASSERT_EQ(result2.size(), result3.size());
    for (size_t i = 0; i < result1.size(); ++i) {
        MultiIndex alpha = level.grid().multiindex(i);
        if(alpha[0] == 0 && !DirBC_Interior) ASSERT_NEAR(result1[i], result4[i], 1e-8);
        else ASSERT_NEAR(result1[i], result4[i], 1e-10);
    }
}


TEST(OperatorATest, applyA_AcrossOrigin) {
    std::vector<double> radii = {1e-5, 0.2, 0.25, 0.5, 0.8, 0.9, 0.95, 1.2, 1.3};
    std::vector<double> angles = {0, M_PI/16, M_PI/8, M_PI/2, M_PI, M_PI+M_PI/16, M_PI+M_PI/8, M_PI+M_PI/2, M_PI+M_PI};

    auto grid = std::make_unique<PolarGrid>(radii, angles);

    auto levelCache = std::make_unique<LevelCache>(*grid);

    int extrapolation = 0;
    Level level(0, std::move(grid), std::move(levelCache), extrapolation);

    double Rmax = radii.back();
    double kappa_eps=0.3;
    double delta_e=1.4;

    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    double alpha_jump = 0.7081;
    std::unique_ptr<DensityProfileCoefficients> coefficients = std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
    std::unique_ptr<BoundaryConditions> boundary_conditions = std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
    std::unique_ptr<SourceTerm> source_term = std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);

    SystemParameters system_parameters(std::move(coefficients), std::move(boundary_conditions), std::move(source_term));

    bool DirBC_Interior = false;
    int maxOpenMPThreads = 1;
    int openMPTaskThreads = 1;
    
    Residual residual_operator(level.grid(), level.levelCache(), domain_geometry, system_parameters, DirBC_Interior, maxOpenMPThreads, openMPTaskThreads);

    unsigned int seed = 42;
    Vector<double> x = generate_random_sample_data_Residual(level.grid(), seed);
    Vector<double> rhs = generate_random_sample_data_Residual(level.grid(), seed);

    Vector<double> result1(level.grid().number_of_nodes());
    Vector<double> result2(level.grid().number_of_nodes());
    Vector<double> result3(level.grid().number_of_nodes());

    residual_operator.computeResidual_V1(result1, rhs, x);
    residual_operator.computeResidual_V2(result2, rhs, x);
    residual_operator.computeResidual_V3(result3, rhs, x);

    ASSERT_EQ(result1.size(), result2.size());
    for (size_t i = 0; i < result1.size(); ++i) {
        ASSERT_NEAR(result1[i], result2[i], 1e-12);
    }

    ASSERT_EQ(result2.size(), result3.size());
    for (size_t i = 0; i < result2.size(); ++i) {
        ASSERT_NEAR(result2[i], result3[i], 1e-12);
    }

    Vector<double> result4(level.grid().number_of_nodes());
    result4 = rhs;
    residual_operator.applyATake0(result4, x, -1.0);

    ASSERT_EQ(result2.size(), result3.size());
    for (size_t i = 0; i < result1.size(); ++i) {
        MultiIndex alpha = level.grid().multiindex(i);
        if(alpha[0] == 0 && !DirBC_Interior) ASSERT_NEAR(result1[i], result4[i], 1e-8);
        else ASSERT_NEAR(result1[i], result4[i], 1e-11);
    }
}

int main(int argc, char* argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}