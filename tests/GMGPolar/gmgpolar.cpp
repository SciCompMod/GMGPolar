#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>

#include "../../include/GMGPolar/gmgpolar.h"

double constexpr Rmax = 1.3;

template <class T>
class GMGPolarPaperTestCase;

template <class DensityProfileCoefficientsType, class BoundaryConditionsType, class SourceTermType,
          class ExactSolutionType>
class GMGPolarPaperTestCase<
    std::tuple<DensityProfileCoefficientsType, BoundaryConditionsType, SourceTermType, ExactSolutionType>>
    : public testing::Test
{
public:
    using DensityProfileCoefficients = DensityProfileCoefficientsType;
    using BoundaryConditions         = BoundaryConditionsType;
    using SourceTerm                 = SourceTermType;
    using ExactSolution              = ExactSolutionType;
};

using gmg_paper_types =
    testing::Types<std::tuple<ZoniShiftedGyroCoefficients, PolarR6_Boundary_CzarnyGeometry,
                              PolarR6_ZoniShiftedGyro_CzarnyGeometry, PolarR6_CzarnyGeometry>,
                   std::tuple<ZoniShiftedGyroCoefficients, CartesianR6_Boundary_CzarnyGeometry,
                              CartesianR6_ZoniShiftedGyro_CzarnyGeometry, CartesianR6_CzarnyGeometry>,
                   std::tuple<ZoniShiftedGyroCoefficients, Refined_Boundary_CzarnyGeometry,
                              Refined_ZoniShiftedGyro_CzarnyGeometry, Refined_CzarnyGeometry>>;

TYPED_TEST_SUITE(GMGPolarPaperTestCase, gmg_paper_types);

std::vector<double> get_non_uniform_points(double min, double max, int n_pts, double non_uniformity = 0.1)
{
    std::vector<double> points(n_pts);

    std::srand(std::time(nullptr)); // Seed with random value (the time)

    int n_cells = n_pts - 1;
    double const delta((max - min) / double(n_cells));

    points[0] = min;
    for (int i(1); i < n_cells; ++i) {
        double const random_perturbation = double(rand()) / RAND_MAX - 0.5;
        points[i]                        = min + (i + random_perturbation * non_uniformity) * delta;
    }
    points[n_cells] = max;

    return points;
}

std::tuple<double, double> get_gmgpolar_error(int n_r, int n_angles, double non_uniformity,
                                              CzarnyGeometry const& domain_geometry,
                                              DensityProfileCoefficients const& coefficients,
                                              BoundaryConditions const& boundary_conditions,
                                              SourceTerm const& source_term, ExactSolution const& solution)
{
    std::vector<double> radii  = get_non_uniform_points(0.05, Rmax, n_r - 1, non_uniformity); // remove central point
    std::vector<double> angles = get_non_uniform_points(0.0, M_PI, n_angles / 2 + 1, non_uniformity);
    // Every node in the interior ring needs to have an opposite neighboring node
    for (int i(1); i < n_angles / 2 + 1; ++i) {
        angles.push_back(angles[i] + M_PI);
    }
    PolarGrid grid(radii, angles);

    GMGPolar gmgpolar(grid, domain_geometry, coefficients);
    gmgpolar.setSolution(&solution);

    gmgpolar.DirBC_Interior(false); // Use across-origin calculation

    // ----------------------------------------------------------------
    // Parameters to be identified and values confirmed
    gmgpolar.FMG(true);
    gmgpolar.FMG_iterations(3);
    gmgpolar.FMG_cycle(MultigridCycleType::F_CYCLE);

    gmgpolar.extrapolation(ExtrapolationType::IMPLICIT_EXTRAPOLATION);
    gmgpolar.maxLevels(2);
    gmgpolar.preSmoothingSteps(1);
    gmgpolar.postSmoothingSteps(1);
    gmgpolar.multigridCycle(MultigridCycleType::F_CYCLE);

    gmgpolar.maxIterations(150);
    gmgpolar.residualNormType(ResidualNormType::EUCLIDEAN);
    gmgpolar.absoluteTolerance(1e-50);
    gmgpolar.relativeTolerance(1e-50);
    // ----------------------------------------------------------------

    // Perform setup and solve
    gmgpolar.setup();
    gmgpolar.solve(boundary_conditions, source_term);

    double euclid_error = *gmgpolar.exactErrorWeightedEuclidean();
    double inf_error    = *gmgpolar.exactErrorInfinity();
    return std::make_tuple(euclid_error, inf_error);
}

TYPED_TEST(GMGPolarPaperTestCase, TestUniform)
{
    int n_r      = 32;
    int n_angles = 64;

    double kappa_eps = 0.3;
    double delta_e   = 1.4;
    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    const double alpha_jump = 0.0; // Unused value
    typename TestFixture::DensityProfileCoefficients coefficients(Rmax, alpha_jump);
    typename TestFixture::BoundaryConditions boundary_conditions(Rmax, kappa_eps, delta_e);
    typename TestFixture::SourceTerm source_term(Rmax, kappa_eps, delta_e);

    typename TestFixture::ExactSolution solution(Rmax, kappa_eps, delta_e);

    auto [euclid_error, inf_error] = get_gmgpolar_error(n_r, n_angles, 0.0, domain_geometry, coefficients,
                                                        boundary_conditions, source_term, solution);
    auto [euclid_error_refined, inf_error_refined] = get_gmgpolar_error(
        n_r * 2, n_angles * 2, 0.0, domain_geometry, coefficients, boundary_conditions, source_term, solution);

    double euclid_order = log(euclid_error / euclid_error_refined) / log(2);
    double inf_order    = log(inf_error / inf_error_refined) / log(2);

    ASSERT_GT(euclid_order, 3);
    ASSERT_GT(inf_order, 3);
}

TYPED_TEST(GMGPolarPaperTestCase, TestNonUniform)
{
    int n_r      = 32;
    int n_angles = 64;

    double kappa_eps = 0.3;
    double delta_e   = 1.4;
    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    const double alpha_jump = 0.0; // Unused value
    typename TestFixture::DensityProfileCoefficients coefficients(Rmax, alpha_jump);
    typename TestFixture::BoundaryConditions boundary_conditions(Rmax, kappa_eps, delta_e);
    typename TestFixture::SourceTerm source_term(Rmax, kappa_eps, delta_e);

    typename TestFixture::ExactSolution solution(Rmax, kappa_eps, delta_e);

    double non_uniformity = 0.1;

    auto [euclid_error, inf_error] = get_gmgpolar_error(n_r, n_angles, non_uniformity, domain_geometry, coefficients,
                                                        boundary_conditions, source_term, solution);
    auto [euclid_error_refined, inf_error_refined] =
        get_gmgpolar_error(n_r * 2, n_angles * 2, non_uniformity, domain_geometry, coefficients, boundary_conditions,
                           source_term, solution);

    double euclid_order = log(euclid_error / euclid_error_refined) / log(2);
    double inf_order    = log(inf_error / inf_error_refined) / log(2);

    ASSERT_GT(euclid_order, 2.5);
    ASSERT_GT(inf_order, 2.5);
}
