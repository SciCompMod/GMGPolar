#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>

#include "../../include/GMGPolar/gmgpolar.h"

double constexpr Rmax = 1.3;

template <class T>
class GMGPolarPaperTestCase;

template <ExtrapolationType extrapolation_, class DensityProfileCoefficientsType, class BoundaryConditionsType,
          class SourceTermType, class ExactSolutionType, double expected_order_>
class GMGPolarPaperTestCase<std::tuple<std::integral_constant<ExtrapolationType, extrapolation_>,
                                       DensityProfileCoefficientsType, BoundaryConditionsType, SourceTermType,
                                       ExactSolutionType, std::integral_constant<double, expected_order_>>>
    : public testing::Test
{
public:
    using DensityProfileCoefficients = DensityProfileCoefficientsType;
    using BoundaryConditions         = BoundaryConditionsType;
    using SourceTerm                 = SourceTermType;
    using ExactSolution              = ExactSolutionType;

    static constexpr ExtrapolationType extrapolation = extrapolation_;
    static constexpr double expected_order           = expected_order_;
};

using gmg_paper_types = testing::Types<
    /* PolarR6 - NONE_EXTRAPOLATION: Order >= 1.9 */
    std::tuple<std::integral_constant<ExtrapolationType, ExtrapolationType::NONE>, ZoniShiftedGyroCoefficients,
               PolarR6_Boundary_CzarnyGeometry, PolarR6_ZoniShiftedGyro_CzarnyGeometry, PolarR6_CzarnyGeometry,
               std::integral_constant<double, 1.9>>,
    /* CartesianR6 - NONE_EXTRAPOLATION: order >= 1.9 */
    std::tuple<std::integral_constant<ExtrapolationType, ExtrapolationType::NONE>, ZoniShiftedGyroCoefficients,
               CartesianR6_Boundary_CzarnyGeometry, CartesianR6_ZoniShiftedGyro_CzarnyGeometry,
               CartesianR6_CzarnyGeometry, std::integral_constant<double, 1.9>>,
    /* PolarR6 - IMPLICIT_EXTRAPOLATION: order >= 3.9 */
    std::tuple<std::integral_constant<ExtrapolationType, ExtrapolationType::IMPLICIT_EXTRAPOLATION>,
               ZoniShiftedGyroCoefficients, PolarR6_Boundary_CzarnyGeometry, PolarR6_ZoniShiftedGyro_CzarnyGeometry,
               PolarR6_CzarnyGeometry, std::integral_constant<double, 2.9>>,
    /* CartesianR6 - IMPLICIT_EXTRAPOLATION: order >= 3.4 */
    std::tuple<std::integral_constant<ExtrapolationType, ExtrapolationType::IMPLICIT_EXTRAPOLATION>,
               ZoniShiftedGyroCoefficients, CartesianR6_Boundary_CzarnyGeometry,
               CartesianR6_ZoniShiftedGyro_CzarnyGeometry, CartesianR6_CzarnyGeometry,
               std::integral_constant<double, 2.4>>,
    /* PolarR6 - COMBINED_EXTRAPOLATION: order >= 3.9 */
    std::tuple<std::integral_constant<ExtrapolationType, ExtrapolationType::COMBINED>, ZoniShiftedGyroCoefficients,
               PolarR6_Boundary_CzarnyGeometry, PolarR6_ZoniShiftedGyro_CzarnyGeometry, PolarR6_CzarnyGeometry,
               std::integral_constant<double, 2.9>>,
    /* CartesianR6 - COMBINED_EXTRAPOLATION: order >= 3.5 */
    std::tuple<std::integral_constant<ExtrapolationType, ExtrapolationType::COMBINED>, ZoniShiftedGyroCoefficients,
               CartesianR6_Boundary_CzarnyGeometry, CartesianR6_ZoniShiftedGyro_CzarnyGeometry,
               CartesianR6_CzarnyGeometry, std::integral_constant<double, 2.4>>,
    /* PolarR6 - IMPLICIT_EXTRAPOLATION_FULL_GRID_SMOOTHING: order >= 3.4 */
    std::tuple<std::integral_constant<ExtrapolationType, ExtrapolationType::IMPLICIT_FULL_GRID_SMOOTHING>,
               ZoniShiftedGyroCoefficients, PolarR6_Boundary_CzarnyGeometry, PolarR6_ZoniShiftedGyro_CzarnyGeometry,
               PolarR6_CzarnyGeometry, std::integral_constant<double, 2.4>>,
    /* CartesianR6 - IMPLICIT_EXTRAPOLATION_FULL_GRID_SMOOTHING: order >= 2.9 */
    std::tuple<std::integral_constant<ExtrapolationType, ExtrapolationType::IMPLICIT_FULL_GRID_SMOOTHING>,
               ZoniShiftedGyroCoefficients, CartesianR6_Boundary_CzarnyGeometry,
               CartesianR6_ZoniShiftedGyro_CzarnyGeometry, CartesianR6_CzarnyGeometry,
               std::integral_constant<double, 2.4>>>;

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

std::vector<double> refine(std::vector<double> const& original_points, double non_uniformity)
{
    std::vector<double> refined(2 * original_points.size() - 1);
    refined[0] = original_points[0];
    for (std::size_t i(1); i < original_points.size(); ++i) {
        double const random_perturbation = double(rand()) / RAND_MAX - 0.5;

        refined[2 * i - 1] = original_points[i - 1] + 0.5 * (1 + random_perturbation * non_uniformity) *
                                                          (original_points[i] - original_points[i - 1]);
        refined[2 * i] = original_points[i];
    }
    return refined;
}

std::tuple<double, double> get_gmgpolar_error(PolarGrid const& grid, CzarnyGeometry const& domain_geometry,
                                              DensityProfileCoefficients const& coefficients,
                                              BoundaryConditions const& boundary_conditions,
                                              SourceTerm const& source_term, ExactSolution const& solution,
                                              ExtrapolationType extrapolation)
{
    GMGPolar gmgpolar(grid, domain_geometry, coefficients);
    gmgpolar.setSolution(&solution);

    // --- General solver output and visualization settings --- //
    gmgpolar.verbose(0);
    gmgpolar.paraview(false);

    // --- Parallelization and threading settings --- //
    gmgpolar.maxOpenMPThreads(1);
    gmgpolar.threadReductionFactor(1.0);

    // --- Discretization and method settings --- //
    gmgpolar.DirBC_Interior(false); // Use across-origin calculation
    gmgpolar.stencilDistributionMethod(StencilDistributionMethod::CPU_TAKE);
    gmgpolar.cacheDensityProfileCoefficients(true);
    gmgpolar.cacheDomainGeometry(true);

    // --- Multigrid settings --- //
    gmgpolar.FMG(true);
    gmgpolar.FMG_iterations(3);
    gmgpolar.FMG_cycle(MultigridCycleType::F_CYCLE);
    gmgpolar.extrapolation(extrapolation);
    gmgpolar.maxLevels(-1);
    gmgpolar.preSmoothingSteps(1);
    gmgpolar.postSmoothingSteps(1);
    gmgpolar.multigridCycle(MultigridCycleType::V_CYCLE);

    // --- Iterative gmgpolar controls --- //
    gmgpolar.maxIterations(25);
    gmgpolar.residualNormType(ResidualNormType::EUCLIDEAN);
    gmgpolar.absoluteTolerance(1e-7);
    gmgpolar.relativeTolerance(1e-6);
    // ----------------------------------------------------------------

    // Perform setup and solve
    gmgpolar.setup();
    gmgpolar.solve(boundary_conditions, source_term);

    double euclid_error = *gmgpolar.exactErrorWeightedEuclidean();
    double inf_error    = *gmgpolar.exactErrorInfinity();
    return std::make_tuple(euclid_error, inf_error);
}

template <class TestFixture>
void test_convergence(double non_uniformity)
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

    std::vector<double> radii  = get_non_uniform_points(1e-8, Rmax, n_r + 1, non_uniformity); // remove central point
    std::vector<double> angles = get_non_uniform_points(0.0, M_PI, n_angles / 2 + 1, non_uniformity);
    // Every node in the interior ring needs to have an opposite neighboring node
    for (int i(1); i < n_angles / 2 + 1; ++i) {
        angles.push_back(angles[i] + M_PI);
    }
    std::vector<double> radii_refined  = refine(radii, non_uniformity);
    std::vector<double> angles_refined = refine(angles, non_uniformity);

    auto [euclid_error, inf_error] =
        get_gmgpolar_error(PolarGrid(radii, angles), domain_geometry, coefficients, boundary_conditions, source_term,
                           solution, TestFixture::extrapolation);
    auto [euclid_error_refined, inf_error_refined] =
        get_gmgpolar_error(PolarGrid(radii_refined, angles_refined), domain_geometry, coefficients, boundary_conditions,
                           source_term, solution, TestFixture::extrapolation);

    double euclid_order = log(euclid_error / euclid_error_refined) / log(2);
    double inf_order    = log(inf_error / inf_error_refined) / log(2);

    double expected_order = TestFixture::expected_order;

    ASSERT_GT(euclid_order, expected_order);
    ASSERT_GT(inf_order, expected_order);
}

TYPED_TEST(GMGPolarPaperTestCase, TestUniform)
{
    test_convergence<TestFixture>(0.0);
}

TYPED_TEST(GMGPolarPaperTestCase, TestNonUniform)
{
    test_convergence<TestFixture>(0.1);
}
