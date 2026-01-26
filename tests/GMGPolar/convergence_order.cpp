#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <vector>
#include <random>

#include "../../include/GMGPolar/gmgpolar.h"

template <class T>
class GMGPolarPaperTestCase;

// clang-format off
template <
    class DensityProfileCoefficientsType,
    class BoundaryConditionsType,
    class SourceTermType,
    class ExactSolutionType,
    ExtrapolationType extrapolation_,
    double expected_l2_order_,
    double expected_inf_order_
>
class GMGPolarPaperTestCase<
    std::tuple<
        DensityProfileCoefficientsType,
        BoundaryConditionsType,
        SourceTermType,
        ExactSolutionType,
        std::integral_constant<ExtrapolationType, extrapolation_>,
        std::integral_constant<double, expected_l2_order_>,
        std::integral_constant<double, expected_inf_order_>
    >
> : public testing::Test
// clang-format on

{
public:
    using DensityProfileCoefficients = DensityProfileCoefficientsType;
    using BoundaryConditions         = BoundaryConditionsType;
    using SourceTerm                 = SourceTermType;
    using ExactSolution              = ExactSolutionType;

    static constexpr ExtrapolationType extrapolation = extrapolation_;
    static constexpr double expected_l2_order        = expected_l2_order_;
    static constexpr double expected_inf_order       = expected_inf_order_;
};

// clang-format off
using gmg_paper_types = testing::Types<
    /* ----------------- */
    /* Solution: PolarR6 */
    /* ----------------- */
    /* No Extrapolation */
    /* Expected order: L2 = 2.0, Infinity = 2.0 */
    std::tuple<
        ZoniShiftedGyroCoefficients, PolarR6_Boundary_CzarnyGeometry,
        PolarR6_ZoniShiftedGyro_CzarnyGeometry, PolarR6_CzarnyGeometry,
        std::integral_constant<ExtrapolationType, ExtrapolationType::NONE>,
        std::integral_constant<double, 1.9>, std::integral_constant<double, 1.9>
    >,
    /* Implicit Extrapolation. */
    /* Expected rder: L2 = 4.0, Infinity = 3.0 */
    std::tuple<
        ZoniShiftedGyroCoefficients, PolarR6_Boundary_CzarnyGeometry,
        PolarR6_ZoniShiftedGyro_CzarnyGeometry, PolarR6_CzarnyGeometry,
        std::integral_constant<ExtrapolationType, ExtrapolationType::IMPLICIT_EXTRAPOLATION>,
        std::integral_constant<double, 3.9>, std::integral_constant<double, 2.9>
    >,
    /* Combined Extrapolation */
    /* Expected order: L2 = 4.0, Infinity = 3.0 */
    std::tuple<ZoniShiftedGyroCoefficients, PolarR6_Boundary_CzarnyGeometry,
        PolarR6_ZoniShiftedGyro_CzarnyGeometry, PolarR6_CzarnyGeometry,
        std::integral_constant<ExtrapolationType, ExtrapolationType::COMBINED>,
        std::integral_constant<double, 3.9>, std::integral_constant<double, 2.9>
    >,
    /* Implicit Extrapolation with Fullgrid-Smoothing */
    /* Expected order: L2 = 4.0, Infinity = 3.0 */
    std::tuple<
        ZoniShiftedGyroCoefficients, PolarR6_Boundary_CzarnyGeometry,
        PolarR6_ZoniShiftedGyro_CzarnyGeometry, PolarR6_CzarnyGeometry,
        std::integral_constant<ExtrapolationType, ExtrapolationType::IMPLICIT_FULL_GRID_SMOOTHING>,
        std::integral_constant<double, 3.0>, std::integral_constant<double, 2.5>
    >,

    /* --------------------- */
    /* Solution: CartesianR6 */
    /* --------------------- */
    /* No Extrapolation */
    /* Expected order: L2 = 2.0, Infinity = 2.0 */
    std::tuple<
        ZoniShiftedGyroCoefficients, CartesianR6_Boundary_CzarnyGeometry,
        CartesianR6_ZoniShiftedGyro_CzarnyGeometry, CartesianR6_CzarnyGeometry,
        std::integral_constant<ExtrapolationType, ExtrapolationType::NONE>,
        std::integral_constant<double, 1.9>, std::integral_constant<double, 1.9>
    >,
    /* Implicit Extrapolation. */
    /* Expected rder: L2 = 4.0, Infinity = 3.0 */
    std::tuple<
        ZoniShiftedGyroCoefficients, CartesianR6_Boundary_CzarnyGeometry,
        CartesianR6_ZoniShiftedGyro_CzarnyGeometry, CartesianR6_CzarnyGeometry,
        std::integral_constant<ExtrapolationType, ExtrapolationType::IMPLICIT_EXTRAPOLATION>,
        std::integral_constant<double, 3.4>, std::integral_constant<double, 2.9>
    >,
    /* Combined Extrapolation */
    /* Expected order: L2 = 4.0, Infinity = 3.0 */
    std::tuple<
        ZoniShiftedGyroCoefficients, CartesianR6_Boundary_CzarnyGeometry,
        CartesianR6_ZoniShiftedGyro_CzarnyGeometry, CartesianR6_CzarnyGeometry,
        std::integral_constant<ExtrapolationType, ExtrapolationType::COMBINED>,
        std::integral_constant<double, 3.4>, std::integral_constant<double, 2.9>
    >,
    /* Implicit Extrapolation with Fullgrid-Smoothing */
    /* Expected order: L2 = 4.0, Infinity = 3.0 */
    std::tuple<
        ZoniShiftedGyroCoefficients, CartesianR6_Boundary_CzarnyGeometry,
        CartesianR6_ZoniShiftedGyro_CzarnyGeometry, CartesianR6_CzarnyGeometry,
        std::integral_constant<ExtrapolationType, ExtrapolationType::IMPLICIT_FULL_GRID_SMOOTHING>,
        std::integral_constant<double, 3.0>, std::integral_constant<double, 2.5>
    >
>;
// clang-format on

TYPED_TEST_SUITE(GMGPolarPaperTestCase, gmg_paper_types);

std::vector<double> get_non_uniform_points(double min, double max, int n_pts, double non_uniformity = 0.1)
{
    // Fixed seed for reproducibility in tests
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(-0.5, 0.5);

    std::vector<double> points(n_pts);
    int n_cells        = n_pts - 1;
    double const delta = (max - min) / double(n_cells);

    points[0] = min;
    for (int i = 1; i < n_cells; ++i) {
        points[i] = min + (i + dist(gen) * non_uniformity) * delta;
    }
    points[n_cells] = max;

    return points;
}

std::vector<double> refine(std::vector<double> const& original_points)
{
    std::vector<double> refined(2 * original_points.size() - 1);
    refined[0] = original_points[0];
    for (std::size_t i(1); i < original_points.size(); ++i) {
        refined[2 * i - 1] = 0.5 * (original_points[i] + original_points[i - 1]);
        refined[2 * i]     = original_points[i];
    }
    return refined;
}

template <class DensityProfileCoefficients>
std::tuple<double, double>
get_gmgpolar_error(PolarGrid const& grid, CzarnyGeometry const& domain_geometry,
                   DensityProfileCoefficients const& coefficients, BoundaryConditions const& boundary_conditions,
                   SourceTerm const& source_term, ExactSolution const& solution, ExtrapolationType extrapolation)
{
    GMGPolar gmgpolar(grid, domain_geometry, coefficients);
    gmgpolar.setSolution(&solution);

    // --- General solver output and visualization settings --- //
    gmgpolar.verbose(0);
    gmgpolar.paraview(false);

    // --- Parallelization settings --- //
    gmgpolar.maxOpenMPThreads(1);

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

    double Rmax      = 1.3;
    double kappa_eps = 0.3;
    double delta_e   = 1.4;
    CzarnyGeometry domain_geometry(Rmax, kappa_eps, delta_e);

    const double alpha_jump = 0.0; // Unused value
    typename TestFixture::DensityProfileCoefficients coefficients(Rmax, alpha_jump);
    typename TestFixture::BoundaryConditions boundary_conditions(Rmax, kappa_eps, delta_e);
    typename TestFixture::SourceTerm source_term(Rmax, kappa_eps, delta_e);
    typename TestFixture::ExactSolution solution(Rmax, kappa_eps, delta_e);

    // For extrapolation, we need a uniform refinement
    std::vector<double> non_uniform_radii  = get_non_uniform_points(1e-8, Rmax, n_r / 2 + 1, non_uniformity);
    std::vector<double> non_uniform_angles = get_non_uniform_points(0.0, 2 * M_PI, n_angles / 2 + 1, non_uniformity);
    // Extrapolation requires at least one uniform refinement
    std::vector<double> radii  = refine(non_uniform_radii);
    std::vector<double> angles = refine(non_uniform_angles);
    // Get refined points
    std::vector<double> radii_refined  = refine(radii);
    std::vector<double> angles_refined = refine(angles);

    auto [euclid_error, inf_error] =
        get_gmgpolar_error(PolarGrid(radii, angles), domain_geometry, coefficients, boundary_conditions, source_term,
                           solution, TestFixture::extrapolation);
    auto [euclid_error_refined, inf_error_refined] =
        get_gmgpolar_error(PolarGrid(radii_refined, angles_refined), domain_geometry, coefficients, boundary_conditions,
                           source_term, solution, TestFixture::extrapolation);

    double euclid_order = log(euclid_error / euclid_error_refined) / log(2);
    double inf_order    = log(inf_error / inf_error_refined) / log(2);

    ASSERT_GT(euclid_order, TestFixture::expected_l2_order);
    ASSERT_GT(inf_order, TestFixture::expected_inf_order);
}

TYPED_TEST(GMGPolarPaperTestCase, TestUniform)
{
    test_convergence<TestFixture>(0.0);
}

TYPED_TEST(GMGPolarPaperTestCase, TestNonUniform)
{
    test_convergence<TestFixture>(0.5);
}
