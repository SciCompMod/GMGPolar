#include <gtest/gtest.h>
#include "../../include/ConfigParser/config_parser.h"

struct TestParams {
    int geometry;
    int problem;
    int alpha_coeff;
    int beta_coeff;
    double kappa_eps;
    double delta_e;
    int case_id;
};

class ConfigParserTest : public testing::TestWithParam<TestParams>
{
protected:
    void SetUp() override
    {
        params = GetParam();
    }

    TestParams params;
};

std::vector<char*> make_argv(const std::vector<std::string>& args)
{
    std::vector<char*> argv;
    for (const auto& s : args) {
        argv.push_back(const_cast<char*>(s.c_str()));
    }
    return argv;
}

TEST_P(ConfigParserTest, ParseAllGeometryAndProblemCombinations)
{
    const int verbose                          = (params.case_id > 0 ? 0 : 1);
    const bool paraview                        = false;
    const int maxOpenMPThreads                 = 2;
    const bool DirBC_Interior                  = false;
    const int stencilDistributionMethod        = params.case_id % 2;
    const bool cacheDensityProfileCoefficients = true;
    const bool cacheDomainGeometry             = true;
    const double R0                            = 1e-8;
    const double Rmax                          = 1.3;
    const int nr_exp                           = 3;
    const int ntheta_exp                       = 4;
    const int anisotropic_factor               = 1;
    const int divideBy2                        = 1;
    const bool FMG                             = true;
    const int FMG_iterations                   = 1;
    const int FMG_cycle                        = params.case_id % 3;
    const int extrapolation                    = params.case_id % 4;
    const int maxLevels                        = 7;
    const int preSmoothingSteps                = 1;
    const int postSmoothingSteps               = 1;
    const int multigridCycle                   = params.case_id % 3;
    const int maxIterations                    = 150;
    const int residualNormType                 = params.case_id % 3;
    const double absoluteTolerance             = 1e-8;
    const double relativeTolerance             = 1e-5;

    // Calculate alpha_jump based on alpha_coeff
    double alpha_jump;
    switch (params.alpha_coeff) {
    case 0:
        alpha_jump = 0.5 * Rmax;
        break;
    case 1:
        alpha_jump = 0.66 * Rmax;
        break;
    case 2:
        alpha_jump = 0.4837 * Rmax;
        break;
    case 3:
        alpha_jump = 0.678 * Rmax;
        break;
    default:
        FAIL() << "Invalid alpha_coeff value";
    }

    auto double_to_string = [](double x) {
        std::ostringstream ss;
        ss << std::setprecision(17) << x;
        return ss.str();
    };

    // Construct command-line arguments
    std::vector<std::string> args = {"program_name",
                                     "--verbose",
                                     std::to_string(verbose),
                                     "--paraview",
                                     paraview ? "1" : "0",
                                     "--maxOpenMPThreads",
                                     std::to_string(maxOpenMPThreads),
                                     "--DirBC_Interior",
                                     DirBC_Interior ? "1" : "0",
                                     "--stencilDistributionMethod",
                                     std::to_string(stencilDistributionMethod),
                                     "--cacheDensityProfileCoefficients",
                                     cacheDensityProfileCoefficients ? "1" : "0",
                                     "--cacheDomainGeometry",
                                     cacheDomainGeometry ? "1" : "0",
                                     "--R0",
                                     double_to_string(R0),
                                     "--Rmax",
                                     double_to_string(Rmax),
                                     "--nr_exp",
                                     std::to_string(nr_exp),
                                     "--ntheta_exp",
                                     std::to_string(ntheta_exp),
                                     "--anisotropic_factor",
                                     std::to_string(anisotropic_factor),
                                     "--divideBy2",
                                     std::to_string(divideBy2),
                                     "--geometry",
                                     std::to_string(params.geometry),
                                     "--kappa_eps",
                                     double_to_string(params.kappa_eps),
                                     "--delta_e",
                                     double_to_string(params.delta_e),
                                     "--problem",
                                     std::to_string(params.problem),
                                     "--alpha_coeff",
                                     std::to_string(params.alpha_coeff),
                                     "--alpha_jump",
                                     double_to_string(alpha_jump),
                                     "--beta_coeff",
                                     std::to_string(params.beta_coeff),
                                     "--FMG",
                                     FMG ? "1" : "0",
                                     "--FMG_iterations",
                                     std::to_string(FMG_iterations),
                                     "--FMG_cycle",
                                     std::to_string(FMG_cycle),
                                     "--extrapolation",
                                     std::to_string(extrapolation),
                                     "--maxLevels",
                                     std::to_string(maxLevels),
                                     "--preSmoothingSteps",
                                     std::to_string(preSmoothingSteps),
                                     "--postSmoothingSteps",
                                     std::to_string(postSmoothingSteps),
                                     "--multigridCycle",
                                     std::to_string(multigridCycle),
                                     "--maxIterations",
                                     std::to_string(maxIterations),
                                     "--residualNormType",
                                     std::to_string(residualNormType),
                                     "--absoluteTolerance",
                                     double_to_string(absoluteTolerance),
                                     "--relativeTolerance",
                                     double_to_string(relativeTolerance)};

    std::vector<char*> argv = make_argv(args);
    int argc                = argv.size();

    ConfigParser parser;
    ASSERT_TRUE(parser.parse(argc, argv.data()));

    // Test case objects are initialized
    EXPECT_NE(&parser.domainGeometry(), nullptr);
    EXPECT_NE(&parser.densityProfileCoefficients(), nullptr);
    EXPECT_NE(&parser.exactSolution(), nullptr);

    // Control parameters
    EXPECT_EQ(parser.verbose(), verbose);
    EXPECT_EQ(parser.paraview(), paraview);
    EXPECT_EQ(parser.maxOpenMPThreads(), maxOpenMPThreads);
    EXPECT_EQ(parser.DirBC_Interior(), DirBC_Interior);
    EXPECT_EQ(parser.stencilDistributionMethod(), static_cast<StencilDistributionMethod>(stencilDistributionMethod));
    EXPECT_EQ(parser.cacheDensityProfileCoefficients(), cacheDensityProfileCoefficients);
    EXPECT_EQ(parser.cacheDomainGeometry(), cacheDomainGeometry);

    // Grid
    const PolarGrid& grid = parser.grid();
    EXPECT_NE(&grid, nullptr);
    EXPECT_DOUBLE_EQ(grid.radius(0), R0);
    EXPECT_DOUBLE_EQ(grid.radius(grid.nr() - 1), Rmax);

    // Multigrid parameters
    EXPECT_EQ(parser.FMG(), FMG);
    EXPECT_EQ(parser.FMG_iterations(), FMG_iterations);
    EXPECT_EQ(parser.FMG_cycle(), static_cast<MultigridCycleType>(FMG_cycle));
    EXPECT_EQ(parser.extrapolation(), static_cast<ExtrapolationType>(extrapolation));
    EXPECT_EQ(parser.maxLevels(), maxLevels);
    EXPECT_EQ(parser.preSmoothingSteps(), preSmoothingSteps);
    EXPECT_EQ(parser.postSmoothingSteps(), postSmoothingSteps);
    EXPECT_EQ(parser.multigridCycle(), static_cast<MultigridCycleType>(multigridCycle));
    EXPECT_EQ(parser.maxIterations(), maxIterations);
    EXPECT_EQ(parser.residualNormType(), static_cast<ResidualNormType>(residualNormType));
    ASSERT_TRUE(parser.absoluteTolerance().has_value());
    EXPECT_DOUBLE_EQ(parser.absoluteTolerance().value(), absoluteTolerance);
    ASSERT_TRUE(parser.relativeTolerance().has_value());
    EXPECT_DOUBLE_EQ(parser.relativeTolerance().value(), relativeTolerance);

    // Solver
    std::visit(
        [&parser](auto const& domain_geometry, auto const& density_profile_coeffs) {
            // Get the types of the domain geometry and the density profile coefficients
            using DG = std::decay_t<decltype(domain_geometry)>;
            using DC = std::decay_t<decltype(density_profile_coeffs)>;
            // Create GMGPolar solver for the selected geometry and coefficient types
            GMGPolar<DG, DC> solver(parser.grid(), domain_geometry, density_profile_coeffs);

            // --- General solver output and visualization settings --- //
            solver.verbose(parser.verbose()); // Enable/disable verbose output
            solver.paraview(parser.paraview()); // Enable/disable ParaView output

            // --- Parallelization and threading settings --- //
            solver.maxOpenMPThreads(parser.maxOpenMPThreads()); // Maximum OpenMP threads to use
            omp_set_num_threads(parser.maxOpenMPThreads()); // Global OpenMP thread limit

            // --- Numerical method setup --- //
            solver.DirBC_Interior(
                parser.DirBC_Interior()); // Interior boundary conditions: Dirichlet, Across-the-origin
            solver.stencilDistributionMethod(
                parser.stencilDistributionMethod()); // Stencil distribution strategy: Take, Give
            solver.cacheDensityProfileCoefficients(
                parser.cacheDensityProfileCoefficients()); // Cache density profile coefficients: alpha, beta
            solver.cacheDomainGeometry(
                parser.cacheDomainGeometry()); // Cache domain geometry data: arr, att, art, detDF

            // --- Multigrid settings --- //
            solver.extrapolation(parser.extrapolation()); // Enable/disable extrapolation
            solver.maxLevels(parser.maxLevels()); // Max multigrid levels (-1 = use deepest possible)
            solver.preSmoothingSteps(parser.preSmoothingSteps()); // Smoothing before coarse-grid correction
            solver.postSmoothingSteps(parser.postSmoothingSteps()); // Smoothing after coarse-grid correction
            solver.multigridCycle(parser.multigridCycle()); // Multigrid cycle type
            solver.FMG(parser.FMG()); // Full Multigrid mode on/off
            solver.FMG_iterations(parser.FMG_iterations()); // FMG iteration count
            solver.FMG_cycle(parser.FMG_cycle()); // FMG cycle type

            // --- Preconditioned Conjugate Gradient settings --- //
            solver.PCG(parser.PCG()); // Preconditioned Conjugate Gradient mode on/off
            solver.PCG_FMG(parser.PCG_FMG()); // Use FMG as preconditioner for PCG
            solver.PCG_FMG_iterations(parser.PCG_FMG_iterations()); // FMG iterations for PCG preconditioner
            solver.PCG_FMG_cycle(parser.PCG_FMG_cycle()); // FMG cycle type for PCG preconditioner
            solver.PCG_MG_iterations(parser.PCG_MG_iterations()); // Multigrid iterations for PCG preconditioner
            solver.PCG_MG_cycle(parser.PCG_MG_cycle()); // Multigrid cycle type for PCG preconditioner

            // --- Iterative solver controls --- //
            solver.maxIterations(parser.maxIterations()); // Max number of iterations
            solver.residualNormType(parser.residualNormType()); // Residual norm type (L2, weighted-L2, L∞)
            solver.absoluteTolerance(parser.absoluteTolerance()); // Absolute residual tolerance
            solver.relativeTolerance(parser.relativeTolerance()); // Relative residual tolerance

            // --- Finalize solver setup --- //
            solver.setup(); // (allocates internal data, prepares operators, etc.)

            // --- Provide optional exact solution --- //
            solver.setSolution(&parser.exactSolution());

            // --- Solve Phase --- //
            parser.solve(solver);

            // --- Retrieve solution and associated grid --- //
            Vector<double> solution = solver.solution();
            const PolarGrid& grid   = solver.grid();

            // --- Verify results ---
            int number_of_iterations                             = solver.numberOfIterations();
            std::optional<double> exact_error_weighted_euclidean = solver.exactErrorWeightedEuclidean();
            std::optional<double> exact_infinity_error           = solver.exactErrorInfinity();
            double reductionFactor                               = solver.meanResidualReductionFactor();

            ASSERT_TRUE(exact_error_weighted_euclidean.has_value());
            ASSERT_TRUE(exact_infinity_error.has_value());

            EXPECT_LE(exact_error_weighted_euclidean.value(), 1e-3);
            EXPECT_LE(exact_infinity_error.value(), 1e-3);
            EXPECT_GE(number_of_iterations, 1);
        },
        parser.domainGeometry(), parser.densityProfileCoefficients());
}

// Define test cases covering all combinations
std::vector<TestParams> generate_all_combinations()
{
    std::vector<TestParams> cases;
    int case_id = 0;

    // Non-Culham geometries
    for (int geometry = 0; geometry <= 2; ++geometry) {
        for (int problem = 0; problem <= 3; ++problem) {
            if (problem == 3) {
                double kappa_eps = (geometry == 1 || geometry == 2) ? 0.3 : 0.0;
                double delta_e   = (geometry == 1) ? 0.2 : (geometry == 2) ? 1.4 : 0.0;
                cases.push_back(TestParams{geometry, 3, 3, 1, kappa_eps, delta_e, case_id++});
            }
            else {
                for (int alpha_coeff = 0; alpha_coeff <= 3; ++alpha_coeff) {
                    for (int beta_coeff = 0; beta_coeff <= 1; ++beta_coeff) {
                        double kappa_eps = (geometry == 1 || geometry == 2) ? 0.3 : 0.0;
                        double delta_e   = (geometry == 1) ? 0.2 : (geometry == 2) ? 1.4 : 0.0;
                        cases.push_back(
                            TestParams{geometry, problem, alpha_coeff, beta_coeff, kappa_eps, delta_e, case_id++});
                    }
                }
            }
        }
    }

    // Culham geometry
    for (int problem : {2, 3}) {
        cases.push_back(TestParams{3, problem, 3, 1, 0.0, 0.0, case_id++});
    }

    return cases;
}
INSTANTIATE_TEST_SUITE_P(AllCombinations, ConfigParserTest, ::testing::ValuesIn(generate_all_combinations()));
