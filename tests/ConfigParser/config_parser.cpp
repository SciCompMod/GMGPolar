#include <gtest/gtest.h>

#include "../../include/ConfigParser/config_parser.h"

// Helper to convert string to char* array for argv
std::vector<char*> make_argv(const std::vector<std::string>& args)
{
    std::vector<char*> argv;
    for (const auto& s : args) {
        argv.push_back(const_cast<char*>(s.c_str()));
    }
    return argv;
}

TEST(ConfigParserTest, ParseCommandLineArguments)
{
    // Example values
    int verbose                          = 0;
    bool paraview                        = true;
    int maxOpenMPThreads                 = 8;
    double threadReductionFactor         = 0.5;
    bool DirBC_Interior                  = true;
    int stencilDistributionMethod        = 1;
    bool cacheDensityProfileCoefficients = true;
    bool cacheDomainGeometry             = false;
    double R0                            = 1e-9;
    double Rmax                          = 2.0;
    int nr_exp                           = 3;
    int ntheta_exp                       = 5;
    int anisotropic_factor               = 2;
    int divideBy2                        = 1;
    int geometry                         = 2;
    double kappa_eps                     = 0.5;
    double delta_e                       = 1.5;
    int problem                          = 1;
    int alpha_coeff                      = 2;
    int beta_coeff                       = 1;
    double alpha_jump                    = 1.33;
    bool FMG                             = true;
    int FMG_iterations                   = 2;
    int FMG_cycle                        = 0;
    int extrapolation                    = 1;
    int maxLevels                        = 5;
    int preSmoothingSteps                = 3;
    int postSmoothingSteps               = 5;
    int multigridCycle                   = 1;
    int maxIterations                    = 400;
    int residualNormType                 = 0;
    double absoluteTolerance             = 1e-7;
    double relativeTolerance             = 1e-10;

    auto double_to_string = [](double x) {
        std::ostringstream ss;
        ss << std::setprecision(17) << x; // maximum double precision
        return ss.str();
    };

    // Construct simulated command-line arguments
    // clang-format off
    std::vector<std::string> args = {
        "program_name",
        "--verbose", std::to_string(verbose),
        "--paraview", paraview ? "1" : "0",
        "--maxOpenMPThreads", std::to_string(maxOpenMPThreads),
        "--threadReductionFactor", double_to_string(threadReductionFactor),
        "--DirBC_Interior", DirBC_Interior ? "1" : "0",
        "--stencilDistributionMethod", std::to_string(stencilDistributionMethod),
        "--cacheDensityProfileCoefficients", cacheDensityProfileCoefficients ? "1" : "0",
        "--cacheDomainGeometry", cacheDomainGeometry ? "1" : "0",
        "--R0", double_to_string(R0),
        "--Rmax", double_to_string(Rmax),
        "--nr_exp", std::to_string(nr_exp),
        "--ntheta_exp", std::to_string(ntheta_exp),
        "--anisotropic_factor", std::to_string(anisotropic_factor),
        "--divideBy2", std::to_string(divideBy2),
        "--geometry", std::to_string(geometry),
        "--kappa_eps", double_to_string(kappa_eps),
        "--delta_e", double_to_string(delta_e),
        "--problem", std::to_string(problem),
        "--alpha_coeff", std::to_string(alpha_coeff),
        "--alpha_jump", double_to_string(alpha_jump),
        "--beta_coeff", std::to_string(beta_coeff),
        "--FMG", FMG ? "1" : "0",
        "--FMG_iterations", std::to_string(FMG_iterations),
        "--FMG_cycle", std::to_string(FMG_cycle),
        "--extrapolation", std::to_string(extrapolation),
        "--maxLevels", std::to_string(maxLevels),
        "--preSmoothingSteps", std::to_string(preSmoothingSteps),
        "--postSmoothingSteps", std::to_string(postSmoothingSteps),
        "--multigridCycle", std::to_string(multigridCycle),
        "--maxIterations", std::to_string(maxIterations),
        "--residualNormType", std::to_string(residualNormType),
        "--absoluteTolerance", double_to_string(absoluteTolerance),
        "--relativeTolerance", double_to_string(relativeTolerance)
    };
    // clang-format on

    std::vector<char*> argv = make_argv(args);
    int argc                = static_cast<int>(argv.size());

    ConfigParser parser;
    ASSERT_TRUE(parser.parse(argc, argv.data()));

    // Test case objects are initialized
    EXPECT_NE(&parser.domainGeometry(), nullptr);
    EXPECT_NE(&parser.densityProfileCoefficients(), nullptr);
    EXPECT_NE(&parser.boundaryConditions(), nullptr);
    EXPECT_NE(&parser.sourceTerm(), nullptr);
    EXPECT_NE(&parser.exactSolution(), nullptr);

    // Control parameters
    EXPECT_EQ(parser.verbose(), verbose);
    EXPECT_EQ(parser.paraview(), paraview);
    EXPECT_EQ(parser.maxOpenMPThreads(), maxOpenMPThreads);
    EXPECT_DOUBLE_EQ(parser.threadReductionFactor(), threadReductionFactor);
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
}