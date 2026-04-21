#pragma once

#include <optional>
#include <string>
#include <stdexcept>
#include <memory>

#include "../../include/ConfigParser/cmdline.h"
#include "../../include/Definitions/global_definitions.h"
#include "../../include/PolarGrid/polargrid.h"
#include "../../include/GMGPolar/test_cases.h"
#include "../../include/GMGPolar/igmgpolar.h"
#include "../../include/GMGPolar/gmgpolar.h"
#include "test_selection.h"

namespace gmgpolar
{

class ConfigParser
{
public:
    /* Specifies whether user input is required */
    static constexpr bool OPTIONAL = false;
    static constexpr bool REQUIRED = true;

    ConfigParser();

    bool parse(int argc, char* argv[]);

    // Test Case
    const DomainGeometryVariant& domainGeometry() const;
    const DensityProfileCoefficientsVariant& densityProfileCoefficients() const;
    const ExactSolution& exactSolution() const;
    std::unique_ptr<IGMGPolar> solver() const;

    template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
    void solve(GMGPolar<DomainGeometry, DensityProfileCoefficients>& solver) const;

    // Control Parameters
    int verbose() const;
    bool paraview() const;
    int maxOpenMPThreads() const;
    bool DirBC_Interior() const;
    StencilDistributionMethod stencilDistributionMethod() const;
    bool cacheDensityProfileCoefficients() const;
    bool cacheDomainGeometry() const;

    const PolarGrid& grid() const;

    // Full Multigrid Method
    bool FMG() const;
    int FMG_iterations() const;
    MultigridCycleType FMG_cycle() const;

    // Preconditioned Conjugate Gradient Method
    bool PCG() const;
    bool PCG_FMG() const;
    int PCG_FMG_iterations() const;
    MultigridCycleType PCG_FMG_cycle() const;
    int PCG_MG_iterations() const;
    MultigridCycleType PCG_MG_cycle() const;

    // Multigrid Parameters
    ExtrapolationType extrapolation() const;
    int maxLevels() const;
    int preSmoothingSteps() const;
    int postSmoothingSteps() const;
    MultigridCycleType multigridCycle() const;
    int maxIterations() const;
    ResidualNormType residualNormType() const;
    std::optional<double> absoluteTolerance() const;
    std::optional<double> relativeTolerance() const;

private:
    // Parse command-line arguments to extract problem configuration
    cmdline::parser parser_;
    // Input Functions
    std::unique_ptr<const DomainGeometryVariant> domain_geometry_;
    std::unique_ptr<const DensityProfileCoefficientsVariant> density_profile_coefficients_;
    std::unique_ptr<const ExactSolution> exact_solution_;
    GeometryType geometry_type_;
    ProblemType problem_type_;
    AlphaCoeff alpha_type_;
    BetaCoeff beta_type_;
    double Rmax_;
    double kappa_eps_;
    double delta_e_;
    double alpha_jump_;
    // General solver output and visualization settings
    int verbose_;
    bool paraview_;
    // Parallelization and threading settings
    int max_omp_threads_;
    // Numerical method setup
    bool DirBC_Interior_;
    StencilDistributionMethod stencil_distribution_method_;
    bool cache_density_profile_coefficients_;
    bool cache_domain_geometry_;
    // Grid configuration
    PolarGrid grid_;
    // Multigrid settings
    ExtrapolationType extrapolation_;
    int max_levels_;
    int pre_smoothing_steps_;
    int post_smoothing_steps_;
    MultigridCycleType multigrid_cycle_;
    bool FMG_;
    int FMG_iterations_;
    MultigridCycleType FMG_cycle_;
    bool PCG_;
    bool PCG_FMG_;
    int PCG_FMG_iterations_;
    MultigridCycleType PCG_FMG_cycle_;
    int PCG_MG_iterations_;
    MultigridCycleType PCG_MG_cycle_;
    // Iterative solver controls
    int max_iterations_;
    ResidualNormType residual_norm_type_;
    std::optional<double> absolute_tolerance_;
    std::optional<double> relative_tolerance_;

    void selectTestCase(GeometryType geometry_type, ProblemType problem_type, AlphaCoeff alpha_type,
                        BetaCoeff beta_type, double Rmax, double kappa_eps, double delta_e, double alpha_jump);
};
} // namespace gmgpolar
