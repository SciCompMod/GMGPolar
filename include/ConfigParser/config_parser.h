#pragma once

#include <optional>
#include <string>
#include <stdexcept>
#include <memory>

#include "../../include/ConfigParser/cmdline.h"
#include "../../include/common/global_definitions.h"
#include "../../include/PolarGrid/polargrid.h"
#include "../../include/GMGPolar/test_cases.h"
#include "../../include/GMGPolar/igmgpolar.h"
#include "test_selection.h"

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
    const BoundaryConditions& boundaryConditions() const;
    const SourceTerm& sourceTerm() const;
    const ExactSolution& exactSolution() const;
    std::unique_ptr<IGMGPolar> solver() const;

    // Control Parameters
    int verbose() const;
    bool paraview() const;
    int maxOpenMPThreads() const;
    bool DirBC_Interior() const;
    StencilDistributionMethod stencilDistributionMethod() const;
    bool cacheDensityProfileCoefficients() const;
    bool cacheDomainGeometry() const;

    const PolarGrid& grid() const;

    // Multigrid Parameters
    bool FMG() const;
    int FMG_iterations() const;
    MultigridCycleType FMG_cycle() const;
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
    std::unique_ptr<const BoundaryConditions> boundary_conditions_;
    std::unique_ptr<const SourceTerm> source_term_;
    std::unique_ptr<const ExactSolution> exact_solution_;
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
    // Iterative solver controls
    int max_iterations_;
    ResidualNormType residual_norm_type_;
    std::optional<double> absolute_tolerance_;
    std::optional<double> relative_tolerance_;

    void selectTestCase(GeometryType geometry_type, ProblemType problem_type, AlphaCoeff alpha_type,
                        BetaCoeff beta_type, double Rmax, double kappa_eps, double delta_e, double alpha_jump);
};
