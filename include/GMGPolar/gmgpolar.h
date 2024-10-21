#pragma once

#include <iostream>
#include <memory>
#include <filesystem>
#include <optional>
#include <utility>
#include <chrono>
#include <omp.h>

class Level;
class LevelCache;

#include "../Level/level.h"
#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/densityProfileCoefficients.h"
#include "../InputFunctions/boundaryConditions.h"
#include "../InputFunctions/sourceTerm.h"
#include "../InputFunctions/exactSolution.h"
#include "../common/constants.h"
#include "../PolarGrid/polargrid.h"
#include "../Utilities/cmdline.h"
#include "../LinearAlgebra/vector.h"
#include "../LinearAlgebra/matrix.h"
#include "../LinearAlgebra/vector_operations.h"
#include "../Interpolation/interpolation.h"
#include "test_cases.h"

class GMGPolar {
public:
    /* ------------------------ */
    /* GMGPoloar initialization */
    GMGPolar();
    GMGPolar(std::unique_ptr<const DomainGeometry> domain_geometry, 
             std::unique_ptr<const DensityProfileCoefficients> density_profile_coefficients,
             std::unique_ptr<const BoundaryConditions> boundary_conditions,
             std::unique_ptr<const SourceTerm> source_term);

    void setParameters(int argc, char* argv[]);
    void setSolution(std::unique_ptr<const ExactSolution> exact_solution);

    /* ---------------------- */
    /* GMGPolar Setup & Solve */
    void setup();
    void solve();

    /* ----------------- */
    /* GMGPolar Solution */
    Vector<double>& solution();
    const Vector<double>& solution() const;
    const PolarGrid& grid() const;

    /* --------------- */
    /* Grid Parameters */
    double R0() const;
    void R0(double R0);
    double Rmax() const;
    void Rmax(double Rmax);
    
    int nr_exp() const;
    void nr_exp(int nr_exp);
    int ntheta_exp() const;
    void ntheta_exp(int ntheta_exp);
    
    int anisotropic_factor() const;
    void anisotropic_factor(int anisotropic_factor);
    int divideBy2() const;
    void divideBy2(int divideBy2);
    
    bool write_grid_file() const;
    void write_grid_file(bool write_grid_file);
    bool load_grid_file() const;
    void load_grid_file(bool load_grid_file);

    std::string file_grid_radii() const;
    void file_grid_radii(const std::string& file_name);
    std::string file_grid_angles() const;
    void file_grid_angles(const std::string& file_name);

    /* ------------------- */
    /* Geometry Parameters */
    bool DirBC_Interior() const;
    void DirBC_Interior(bool DirBC_Interior);

    /* -------------------- */
    /* Multigrid Parameters */
    bool FMG() const;
    void FMG(bool FMG);
    ExtrapolationType extrapolation() const;
    void extrapolation(ExtrapolationType extrapolation);

    int maxLevels() const;
    void maxLevels(int max_levels);
    MultigridCycleType multigridCycle() const;
    void multigridCycle(MultigridCycleType multigrid_cycle);
    
    int preSmoothingSteps() const;
    void preSmoothingSteps(int pre_smoothing_steps);
    int postSmoothingSteps() const;
    void postSmoothingSteps(int post_smoothing_steps);

    int maxIterations() const;
    void maxIterations(int max_iterations);
    ResidualNormType residualNormType() const;
    void residualNormType(ResidualNormType residual_norm_type);
    double absoluteTolerance() const;
    void absoluteTolerance(double absolute_tolerance);
    double relativeTolerance() const;
    void relativeTolerance(double relative_tolerance);

    /* ------------------ */
    /* Control Parameters */
    int verbose() const;
    void verbose(int verbose);
    bool paraview() const;
    void paraview(bool paraview);
    int maxOpenMPThreads() const;
    void maxOpenMPThreads(int max_omp_threads);
    double threadReductionFactor() const;
    void threadReductionFactor(double thread_reduction_factor);
    ImplementationType implementationType() const;
    void implementationType(ImplementationType implementation_type);
    bool cacheDensityProfileCoefficients() const;
    void cacheDensityProfileCoefficients(bool cache_density_profile_coefficients);
    bool cacheDomainGeometry() const;
    void cacheDomainGeometry(bool cache_domain_geometry);

    /* --------*/
    /* Timings */
    void printTimings() const;
    void resetTimings();

    double t_setup_total;
    double t_setup_createLevels;
    double t_setup_rhs;
    double t_setup_smoother;
    double t_setup_directSolver;

    double t_solve_total;
    double t_solve_initial_approximation;
    double t_solve_multigrid_iterations;
    double t_check_convergence;
    double t_check_exact_error;

    double t_avg_MGC_total;
    double t_avg_MGC_preSmoothing;
    double t_avg_MGC_postSmoothing;
    double t_avg_MGC_residual;
    double t_avg_MGC_directSolver;

private:
    /* --------------- */
    /* Grid Parameters */
    double R0_;
    double Rmax_;
    int nr_exp_;
    int ntheta_exp_;
    int anisotropic_factor_;
    int divideBy2_;
    bool write_grid_file_;
    bool load_grid_file_;
    std::string file_grid_radii_;
    std::string file_grid_angles_;
    /* ------------------- */
    /* Geometry Parameters */
    bool DirBC_Interior_;
    /* ---------- */
    /* Test Cases */
    GeometryType geometry_;
    double kappa_eps_;
    double delta_e_;
    ProblemType problem_;
    AlphaCoeff alpha_;
    double alpha_jump_;
    BetaCoeff beta_;
    /* -------------------- */
    /* Multigrid Parameters */
    bool FMG_;
    int FMG_iterations_;
    MultigridCycleType FMG_cycle_;
    ExtrapolationType extrapolation_;
    int max_levels_;
    int pre_smoothing_steps_;
    int post_smoothing_steps_;
    MultigridCycleType multigrid_cycle_;
    int max_iterations_;
    ResidualNormType residual_norm_type_;
    std::optional<double> absolute_tolerance_;
    std::optional<double> relative_tolerance_;
    /* ------------------ */
    /* Control Parameters */
    int verbose_;
    bool paraview_;
    int max_omp_threads_;
    double thread_reduction_factor_;
    ImplementationType implementation_type_;
    bool cache_density_profile_coefficients_;
    bool cache_domain_geometry_;

    /* --------------- */
    /* Input Functions */
    std::unique_ptr<const DomainGeometry> domain_geometry_;
    std::unique_ptr<const DensityProfileCoefficients> density_profile_coefficients_;
    std::unique_ptr<const BoundaryConditions> boundary_conditions_;
    std::unique_ptr<const SourceTerm> source_term_;
    std::unique_ptr<const ExactSolution> exact_solution_ = nullptr; // Optional exact solution for validation

    /* ------------------------------ */
    /* Parser for GMGPolar parameters */
    cmdline::parser parser_;

    /* ---------------- */
    /* Multigrid levels */
    int number_of_levels_;
    std::vector<Level> levels_;
    std::vector<int> threads_per_level_;

    std::unique_ptr<Interpolation> interpolation_;

    /* Chooses if full grid smoothing is active on level 0 for extrapolation > 0. */
    bool full_grid_smoothing_ = false;

    /* -------------------- */
    /* Convergence criteria */
    int number_of_iterations_;
    std::vector<double> residual_norms_;
    double mean_residual_reduction_factor_;
    bool converged(const double& current_residual_norm, const double& first_residual_norm);

    std::vector<std::pair<double,double>> exact_errors_; // Only when exact solution provided
    std::pair<double, double> computeExactError(Level& level, const Vector<double>& solution, Vector<double>& error);

    /* ---------------------------------------- */
    /* Parser Functions for GMGPolar Parameters */
    void initializeGrid();
    void initializeGeometry();
    void initializeMultigrid();
    void initializeGeneral();

    void parseGrid();
    void parseGeometry();
    void parseMultigrid();
    void parseGeneral();

    void selectTestCase();

    /* --------------- */
    /* Setup Functions */
    PolarGrid createFinestGrid();
    int chooseNumberOfLevels(const PolarGrid& finest_grid);

    void build_rhs_f(const Level& level, Vector<double>& rhs_f);
    void discretize_rhs_f(const Level& level, Vector<double>& rhs_f);

    /* ------------------- */
    /* Multigrid Functions */
    void multigrid_V_Cycle(const int level_depth, Vector<double>& solution, Vector<double>& rhs, Vector<double>& residual);
    void multigrid_W_Cycle(const int level_depth, Vector<double>& solution, Vector<double>& rhs, Vector<double>& residual);
    void multigrid_F_Cycle(const int level_depth, Vector<double>& solution, Vector<double>& rhs, Vector<double>& residual);
    void implicitlyExtrapolatedMultigrid_V_Cycle(const int level_depth, Vector<double>& solution, Vector<double>& rhs, Vector<double>& residual);
    void implicitlyExtrapolatedMultigrid_W_Cycle(const int level_depth, Vector<double>& solution, Vector<double>& rhs, Vector<double>& residual);
    void implicitlyExtrapolatedMultigrid_F_Cycle(const int level_depth, Vector<double>& solution, Vector<double>& rhs, Vector<double>& residual);

    void prolongation(const int current_level, Vector<double>& result, const Vector<double>& x) const;
    void restriction(const int current_level, Vector<double>& result, const Vector<double>& x) const;
    void injection(const int current_level, Vector<double>& result, const Vector<double>& x) const;
    void extrapolatedProlongation(const int current_level, Vector<double>& result, const Vector<double>& x) const;
    void extrapolatedRestriction(const int current_level, Vector<double>& result, const Vector<double>& x) const;
    void FMGInterpolation(const int current_level, Vector<double>& result, const Vector<double>& x) const;

    void extrapolatedResidual(const int current_level, Vector<double>& residual, const Vector<double>& residual_next_level);

    /* ------------- */
    /* Visualization */
    void writeToVTK(const std::filesystem::path& file_path, const PolarGrid& grid);
    void writeToVTK(const std::filesystem::path& file_path, const Level& level, const Vector<double>& grid_function);
};