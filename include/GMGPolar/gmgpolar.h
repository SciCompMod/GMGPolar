#pragma once

class Level;
class LevelCache;

#include "../Level/level.h"

#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/systemParameters.h"
#include "../InputFunctions/exactSolution.h"

#include "../common/constants.h"

#include "../PolarGrid/polargrid.h"

#include "../Utilities/cmdline.h"

#include "../Level/level.h"

#include "../LinearAlgebra/vector.h"
#include "../LinearAlgebra/matrix.h"
#include "../LinearAlgebra/operations.h"

#include "../Interpolation/interpolation.h"

#include <iostream>
#include <memory>
#include <filesystem>
#include <optional>
#include <utility>
#include <chrono>

#include <omp.h>

#include "test_cases.h"

class GMGPolar {
public:
    /* ------------------------ */
    /* GMGPoloar initialization */
    GMGPolar();
    GMGPolar(std::unique_ptr<const DomainGeometry> domain_geometry, 
             std::unique_ptr<const SystemParameters> system_parameters);

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
    int extrapolation() const;
    void extrapolation(int extrapolation);

    int maxLevels() const;
    void maxLevels(int maxLevels);
    MultigridCycleType multigrid_cycle() const;
    void multigrid_cycle(MultigridCycleType multigrid_cycle);
    
    int preSmoothingSteps() const;
    void preSmoothingSteps(int preSmoothingSteps);
    int postSmoothingSteps() const;
    void postSmoothingSteps(int postSmoothingSteps);

    int maxIterations() const;
    void maxIterations(int maxIterations);
    ResidualNormType residualNormType() const;
    void residualNormType(ResidualNormType residualNormType);
    double absoluteTolerance() const;
    void absoluteTolerance(double absoluteTolerance);
    double relativeTolerance() const;
    void relativeTolerance(double relativeTolerance);

    /* ------------------ */
    /* Control Parameters */
    int maxOpenMPThreads() const;
    void maxOpenMPThreads(int maxOpenMPThreads);
    int finestLevelThreads() const;
    void finestLevelThreads(int finestLevelThreads);
    double threadReductionFactor() const;
    void threadReductionFactor(double threadReductionFactor);

    /* --------*/
    /* Timings */
    void printTimings() const;

    double t_setup_total = 0.0;
    double t_setup_createLevels = 0.0;
    double t_setup_rhs = 0.0;
    double t_setup_smoother = 0.0;
    double t_setup_directSolver = 0.0;

    double t_solve_total = 0.0;
    double t_solve_multigrid_iterations = 0.0;
    double t_check_convergence = 0.0;
    double t_check_exact_error = 0.0;

    double t_avg_MGC_total = 0.0;
    double t_avg_MGC_preSmoothing = 0.0;
    double t_avg_MGC_postSmoothing = 0.0;
    double t_avg_MGC_residual = 0.0;
    double t_avg_MGC_directSolver = 0.0;

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
    int extrapolation_;
    int maxLevels_;
    int preSmoothingSteps_;
    int postSmoothingSteps_;
    MultigridCycleType multigrid_cycle_;
    int max_iterations_;
    ResidualNormType residual_norm_type_;
    std::optional<double> absolute_tolerance_;
    std::optional<double> relative_tolerance_;
    /* ------------------ */
    /* Control Parameters */
    int maxOpenMPThreads_;
    int finestLevelThreads_;
    double threadReductionFactor_;

    /* ------------------------------ */
    /* Parser for GMGPolar parameters */
    cmdline::parser parser_;

    /* --------------- */
    /* Input Functions */
    std::unique_ptr<const DomainGeometry> domain_geometry_;
    std::unique_ptr<const SystemParameters> system_parameters_;
    std::unique_ptr<const ExactSolution> exact_solution_ = nullptr; // Optional exact solution for validation

    /* ---------------- */
    /* Multigrid levels */
    int numberOflevels_;
    std::vector<Level> levels_;
    std::vector<int> taskingThreads_; // Optimal number of threads for OpenMP tasks on each level

    std::unique_ptr<Interpolation> interpolation_;

    /* -------------------- */
    /* Convergence criteria */
    int number_of_iterations_;
    std::vector<double> residual_norms_;
    double mean_residual_reduction_factor_rho_;
    bool converged(const double& current_residual_norm, const double& first_residual_norm);

    std::vector<std::pair<double,double>> exact_errors_; // Only when exact solution provided
    std::pair<double, double> compute_exact_error(Level& level, const Vector<double>& solution, Vector<double>& error);

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
    void implicitly_extrapolated_multigrid_V_Cycle(const int level_depth, Vector<double>& solution, Vector<double>& rhs, Vector<double>& residual);
    void implicitly_extrapolated_multigrid_W_Cycle(const int level_depth, Vector<double>& solution, Vector<double>& rhs, Vector<double>& residual);
    void implicitly_extrapolated_multigrid_F_Cycle(const int level_depth, Vector<double>& solution, Vector<double>& rhs, Vector<double>& residual);

    void prolongateToUpperLevel(const int current_level, Vector<double>& result, const Vector<double>& x) const;
    void restrictToLowerLevel(const int current_level, Vector<double>& result, const Vector<double>& x) const;
    void injectToLowerLevel(const int current_level, Vector<double>& result, const Vector<double>& x) const;
    void extrapolation_prolongateToUpperLevel(const int current_level, Vector<double>& result, const Vector<double>& x) const;
    void extrapolation_restrictToLowerLevel(const int current_level, Vector<double>& result, const Vector<double>& x) const;

    void extrapolated_residual(const int current_level, Vector<double>& residual, const Vector<double>& residual_next_level);

    /* ------------- */
    /* Visualization */
    void write_to_vtk(const std::filesystem::path& file_path, const PolarGrid& grid);
    void write_to_vtk(const std::filesystem::path& file_path, const Level& level, const Vector<double>& grid_function);

    void resetTimings();
};