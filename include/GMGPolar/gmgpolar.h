#pragma once

class Level;
class LevelCache;

#include "../Level/level.h"

#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/systemParameters.h"
#include "../InputFunctions/exactSolution.h"

#include "../PolarGrid/polargrid.h"

#include "../Utilities/cmdline.h"

#include "../Level/level.h"

#include "../LinearAlgebra/vector.h"
#include "../LinearAlgebra/matrix.h"
#include "../LinearAlgebra/operations.h"

#include <omp.h>
#include <optional>
#include <chrono>
#include <filesystem>

class GMGPolar {
public:
    explicit GMGPolar(const DomainGeometry& domain_geometry, const SystemParameters& system_parameters);

    void setSolution(const ExactSolution& exact_solution);

    void setParameters(int argc, char* argv[]);

    void setup();
    void solve();

    // Grid Parameters
    double R0;
    double Rmax;
    int nr_exp;
    int ntheta_exp;
    int anisotropic_factor;
    int divideBy2;
    bool write_grid_file;
    bool load_grid_file;
    std::string file_grid_r;
    std::string file_grid_theta;
    // Geometry Parameters
    bool DirBC_Interior;
    // Multigrid Parameters
    int extrapolation;
    int maxLevels;
    int v1;
    int v2;
    int cycle;
    // Control Parameters
    int maxOpenMPThreads;
    int finestLevelThreads;
    double threadReductionFactor;

private:
    const DomainGeometry& domain_geometry_;
    const SystemParameters& system_parameters_;
    std::shared_ptr<ExactSolution> exact_solution_ = nullptr; // Optional exact solution for validation

    // Parser for GMGPolar parameters
    cmdline::parser parser_;

    // Multigrid levels
    int numberOflevels_;
    std::vector<Level> levels_;

    // Setup functions
    PolarGrid createFinestGrid();
    int chooseNumberOfLevels(const PolarGrid& finest_grid);

    
    void write_to_vtk(const std::filesystem::path& file_path, const PolarGrid& grid, const Vector<double>& grid_function, const LevelCache& leveldata);

    // Parser functions for GMGPolar parameters
    // Implementation in GMGPolar/parser.cpp
    void initializeGrid();
    void initializeGeometry();
    void initializeMultigrid();
    void initializeGeneral();
    void parseGrid();
    void parseGeometry();
    void parseMultigrid();
    void parseGeneral();
};