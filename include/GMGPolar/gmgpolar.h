#pragma once

class Level;
class Interpolation;
class Operator;

#include <limits>
#include <cassert>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <stdexcept>

#include "../common/scalar.h"
#include "../common/constants.h"
#include "../utilities/display_util.h"
#include "../utilities/cmdline.h"
#include "../linear_algebra/vector.h"
#include "../PolarGrid/polargrid.h"
#include "../Level/level.h"
#include "../Interpolation/interpolation.h"

#include "../ExactFunctions/exactfunctions.h"

#include "../../include/test_cases/CartesianR2GyroSonnendruckerCircular.h"
#include "../../include/test_cases/CartesianR2GyroSonnendruckerShafranov.h"
#include "../../include/test_cases/CartesianR2GyroSonnendruckerTriangular.h"
#include "../../include/test_cases/PolarR6GyroSonnendruckerCircular.h"
#include "../../include/test_cases/PolarR6GyroSonnendruckerShafranov.h"
#include "../../include/test_cases/PolarR6GyroSonnendruckerTriangular.h"
#include "../../include/test_cases/CartesianR6GyroSonnendruckerCircular.h"
#include "../../include/test_cases/CartesianR6GyroSonnendruckerShafranov.h"
#include "../../include/test_cases/CartesianR6GyroSonnendruckerTriangular.h"
#include "../../include/test_cases/CartesianR2GyroZoniCircular.h"
#include "../../include/test_cases/CartesianR2GyroZoniShafranov.h"
#include "../../include/test_cases/CartesianR2GyroZoniTriangular.h"
#include "../../include/test_cases/PolarR6GyroZoniCircular.h"
#include "../../include/test_cases/PolarR6GyroZoniShafranov.h"
#include "../../include/test_cases/PolarR6GyroZoniTriangular.h"
#include "../../include/test_cases/CartesianR6GyroZoniCircular.h"
#include "../../include/test_cases/CartesianR6GyroZoniShafranov.h"
#include "../../include/test_cases/CartesianR6GyroZoniTriangular.h"
#include "../../include/test_cases/CartesianR2GyroZoniShiftedCircular.h"
#include "../../include/test_cases/CartesianR2GyroZoniShiftedShafranov.h"
#include "../../include/test_cases/CartesianR2GyroZoniShiftedTriangular.h"
#include "../../include/test_cases/PolarR6GyroZoniShiftedCircular.h"
#include "../../include/test_cases/PolarR6GyroZoniShiftedShafranov.h"
#include "../../include/test_cases/PolarR6GyroZoniShiftedTriangular.h"
#include "../../include/test_cases/PolarR6GyroZoniShiftedCulham.h"
#include "../../include/test_cases/CartesianR6GyroZoniShiftedCircular.h"
#include "../../include/test_cases/CartesianR6GyroZoniShiftedShafranov.h"
#include "../../include/test_cases/CartesianR6GyroZoniShiftedTriangular.h"
#include "../../include/test_cases/CartesianR2SonnendruckerCircular.h"
#include "../../include/test_cases/CartesianR2SonnendruckerShafranov.h"
#include "../../include/test_cases/CartesianR2SonnendruckerTriangular.h"
#include "../../include/test_cases/PolarR6SonnendruckerCircular.h"
#include "../../include/test_cases/PolarR6SonnendruckerShafranov.h"
#include "../../include/test_cases/PolarR6SonnendruckerTriangular.h"
#include "../../include/test_cases/CartesianR6SonnendruckerCircular.h"
#include "../../include/test_cases/CartesianR6SonnendruckerShafranov.h"
#include "../../include/test_cases/CartesianR6SonnendruckerTriangular.h"
#include "../../include/test_cases/CartesianR2ZoniCircular.h"
#include "../../include/test_cases/CartesianR2ZoniShafranov.h"
#include "../../include/test_cases/CartesianR2ZoniTriangular.h"
#include "../../include/test_cases/PolarR6ZoniCircular.h"
#include "../../include/test_cases/PolarR6ZoniShafranov.h"
#include "../../include/test_cases/PolarR6ZoniTriangular.h"
#include "../../include/test_cases/CartesianR6ZoniCircular.h"
#include "../../include/test_cases/CartesianR6ZoniShafranov.h"
#include "../../include/test_cases/CartesianR6ZoniTriangular.h"
#include "../../include/test_cases/CartesianR2ZoniShiftedCircular.h"
#include "../../include/test_cases/CartesianR2ZoniShiftedShafranov.h"
#include "../../include/test_cases/CartesianR2ZoniShiftedTriangular.h"
#include "../../include/test_cases/PolarR6ZoniShiftedCircular.h"
#include "../../include/test_cases/PolarR6ZoniShiftedShafranov.h"
#include "../../include/test_cases/PolarR6ZoniShiftedTriangular.h"
#include "../../include/test_cases/CartesianR6ZoniShiftedCircular.h"
#include "../../include/test_cases/CartesianR6ZoniShiftedShafranov.h"
#include "../../include/test_cases/CartesianR6ZoniShiftedTriangular.h"
#include "../../include/test_cases/CartesianR2PoissonCircular.h"
#include "../../include/test_cases/CartesianR2PoissonShafranov.h"
#include "../../include/test_cases/CartesianR2PoissonTriangular.h"
#include "../../include/test_cases/PolarR6PoissonCircular.h"
#include "../../include/test_cases/PolarR6PoissonShafranov.h"
#include "../../include/test_cases/PolarR6PoissonTriangular.h"
#include "../../include/test_cases/CartesianR6PoissonCircular.h"
#include "../../include/test_cases/CartesianR6PoissonShafranov.h"
#include "../../include/test_cases/CartesianR6PoissonTriangular.h"
#include "../../include/test_cases/RefinedGyroZoniShiftedCircular.h"
#include "../../include/test_cases/RefinedGyroZoniShiftedShafranov.h"
#include "../../include/test_cases/RefinedGyroZoniShiftedTriangular.h"
#include "../../include/test_cases/RefinedGyroZoniShiftedCulham.h"

#include <random>
#include <chrono>

class GMGPolar {
public:
    GMGPolar();
    void setParameters(int argc, char* argv[]);
    void setup();
    void solve();

    // Grid Parameters
    scalar_t R0;
    scalar_t Rmax;
    int nr_exp;
    int ntheta_exp;
    int anisotropic_factor;
    int divideBy2;
    bool write_grid_file;
    bool load_grid_file;
    std::string file_grid_r;
    std::string file_grid_theta;
    // Geometry Parameters
    alpha_coeff alpha;
    beta_coeff beta;
    problem_type problem;
    geometry_type geometry;
    scalar_t kappa_eps;
    scalar_t delta_e;
    // Multigrid Parameters
    int extrapolation;
    int maxLevels;
    int v1;
    int v2;
    int cycle;
    // Control Parameters
    int verbose;
private:
    // GMGPolar Parameters
    cmdline::parser parser_;

    std::vector<Level> levels_;
    int numberOflevels_;

    std::unique_ptr<Interpolation> interpolation_;

    const Level& getLevel(const int index) const;

    PolarGrid createFinestGrid();
    int numberOfLevels(const PolarGrid& finest_grid);
    void initializeLevels(const int levels, PolarGrid& finest_grid);
    std::shared_ptr<ExactFunctions> selectExactFunctionsClass();

    void prolongateToNextLevel(const int current_level, Vector<scalar_t> &result, const Vector<scalar_t> &x) const;
    void restrictToLowerLevel(const int current_level, Vector<scalar_t> &result, const Vector<scalar_t> &x) const;

    // Implementation in GMGPolar/gmgpolar_parser.cpp
    enum {
        OPTIONAL = 0,
        REQUIRED = 1
    };
    void initializeGrid();
    void initializeGeometry();
    void initializeMultigrid();
    void initializeGeneral();
    void parseGrid();
    void parseGeometry();
    void parseMultigrid();
    void parseGeneral();
};



