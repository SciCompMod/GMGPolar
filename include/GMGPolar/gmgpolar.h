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

#include "../InputFunctions/domain_geometry.h"
#include "../InputFunctions/system_parameters.h"
#include "../InputFunctions/exact_solution.h"

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

#include <functional>


class GMGPolar {
public:
    explicit GMGPolar();

    void setRadialRefinement(double r_jump);

    void setGeometry(
        const dFx_dr_Functor& dFx_dr, 
        const dFy_dr_Functor& dFy_dr, 
        const dFx_dt_Functor& dFx_dt, 
        const dFy_dt_Functor& dFy_dt
    );
    void setParameters(
        const alpha_Functor& alpha, 
        const beta_Functor& beta, 
        const rhs_f_Functor& rhs_f, 
        const u_D_Functor& u_D
    );
    void setSystemParameters(const exact_solution_Functor& exact_solution);

    void setParameters(int argc, char* argv[]);
    void setup();
    void solve();

    std::shared_ptr<dFx_dr_Functor> dFx_dr_;
    std::shared_ptr<dFy_dr_Functor> dFy_dr_;
    std::shared_ptr<dFx_dt_Functor> dFx_dt_;
    std::shared_ptr<dFy_dt_Functor> dFy_dt_;

    std::shared_ptr<alpha_Functor> alpha_;
    std::shared_ptr<beta_Functor> beta_;
    std::shared_ptr<rhs_f_Functor> rhs_f_;
    std::shared_ptr<u_D_Functor> u_D_;

    std::shared_ptr<exact_solution_Functor> exact_solution_;


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
    double r_jump_;
    // Geometry Parameters
    alpha_coeff alpha; // UNUSED
    beta_coeff beta; // UNUSED
    problem_type problem; // UNUSED
    geometry_type geometry; // UNUSED
    scalar_t kappa_eps; // UNUSED
    scalar_t delta_e; // UNUSED
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



