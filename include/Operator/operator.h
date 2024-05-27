#pragma once

class Level;
class GMGPolar;

#include <vector>
#include <cassert>
#include <cmath>
#include <thread>
#include <mutex>

#include <omp.h>

#include "../include/Level/level.h"
#include "../include/PolarGrid/polargrid.h"
#include "../common/scalar.h"
#include "../include/TaskDistribution/TaskDistribution.h"
#include "../include/linear_algebra/vector.h"
#include "../include/linear_algebra/matrix.h"
#include "../include/linear_algebra/operations.h"
#include "../include/GMGPolar/gmgpolar.h"
#include "../include/ExactFunctions/exactfunctions.h"
#include "../include/InputFunctions/domain_geometry.h"

#include "../include/InputFunctions/domain_geometry.h"
#include "../include/InputFunctions/system_parameters.h"

class Operator {
public:    
    Operator(const GMGPolar& gmgpolar, const PolarGrid& grid);

    void applyATake0(const Level& onLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x, const scalar_t& scaleAx) const;

    void applyAGive(const Level& onLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x, const scalar_t& scaleAx) const;
    void applyAGiveTasks(const Level& onLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x, const scalar_t& scaleAx) const;
    void applyAGiveMutex(Level& onLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x, const scalar_t& scaleAx);

    void arr_att_art(
        const double& r, const double& theta, 
        const double& sin_theta, const double& cos_theta, 
        const double& coeff_alpha, 
        double& arr, double& att, double& art, double& detDF
    ) const;

    void applySmoothing(const Level& onLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const;

    std::pair<SparseMatrix<scalar_t>, Vector<scalar_t>> build_system(const Level& onLevel) const;

private:
    bool DirBC_Interior_;
    double kappa_eps_;
    double delta_e_;
    double Rmax_;
    geometry_type geometry_;
    std::vector<scalar_t> sin_theta_;
    std::vector<scalar_t> cos_theta_;

    std::shared_ptr<dFx_dr_Functor> dFx_dr_;
    std::shared_ptr<dFy_dr_Functor> dFy_dr_;
    std::shared_ptr<dFx_dt_Functor> dFx_dt_;
    std::shared_ptr<dFy_dt_Functor> dFy_dt_;

    std::shared_ptr<alpha_Functor> alpha_;
    std::shared_ptr<beta_Functor> beta_;

    // Directed graph dependencies for applyAGiveMutex

    int numberMod0Circles;
    int numberMod1Circles;
    int numberMod2Circles;

    int numberDiv3Radials;

    int additional_circle_task;
    int additional_radial_task;

    std::vector<std::mutex> Mod1CircleMutexes;
    std::vector<int> Mod1CircleDepCounter;

    std::vector<std::mutex> Mod2CircleMutexes;
    std::vector<int> Mod2CircleDepCounter;

    std::vector<std::mutex> Mod1RadialMutexes;
    std::vector<int> Mod1RadialDepCounter;

    std::vector<std::mutex> Mod2RadialMutexes;
    std::vector<int> Mod2RadialDepCounter;


//     int Circles = 300;
//     int BlackCircles = 150;
//     int WhiteCircles = 150;
//     // BlackCirclesTaskIndex [8,6,4,2,0]
//     // WhiteCircleTaskIndex [9,7,5,3,1]

//     int Radials = 200;
//     int BlackRadials = 100;
//     int WhiteRadials = 100;
//     // BlackRadialTaskIndex [10,12,14,16,18]
//     // WhiteRadialTaskIndex [11,13,15,17,19]

//     assert(BlackCircles == WhiteCircles || BlackCircles == WhiteCircles + 1);
//     std::vector<std::mutex> WhiteCircleMutexes(WhiteCircles);
//     std::vector<int> WhiteCircleDepCounter(WhiteCircles, 2);
//     if(BlackCircles == WhiteCircles) WhiteCircleDepCounter.back() = 1;


//     assert(BlackRadials == WhiteRadials);
//     std::vector<std::mutex> WhiteRadialMutexes(WhiteRadials);
//     std::vector<int> WhiteRadialDepCounter(WhiteRadials, 2);


};

#include "operator.inl"
