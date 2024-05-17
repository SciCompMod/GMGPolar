#pragma once

class Level;
class GMGPolar;

#include <vector>
#include <cassert>
#include <cmath>

#include <omp.h>

#include "../include/Level/level.h"
#include "../include/PolarGrid/polargrid.h"
#include "../common/scalar.h"
#include "../include/TaskDistribution/TaskDistribution.h"
#include "../include/linear_algebra/vector.h"
#include "../include/linear_algebra/operations.h"
#include "../include/GMGPolar/gmgpolar.h"
#include "../include/ExactFunctions/exactfunctions.h"

class Operator {
public:    
    Operator(const GMGPolar& gmgpolar, const PolarGrid& grid);

    void applyATake0(const Level& onLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const;
    void applyATake(const Level& onLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const;
    void applyAGive(const Level& onLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const;
    void applyAGiveTasks(const Level& onLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const;

    void arr_att_art(const ExactFunctions& exactFuncs, double r, double theta, int i_theta, double coeff_alpha, double& arr, double& att, double& art, double& detDFinv) const;

    // void applyA_orthogonal_Give(const Level& onLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x, const int smoother) const;

private:
    double kappa_eps_;
    double delta_e_;
    double Rmax_;
    geometry_type geometry_;
    std::vector<scalar_t> sin_theta_;
    std::vector<scalar_t> cos_theta_;
};


