#pragma once

#include <cmath>

#include "../sourceTerm.h"

#include "../../PolarGrid/polargrid.h"

class CartesianR2_ZoniShiftedGyro_CircularGeometry
{
public:
    explicit CartesianR2_ZoniShiftedGyro_CircularGeometry(PolarGrid const& grid, double Rmax);

    double operator()(std::size_t i_r, std::size_t i_theta) const;

private:
    PolarGrid const& grid_;
    const double Rmax = 1.3;
};
