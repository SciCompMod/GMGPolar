#pragma once

#include <cmath>

#include "../sourceTerm.h"

#include "../../PolarGrid/polargrid.h"

class PolarR6_ZoniShifted_CircularGeometry
{
public:
    explicit PolarR6_ZoniShifted_CircularGeometry(PolarGrid const& grid, double Rmax);

    double operator()(std::size_t i_r, std::size_t i_theta) const;

private:
    PolarGrid const& grid_;
    const double Rmax = 1.3;
};
