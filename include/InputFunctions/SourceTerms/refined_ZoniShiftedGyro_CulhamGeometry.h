#pragma once

#include <cmath>

#include "../sourceTerm.h"

#include "../../PolarGrid/polargrid.h"

namespace gmgpolar
{

class Refined_ZoniShiftedGyro_CulhamGeometry
{
public:
    explicit Refined_ZoniShiftedGyro_CulhamGeometry(PolarGrid const& grid, double Rmax);

    double operator()(std::size_t i_r, std::size_t i_theta) const;

private:
    PolarGrid const& grid_;
    const double Rmax = 1.3;
};
} // namespace gmgpolar
