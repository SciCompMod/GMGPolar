#pragma once

#include <cmath>

#include "../sourceTerm.h"

#include "../../PolarGrid/polargrid.h"

namespace gmgpolar
{

class CartesianR6_ZoniGyro_ShafranovGeometry
{
public:
    explicit CartesianR6_ZoniGyro_ShafranovGeometry(PolarGrid const& grid, double Rmax, double elongation_kappa,
                                                    double shift_delta);

    double operator()(std::size_t i_r, std::size_t i_theta) const;

private:
    PolarGrid const& grid_;
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};
} // namespace gmgpolar
