#pragma once

#include <cmath>

#include "../sourceTerm.h"

#include "../../PolarGrid/polargrid.h"

namespace gmgpolar
{

class CartesianR2_SonnendruckerGyro_CircularGeometry
{
public:
    explicit CartesianR2_SonnendruckerGyro_CircularGeometry(PolarGrid const& grid, double Rmax);

    KOKKOS_FUNCTION double operator()(std::size_t i_r, std::size_t i_theta) const;

private:
    PolarGrid grid_;
    const double Rmax = 1.3;
};
} // namespace gmgpolar
