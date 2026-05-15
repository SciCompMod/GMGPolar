#pragma once

#include <cmath>

#include "../sourceTerm.h"

#include "../../PolarGrid/polargrid.h"

namespace gmgpolar
{

class CartesianR2_Zoni_CircularGeometry
{
public:
    explicit CartesianR2_Zoni_CircularGeometry(PolarGrid const& grid, double Rmax);
    KOKKOS_DEFAULTED_FUNCTION CartesianR2_Zoni_CircularGeometry(const CartesianR2_Zoni_CircularGeometry&) = default;

    KOKKOS_FUNCTION double operator()(std::size_t i_r, std::size_t i_theta) const;

private:
    PolarGrid grid_;
    const double Rmax = 1.3;
};
} // namespace gmgpolar
