#pragma once

#include <cmath>

#include "../sourceTerm.h"

#include "../../PolarGrid/polargrid.h"

namespace gmgpolar
{

class CartesianR6_Zoni_CircularGeometry
{
public:
    explicit CartesianR6_Zoni_CircularGeometry(PolarGrid<Kokkos::HostSpace> const& grid, double Rmax);
    KOKKOS_DEFAULTED_FUNCTION CartesianR6_Zoni_CircularGeometry(const CartesianR6_Zoni_CircularGeometry&) = default;

    KOKKOS_FUNCTION double operator()(std::size_t i_r, std::size_t i_theta) const;

private:
    PolarGrid<Kokkos::HostSpace> grid_;
    const double Rmax = 1.3;
};
} // namespace gmgpolar
