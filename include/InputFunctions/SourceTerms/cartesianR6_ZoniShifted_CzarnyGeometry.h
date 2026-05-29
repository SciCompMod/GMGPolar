#pragma once

#include <cmath>

#include "../sourceTerm.h"

#include "../../PolarGrid/polargrid.h"

namespace gmgpolar
{

class CartesianR6_ZoniShifted_CzarnyGeometry
{
public:
    explicit CartesianR6_ZoniShifted_CzarnyGeometry(PolarGrid<Kokkos::HostSpace> const& grid, double Rmax,
                                                    double inverse_aspect_ratio_epsilon, double ellipticity_e);
    KOKKOS_DEFAULTED_FUNCTION
    CartesianR6_ZoniShifted_CzarnyGeometry(const CartesianR6_ZoniShifted_CzarnyGeometry&) = default;

    KOKKOS_FUNCTION double operator()(std::size_t i_r, std::size_t i_theta) const;

private:
    PolarGrid<Kokkos::HostSpace> grid_;
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;

    void initializeGeometry();
    double factor_xi;
};
} // namespace gmgpolar
