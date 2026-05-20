#include "../include/InputFunctions/SourceTerms/polarR6_ZoniShiftedGyro_CulhamGeometry.h"
using namespace gmgpolar;

PolarR6_ZoniShiftedGyro_CulhamGeometry::PolarR6_ZoniShiftedGyro_CulhamGeometry(PolarGrid<Kokkos::HostSpace> const& grid,
                                                                               double Rmax)
    : grid_(grid)
    , Rmax(Rmax)
{
}

KOKKOS_FUNCTION double PolarR6_ZoniShiftedGyro_CulhamGeometry::operator()(std::size_t i_r, std::size_t i_theta) const
{
    double r     = grid_.radius(i_r);
    double theta = grid_.theta(i_theta);
    return 0.4096 * pow((r / Rmax), 6.0) * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}
