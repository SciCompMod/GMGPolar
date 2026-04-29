#include <InputFunctions/DomainGeometry/shafranovGeometry.h>
using namespace gmgpolar;

ShafranovGeometry::ShafranovGeometry(double Rmax, double elongation_kappa, double shift_delta)
    : Rmax(Rmax)
    , elongation_kappa(elongation_kappa)
    , shift_delta(shift_delta)
{
}