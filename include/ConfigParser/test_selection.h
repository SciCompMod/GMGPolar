#pragma once

#include <variant>

#include "../../include/GMGPolar/test_cases.h"

namespace gmgpolar
{

using DomainGeometryVariant = std::variant<CircularGeometry, ShafranovGeometry, CzarnyGeometry, CulhamGeometry>;

using DensityProfileCoefficientsVariant =
    std::variant<PoissonCoefficients, SonnendruckerCoefficients, SonnendruckerGyroCoefficients, ZoniCoefficients,
                 ZoniGyroCoefficients, ZoniShiftedCoefficients, ZoniShiftedGyroCoefficients>;

using ExactSolutionVariant =
    std::variant<CartesianR2_CircularGeometry, CartesianR2_CzarnyGeometry, CartesianR2_ShafranovGeometry,
                 CartesianR6_CircularGeometry, CartesianR6_CzarnyGeometry, CartesianR6_ShafranovGeometry,
                 PolarR6_CircularGeometry, PolarR6_CulhamGeometry, PolarR6_CzarnyGeometry, PolarR6_ShafranovGeometry,
                 Refined_CircularGeometry, Refined_CulhamGeometry, Refined_CzarnyGeometry, Refined_ShafranovGeometry>;

} // namespace gmgpolar
