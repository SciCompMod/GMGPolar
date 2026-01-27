#pragma once
#include <variant>
#include "../../include/GMGPolar/test_cases.h"
#include "../../include/InputFunctions/BoundaryConditions/cartesianR2_Boundary_CircularGeometry.h"
#include "../../include/InputFunctions/BoundaryConditions/cartesianR6_Boundary_CzarnyGeometry.h"
#include "../../include/InputFunctions/BoundaryConditions/polarR6_Boundary_CzarnyGeometry.h"
#include "../../include/InputFunctions/BoundaryConditions/refined_Boundary_CzarnyGeometry.h"
#include "../../include/InputFunctions/BoundaryConditions/cartesianR2_Boundary_CzarnyGeometry.h"
#include "../../include/InputFunctions/BoundaryConditions/cartesianR6_Boundary_ShafranovGeometry.h"
#include "../../include/InputFunctions/BoundaryConditions/polarR6_Boundary_ShafranovGeometry.h"
#include "../../include/InputFunctions/BoundaryConditions/refined_Boundary_ShafranovGeometry.h"
#include "../../include/InputFunctions/BoundaryConditions/cartesianR2_Boundary_ShafranovGeometry.h"
#include "../../include/InputFunctions/BoundaryConditions/polarR6_Boundary_CircularGeometry.h"
#include "../../include/InputFunctions/BoundaryConditions/refined_Boundary_CircularGeometry.h"
#include "../../include/InputFunctions/BoundaryConditions/cartesianR6_Boundary_CircularGeometry.h"
#include "../../include/InputFunctions/BoundaryConditions/polarR6_Boundary_CulhamGeometry.h"
#include "../../include/InputFunctions/BoundaryConditions/refined_Boundary_CulhamGeometry.h"

using DomainGeometryVariant = std::variant<CircularGeometry, ShafranovGeometry, CzarnyGeometry, CulhamGeometry>;

using DensityProfileCoefficientsVariant =
    std::variant<PoissonCoefficients, SonnendruckerCoefficients, SonnendruckerGyroCoefficients, ZoniCoefficients,
                 ZoniGyroCoefficients, ZoniShiftedCoefficients, ZoniShiftedGyroCoefficients>;

using BoundaryConditionsVariant = std::variant<
    CartesianR2_Boundary_CircularGeometry, CartesianR6_Boundary_CzarnyGeometry, PolarR6_Boundary_CzarnyGeometry,
    Refined_Boundary_CzarnyGeometry, CartesianR2_Boundary_CzarnyGeometry, CartesianR6_Boundary_ShafranovGeometry,
    PolarR6_Boundary_ShafranovGeometry, Refined_Boundary_ShafranovGeometry, CartesianR2_Boundary_ShafranovGeometry,
    PolarR6_Boundary_CircularGeometry, Refined_Boundary_CircularGeometry, CartesianR6_Boundary_CircularGeometry,
    PolarR6_Boundary_CulhamGeometry, Refined_Boundary_CulhamGeometry>;
