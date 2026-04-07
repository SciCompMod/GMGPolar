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

#include "../../include/InputFunctions/SourceTerms/cartesianR2_Poisson_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_Poisson_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_Poisson_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_SonnendruckerGyro_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_SonnendruckerGyro_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_SonnendruckerGyro_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_Sonnendrucker_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_Sonnendrucker_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_Sonnendrucker_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_ZoniGyro_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_ZoniGyro_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_ZoniGyro_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_ZoniShiftedGyro_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_ZoniShiftedGyro_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_ZoniShiftedGyro_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_ZoniShifted_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_ZoniShifted_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_ZoniShifted_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_Zoni_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_Zoni_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR2_Zoni_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_Poisson_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_Poisson_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_Poisson_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_SonnendruckerGyro_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_SonnendruckerGyro_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_SonnendruckerGyro_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_Sonnendrucker_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_Sonnendrucker_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_Sonnendrucker_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_ZoniGyro_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_ZoniGyro_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_ZoniGyro_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_ZoniShiftedGyro_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_ZoniShiftedGyro_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_ZoniShiftedGyro_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_ZoniShifted_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_ZoniShifted_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_ZoniShifted_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_Zoni_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_Zoni_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/cartesianR6_Zoni_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_Poisson_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_Poisson_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_Poisson_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_SonnendruckerGyro_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_SonnendruckerGyro_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_SonnendruckerGyro_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_Sonnendrucker_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_Sonnendrucker_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_Sonnendrucker_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_ZoniGyro_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_ZoniGyro_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_ZoniGyro_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_ZoniShiftedGyro_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_ZoniShiftedGyro_CulhamGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_ZoniShiftedGyro_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_ZoniShiftedGyro_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_ZoniShifted_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_ZoniShifted_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_ZoniShifted_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_Zoni_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_Zoni_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/polarR6_Zoni_ShafranovGeometry.h"
#include "../../include/InputFunctions/SourceTerms/refined_ZoniShiftedGyro_CircularGeometry.h"
#include "../../include/InputFunctions/SourceTerms/refined_ZoniShiftedGyro_CulhamGeometry.h"
#include "../../include/InputFunctions/SourceTerms/refined_ZoniShiftedGyro_CzarnyGeometry.h"
#include "../../include/InputFunctions/SourceTerms/refined_ZoniShiftedGyro_ShafranovGeometry.h"

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

using SourceTermVariant = std::variant<
    CartesianR2_Poisson_CircularGeometry, CartesianR2_Poisson_CzarnyGeometry, CartesianR2_Poisson_ShafranovGeometry,
    CartesianR2_SonnendruckerGyro_CircularGeometry, CartesianR2_SonnendruckerGyro_CzarnyGeometry,
    CartesianR2_SonnendruckerGyro_ShafranovGeometry, CartesianR2_Sonnendrucker_CircularGeometry,
    CartesianR2_Sonnendrucker_CzarnyGeometry, CartesianR2_Sonnendrucker_ShafranovGeometry,
    CartesianR2_ZoniGyro_CircularGeometry, CartesianR2_ZoniGyro_CzarnyGeometry, CartesianR2_ZoniGyro_ShafranovGeometry,
    CartesianR2_ZoniShiftedGyro_CircularGeometry, CartesianR2_ZoniShiftedGyro_CzarnyGeometry,
    CartesianR2_ZoniShiftedGyro_ShafranovGeometry, CartesianR2_ZoniShifted_CircularGeometry,
    CartesianR2_ZoniShifted_CzarnyGeometry, CartesianR2_ZoniShifted_ShafranovGeometry,
    CartesianR2_Zoni_CircularGeometry, CartesianR2_Zoni_CzarnyGeometry, CartesianR2_Zoni_ShafranovGeometry,
    CartesianR6_Poisson_CircularGeometry, CartesianR6_Poisson_CzarnyGeometry, CartesianR6_Poisson_ShafranovGeometry,
    CartesianR6_SonnendruckerGyro_CircularGeometry, CartesianR6_SonnendruckerGyro_CzarnyGeometry,
    CartesianR6_SonnendruckerGyro_ShafranovGeometry, CartesianR6_Sonnendrucker_CircularGeometry,
    CartesianR6_Sonnendrucker_CzarnyGeometry, CartesianR6_Sonnendrucker_ShafranovGeometry,
    CartesianR6_ZoniGyro_CircularGeometry, CartesianR6_ZoniGyro_CzarnyGeometry, CartesianR6_ZoniGyro_ShafranovGeometry,
    CartesianR6_ZoniShiftedGyro_CircularGeometry, CartesianR6_ZoniShiftedGyro_CzarnyGeometry,
    CartesianR6_ZoniShiftedGyro_ShafranovGeometry, CartesianR6_ZoniShifted_CircularGeometry,
    CartesianR6_ZoniShifted_CzarnyGeometry, CartesianR6_ZoniShifted_ShafranovGeometry,
    CartesianR6_Zoni_CircularGeometry, CartesianR6_Zoni_CzarnyGeometry, CartesianR6_Zoni_ShafranovGeometry,
    PolarR6_Poisson_CircularGeometry, PolarR6_Poisson_CzarnyGeometry, PolarR6_Poisson_ShafranovGeometry,
    PolarR6_SonnendruckerGyro_CircularGeometry, PolarR6_SonnendruckerGyro_CzarnyGeometry,
    PolarR6_SonnendruckerGyro_ShafranovGeometry, PolarR6_Sonnendrucker_CircularGeometry,
    PolarR6_Sonnendrucker_CzarnyGeometry, PolarR6_Sonnendrucker_ShafranovGeometry, PolarR6_ZoniGyro_CircularGeometry,
    PolarR6_ZoniGyro_CzarnyGeometry, PolarR6_ZoniGyro_ShafranovGeometry, PolarR6_ZoniShiftedGyro_CircularGeometry,
    PolarR6_ZoniShiftedGyro_CulhamGeometry, PolarR6_ZoniShiftedGyro_CzarnyGeometry,
    PolarR6_ZoniShiftedGyro_ShafranovGeometry, PolarR6_ZoniShifted_CircularGeometry, PolarR6_ZoniShifted_CzarnyGeometry,
    PolarR6_ZoniShifted_ShafranovGeometry, PolarR6_Zoni_CircularGeometry, PolarR6_Zoni_CzarnyGeometry,
    PolarR6_Zoni_ShafranovGeometry, Refined_ZoniShiftedGyro_CircularGeometry, Refined_ZoniShiftedGyro_CulhamGeometry,
    Refined_ZoniShiftedGyro_CzarnyGeometry, Refined_ZoniShiftedGyro_ShafranovGeometry>;