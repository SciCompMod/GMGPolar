#pragma once

#include "../InputFunctions/boundaryConditions.h"
#include "../InputFunctions/densityProfileCoefficients.h"
#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/exactSolution.h"
#include "../InputFunctions/sourceTerm.h"

/* ---------- */
/* Test Cases */
/* ---------- */

/* --------------- */
/* Domain Geometry */
#include <InputFunctions/DomainGeometry/circularGeometry.h>
#include <InputFunctions/DomainGeometry/culhamGeometry.h>
#include <InputFunctions/DomainGeometry/czarnyGeometry.h>
#include <InputFunctions/DomainGeometry/shafranovGeometry.h>

/* --------------- */
/* Exact Solutions */
#include <InputFunctions/ExactSolution/cartesianR2_CircularGeometry.h>
#include <InputFunctions/ExactSolution/cartesianR2_CzarnyGeometry.h>
#include <InputFunctions/ExactSolution/cartesianR2_ShafranovGeometry.h>
#include <InputFunctions/ExactSolution/cartesianR6_CircularGeometry.h>
#include <InputFunctions/ExactSolution/cartesianR6_CzarnyGeometry.h>
#include <InputFunctions/ExactSolution/cartesianR6_ShafranovGeometry.h>
#include <InputFunctions/ExactSolution/polarR6_CircularGeometry.h>
#include <InputFunctions/ExactSolution/polarR6_CulhamGeometry.h>
#include <InputFunctions/ExactSolution/polarR6_CzarnyGeometry.h>
#include <InputFunctions/ExactSolution/polarR6_ShafranovGeometry.h>
#include <InputFunctions/ExactSolution/refined_CircularGeometry.h>
#include <InputFunctions/ExactSolution/refined_CulhamGeometry.h>
#include <InputFunctions/ExactSolution/refined_CzarnyGeometry.h>
#include <InputFunctions/ExactSolution/refined_ShafranovGeometry.h>

/* ------------------- */
/* Boundary Conditions */
#include <InputFunctions/BoundaryConditions/cartesianR2_Boundary_CircularGeometry.h>
#include <InputFunctions/BoundaryConditions/cartesianR2_Boundary_CzarnyGeometry.h>
#include <InputFunctions/BoundaryConditions/cartesianR2_Boundary_ShafranovGeometry.h>
#include <InputFunctions/BoundaryConditions/cartesianR6_Boundary_CircularGeometry.h>
#include <InputFunctions/BoundaryConditions/cartesianR6_Boundary_CzarnyGeometry.h>
#include <InputFunctions/BoundaryConditions/cartesianR6_Boundary_ShafranovGeometry.h>
#include <InputFunctions/BoundaryConditions/polarR6_Boundary_CircularGeometry.h>
#include <InputFunctions/BoundaryConditions/polarR6_Boundary_CulhamGeometry.h>
#include <InputFunctions/BoundaryConditions/polarR6_Boundary_CzarnyGeometry.h>
#include <InputFunctions/BoundaryConditions/polarR6_Boundary_ShafranovGeometry.h>
#include <InputFunctions/BoundaryConditions/refined_Boundary_CircularGeometry.h>
#include <InputFunctions/BoundaryConditions/refined_Boundary_CulhamGeometry.h>
#include <InputFunctions/BoundaryConditions/refined_Boundary_CzarnyGeometry.h>
#include <InputFunctions/BoundaryConditions/refined_Boundary_ShafranovGeometry.h>

/* -----------------------------*/
/* Density Profile Coefficients */
#include <InputFunctions/DensityProfileCoefficients/poissonCoefficients.h>
#include <InputFunctions/DensityProfileCoefficients/sonnendruckerCoefficients.h>
#include <InputFunctions/DensityProfileCoefficients/sonnendruckerGyroCoefficients.h>
#include <InputFunctions/DensityProfileCoefficients/zoniCoefficients.h>
#include <InputFunctions/DensityProfileCoefficients/zoniGyroCoefficients.h>
#include <InputFunctions/DensityProfileCoefficients/zoniShiftedCoefficients.h>
#include <InputFunctions/DensityProfileCoefficients/zoniShiftedGyroCoefficients.h>

/* ------------------------- */
/* Source Terms: CartesianR2 */
#include <InputFunctions/SourceTerms/cartesianR2_Poisson_CircularGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_Poisson_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_Poisson_ShafranovGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_SonnendruckerGyro_CircularGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_SonnendruckerGyro_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_SonnendruckerGyro_ShafranovGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_Sonnendrucker_CircularGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_Sonnendrucker_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_Sonnendrucker_ShafranovGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_ZoniGyro_CircularGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_ZoniGyro_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_ZoniGyro_ShafranovGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_ZoniShiftedGyro_CircularGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_ZoniShiftedGyro_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_ZoniShiftedGyro_ShafranovGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_ZoniShifted_CircularGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_ZoniShifted_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_ZoniShifted_ShafranovGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_Zoni_CircularGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_Zoni_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR2_Zoni_ShafranovGeometry.h>
/* ------------------------- */
/* Source Terms: CartesianR6 */
#include <InputFunctions/SourceTerms/cartesianR6_Poisson_CircularGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_Poisson_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_Poisson_ShafranovGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_SonnendruckerGyro_CircularGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_SonnendruckerGyro_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_SonnendruckerGyro_ShafranovGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_Sonnendrucker_CircularGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_Sonnendrucker_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_Sonnendrucker_ShafranovGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_ZoniGyro_CircularGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_ZoniGyro_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_ZoniGyro_ShafranovGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_ZoniShiftedGyro_CircularGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_ZoniShiftedGyro_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_ZoniShiftedGyro_ShafranovGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_ZoniShifted_CircularGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_ZoniShifted_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_ZoniShifted_ShafranovGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_Zoni_CircularGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_Zoni_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/cartesianR6_Zoni_ShafranovGeometry.h>
/* --------------------- */
/* Source Terms: PolarR6 */
#include <InputFunctions/SourceTerms/polarR6_Poisson_CircularGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_Poisson_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_Poisson_ShafranovGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_SonnendruckerGyro_CircularGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_SonnendruckerGyro_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_SonnendruckerGyro_ShafranovGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_Sonnendrucker_CircularGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_Sonnendrucker_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_Sonnendrucker_ShafranovGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_ZoniGyro_CircularGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_ZoniGyro_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_ZoniGyro_ShafranovGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_ZoniShiftedGyro_CircularGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_ZoniShiftedGyro_CulhamGeometry.h> /* Culham */
#include <InputFunctions/SourceTerms/polarR6_ZoniShiftedGyro_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_ZoniShiftedGyro_ShafranovGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_ZoniShifted_CircularGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_ZoniShifted_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_ZoniShifted_ShafranovGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_Zoni_CircularGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_Zoni_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/polarR6_Zoni_ShafranovGeometry.h>

/* --------------------- */
/* Source Terms: Refined */
#include <InputFunctions/SourceTerms/refined_ZoniShiftedGyro_CircularGeometry.h>
#include <InputFunctions/SourceTerms/refined_ZoniShiftedGyro_CulhamGeometry.h> /* Culham */
#include <InputFunctions/SourceTerms/refined_ZoniShiftedGyro_CzarnyGeometry.h>
#include <InputFunctions/SourceTerms/refined_ZoniShiftedGyro_ShafranovGeometry.h>