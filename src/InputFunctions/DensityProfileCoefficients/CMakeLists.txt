set(DENSITY_PROFILE_COEFFICIENTS_SOURCES
    poissonCoefficients.cpp
    sonnendruckerCoefficients.cpp
    sonnendruckerGyroCoefficients.cpp
    zoniCoefficients.cpp
    zoniGyroCoefficients.cpp
    zoniShiftedCoefficients.cpp
    zoniShiftedGyroCoefficients.cpp
)

add_library(InputFunctions_DensityProfileCoefficients STATIC ${DENSITY_PROFILE_COEFFICIENTS_SOURCES})
target_include_directories(InputFunctions_DensityProfileCoefficients PUBLIC
    ${GMGPOLAR_INCLUDE_DIRS}
)
