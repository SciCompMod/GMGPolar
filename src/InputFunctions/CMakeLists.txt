add_subdirectory(BoundaryConditions)
add_subdirectory(DensityProfileCoefficients)
add_subdirectory(DomainGeometry)
add_subdirectory(ExactSolution)
add_subdirectory(SourceTerms)

add_library(InputFunctions INTERFACE)

target_link_libraries(InputFunctions INTERFACE
    InputFunctions_BoundaryConditions
    InputFunctions_DensityProfileCoefficients
    InputFunctions_DomainGeometry
    InputFunctions_ExactSolution
    InputFunctions_SourceTerms
)

install(TARGETS
    InputFunctions
    InputFunctions_BoundaryConditions
    InputFunctions_DensityProfileCoefficients
    InputFunctions_DomainGeometry
    InputFunctions_ExactSolution
    InputFunctions_SourceTerms
    EXPORT GMGPolarTargets
    ARCHIVE DESTINATION lib
)
