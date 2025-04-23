set(BOUNDARY_CONDITIONS_SOURCES
    cartesianR2_Boundary_CircularGeometry.cpp
    cartesianR2_Boundary_CzarnyGeometry.cpp
    cartesianR2_Boundary_ShafranovGeometry.cpp
    cartesianR6_Boundary_CircularGeometry.cpp
    cartesianR6_Boundary_CzarnyGeometry.cpp
    cartesianR6_Boundary_ShafranovGeometry.cpp
    polarR6_Boundary_CircularGeometry.cpp
    polarR6_Boundary_CulhamGeometry.cpp
    polarR6_Boundary_CzarnyGeometry.cpp
    polarR6_Boundary_ShafranovGeometry.cpp
    refined_Boundary_CircularGeometry.cpp
    refined_Boundary_CulhamGeometry.cpp
    refined_Boundary_CzarnyGeometry.cpp
    refined_Boundary_ShafranovGeometry.cpp
)

add_library(InputFunctions_BoundaryConditions STATIC ${BOUNDARY_CONDITIONS_SOURCES})

target_include_directories(InputFunctions_BoundaryConditions PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/include/InputFunctions/BoundaryConditions
)
