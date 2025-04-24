set(EXACT_SOLUTION_SOURCES
    cartesianR2_CircularGeometry.cpp
    cartesianR2_CzarnyGeometry.cpp
    cartesianR2_ShafranovGeometry.cpp
    cartesianR6_CircularGeometry.cpp
    cartesianR6_CzarnyGeometry.cpp
    cartesianR6_ShafranovGeometry.cpp
    polarR6_CircularGeometry.cpp
    polarR6_CulhamGeometry.cpp
    polarR6_CzarnyGeometry.cpp
    polarR6_ShafranovGeometry.cpp
    refined_CircularGeometry.cpp
    refined_CulhamGeometry.cpp
    refined_CzarnyGeometry.cpp
    refined_ShafranovGeometry.cpp
)

add_library(InputFunctions_ExactSolution STATIC ${EXACT_SOLUTION_SOURCES})
target_include_directories(InputFunctions_ExactSolution PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/include/InputFunctions/ExactSolution
)