
# Add subdirectories for components
add_subdirectory(InputFunctions)

# Gather all source files
# file(GLOB_RECURSE POLAR_GRID_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/PolarGrid/*.cpp)
set(POLAR_GRID_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/PolarGrid/polargrid.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PolarGrid/anisotropic_division.cpp
)

# Gather all source files
# file(GLOB_RECURSE LINEAR_ALGEBRA_SOURCES 
# ${CMAKE_CURRENT_SOURCE_DIR}/LinearAlgebra/Solvers/*.cpp)
set(LINEAR_ALGEBRA_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/LinearAlgebra/Solvers/coo_mumps_solver.cpp
)

# file(GLOB_RECURSE GMG_POLAR_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/*.cpp)
set(GMG_POLAR_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/gmgpolar.cpp
)

# file(GLOB_RECURSE STENCIL_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Stencil/*.cpp)
set(STENCIL_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/Stencil/stencil.cpp
)

# file(GLOB_RECURSE INTERPOLATION_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Interpolation/*.cpp)
set(INTERPOLATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/Interpolation/extrapolated_prolongation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Interpolation/extrapolated_restriction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Interpolation/fmg_interpolation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Interpolation/injection.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Interpolation/interpolation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Interpolation/prolongation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Interpolation/restriction.cpp
)

# Gather all source files
# file(GLOB_RECURSE CONFIG_PARSER_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/ConfigParser/*.cpp)
set(CONFIG_PARSER_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/ConfigParser/config_parser.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ConfigParser/select_test_case.cpp
)

# Create the main library
add_library(GMGPolarLib STATIC
    ${POLAR_GRID_SOURCES}
    ${LINEAR_ALGEBRA_SOURCES}
    ${GMG_POLAR_SOURCES}
    ${STENCIL_SOURCES}
    ${INTERPOLATION_SOURCES}
    ${CONFIG_PARSER_SOURCES}
)

# Basic library configuration
target_include_directories(GMGPolarLib PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

# Handle OpenMP
target_link_libraries(GMGPolarLib
    PUBLIC
    OpenMP::OpenMP_CXX
    Kokkos::kokkos
)


install(TARGETS GMGPolarLib
    ARCHIVE DESTINATION lib
)

# Handle MUMPS configuration
if(GMGPOLAR_USE_MUMPS)
    target_link_libraries(GMGPolarLib PUBLIC MUMPS::MUMPS)
    target_compile_definitions(GMGPolarLib PUBLIC GMGPOLAR_USE_MUMPS)
endif()

# Handle LIKWID configuration
if(GMGPOLAR_USE_LIKWID)
    target_link_libraries(GMGPolarLib PUBLIC LIKWID::LIKWID)
    target_compile_definitions(GMGPolarLib PUBLIC GMGPOLAR_USE_LIKWID)
endif()
