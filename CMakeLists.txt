cmake_minimum_required(VERSION 3.12)
project(GMGPolar VERSION 2.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

set(MUMPS_PREFIX_PATH "~/spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/mumps-5.6.2-m4xrhc4mshzrxmgptzbpult3nbf6qrzk")
include_directories(${MUMPS_PREFIX_PATH}/include)
link_directories(${MUMPS_PREFIX_PATH}/lib)

set(METIS_PREFIX_PATH "~/spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/metis-5.1.0-bgoncx22w55soviybggl5ydjakvkm34v")
include_directories(${METIS_PREFIX_PATH}/include)
link_directories(${METIS_PREFIX_PATH}/lib)

set(LIKWID_PREFIX_PATH "~/spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/likwid-5.3.0-6mvx2snsqnamuyhaqspd6gxkfuaso36g")
include_directories(${LIKWID_PREFIX_PATH}/include)
link_directories(${LIKWID_PREFIX_PATH}/lib)

set(MUMPS_LIBRARIES
    mumps_common
    smumps
    dmumps
    mpiseq
    metis
)

include_directories(include)

# Add source files
file(GLOB_RECURSE GMG_POLAR_SOURCES "src/GMGPolar/*.cpp")
file(GLOB_RECURSE POLAR_GRID_SOURCES "src/PolarGrid/*.cpp")
file(GLOB_RECURSE LEVEL_SOURCES "src/Level/*.cpp")

file(GLOB_RECURSE TASK_DISTRIBUTION_SOURCES "src/TaskDistribution/*.cpp")
file(GLOB_RECURSE STENCIL_SOURCES "src/Stencil/*.cpp")
file(GLOB_RECURSE INTERPOLATION_SOURCES "src/Interpolation/*.cpp")
file(GLOB_RECURSE DIRECT_SOLVER_SOURCES "src/DirectSolver/*.cpp")
file(GLOB_RECURSE RESIDUAL_SOURCES "src/Residual/*.cpp")
file(GLOB_RECURSE SMOOTHER_SOURCES "src/Smoother/*.cpp")

file(GLOB_RECURSE DENSITY_PROFILE_COEFFICIENTS "src/InputFunctions/DensityProfileCoefficients/*.cpp")
file(GLOB_RECURSE DOMAIN_GEOMETRY_SOURCES "src/InputFunctions/DomainGeometry/*.cpp")
file(GLOB_RECURSE EXACT_SOLUTION_SOURCES "src/InputFunctions/ExactSolution/*.cpp")
file(GLOB_RECURSE BOUNDARY_CONDITION_SOURCES "src/InputFunctions/BoundaryConditions/*.cpp")
file(GLOB_RECURSE SOURCE_TERMS_SOURCES "src/InputFunctions/SourceTerms/*.cpp")

# file(GLOB_RECURSE GMG_POLAR_SOURCES "src/GMGPolar/*.cpp")
# file(GLOB_RECURSE INTERPOLATION_SOURCES "src/Interpolation/*.cpp")
# file(GLOB_RECURSE OPERATOR_SOURCES "src/Operator/*.cpp")
# file(GLOB_RECURSE POLAR_GRID_SOURCES "src/PolarGrid/*.cpp")
# file(GLOB_RECURSE TASK_DISTRIBUTION_SOURCES "src/TaskDistribution/*.cpp")
# file(GLOB_RECURSE TEST_SOURCES "src/test_cases/*.cpp")
# file(GLOB_RECURSE INPUT_FUNCTIONS_SOURCES "src/InputFunctions/*.cpp")

# Add main.cpp to the list of source files
set(MAIN_SOURCE "src/main.cpp")

# Add executable target
add_executable(gmgpolar 
    ${MAIN_SOURCE}
    ${GMG_POLAR_SOURCES}
    ${POLAR_GRID_SOURCES}
    ${LEVEL_SOURCES}
    ${TASK_DISTRIBUTION_SOURCES}
    ${STENCIL_SOURCES}
    ${INTERPOLATION_SOURCES}
    ${DIRECT_SOLVER_SOURCES}
    ${RESIDUAL_SOURCES}
    ${SMOOTHER_SOURCES}

    ${DENSITY_PROFILE_COEFFICIENTS}
    ${DOMAIN_GEOMETRY_SOURCES}
    ${EXACT_SOLUTION_SOURCES}
    ${BOUNDARY_CONDITION_SOURCES}
    ${SOURCE_TERMS_SOURCES}

    # ${GMG_POLAR_SOURCES}
    # ${INTERPOLATION_SOURCES}
    # ${OPERATOR_SOURCES}
    # ${POLAR_GRID_SOURCES}
    # ${TASK_DISTRIBUTION_SOURCES}
    # ${TEST_SOURCES}
)

# Add OpenMP flags if available
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
target_link_libraries(gmgpolar PUBLIC ${MUMPS_LIBRARIES} OpenMP::OpenMP_CXX)

# enable_testing()
# add_subdirectory(third-party)
# add_subdirectory(tests)