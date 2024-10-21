cmake_minimum_required(VERSION 3.12)
project(GMGPolar VERSION 2.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic -Wno-unused")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -mtune=generic")

# External library paths
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

# Include directories for the project
include_directories(include)

# 1. Create a library target for the PolarGrid module
file(GLOB_RECURSE POLAR_GRID_SOURCES "src/PolarGrid/*.cpp")
add_library(PolarGrid STATIC ${POLAR_GRID_SOURCES})
target_include_directories(PolarGrid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/PolarGrid)

# 2. Create a library target for the InputFunctions module
file(GLOB_RECURSE INPUT_FUNCTIONS_SOURCES
    "src/InputFunctions/DensityProfileCoefficients/*.cpp"
    "src/InputFunctions/DomainGeometry/*.cpp"
    "src/InputFunctions/ExactSolution/*.cpp"
    "src/InputFunctions/BoundaryConditions/*.cpp"
    "src/InputFunctions/SourceTerms/*.cpp"
)
add_library(InputFunctions STATIC ${INPUT_FUNCTIONS_SOURCES})
target_include_directories(InputFunctions PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/InputFunctions)

# 3. Add source files for the rest of the project
file(GLOB_RECURSE GMG_POLAR_SOURCES "src/GMGPolar/*.cpp")
file(GLOB_RECURSE MULTIGRID_METHODS_SOURCES "src/GMGPolar/MultigridMethods/*.cpp")
file(GLOB_RECURSE LEVEL_SOURCES "src/Level/*.cpp")
file(GLOB_RECURSE STENCIL_SOURCES "src/Stencil/*.cpp")
file(GLOB_RECURSE INTERPOLATION_SOURCES "src/Interpolation/*.cpp")
file(GLOB_RECURSE DIRECT_SOLVER_SOURCES "src/DirectSolver/*.cpp" "src/DirectSolverGive/*.cpp" "src/DirectSolverTake/*.cpp")
file(GLOB_RECURSE RESIDUAL_SOURCES "src/Residual/*.cpp" "src/Residual/ResidualGive/*.cpp" "src/Residual/ResidualTake/*.cpp")
file(GLOB_RECURSE SMOOTHER_SOURCES "src/Smoother/*.cpp" "src/SmootherGive/*.cpp" "src/SmootherTake/*.cpp")
file(GLOB_RECURSE EXTRAPOLATED_SMOOTHER_SOURCES "src/ExtrapolatedSmoother/*.cpp" "src/ExtrapolatedSmoother/ExtrapolatedSmootherGive/*.cpp" "src/ExtrapolatedSmoother/ExtrapolatedSmootherTake/*.cpp")

# Add main.cpp to the list of source files
set(MAIN_SOURCE "src/main.cpp")

# 4. Add the main executable target
add_executable(gmgpolar 
    ${MAIN_SOURCE}
    ${GMG_POLAR_SOURCES}
    ${MULTIGRID_METHODS_SOURCES}
    ${LEVEL_SOURCES}
    ${STENCIL_SOURCES}
    ${INTERPOLATION_SOURCES}
    ${DIRECT_SOLVER_SOURCES}
    ${RESIDUAL_SOURCES}
    ${SMOOTHER_SOURCES}
    ${EXTRAPOLATED_SMOOTHER_SOURCES}
)

# 5. Link the PolarGrid and InputFunctions libraries to the main executable
target_link_libraries(gmgpolar PRIVATE PolarGrid InputFunctions ${MUMPS_LIBRARIES} OpenMP::OpenMP_CXX)

# 6. Add OpenMP flags if available
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# 7. Enable testing and include other directories
enable_testing()
add_subdirectory(third-party)
add_subdirectory(tests)
