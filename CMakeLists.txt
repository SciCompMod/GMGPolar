cmake_minimum_required(VERSION 3.12)
project(GMGPolar VERSION 2.0.0 LANGUAGES CXX)

option(GMGPOLAR_BUILD_TESTS "Build GMGPolar unit tests." ON)
option(GMGPOLAR_USE_LIKWID "Use LIKWID to measure code (regions)." ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic -Wno-unused")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -mtune=generic")

# Mumps: Sparse Matrix Solver
set(MUMPS_PREFIX_PATH "~/spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/mumps-5.6.2-m4xrhc4mshzrxmgptzbpult3nbf6qrzk")
include_directories(${MUMPS_PREFIX_PATH}/include)
link_directories(${MUMPS_PREFIX_PATH}/lib)

# Metis: Matrix reordering for Mumps
set(METIS_PREFIX_PATH "~/spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/metis-5.1.0-bgoncx22w55soviybggl5ydjakvkm34v")
include_directories(${METIS_PREFIX_PATH}/include)
link_directories(${METIS_PREFIX_PATH}/lib)

# Likwid: Performance monitoring
set(LIKWID_PREFIX_PATH "~/spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/likwid-5.3.0-6mvx2snsqnamuyhaqspd6gxkfuaso36g")
include_directories(${LIKWID_PREFIX_PATH}/include)
link_directories(${LIKWID_PREFIX_PATH}/lib)

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

# 3. Collect all the common source files for GMGPolar into a library
file(GLOB_RECURSE GMG_POLAR_SOURCES "src/GMGPolar/*.cpp")
file(GLOB_RECURSE MULTIGRID_METHODS_SOURCES "src/GMGPolar/MultigridMethods/*.cpp")
file(GLOB_RECURSE LEVEL_SOURCES "src/Level/*.cpp")
file(GLOB_RECURSE STENCIL_SOURCES "src/Stencil/*.cpp")
file(GLOB_RECURSE INTERPOLATION_SOURCES "src/Interpolation/*.cpp")
file(GLOB_RECURSE DIRECT_SOLVER_SOURCES "src/DirectSolver/*.cpp" "src/DirectSolverGive/*.cpp" "src/DirectSolverTake/*.cpp" "src/DirectSolverTakeNoMumps/*.cpp")
file(GLOB_RECURSE RESIDUAL_SOURCES "src/Residual/*.cpp" "src/Residual/ResidualGive/*.cpp" "src/Residual/ResidualTake/*.cpp")
file(GLOB_RECURSE SMOOTHER_SOURCES "src/Smoother/*.cpp" "src/SmootherGive/*.cpp" "src/SmootherTake/*.cpp")
file(GLOB_RECURSE EXTRAPOLATED_SMOOTHER_SOURCES "src/ExtrapolatedSmoother/*.cpp" "src/ExtrapolatedSmoother/ExtrapolatedSmootherGive/*.cpp" "src/ExtrapolatedSmoother/ExtrapolatedSmootherTake/*.cpp")

# 4. Create the GMGPolarLib library and link PolarGrid and InputFunctions
add_library(GMGPolarLib STATIC
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
# Link PolarGrid and InputFunctions to GMGPolarLib
target_link_libraries(GMGPolarLib PUBLIC PolarGrid InputFunctions)

# Link Likwid to GMGPolarLib
if(GMGPOLAR_USE_LIKWID)
    target_link_libraries(GMGPolarLib PUBLIC likwid)
    target_compile_definitions(GMGPolarLib PUBLIC "-DLIKWID_PERFMON")
    add_compile_definitions(GMGPOLAR_USE_LIKWID)
endif()

# Mumps options
set(MUMPS_LIBRARIES
    mumps_common
    smumps
    dmumps
    mpiseq
    metis
)

# 5. Add the main executable target (gmgpolar)
set(MAIN_SOURCE "src/main.cpp")
add_executable(gmgpolar ${MAIN_SOURCE})

# 6. Add another executable target (convergence_order)
set(CONVERGENCE_ORDER_SOURCE "src/convergence_order.cpp")
add_executable(convergence_order ${CONVERGENCE_ORDER_SOURCE})

set(WEAK_SCALING_SOURCE "src/weak_scaling.cpp")
add_executable(weak_scaling ${WEAK_SCALING_SOURCE})

set(STRONG_SCALING_SOURCE "src/strong_scaling.cpp")
add_executable(strong_scaling ${STRONG_SCALING_SOURCE})

# 7. Link GMGPolarLib (which now includes PolarGrid and InputFunctions) to both executables
target_link_libraries(gmgpolar PRIVATE GMGPolarLib ${MUMPS_LIBRARIES} OpenMP::OpenMP_CXX)
target_link_libraries(convergence_order PRIVATE GMGPolarLib ${MUMPS_LIBRARIES} OpenMP::OpenMP_CXX)
target_link_libraries(weak_scaling PRIVATE GMGPolarLib ${MUMPS_LIBRARIES} OpenMP::OpenMP_CXX)
target_link_libraries(strong_scaling PRIVATE GMGPolarLib ${MUMPS_LIBRARIES} OpenMP::OpenMP_CXX)

# 8. Add OpenMP flags if available
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# 9. Enable testing and include other directories
if(GMGPOLAR_BUILD_TESTS)
    enable_testing()
    add_subdirectory(third-party)
    add_subdirectory(tests)
endif()