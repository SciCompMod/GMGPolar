cmake_minimum_required(VERSION 3.12)

# Set the project name
project(GMGPolarTests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

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

include_directories(../include)

# Add source files
file(GLOB_RECURSE GMG_POLAR_SOURCES "../src/GMGPolar/*.cpp")
file(GLOB_RECURSE MULTIGRID_METHODS_SOURCES "../src/GMGPolar/MultigridMethods/*.cpp")
file(GLOB_RECURSE POLAR_GRID_SOURCES "../src/PolarGrid/*.cpp")
file(GLOB_RECURSE LEVEL_SOURCES "../src/Level/*.cpp")

file(GLOB_RECURSE TASK_DISTRIBUTION_SOURCES "../src/TaskDistribution/*.cpp")
file(GLOB_RECURSE STENCIL_SOURCES "../src/Stencil/*.cpp")
file(GLOB_RECURSE INTERPOLATION_SOURCES "../src/Interpolation/*.cpp")
file(GLOB_RECURSE DIRECT_SOLVER_SOURCES "../src/DirectSolver/*.cpp")
file(GLOB_RECURSE RESIDUAL_SOURCES "../src/Residual/*.cpp")
file(GLOB_RECURSE SMOOTHER_SOURCES "../src/Smoother/*.cpp")
file(GLOB_RECURSE EXTRAPOLATED_SMOOTHER_SOURCES "../src/ExtrapolatedSmoother/*.cpp")

file(GLOB_RECURSE DENSITY_PROFILE_COEFFICIENTS "../src/InputFunctions/DensityProfileCoefficients/*.cpp")
file(GLOB_RECURSE DOMAIN_GEOMETRY_SOURCES "../src/InputFunctions/DomainGeometry/*.cpp")
file(GLOB_RECURSE EXACT_SOLUTION_SOURCES "../src/InputFunctions/ExactSolution/*.cpp")
file(GLOB_RECURSE BOUNDARY_CONDITION_SOURCES "../src/InputFunctions/BoundaryConditions/*.cpp")
file(GLOB_RECURSE SOURCE_TERMS_SOURCES "../src/InputFunctions/SourceTerms/*.cpp")

# Add test files
file(GLOB_RECURSE TEST_SOURCES "../tests/*.cpp")

# Set the list of test sources
set(TESTSOURCES
    LinearAlgebra/vector.cpp
    LinearAlgebra/vector_operations.cpp
    LinearAlgebra/matrix.cpp
    LinearAlgebra/tridiagonal_solver.cpp
    LinearAlgebra/cyclic_tridiagonal_solver.cpp
    PolarGrid/polargrid.cpp
    Interpolation/prolongation.cpp
    Interpolation/restriction.cpp
    Interpolation/extrapolated_prolongation.cpp
    Interpolation/extrapolated_restriction.cpp
    Residual/residual.cpp
    DirectSolver/directSolver.cpp
    Smoother/smoother.cpp
    ExtrapolatedSmoother/extrapolated_smoother.cpp
)

# Create static libraries from source files
add_library(GMGPolarLib STATIC
    ${GMG_POLAR_SOURCES}
    ${MULTIGRID_METHODS_SOURCES}
    ${POLAR_GRID_SOURCES}
    ${LEVEL_SOURCES}
    ${TASK_DISTRIBUTION_SOURCES}
    ${STENCIL_SOURCES}
    ${INTERPOLATION_SOURCES}
    ${DIRECT_SOLVER_SOURCES}
    ${RESIDUAL_SOURCES}
    ${SMOOTHER_SOURCES}
    ${EXTRAPOLATED_SMOOTHER_SOURCES}
    ${DENSITY_PROFILE_COEFFICIENTS}
    ${DOMAIN_GEOMETRY_SOURCES}
    ${EXACT_SOLUTION_SOURCES}
    ${BOUNDARY_CONDITION_SOURCES}
    ${SOURCE_TERMS_SOURCES}
)

# Add OpenMP flags if available
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(GMGPolarLib PUBLIC OpenMP::OpenMP_CXX)
endif()

# Add executables for each test source
foreach(TESTSOURCE ${TEST_SOURCES})
    get_filename_component(TESTNAME ${TESTSOURCE} NAME_WE)
    add_executable(${TESTNAME} ${TESTSOURCE})
    target_link_libraries(${TESTNAME} GMGPolarLib GTest::gtest_main ${MUMPS_LIBRARIES})
    target_compile_features(${TESTNAME} PRIVATE cxx_std_20)
    list(APPEND TESTS ${TESTNAME})
endforeach()

# Enable testing
include(GoogleTest)
gtest_discover_tests(${TESTS})