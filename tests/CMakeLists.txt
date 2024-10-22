cmake_minimum_required(VERSION 3.12)

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
file(GLOB_RECURSE POLAR_GRID_SOURCES "../src/PolarGrid/*.cpp")

file(GLOB_RECURSE INPUT_FUNCTIONS_SOURCES
    "../src/InputFunctions/DensityProfileCoefficients/*.cpp"
    "../src/InputFunctions/DomainGeometry/*.cpp"
    "../src/InputFunctions/ExactSolution/*.cpp"
    "../src/InputFunctions/BoundaryConditions/*.cpp"
    "../src/InputFunctions/SourceTerms/*.cpp"
)

file(GLOB_RECURSE GMG_POLAR_SOURCES "../src/GMGPolar/*.cpp")
file(GLOB_RECURSE MULTIGRID_METHODS_SOURCES "../src/GMGPolar/MultigridMethods/*.cpp")
file(GLOB_RECURSE LEVEL_SOURCES "../src/Level/*.cpp")
file(GLOB_RECURSE STENCIL_SOURCES "../src/Stencil/*.cpp")
file(GLOB_RECURSE INTERPOLATION_SOURCES "../src/Interpolation/*.cpp")

file(GLOB_RECURSE DIRECT_SOLVER_SOURCES 
    "../src/DirectSolver/*.cpp" 
    "../src/DirectSolverGive/*.cpp" 
    "../src/DirectSolverTake/*.cpp"
)
file(GLOB_RECURSE RESIDUAL_SOURCES 
    "../src/Residual/*.cpp" 
    "../src/Residual/ResidualGive/*.cpp" 
    "../src/Residual/ResidualTake/*.cpp"
)
file(GLOB_RECURSE SMOOTHER_SOURCES
    "../src/Smoother/*.cpp"
    "../src/SmootherGive/*.cpp"
    "../src/SmootherTake/*.cpp"
)
file(GLOB_RECURSE EXTRAPOLATED_SMOOTHER_SOURCES 
    "../src/ExtrapolatedSmoother/*.cpp" 
    "../src/ExtrapolatedSmoother/ExtrapolatedSmootherGive/*.cpp"
    "../src/ExtrapolatedSmoother/ExtrapolatedSmootherTake/*.cpp"
)

add_library(GMGPolarLib STATIC
    ${GMG_POLAR_SOURCES}
    ${MULTIGRID_METHODS_SOURCES}
    ${POLAR_GRID_SOURCES}
    ${LEVEL_SOURCES}
    ${STENCIL_SOURCES}
    ${INTERPOLATION_SOURCES}
    ${DIRECT_SOLVER_SOURCES}
    ${RESIDUAL_SOURCES}
    ${SMOOTHER_SOURCES}
    ${EXTRAPOLATED_SMOOTHER_SOURCES}
    ${INPUT_FUNCTIONS_SOURCES}
)

add_executable(test_vector LinearAlgebra/vector.cpp)
target_compile_features(test_vector PRIVATE cxx_std_17)
target_link_libraries(test_vector GMGPolarLib GTest::gtest_main ${MUMPS_LIBRARIES})
add_test(NAME test_vector COMMAND test_vector)

add_executable(test_vector_operations LinearAlgebra/vector_operations.cpp)
target_compile_features(test_vector_operations PRIVATE cxx_std_17)
target_link_libraries(test_vector_operations GMGPolarLib GTest::gtest_main ${MUMPS_LIBRARIES})
add_test(NAME test_vector_operations COMMAND test_vector_operations)

add_executable(test_matrix LinearAlgebra/matrix.cpp)
target_compile_features(test_matrix PRIVATE cxx_std_17)
target_link_libraries(test_matrix GMGPolarLib GTest::gtest_main ${MUMPS_LIBRARIES})
add_test(NAME test_matrix COMMAND test_matrix)

add_executable(test_tridiagonal_solver LinearAlgebra/tridiagonal_solver.cpp)
target_compile_features(test_tridiagonal_solver PRIVATE cxx_std_17)
target_link_libraries(test_tridiagonal_solver GMGPolarLib GTest::gtest_main ${MUMPS_LIBRARIES})
add_test(NAME test_tridiagonal_solver COMMAND test_tridiagonal_solver)

add_executable(test_cyclic_tridiagonal_solver LinearAlgebra/cyclic_tridiagonal_solver.cpp)
target_compile_features(test_cyclic_tridiagonal_solver PRIVATE cxx_std_17)
target_link_libraries(test_cyclic_tridiagonal_solver GMGPolarLib GTest::gtest_main ${MUMPS_LIBRARIES})
add_test(NAME test_cyclic_tridiagonal_solver COMMAND test_cyclic_tridiagonal_solver)

add_executable(test_polargrid PolarGrid/polargrid.cpp)
target_compile_features(test_polargrid PRIVATE cxx_std_17)
target_link_libraries(test_polargrid GMGPolarLib GTest::gtest_main ${MUMPS_LIBRARIES})
add_test(NAME test_polargrid COMMAND test_polargrid)

add_executable(test_prolongation Interpolation/prolongation.cpp)
target_compile_features(test_prolongation PRIVATE cxx_std_17)
target_link_libraries(test_prolongation GMGPolarLib GTest::gtest_main ${MUMPS_LIBRARIES})
add_test(NAME test_prolongation COMMAND test_prolongation)

add_executable(test_restriction Interpolation/restriction.cpp)
target_compile_features(test_restriction PRIVATE cxx_std_17)
target_link_libraries(test_restriction GMGPolarLib GTest::gtest_main ${MUMPS_LIBRARIES})
add_test(NAME test_restriction COMMAND test_restriction)

add_executable(test_extrapolated_prolongation Interpolation/extrapolated_prolongation.cpp)
target_compile_features(test_extrapolated_prolongation PRIVATE cxx_std_17)
target_link_libraries(test_extrapolated_prolongation GMGPolarLib GTest::gtest_main ${MUMPS_LIBRARIES})
add_test(NAME test_extrapolated_prolongation COMMAND test_extrapolated_prolongation)

add_executable(test_extrapolated_restriction Interpolation/extrapolated_restriction.cpp)
target_compile_features(test_extrapolated_restriction PRIVATE cxx_std_17)
target_link_libraries(test_extrapolated_restriction GMGPolarLib GTest::gtest_main ${MUMPS_LIBRARIES})
add_test(NAME test_extrapolated_restriction COMMAND test_extrapolated_restriction)

add_executable(test_residual Residual/residual.cpp)
target_compile_features(test_residual PRIVATE cxx_std_17)
target_link_libraries(test_residual GMGPolarLib GTest::gtest_main ${MUMPS_LIBRARIES})
add_test(NAME test_residual COMMAND test_residual)

add_executable(test_directSolver DirectSolver/directSolver.cpp)
target_compile_features(test_directSolver PRIVATE cxx_std_17)
target_link_libraries(test_directSolver GMGPolarLib GTest::gtest_main ${MUMPS_LIBRARIES})
add_test(NAME test_directSolver COMMAND test_directSolver)

add_executable(test_smoother Smoother/smoother.cpp)
target_compile_features(test_smoother PRIVATE cxx_std_17)
target_link_libraries(test_smoother GMGPolarLib GTest::gtest_main ${MUMPS_LIBRARIES})
add_test(NAME test_smoother COMMAND test_smoother)

add_executable(test_extrapolated_smoother ExtrapolatedSmoother/extrapolated_smoother.cpp)
target_compile_features(test_extrapolated_smoother PRIVATE cxx_std_17)
target_link_libraries(test_extrapolated_smoother GMGPolarLib GTest::gtest_main ${MUMPS_LIBRARIES})
add_test(NAME test_extrapolated_smoother COMMAND test_extrapolated_smoother)

include(GoogleTest)
