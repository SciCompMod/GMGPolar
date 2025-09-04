cmake_minimum_required(VERSION 3.12)

project(GMGPolarTests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Create a single test executable that includes all test sources
add_executable(gmgpolar_tests
    gmgpolar_tests.cpp
    LinearAlgebra/vector.cpp
    LinearAlgebra/vector_operations.cpp
    LinearAlgebra/matrix.cpp
    LinearAlgebra/csr_matrix.cpp
    LinearAlgebra/csr_solver.cpp
    LinearAlgebra/tridiagonal_solver.cpp
    LinearAlgebra/cyclic_tridiagonal_solver.cpp
    PolarGrid/polargrid.cpp
    Interpolation/prolongation.cpp
    Interpolation/restriction.cpp
    Interpolation/extrapolated_prolongation.cpp
    Interpolation/extrapolated_restriction.cpp
    Residual/residual.cpp
    DirectSolver/directSolver.cpp
    DirectSolver/directSolverNoMumps.cpp
    Smoother/smoother.cpp
    ExtrapolatedSmoother/extrapolated_smoother.cpp
    ConfigParser/config_parser.cpp
    GMGPolar/polarR6.cpp
    GMGPolar/cartesianR2.cpp
    GMGPolar/non_uniform.cpp
    GMGPolar/stencil_distribution.cpp
    GMGPolar/extrapolation.cpp
    GMGPolar/cycle_types.cpp
)

# Set the compile features and link libraries
target_compile_features(gmgpolar_tests PRIVATE cxx_std_17)
target_link_libraries(gmgpolar_tests GMGPolarLib GTest::gtest_main)

# Register the test
add_test(NAME gmgpolar_tests COMMAND gmgpolar_tests)

include(GoogleTest)
gtest_discover_tests(gmgpolar_tests)
