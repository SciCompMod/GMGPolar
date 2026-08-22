
set(gmgpolar_test_sources
    gmgpolar_tests.cpp
    LinearAlgebra/Vector/vector_operations.cpp
    LinearAlgebra/Matrix/coo_matrix.cpp
    LinearAlgebra/Matrix/csr_matrix.cpp
    LinearAlgebra/Solvers/csr_lu_solver.cpp
    LinearAlgebra/Solvers/coo_mumps_solver.cpp
    LinearAlgebra/Solvers/tridiagonal_solver.cpp
    PolarGrid/polargrid.cpp
    Interpolation/prolongation.cpp
    Interpolation/restriction.cpp
    Interpolation/extrapolated_prolongation.cpp
    Interpolation/extrapolated_restriction.cpp
    Residual/residual.cpp
    DirectSolver/directSolver.cpp
    Smoother/smoother.cpp
    ExtrapolatedSmoother/extrapolated_smoother.cpp
    ConfigParser/config_parser.cpp
    GMGPolar/solve_tests.cpp
    GMGPolar/pcg_tests.cpp
    GMGPolar/convergence_order.cpp
)

include(GoogleTest)

add_executable(gmgpolar_tests ${gmgpolar_test_sources})
target_compile_features(gmgpolar_tests PRIVATE cxx_std_20)
target_link_libraries(gmgpolar_tests GMGPolarLib GMGPolarInterface GTest::gtest)

gtest_discover_tests(gmgpolar_tests DISCOVERY_MODE PRE_TEST)

add_executable(convergence_order ConvergenceOrder/convergence_order.cpp)
target_link_libraries(convergence_order PRIVATE GMGPolarLib GMGPolarInterface)
add_test(NAME convergence_order COMMAND convergence_order)
