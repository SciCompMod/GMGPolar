
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

# Build the test suite twice: once for serial (1 thread) and once for parallel (16 threads).
# GMGPOLAR_TEST_THREADS is passed to Kokkos::InitializationSettings at startup.
foreach(label_threads "Serial;1" "Parallel;16")
    list(GET label_threads 0 label)
    list(GET label_threads 1 num_threads)
    set(test_target gmgpolar_tests_${label})

    add_executable(${test_target} ${gmgpolar_test_sources})
    target_compile_features(${test_target} PRIVATE cxx_std_20)
    target_compile_definitions(${test_target} PRIVATE GMGPOLAR_TEST_THREADS=${num_threads})
    target_link_libraries(gmgpolar_tests GMGPolarLib GMGPolarInterface GTest::gtest_main)

    gtest_discover_tests(${test_target}
        TEST_PREFIX ${label}
        DISCOVERY_MODE PRE_TEST)
endforeach()
