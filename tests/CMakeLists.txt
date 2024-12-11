cmake_minimum_required(VERSION 3.12)

project(GMGPolarTests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

### ------------- ###
### LinearAlgebra ###
### ------------- ###

add_executable(test_vector LinearAlgebra/vector.cpp)
target_compile_features(test_vector PRIVATE cxx_std_20)
target_link_libraries(test_vector GTest::gtest_main GMGPolarLib)
add_test(NAME test_vector COMMAND test_vector)

add_executable(test_vector_operations LinearAlgebra/vector_operations.cpp)
target_compile_features(test_vector_operations PRIVATE cxx_std_20)
target_link_libraries(test_vector_operations GTest::gtest_main GMGPolarLib)
add_test(NAME test_vector_operations COMMAND test_vector_operations)

add_executable(test_gpu_vector LinearAlgebra/gpu_vector.cu)
target_compile_features(test_gpu_vector PRIVATE cxx_std_20)
target_link_libraries(test_gpu_vector GTest::gtest_main GMGPolarLib)
add_test(NAME test_gpu_vector COMMAND test_gpu_vector)

add_executable(test_gpu_vector_operations LinearAlgebra/gpu_vector_operations.cu)
target_compile_features(test_gpu_vector_operations PRIVATE cxx_std_20)
target_link_libraries(test_gpu_vector_operations GTest::gtest_main GMGPolarLib)
add_test(NAME test_gpu_vector_operations COMMAND test_gpu_vector_operations)

add_executable(test_matrix LinearAlgebra/matrix.cpp)
target_compile_features(test_matrix PRIVATE cxx_std_20)
target_link_libraries(test_matrix GTest::gtest_main GMGPolarLib)
add_test(NAME test_matrix COMMAND test_matrix)

# add_executable(test_tridiagonal_solver LinearAlgebra/tridiagonal_solver.cpp)
# target_compile_features(test_tridiagonal_solver PRIVATE cxx_std_17)
# target_link_libraries(test_tridiagonal_solver GMGPolarLib GTest::gtest_main ${MUMPS_LIBRARIES})
# add_test(NAME test_tridiagonal_solver COMMAND test_tridiagonal_solver)

# add_executable(test_cyclic_tridiagonal_solver LinearAlgebra/cyclic_tridiagonal_solver.cpp)
# target_compile_features(test_cyclic_tridiagonal_solver PRIVATE cxx_std_17)
# target_link_libraries(test_cyclic_tridiagonal_solver GMGPolarLib GTest::gtest_main ${MUMPS_LIBRARIES})
# add_test(NAME test_cyclic_tridiagonal_solver COMMAND test_cyclic_tridiagonal_solver)

add_executable(test_polargrid PolarGrid/polargrid.cpp)
target_compile_features(test_polargrid PRIVATE cxx_std_20)
target_link_libraries(test_polargrid GTest::gtest_main GMGPolarLib)
add_test(NAME test_polargrid COMMAND test_polargrid)


add_executable(test_prolongation Interpolation/prolongation.cu)
target_compile_features(test_prolongation PRIVATE cxx_std_20)
target_link_libraries(test_prolongation GTest::gtest_main GMGPolarLib)
add_test(NAME test_prolongation COMMAND test_prolongation)

add_executable(test_restriction Interpolation/restriction.cu)
target_compile_features(test_restriction PRIVATE cxx_std_20)
target_link_libraries(test_restriction GTest::gtest_main GMGPolarLib)
add_test(NAME test_restriction COMMAND test_restriction)

add_executable(test_extrapolated_prolongation Interpolation/extrapolated_prolongation.cu)
target_compile_features(test_extrapolated_prolongation PRIVATE cxx_std_17)
target_link_libraries(test_extrapolated_prolongation GTest::gtest_main GMGPolarLib)
add_test(NAME test_extrapolated_prolongation COMMAND test_extrapolated_prolongation)

add_executable(test_extrapolated_restriction Interpolation/extrapolated_restriction.cu)
target_compile_features(test_extrapolated_restriction PRIVATE cxx_std_20)
target_link_libraries(test_extrapolated_restriction GTest::gtest_main GMGPolarLib)
add_test(NAME test_extrapolated_restriction COMMAND test_extrapolated_restriction)

add_executable(test_fmg_interolation Interpolation/fmg_interpolation.cu)
target_compile_features(test_fmg_interolation PRIVATE cxx_std_20)
target_link_libraries(test_fmg_interolation GTest::gtest_main GMGPolarLib)
add_test(NAME test_fmg_interolation COMMAND test_fmg_interolation)

add_executable(test_injection Interpolation/injection.cu)
target_compile_features(test_injection PRIVATE cxx_std_20)
target_link_libraries(test_injection GTest::gtest_main GMGPolarLib)
add_test(NAME test_injection COMMAND test_injection)

add_executable(test_residual Residual/residual.cu)
target_compile_features(test_residual PRIVATE cxx_std_20)
target_link_libraries(test_residual GTest::gtest_main GMGPolarLib)
add_test(NAME test_residual COMMAND test_residual)

# add_executable(test_directSolver DirectSolver/directSolver.cpp)
# target_compile_features(test_directSolver PRIVATE cxx_std_17)
# target_link_libraries(test_directSolver GMGPolarLib GTest::gtest_main ${MUMPS_LIBRARIES})
# add_test(NAME test_directSolver COMMAND test_directSolver)

add_executable(test_smoother Smoother/smoother.cu)
target_compile_features(test_smoother PRIVATE cxx_std_20)
target_link_libraries(test_smoother GTest::gtest_main GMGPolarLib ${MUMPS_LIBRARIES})
add_test(NAME test_smoother COMMAND test_smoother)

add_executable(test_extrapolated_smoother ExtrapolatedSmoother/extrapolated_smoother.cu)
target_compile_features(test_extrapolated_smoother PRIVATE cxx_std_20)
target_link_libraries(test_extrapolated_smoother GTest::gtest_main GMGPolarLib ${MUMPS_LIBRARIES})
add_test(NAME test_extrapolated_smoother COMMAND test_extrapolated_smoother)








# # Enable CUDA for the project
# enable_language(CUDA)

# # Add the test_gpu_vector executable
# add_executable(test_gpu_vector LinearAlgebra/gpu_vector.cu)

# # Set the required C++ standard
# target_compile_features(test_gpu_vector PRIVATE cxx_std_20)

# # Link with Google Test main library
# target_link_libraries(test_gpu_vector PRIVATE GTest::gtest_main)

# # Add CUDA compile options (optional)
# target_compile_options(test_gpu_vector PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
#     --expt-extended-lambda --expt-relaxed-constexpr>)

# # Link CUDA runtime library
# find_package(CUDAToolkit REQUIRED)
# target_link_libraries(test_gpu_vector PRIVATE CUDA::cudart)

# # Add the test to CTest
# add_test(NAME test_gpu_vector COMMAND test_gpu_vector)





include(GoogleTest)