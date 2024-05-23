
cmake_minimum_required(VERSION 3.12)

# Set the project name
project(GMGPolarTests LANGUAGES CXX)

# Set C++ standard to C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Set the list of test sources
set(TESTSOURCES
    vector.cpp
    operations.cpp
    polargrid.cpp
    # prolongation.cpp
    # restriction.cpp
)

# Add executables for each test source
foreach(TESTSOURCE ${TESTSOURCES})
    get_filename_component(TESTNAME ${TESTSOURCE} NAME_WE)
    add_executable(${TESTNAME} ${TESTSOURCE} ${LEVEL_SOURCES} ${GMG_POLAR_SOURCES} ${INTERPOLATION_SOURCES} ${OPERATOR_SOURCES} ${POLAR_GRID_SOURCES} ${TASK_DISTRIBUTION_SOURCES} ${TEST_SOURCES} ${INPUT_FUNCTIONS_SOURCES})
    target_compile_features(${TESTNAME} PRIVATE cxx_std_20)
    target_link_libraries(${TESTNAME} GTest::gtest_main)
    list(APPEND TESTS ${TESTNAME})
endforeach()

# Enable testing
include(GoogleTest)
gtest_discover_tests(${TESTS})















# # set(TESTSOURCES
# # polargrid.cpp,
# # prolongation.cpp)

# # add_executable(gmgpolar_refac_test${TESTSOURCES})

# add_executable(test_vector vector.cpp)
# target_compile_features(test_vector PRIVATE cxx_std_20)
# target_link_libraries(test_vector GTest::gtest_main)
# list(APPEND test_vector_tests test_vector)

# add_executable(test_operations operations.cpp)
# target_compile_features(test_operations PRIVATE cxx_std_20)
# target_link_libraries(test_operations GTest::gtest_main)
# list(APPEND test_operations_tests test_operations)

# add_executable(test_polargrid polargrid.cpp ${POLAR_GRID_SOURCES})
# target_compile_features(test_polargrid PRIVATE cxx_std_20)
# target_link_libraries(test_polargrid GTest::gtest_main)
# list(APPEND test_polargrid_tests test_polargrid)

# add_executable(test_prolongation prolongation.cpp ${LEVEL_SOURCES} ${GMG_POLAR_SOURCES} ${OPERATOR_SOURCES} ${POLAR_GRID_SOURCES} ${TASK_DISTRIBUTION_SOURCES})
# target_compile_features(test_prolongation PRIVATE cxx_std_20)
# target_link_libraries(test_prolongation GTest::gtest_main)
# list(APPEND test_prolongation_tests test_prolongation)

# # Enable testing
# include(GoogleTest)
# # gtest_discover_tests(test_polargrid)
# # gtest_discover_tests(test_prolongation)













# cmake_minimum_required(VERSION 3.12)
# project(GMGPolarTests)

# set(COMMON_SOURCE
# ../common/equals.h
# ../common/scalar.h
# ../common/space_dimension.h
# ../common/multiindex.h
# ../common/multiindex.cpp
# ../common/point.h
# ../common/point.cpp
# )

# set(POLAR_GRID_SOURCE
# ../PolarGrid/polargrid.hpp
# ../PolarGrid/polargrid.cpp
# ../PolarGrid/anisotropic_division.cpp
# ../PolarGrid/load_write_grid.cpp
# )

# set(HEADER_FILES
# ../include/cmdline.h
# ../include/display_util.hpp
# )

# set(LEVEL_SOURCE
# ../Level/level.hpp
# ../Level/level.cpp
# )

# set(GMG_POLAR_SOURCE
# ../GMGPolar/gmgpolar.hpp
# ../GMGPolar/gmgpolar.cpp
# ../GMGPolar/gmgpolar_parser.cpp
# ../GMGPolar/gmgpolar_mesh.cpp
# )

# set(OPERATOR_SOURCE
# ../Operator/operator.hpp
# ../Operator/operator.cpp
# )

# set(TASKDISTRIBUTION_SOURCE
# ../TaskDistribution/TaskDistribution.hpp
#     ../TaskDistribution/TaskDistribution.cpp
# )

# set(PolarGrid)

# add_executable(test_polargrid polargrid.cpp ${POLAR_GRID_SOURCE} ${COMMON_SOURCE} ${HEADER_FILES} ${LEVEL_SOURCE} ${GMG_POLAR_SOURCE} ${OPERATOR_SOURCE} ${TASKDISTRIBUTION_SOURCE})
# target_compile_features(test_polargrid PRIVATE cxx_std_20)
# target_link_libraries(test_polargrid GTest::gtest_main)
# list(APPEND test_polargrid_tests test_polargrid)


# set(Class)

# add_executable(test_class class.cpp ${POLAR_GRID_SOURCE} ${COMMON_SOURCE} ${HEADER_FILES} ${LEVEL_SOURCE} ${GMG_POLAR_SOURCE} ${OPERATOR_SOURCE} ${TASKDISTRIBUTION_SOURCE})
# target_compile_features(test_class PRIVATE cxx_std_20)
# target_link_libraries(test_class GTest::gtest_main)
# list(APPEND test_class_tests test_class)

