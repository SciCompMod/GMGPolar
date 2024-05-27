cmake_minimum_required(VERSION 3.12)
project(GMGPolar VERSION 2.0.0 LANGUAGES CXX)

# Set C++ standard to C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Ensure CMake build type is set
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Add compiler flags for different build types
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")


set(METIS_PREFIX_PATH "~/spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/metis-5.1.0-bgoncx22w55soviybggl5ydjakvkm34v")
include_directories(${METIS_PREFIX_PATH}/include)
link_directories(${METIS_PREFIX_PATH}/lib)

set(MUMPS_PREFIX_PATH "~/spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/mumps-5.6.2-m4xrhc4mshzrxmgptzbpult3nbf6qrzk")
include_directories(${MUMPS_PREFIX_PATH}/include)
link_directories(${MUMPS_PREFIX_PATH}/lib)

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

# find_package(LIKWID REQUIRED)

# Include directories
include_directories(include)

# Add source files
file(GLOB_RECURSE LEVEL_SOURCES "src/Level/*.cpp")
file(GLOB_RECURSE GMG_POLAR_SOURCES "src/GMGPolar/*.cpp")
file(GLOB_RECURSE INTERPOLATION_SOURCES "src/Interpolation/*.cpp")
file(GLOB_RECURSE OPERATOR_SOURCES "src/Operator/*.cpp")
file(GLOB_RECURSE POLAR_GRID_SOURCES "src/PolarGrid/*.cpp")
file(GLOB_RECURSE TASK_DISTRIBUTION_SOURCES "src/TaskDistribution/*.cpp")
file(GLOB_RECURSE TEST_SOURCES "src/test_cases/*.cpp")
file(GLOB_RECURSE INPUT_FUNCTIONS_SOURCES "src/InputFunctions/*.cpp")

# Add main.cpp to the list of source files
set(MAIN_SOURCE "src/main.cpp")

# Add executable target
add_executable(gmgpolar 
    ${MAIN_SOURCE}
    ${LEVEL_SOURCES}
    ${GMG_POLAR_SOURCES}
    ${INTERPOLATION_SOURCES}
    ${OPERATOR_SOURCES}
    ${POLAR_GRID_SOURCES}
    ${TASK_DISTRIBUTION_SOURCES}
    ${TEST_SOURCES}
    ${INPUT_FUNCTIONS_SOURCES}
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
















# Optionally, if you have libraries to link
# target_link_libraries(MyExecutable ${SOME_LIBRARY})

# Optionally, if you have additional dependencies
# find_package(SomePackage REQUIRED)
# include_directories(${SomePackage_INCLUDE_DIRS})
# target_link_libraries(MyExecutable ${SomePackage_LIBRARIES})




# set(COMMON_SOURCE
#     common/equals.h
#     common/scalar.h
#     common/space_dimension.h
#     common/multiindex.h
#     common/multiindex.cpp
#     common/point.h
#     common/point.cpp
# )

# set(POLAR_GRID_SOURCE
#     PolarGrid/polargrid.hpp
#     PolarGrid/polargrid.cpp
#     PolarGrid/anisotropic_division.cpp
#     PolarGrid/load_write_grid.cpp
# )

# set(HEADER_FILES
#     include/cmdline.h
#     include/display_util.hpp
# )

# set(LEVEL_SOURCE
#     Level/level.hpp
#     Level/level.cpp
# )

# set(GMG_POLAR_SOURCE
#     GMGPolar/gmgpolar.hpp
#     GMGPolar/gmgpolar.cpp
#     GMGPolar/gmgpolar_parser.cpp
#     GMGPolar/gmgpolar_mesh.cpp
# )

# set(OPERATOR_SOURCE
#     Operator/operator.hpp
#     Operator/operator.cpp
# )

# set(TASKDISTRIBUTION_SOURCE
#     TaskDistribution/TaskDistribution.hpp
#     TaskDistribution/TaskDistribution.cpp
# )


# set(MAIN_SOURCE
#     main.cpp
# )

# Create executable target
# add_executable(gmgpolar_simulation ${MAIN_SOURCE} ${POLAR_GRID_SOURCE} ${COMMON_SOURCE} ${HEADER_FILES} ${LEVEL_SOURCE} ${GMG_POLAR_SOURCE} ${OPERATOR_SOURCE} ${TASKDISTRIBUTION_SOURCE})





# enable_testing()
# add_subdirectory(third-party)
# add_subdirectory(tests)





# # Include GoogleTest
# add_subdirectory(Tests)

# # Set up testing
# enable_testing()

# # Add the test executable
# add_executable(polar_grid_test ${GMG_POLAR_SOURCE})

# # Link the test executable with GoogleTest
# target_link_libraries(polar_grid_test gmgpolar_simulation)

# # Register the test executable
# add_test(NAME MyMathFunctionsTests COMMAND test_my_math_functions)


# enable_testing()
# add_subdirectory(third-party)
# add_subdirectory(Tests)