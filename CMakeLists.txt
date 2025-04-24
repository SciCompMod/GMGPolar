cmake_minimum_required(VERSION 3.12)
project(GMGPolar VERSION 2.0.0 LANGUAGES CXX)

option(GMGPOLAR_BUILD_TESTS "Build GMGPolar unit tests." ON)
option(GMGPOLAR_USE_LIKWID "Use LIKWID to measure code (regions)." OFF)
option(GMGPOLAR_USE_MUMPS "Use MUMPS to solve linear systems." OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic -Wno-unused -Wno-psabi")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -mtune=generic -Wno-psabi")

# code coverage analysis
# Note: this only works under linux and with make
# Ninja creates different directory names which do not work together with this scrupt
# as STREQUAL is case-sensitive https://github.com/TriBITSPub/TriBITS/issues/131, also allow DEBUG as accepted input
option(GMGPOLAR_TEST_COVERAGE "Enable GCov coverage analysis (adds a 'coverage' target)" OFF)

if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    if(GMGPOLAR_TEST_COVERAGE)
        message(STATUS "Coverage enabled")
        include(CodeCoverage)
        append_coverage_compiler_flags()
        setup_target_for_coverage_lcov(
            NAME coverage
            EXECUTABLE tests/gmgpolar_tests
            LCOV_ARGS --ignore-errors gcov,mismatch
            EXCLUDE "${CMAKE_SOURCE_DIR}/tests*" "${CMAKE_SOURCE_DIR}/src/InputFunctions*" "${CMAKE_BINARY_DIR}/*" "/usr*"
        )
    endif()
endif()

include_directories(include)
add_subdirectory(src)

add_executable(gmgpolar src/main.cpp)
add_executable(convergence_order src/convergence_order.cpp)
add_executable(weak_scaling src/weak_scaling.cpp)
add_executable(strong_scaling src/strong_scaling.cpp)

target_link_libraries(gmgpolar PRIVATE GMGPolarLib)
target_link_libraries(convergence_order PRIVATE GMGPolarLib)
target_link_libraries(weak_scaling PRIVATE GMGPolarLib)
target_link_libraries(strong_scaling PRIVATE GMGPolarLib)

if(GMGPOLAR_BUILD_TESTS)
    enable_testing()
    add_subdirectory(third-party)
    add_subdirectory(tests)
endif()
