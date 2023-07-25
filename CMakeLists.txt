cmake_minimum_required(VERSION 3.16.3)

project(GMGPolar VERSION 1.0.0)

option(GMGPOLAR_BUILD_TESTS "Build GMGPolar unit tests." ON)
option(GMGPOLAR_TEST_COVERAGE "Create GMGPolar code coverage." OFF)

file(GLOB sources_src "./src/*.cpp" "./src/test_cases/*.cpp")

# code coverage analysis
# Note: this only works under linux and with make
# Ninja creates different directory names which do not work together with this scrupt
# as STREQUAL is case-sensitive https://github.com/TriBITSPub/TriBITS/issues/131, also allow DEBUG as accepted input
option(GMGPOLAR_TEST_COVERAGE "Enable GCov coverage analysis (adds a 'coverage' target)" OFF)
mark_as_advanced(GMGPOLAR_TEST_COVERAGE)

if(GMGPOLAR_TEST_COVERAGE)
    message(STATUS "Coverage enabled")
    include(CodeCoverage)
    append_coverage_compiler_flags()
    setup_target_for_coverage_lcov(
        NAME coverage
        EXECUTABLE gmgpolar_tests
        EXCLUDE "${CMAKE_SOURCE_DIR}/tests*" "${CMAKE_SOURCE_DIR}/src/test_cases*" "${CMAKE_BINARY_DIR}/*" "/usr*"
    )
endif()


add_library(GMGPolar ${sources_src})

add_executable(gmgpolar_simulation ./src/main.cpp)

target_include_directories(gmgpolar_simulation PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/include/test_cases )
target_include_directories(GMGPolar PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/include/test_cases )

find_package(OpenMP)

#currently works on GNU compiler
if((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7))

    set(CMAKE_CXX_FLAGS "-O2 -Wall -MMD -MP -Wwrite-strings")
    set(CMAKE_LINKER_FLAGS "-O2 -Wall -MMD -MP -Wwrite-strings")

    if(OPENMP_FOUND)
        string(APPEND CMAKE_CXX_FLAGS " -fopenmp")
        string(APPEND CMAKE_LINKER_FLAGS " -fopenmp")
        
    else()
        message(FATAL_ERROR "OpenMP needed")
    endif()
else()
    message(FATAL_ERROR "Please use GNU compiler or change CMakeLists manually")
endif()

target_link_libraries(gmgpolar_simulation PRIVATE GMGPolar)

include(thirdparty/CMakeLists.txt)
add_subdirectory(tests)

