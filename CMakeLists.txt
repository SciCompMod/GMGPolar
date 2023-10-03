cmake_minimum_required(VERSION 3.16.3)

project(GMGPolar VERSION 1.0.0)

option(GMGPOLAR_BUILD_TESTS "Build GMGPolar unit tests." ON)
option(GMGPOLAR_USE_MUMPS "Use MUMPS to compute matrix factorizations." OFF)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_CXX_FLAGS "")
set(CMAKE_LINKER_FLAGS "")

add_subdirectory(src)

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
            EXCLUDE "${CMAKE_SOURCE_DIR}/tests*" "${CMAKE_SOURCE_DIR}/src/test_cases*" "${CMAKE_BINARY_DIR}/*" "/usr*"
        )
    endif()
endif()


add_library(GMGPolar ${SOURCES_SRC})

add_executable(gmgpolar_simulation ./src/main.cpp)

configure_file(${CMAKE_SOURCE_DIR}/include/config_internal.h.in ${CMAKE_SOURCE_DIR}/include/config_internal.h)

target_include_directories(gmgpolar_simulation PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/include/test_cases )
target_include_directories(GMGPolar PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/include/test_cases )

if(GMGPOLAR_USE_MUMPS)

    set(INC_DIRS
        /home/kueh_mj/.spack/rev.23.05/install/linux-rocky8-zen2/gcc-10.4.0/mumps-5.4.1-fftqkl/include
        /sw/rev/23.05/linux-rocky8-zen2/gcc-10.4.0/metis-5.1.0-akhgsf/include
    )

    set(LIB_DIRS
        /home/kueh_mj/.spack/rev.23.05/install/linux-rocky8-zen2/gcc-10.4.0/mumps-5.4.1-fftqkl/lib
        /sw/rev/23.05/linux-rocky8-zen2/gcc-10.4.0/metis-5.1.0-akhgsf/lib
    )

    include_directories(
        ${INC_DIRS}
    )

    target_link_directories(
        GMGPolar
        PUBLIC
        ${LIB_DIRS}
    )

    set(LIBS
        mpiseq
        dmumps
        mumps_common
        metis
    )

    target_link_libraries(
        GMGPolar
        PUBLIC
        ${LIBS}
    ) 
endif()


find_package(OpenMP)

#currently works on GNU compiler
if((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7))

    string(APPEND CMAKE_CXX_FLAGS " -O2 -Wall -MMD -MP -Wwrite-strings")
    string(APPEND CMAKE_LINKER_FLAGS " -O2 -Wall -MMD -MP -Wwrite-strings")

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

