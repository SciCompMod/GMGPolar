cmake_minimum_required(VERSION 3.12)

# Ensure custom cmake modules can be found
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Options should be defined before they're used
option(GMGPOLAR_BUILD_TESTS "Build GMGPolar unit tests." ON)
option(GMGPOLAR_USE_LIKWID "Use LIKWID to measure code (regions)." OFF)
option(GMGPOLAR_USE_MUMPS "Use MUMPS to solve linear systems." OFF)
option(GMGPOLAR_ENABLE_COVERAGE "Enable code coverage reporting (requires GCC/Clang)" OFF)

project(GMGPolar VERSION 2.1.0 LANGUAGES CXX)
# MUMPS does not automatically provide Fortran libraries
if (${GMGPOLAR_USE_MUMPS})
    enable_language(Fortran)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Set build type - must come before compiler flags
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Set compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic -Wno-unused -Wno-psabi -Wfloat-conversion")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -mtune=generic -Wno-psabi")

if (CMAKE_COMPILER_IS_GNUCXX AND CMAKE_BUILD_TYPE STREQUAL "Release")
    # If using GNU increase maximum number of instructions considered for inlining
    # in Release mode in this folder and sub-folders
    add_compile_options(--param max-inline-insns-single=1500)
endif()

# Set coverage compiler flags - must come before any targets are defined
if(GMGPOLAR_ENABLE_COVERAGE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        message(STATUS "Enabling code coverage flags")
        # Use generator expressions to apply only to specific configurations
        add_compile_options($<$<CONFIG:Debug>:--coverage>)
        add_link_options($<$<CONFIG:Debug>:--coverage>)
        # Force Debug build when coverage is enabled
        if(CMAKE_BUILD_TYPE STREQUAL "Release")
            message(STATUS "Forcing Debug build for coverage")
            set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
        endif()
    else()
        message(WARNING "Code coverage requires GCC or Clang. Current compiler: ${CMAKE_CXX_COMPILER_ID}")
    endif()
endif()

find_package(OpenMP REQUIRED COMPONENTS CXX)
find_package(Kokkos 4.4.1...<6 QUIET REQUIRED)

if(GMGPOLAR_USE_MUMPS)
    find_package(MUMPS REQUIRED COMPONENTS OpenMP METIS)
endif()

if(GMGPOLAR_USE_LIKWID)
    find_package(LIKWID REQUIRED)
endif()

set(GMGPOLAR_INCLUDE_DIRS
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)

include(CMakePackageConfigHelpers)

configure_package_config_file(
    cmake/GMGPolarConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/GMGPolarConfig.cmake
    INSTALL_DESTINATION lib/cmake/GMGPolar
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/GMGPolarConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/GMGPolarConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/GMGPolarConfigVersion.cmake
    cmake/FindMUMPS.cmake
    cmake/FindLIKWID.cmake
    cmake/FindMetis.cmake
    DESTINATION lib/cmake/GMGPolar
)

install(EXPORT GMGPolarTargets
    NAMESPACE GMGPolar::
    DESTINATION lib/cmake/GMGPolar
)

add_subdirectory(src)

add_executable(gmgpolar src/main.cpp)
add_executable(convergence_order src/convergence_order.cpp)
add_executable(weak_scaling src/weak_scaling.cpp)
add_executable(strong_scaling src/strong_scaling.cpp)

target_link_libraries(gmgpolar PRIVATE GMGPolarLib GMGPolarInterface)
target_link_libraries(convergence_order PRIVATE GMGPolarLib GMGPolarInterface)
target_link_libraries(weak_scaling PRIVATE GMGPolarLib GMGPolarInterface)
target_link_libraries(strong_scaling PRIVATE GMGPolarLib GMGPolarInterface)

if(GMGPOLAR_BUILD_TESTS)
    enable_testing()
    find_package(GTest 1.17 QUIET)
    if (NOT GTest_FOUND)
        add_subdirectory(third-party)
    endif()
    add_subdirectory(tests)

    # Add coverage target - moved after test configuration
    if(GMGPOLAR_ENABLE_COVERAGE)
        find_program(LCOV_PATH lcov)
        find_program(GENHTML_PATH genhtml)

        if(LCOV_PATH AND GENHTML_PATH)
        add_custom_target(coverage
            # Reset counters
            COMMAND ${LCOV_PATH} --directory ${CMAKE_BINARY_DIR} --zerocounters
            # Run tests
            COMMAND ctest --test-dir ${CMAKE_BINARY_DIR} || true
            # Capture coverage data
            COMMAND ${LCOV_PATH} --directory ${CMAKE_BINARY_DIR} --capture 
                --output-file ${CMAKE_BINARY_DIR}/coverage.info
                --ignore-errors mismatch,unused
                --rc geninfo_unexecuted_blocks=1
            # Filter out system and unwanted directories
            COMMAND ${LCOV_PATH} --remove ${CMAKE_BINARY_DIR}/coverage.info
                '/usr/*'
                '*/tests/*'
                '*/third-party/*'
                '*/_deps/googletest-src/*'
                '*/include/InputFunctions/*'
                '*/src/InputFunctions/*'
                '*/include/ConfigParser/cmdline.h'
                --output-file ${CMAKE_BINARY_DIR}/coverage-filtered.info
                --ignore-errors unused
            # Generate HTML report
            COMMAND ${GENHTML_PATH} ${CMAKE_BINARY_DIR}/coverage-filtered.info
                --output-directory ${CMAKE_BINARY_DIR}/coverage-report
                --ignore-errors source
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMENT "Generating code coverage report"
        )
        endif()
    endif()
endif()
