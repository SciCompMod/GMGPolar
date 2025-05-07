cmake_minimum_required(VERSION 3.12)
project(GMGPolar VERSION 2.0.0 LANGUAGES CXX)

# Options should be defined before they're used
option(GMGPOLAR_BUILD_TESTS "Build GMGPolar unit tests." OFF)
option(GMGPOLAR_USE_LIKWID "Use LIKWID to measure code (regions)." ON)
option(GMGPOLAR_USE_MUMPS "Use MUMPS to solve linear systems." ON)
option(GMGPOLAR_ENABLE_COVERAGE "Enable code coverage reporting (requires GCC/Clang)" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Set build type - must come before compiler flags
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Set compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic -Wno-unused -Wno-psabi")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -mtune=generic -Wno-psabi")

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

add_executable(paper_convergence_order src/paper_convergence_order.cpp)
target_link_libraries(paper_convergence_order PRIVATE GMGPolarLib)
add_executable(paper_weak_scaling src/paper_weak_scaling.cpp)
target_link_libraries(paper_weak_scaling PRIVATE GMGPolarLib)
add_executable(paper_strong_scaling src/paper_strong_scaling.cpp)
target_link_libraries(paper_strong_scaling PRIVATE GMGPolarLib)

if(GMGPOLAR_BUILD_TESTS)
    enable_testing()
    add_subdirectory(third-party)
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
