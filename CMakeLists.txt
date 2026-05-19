cmake_minimum_required(VERSION 3.12)

# Options should be defined before they're used
option(GMGPOLAR_BUILD_TESTS "Build GMGPolar unit tests." ON)
option(GMGPOLAR_USE_LIKWID "Use LIKWID to measure code (regions)." OFF)
option(GMGPOLAR_USE_MUMPS "Use MUMPS to solve linear systems." OFF)
option(GMGPOLAR_ENABLE_COVERAGE "Enable code coverage reporting (requires GCC/Clang)" OFF)

if (${GMGPOLAR_USE_MUMPS})
    # MUMPS does not automatically provide Fortran libraries
    project(GMGPolar VERSION 2.0.0 LANGUAGES CXX Fortran)
else()
    project(GMGPolar VERSION 2.0.0 LANGUAGES CXX)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Set build type - must come before compiler flags
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
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

find_package(Kokkos  4.4.1...<5.1 QUIET REQUIRED)
include_directories(include)
add_subdirectory(src)

add_executable(gmgpolar src/main.cpp)
add_executable(convergence_order src/convergence_order.cpp)
target_link_libraries(gmgpolar PRIVATE GMGPolarLib)
target_link_libraries(convergence_order PRIVATE GMGPolarLib)

# Builds a single-run inner executable (NAME_run) from SOURCE and generates a
# bash wrapper script NAME that owns the loop over configurations, calling the
# runner once per Kokkos::ScopeGuard lifetime.  SCRIPT_BODY is shell code that
# uses $RUNNER (absolute path injected at cmake-generate time) and $OUTFILE.
function(add_scaling_wrapper NAME SOURCE SCRIPT_BODY)
    set(runner_target ${NAME}_run)
    add_executable(${runner_target} ${SOURCE})
    target_link_libraries(${runner_target} PRIVATE GMGPolarLib)

    set(script "${CMAKE_BINARY_DIR}/${NAME}")
    file(GENERATE
        OUTPUT "${script}"
        CONTENT "#!/bin/bash\nset -euo pipefail\nRUNNER=$<TARGET_FILE:${runner_target}>\n${SCRIPT_BODY}"
    )
    add_custom_target(${NAME} ALL
        COMMAND chmod 755 "${script}"
        COMMENT "Preparing ${NAME} wrapper script"
    )
endfunction()

set(WEAK_SCALING_BODY
[=[OUTFILE="weak_scaling_results.csv"
echo "Threads,DivideBy2,nr,ntheta,geometry,stencil_method,cacheDensityProfileCoefficients,cacheDomainGeometry,FMG,extrapolation_int,TotalTime,t_setup_total,t_setup_createLevels,t_setup_smoother,t_setup_directSolver,t_solve_total,t_solve_initial_approximation,t_solve_multigrid_iterations,t_check_convergence,t_check_exact_error,t_avg_MGC_total,t_avg_MGC_preSmoothing,t_avg_MGC_postSmoothing,t_avg_MGC_residual,t_avg_MGC_directSolver" > "$OUTFILE"
"$RUNNER" "$OUTFILE" 5 --kokkos-threads=1
"$RUNNER" "$OUTFILE" 6 --kokkos-threads=4
"$RUNNER" "$OUTFILE" 7 --kokkos-threads=16
"$RUNNER" "$OUTFILE" 8 --kokkos-threads=56
echo "Results written to $OUTFILE"
]=])

set(STRONG_SCALING_BODY
[=[OUTFILE="strong_scaling_results.csv"
echo "Threads,DivideBy2,nr,ntheta,geometry,stencil_method,cacheDensityProfileCoefficients,cacheDomainGeometry,FMG,extrapolation_int,TotalTime,t_setup_total,t_setup_createLevels,t_setup_smoother,t_setup_directSolver,t_solve_total,t_solve_initial_approximation,t_solve_multigrid_iterations,t_check_convergence,t_check_exact_error,t_avg_MGC_total,t_avg_MGC_preSmoothing,t_avg_MGC_postSmoothing,t_avg_MGC_residual,t_avg_MGC_directSolver" > "$OUTFILE"
for threads in $(seq 1 56); do
    "$RUNNER" "$OUTFILE" --kokkos-threads=$threads
done
echo "Results written to $OUTFILE"
]=])

add_scaling_wrapper(weak_scaling   src/weak_scaling.cpp   "${WEAK_SCALING_BODY}")
add_scaling_wrapper(strong_scaling src/strong_scaling.cpp "${STRONG_SCALING_BODY}")

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
