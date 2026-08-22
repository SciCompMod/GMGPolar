# Create the executables. These test specific cases and are not very useful
# without access to the source code so they are not installed with the
# library
add_executable(gmgpolar main.cpp)
target_link_libraries(gmgpolar PRIVATE GMGPolarLib GMGPolarInterface)

# Builds a single-run inner executable (NAME_run) from SOURCE and generates a
# bash wrapper script NAME that owns the loop over configurations, calling the
# runner once per Kokkos::ScopeGuard lifetime.  SCRIPT_BODY is shell code that
# uses $RUNNER (absolute path injected at cmake-generate time) and $OUTFILE.
function(add_scaling_wrapper NAME SOURCE SCRIPT_BODY)
    set(runner_target ${NAME}_run)
    add_executable(${runner_target} ${SOURCE})
    target_link_libraries(${runner_target} PRIVATE GMGPolarLib GMGPolarInterface)

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

add_scaling_wrapper(weak_scaling   weak_scaling.cpp   "${WEAK_SCALING_BODY}")
add_scaling_wrapper(strong_scaling strong_scaling.cpp "${STRONG_SCALING_BODY}")
