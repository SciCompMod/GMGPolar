cmake_minimum_required(VERSION 3.12)

# Add subdirectories for components
add_subdirectory(InputFunctions)

# Gather all source files
# file(GLOB_RECURSE POLAR_GRID_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/PolarGrid/*.cpp)
set(POLAR_GRID_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/PolarGrid/anisotropic_division.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PolarGrid/multiindex.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PolarGrid/point.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PolarGrid/polargrid.cpp
)

# file(GLOB_RECURSE GMG_POLAR_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/*.cpp)
set(GMG_POLAR_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/build_rhs_f.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/gmgpolar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/level_interpolation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/setup.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/solver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/writeToVTK.cpp
)

# file(GLOB_RECURSE MULTIGRID_METHODS_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/MultigridMethods/*.cpp)
set(MULTIGRID_METHODS_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/MultigridMethods/implicitly_extrapolated_multigrid_F_Cycle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/MultigridMethods/implicitly_extrapolated_multigrid_V_Cycle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/MultigridMethods/implicitly_extrapolated_multigrid_W_Cycle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/MultigridMethods/multigrid_F_Cycle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/MultigridMethods/multigrid_V_Cycle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/MultigridMethods/multigrid_W_Cycle.cpp
)

#file(GLOB_RECURSE LEVEL_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Level/*.cpp)
set(LEVEL_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/Level/level.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Level/levelCache.cpp
)

# file(GLOB_RECURSE STENCIL_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Stencil/*.cpp)
set(STENCIL_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/Stencil/stencil.cpp
)

# file(GLOB_RECURSE INTERPOLATION_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Interpolation/*.cpp)
set(INTERPOLATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/Interpolation/extrapolated_prolongation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Interpolation/extrapolated_restriction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Interpolation/fmg_interpolation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Interpolation/injection.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Interpolation/interpolation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Interpolation/prolongation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Interpolation/restriction.cpp
)

# file(GLOB_RECURSE DIRECT_SOLVER_SOURCES 
#     ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/*.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-COO-MUMPS-Give/*.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-COO-MUMPS-Take/*.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-CSR-LU-Give/*.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-CSR-LU-Take/*.cpp)
set(DIRECT_SOLVER_SOURCES
    # Main DirectSolver files
    ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/directSolver.cpp
    
    # DirectSolverGive
    ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-COO-MUMPS-Give/applySymmetryShift.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-COO-MUMPS-Give/buildSolverMatrix.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-COO-MUMPS-Give/directSolverGive.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-COO-MUMPS-Give/initializeMumps.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-COO-MUMPS-Give/matrixStencil.cpp
    
    # DirectSolverGiveCustomLU
    ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-CSR-LU-Give/buildSolverMatrix.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-CSR-LU-Give/directSolverGive.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-CSR-LU-Give/matrixStencil.cpp
    
    # DirectSolverTake
    ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-COO-MUMPS-Take/applySymmetryShift.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-COO-MUMPS-Take/buildSolverMatrix.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-COO-MUMPS-Take/directSolverTake.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-COO-MUMPS-Take/initializeMumps.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-COO-MUMPS-Take/matrixStencil.cpp
    
    # DirectSolverTakeCustomLU
    ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-CSR-LU-Take/buildSolverMatrix.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-CSR-LU-Take/directSolverTake.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DirectSolver/DirectSolver-CSR-LU-Take/matrixStencil.cpp
)

# file(GLOB_RECURSE RESIDUAL_SOURCES 
#     ${CMAKE_CURRENT_SOURCE_DIR}/Residual/*.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/Residual/ResidualGive/*.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/Residual/ResidualTake/*.cpp)
set(RESIDUAL_SOURCES
    # Main Residual files
    ${CMAKE_CURRENT_SOURCE_DIR}/Residual/residual.cpp
    
    # ResidualGive
    ${CMAKE_CURRENT_SOURCE_DIR}/Residual/ResidualGive/applyAGive.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Residual/ResidualGive/residualGive.cpp
    
    # ResidualTake
    ${CMAKE_CURRENT_SOURCE_DIR}/Residual/ResidualTake/applyResidualTake.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Residual/ResidualTake/residualTake.cpp
)

# file(GLOB_RECURSE SMOOTHER_SOURCES 
#     ${CMAKE_CURRENT_SOURCE_DIR}/Smoother/*.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/Smoother/SmootherGive/*.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/Smoother/SmootherTake/*.cpp)
set(SMOOTHER_SOURCES
    # Main Smoother files
    ${CMAKE_CURRENT_SOURCE_DIR}/Smoother/smoother.cpp
    
    # SmootherGive
    ${CMAKE_CURRENT_SOURCE_DIR}/Smoother/SmootherGive/buildMatrix.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Smoother/SmootherGive/initializeMumps.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Smoother/SmootherGive/matrixStencil.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Smoother/SmootherGive/smootherGive.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Smoother/SmootherGive/smootherSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Smoother/SmootherGive/task_parallelization.cpp
    
    # SmootherTake
    ${CMAKE_CURRENT_SOURCE_DIR}/Smoother/SmootherTake/buildMatrix.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Smoother/SmootherTake/initializeMumps.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Smoother/SmootherTake/matrixStencil.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Smoother/SmootherTake/smootherSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Smoother/SmootherTake/smootherTake.cpp
)

# file(GLOB_RECURSE EXTRAPOLATED_SMOOTHER_SOURCES 
#     ${CMAKE_CURRENT_SOURCE_DIR}/ExtrapolatedSmoother/*.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/ExtrapolatedSmoother/ExtrapolatedSmootherGive/*.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/ExtrapolatedSmoother/ExtrapolatedSmootherTake/*.cpp)
set(EXTRAPOLATED_SMOOTHER_SOURCES
    # Main ExtrapolatedSmoother files
    ${CMAKE_CURRENT_SOURCE_DIR}/ExtrapolatedSmoother/extrapolatedSmoother.cpp
    
    # ExtrapolatedSmootherGive
    ${CMAKE_CURRENT_SOURCE_DIR}/ExtrapolatedSmoother/ExtrapolatedSmootherGive/buildAscMatrices.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExtrapolatedSmoother/ExtrapolatedSmootherGive/extrapolatedSmootherGive.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExtrapolatedSmoother/ExtrapolatedSmootherGive/initializeMumps.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExtrapolatedSmoother/ExtrapolatedSmootherGive/smootherSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExtrapolatedSmoother/ExtrapolatedSmootherGive/smootherStencil.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExtrapolatedSmoother/ExtrapolatedSmootherGive/task_parallelization.cpp
    
    # ExtrapolatedSmootherTake
    ${CMAKE_CURRENT_SOURCE_DIR}/ExtrapolatedSmoother/ExtrapolatedSmootherTake/buildAscMatrices.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExtrapolatedSmoother/ExtrapolatedSmootherTake/extrapolatedSmootherTake.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExtrapolatedSmoother/ExtrapolatedSmootherTake/initializeMumps.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExtrapolatedSmoother/ExtrapolatedSmootherTake/smootherSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExtrapolatedSmoother/ExtrapolatedSmootherTake/smootherStencil.cpp
)

# Gather all source files
# file(GLOB_RECURSE CONFIG_PARSER_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/ConfigParser/*.cpp)
set(CONFIG_PARSER_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/ConfigParser/config_parser.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ConfigParser/select_test_case.cpp
)

# Create the main library
add_library(GMGPolarLib STATIC
    ${POLAR_GRID_SOURCES}
    ${GMG_POLAR_SOURCES}
    ${MULTIGRID_METHODS_SOURCES}
    ${LEVEL_SOURCES}
    ${STENCIL_SOURCES}
    ${INTERPOLATION_SOURCES}
    ${DIRECT_SOLVER_SOURCES}
    ${RESIDUAL_SOURCES}
    ${SMOOTHER_SOURCES}
    ${EXTRAPOLATED_SMOOTHER_SOURCES}
    ${CONFIG_PARSER_SOURCES}
)

# Basic library configuration
target_include_directories(GMGPolarLib PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(GMGPolarLib PUBLIC 
    InputFunctions
)

# Handle OpenMP
find_package(OpenMP REQUIRED)
#if(OpenMP_CXX_FOUND)
target_link_libraries(GMGPolarLib PUBLIC OpenMP::OpenMP_CXX)
#endif()
find_package(Kokkos 4.4.1 QUIET)
target_link_libraries(GMGPolarLib PUBLIC Kokkos::kokkos)


# Handle MUMPS configuration
if(GMGPOLAR_USE_MUMPS)
    if(DEFINED ENV{MUMPS_DIR})
        set(MUMPS_DIR $ENV{MUMPS_DIR})
        
        if(EXISTS "${MUMPS_DIR}/include" AND EXISTS "${MUMPS_DIR}/lib")
            set(MUMPS_INCLUDE_DIR "${MUMPS_DIR}/include")
            set(MUMPS_LIBRARY_DIR "${MUMPS_DIR}/lib")
            
            target_include_directories(GMGPolarLib PUBLIC ${MUMPS_INCLUDE_DIR})
            target_link_directories(GMGPolarLib PUBLIC ${MUMPS_LIBRARY_DIR})
            set(MUMPS_FOUND TRUE)
        endif()
    endif()

    if(DEFINED ENV{METIS_DIR})
        set(METIS_DIR $ENV{METIS_DIR})

        if(EXISTS "${METIS_DIR}/include" AND EXISTS "${METIS_DIR}/lib")
            set(METIS_INCLUDE_DIR "${METIS_DIR}/include")
            set(METIS_LIBRARY_DIR "${METIS_DIR}/lib")
            
            target_include_directories(GMGPolarLib PUBLIC ${METIS_INCLUDE_DIR})
            target_link_directories(GMGPolarLib PUBLIC ${METIS_LIBRARY_DIR})
            set(METIS_FOUND TRUE)
        endif()
    endif()

    if(MUMPS_FOUND)
        target_link_libraries(GMGPolarLib PUBLIC
            mumps_common
            smumps
            dmumps
            mpiseq
        )
        target_compile_definitions(GMGPolarLib PUBLIC GMGPOLAR_USE_MUMPS)
    endif()

    if(METIS_FOUND)
        target_link_libraries(GMGPolarLib PUBLIC metis)
    endif()
endif()

# Handle LIKWID configuration
if(GMGPOLAR_USE_LIKWID)
    if(NOT DEFINED LIKWID_DIR)
        if(DEFINED ENV{LIKWID_DIR})
            set(LIKWID_DIR $ENV{LIKWID_DIR} CACHE PATH "Path to Likwid installation")
        endif()
    endif()

    if(LIKWID_DIR)
        find_path(LIKWID_INCLUDE_DIR
            NAMES likwid.h
            PATHS "${LIKWID_DIR}/include"
            NO_DEFAULT_PATH
            REQUIRED
        )
        find_library(LIKWID_LIBRARY
            NAMES likwid
            PATHS "${LIKWID_DIR}/lib" "${LIKWID_DIR}/lib64"
            NO_DEFAULT_PATH
            REQUIRED
        )
    else()
        find_package(PkgConfig QUIET)
        if(PKG_CONFIG_FOUND)
            pkg_check_modules(LIKWID QUIET likwid)
        endif()
    endif()

    if(LIKWID_INCLUDE_DIR AND LIKWID_LIBRARY)
        target_include_directories(GMGPolarLib PUBLIC ${LIKWID_INCLUDE_DIR})
        target_link_libraries(GMGPolarLib PUBLIC ${LIKWID_LIBRARY})
        target_compile_definitions(GMGPolarLib PUBLIC
            LIKWID_PERFMON
            GMGPOLAR_USE_LIKWID
        )
    endif()
endif()
