cmake_minimum_required(VERSION 3.12)

if (CMAKE_COMPILER_IS_GNUCXX AND CMAKE_BUILD_TYPE STREQUAL "Release")
    # If using GNU increase maximum number of instructions considered for inlining
    # in Release mode in this folder and sub-folders
    add_compile_options(--param max-inline-insns-single=1500)
endif()

# Add subdirectories for components
add_subdirectory(InputFunctions)

# Gather all source files
# file(GLOB_RECURSE POLAR_GRID_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/PolarGrid/*.cpp)
set(POLAR_GRID_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/PolarGrid/polargrid.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PolarGrid/anisotropic_division.cpp
)

# file(GLOB_RECURSE GMG_POLAR_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/*.cpp)
set(GMG_POLAR_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/gmgpolar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/setup.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GMGPolar/solver.cpp
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
    ${STENCIL_SOURCES}
    ${INTERPOLATION_SOURCES}
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
if(OpenMP_CXX_FOUND)
    target_link_libraries(GMGPolarLib PUBLIC OpenMP::OpenMP_CXX)
endif()
find_package(Kokkos  4.4.1...<5 QUIET REQUIRED)
target_link_libraries(GMGPolarLib PUBLIC Kokkos::kokkos)


# Handle MUMPS configuration
if(GMGPOLAR_USE_MUMPS)
    find_package(MUMPS COMPONENTS OpenMP METIS)
    if (${MUMPS_FOUND})
      target_link_libraries(GMGPolarLib PUBLIC MUMPS::MUMPS)
      target_compile_definitions(GMGPolarLib PUBLIC GMGPOLAR_USE_MUMPS)
    elseif(DEFINED ENV{MUMPS_DIR})
        set(MUMPS_DIR $ENV{MUMPS_DIR})

        if(EXISTS "${MUMPS_DIR}/include" AND EXISTS "${MUMPS_DIR}/lib")
            set(MUMPS_INCLUDE_DIR "${MUMPS_DIR}/include")
            set(MUMPS_LIBRARY_DIR "${MUMPS_DIR}/lib")

            target_include_directories(GMGPolarLib PUBLIC ${MUMPS_INCLUDE_DIR})
            target_link_directories(GMGPolarLib PUBLIC ${MUMPS_LIBRARY_DIR})

            if (EXISTS "${MUMPS_DIR}/libseq")
              target_include_directories(GMGPolarLib PUBLIC "${MUMPS_DIR}/libseq")
              target_link_directories(GMGPolarLib PUBLIC "${MUMPS_DIR}/libseq")
            endif()

            set(MUMPS_FOUND TRUE)
        endif()

        if(MUMPS_FOUND)
            target_link_libraries(GMGPolarLib PUBLIC
                mumps_common
                smumps
                dmumps
                mpiseq
            )
            target_compile_definitions(GMGPolarLib PUBLIC GMGPOLAR_USE_MUMPS)
        else()
            message(SEND_ERROR "MUMPS not found")
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

        if(METIS_FOUND)
            target_link_libraries(GMGPolarLib PUBLIC metis)
        else()
            message(SEND_ERROR "METIS not found")
        endif()
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
