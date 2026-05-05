
# Stage 1: prefer an installed MUMPS CMake config (supports COMPONENTS properly)
find_package(MUMPS CONFIG QUIET COMPONENTS ${MUMPS_FIND_COMPONENTS})
if(MUMPS_FOUND)
    return()
endif()

# Stage 2: manual discovery via MUMPS_DIR hint
if(DEFINED ENV{MUMPS_DIR} AND NOT MUMPS_DIR)
    set(MUMPS_DIR "$ENV{MUMPS_DIR}" CACHE PATH "MUMPS installation directory")
endif()

find_path(MUMPS_INCLUDE_DIR
    NAMES dmumps_c.h
    HINTS ${MUMPS_DIR}/include
)

foreach(_lib dmumps smumps mumps_common)
    find_library(MUMPS_${_lib}_LIBRARY
        NAMES ${_lib}
        HINTS ${MUMPS_DIR}/lib ${MUMPS_DIR}/lib64
    )
    list(APPEND _MUMPS_REQUIRED_VARS MUMPS_${_lib}_LIBRARY)
endforeach()

# mpiseq is the sequential MPI stub — only present in sequential builds
find_library(MUMPS_mpiseq_LIBRARY
    NAMES mpiseq
    HINTS ${MUMPS_DIR}/lib ${MUMPS_DIR}/lib64
          ${MUMPS_DIR}/libseq ${MUMPS_DIR}/lib/SEQ
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS
    REQUIRED_VARS MUMPS_INCLUDE_DIR ${_MUMPS_REQUIRED_VARS}
)

if(MUMPS_FOUND AND NOT TARGET MUMPS::MUMPS)
    find_package(Metis REQUIRED)

    set(_mumps_libs
        ${MUMPS_dmumps_LIBRARY}
        ${MUMPS_smumps_LIBRARY}
        ${MUMPS_mumps_common_LIBRARY}
        metis::metis
    )
    if(MUMPS_mpiseq_LIBRARY)
        list(APPEND _mumps_libs ${MUMPS_mpiseq_LIBRARY})
    endif()

    add_library(MUMPS::MUMPS INTERFACE IMPORTED)
    set_target_properties(MUMPS::MUMPS PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES      "${_mumps_libs}"
    )
endif()
