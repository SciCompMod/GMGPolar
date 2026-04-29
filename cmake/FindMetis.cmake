
if(DEFINED ENV{METIS_DIR} AND NOT METIS_DIR)
    set(METIS_DIR "$ENV{METIS_DIR}" CACHE PATH "METIS installation directory")
endif()

find_path(METIS_INCLUDE_DIR
    NAMES metis.h
    HINTS ${METIS_DIR}/include
)

find_library(METIS_LIBRARY
    NAMES metis
    HINTS ${METIS_DIR}/lib ${METIS_DIR}/lib64
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Metis
    REQUIRED_VARS METIS_INCLUDE_DIR METIS_LIBRARY
)

if(Metis_FOUND AND NOT TARGET metis::metis)
    add_library(metis::metis INTERFACE IMPORTED)
    set_target_properties(metis::metis PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${METIS_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES      "${METIS_LIBRARY}"
    )
endif()
