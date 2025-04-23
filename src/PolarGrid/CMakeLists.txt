set(POLAR_GRID_SOURCES
    anisotropic_division.cpp
    load_write_grid.cpp
    multiindex.cpp
    point.cpp
    polargrid.cpp
)

add_library(PolarGrid STATIC ${POLAR_GRID_SOURCES})
target_include_directories(PolarGrid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})