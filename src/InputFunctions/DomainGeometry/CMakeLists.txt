set(DOMAIN_GEOMETRY_SOURCES
    circularGeometry.cpp
    culhamGeometry.cpp
    czarnyGeometry.cpp
    shafranovGeometry.cpp
)

add_library(InputFunctions_DomainGeometry STATIC ${DOMAIN_GEOMETRY_SOURCES})
target_include_directories(InputFunctions_DomainGeometry PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)
