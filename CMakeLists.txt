cmake_minimum_required(VERSION 3.16.3)

project(GMGPolar VERSION 1.0.0)


file(GLOB sources_src "./src/*.cpp" "./src/test_cases/*.cpp")


add_library(GMGPolar ${sources_src})

add_executable(gmgpolar_simulation ./src/main.cpp)

target_include_directories(gmgpolar_simulation PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/include/test_cases )
target_include_directories(GMGPolar PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/include/test_cases )

find_package(OpenMP)

#currently works on GNU compiler
if((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7))

    set(CMAKE_CXX_FLAGS "-O2 -Wall -MMD -MP -Wwrite-strings")
    set(CMAKE_LINKER_FLAGS "-O2 -Wall -MMD -MP -Wwrite-strings")

    if(OPENMP_FOUND)
        string(APPEND CMAKE_CXX_FLAGS " -fopenmp")
        string(APPEND CMAKE_LINKER_FLAGS " -fopenmp")
        
    else()
        message(FATAL_ERROR "OpenMP needed")
    endif()
else()
    message(FATAL_ERROR "Please use GNU compiler or change CMakeLists manually")
endif()

target_link_libraries(gmgpolar_simulation PRIVATE GMGPolar)

include(thirdparty/CMakeLists.txt)
add_subdirectory(tests)

